#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Aug-9, latest modified on 2021-Aug-9
# Output the final electrostatics energy of a given ATP in protein+DNA+ligand system
# Example in Linux: python this_code.py trajectory.pdb

import simtk.openmm.app
import simtk.openmm
import simtk.unit as unit
import sys
import scipy.spatial.distance as sdist
import argparse
import itertools
import configparser
import numpy as np
import math
import os
import pandas
import pdbfixer

__author__ = 'Shikai Jin'
__version__ = '0.0.1'
__date__ = "Aug-9-2021"

_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor
_df = 1 * unit.angstrom / unit.nanometer  # distance scaling factor
_af = 1 * unit.degree / unit.radian  # angle scaling factor

def get_file_len(filename):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_pdbfile_line_atoms(filename):
    file_len = get_file_len(filename)

    with open(filename, 'r') as fopen:
        lines = fopen.readlines()
    n_atoms = 0
    n_lines = 0
    #print(len(lines))
    for line in lines:
        n_lines += 1
        # print line
        if line.split()[0] == "END" or line.split()[0] == "ENDMDL":
            break
        if line[12:16].strip() == "CA":
            n_atoms += 1

    n_snapshot = file_len / n_lines
    return n_snapshot, n_lines, n_atoms

def getline_range(filename, line1, line2):
    assert (line1 <= line2)
    nline = line2 - line1 + 1

    stdout = []
    with open(filename, "r") as text_file:
        for line in itertools.islice(text_file, line1-1, line2):
            stdout.append(line)
    return stdout

def get_pdbfile_i(pdbfile, i, n_lines):
    line_start = 1 + n_lines * i
    line_end = n_lines * (i + 1)
    pdb_part = getline_range(pdbfile, line_start, line_end)
    return pdb_part

def vector(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]

def vabs(a):
    return math.sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))

def parseConfigTable(config_section):
    """ Parses a section of the configuration file as a table.
        This function is used to parse the 3SPN2.conf file"""

    def readData(config_section, a):
        """Filters comments and returns values as a list"""
        temp = config_section.get(a).split('#')[0].split()
        l = []
        for val in temp:
            val = val.strip()
            try:
                x = int(val)
                l += [x]
            except ValueError:
                try:
                    y = float(val)
                    l += [y]
                except ValueError:
                    l += [val]
        return l

    data = []
    for a in config_section:
        if a == 'name':
            columns = readData(config_section, a)
        elif len(a) > 3 and a[:3] == 'row':
            data += [readData(config_section, a)]
        else:
            print(f'Unexpected row {readData(config_section, a)}')
    return pandas.DataFrame(data, columns=columns)


def fixPDB(pdb_file):
    """Uses the pdbfixer library to fix a pdb file, replacing non standard residues, removing
    hetero-atoms and adding missing hydrogens. The input is a pdb file location,
    the output is a fixer object, which is a pdb in the openawsem format."""
    fixer = pdbfixer.PDBFixer(filename=pdb_file, )
    # fixer.findMissingResidues()
    # chains = list(fixer.topology.chains())
    # keys = fixer.missingResidues.keys()
    # for key in list(keys):
    #     chain_tmp = chains[key[0]]
    #     if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
    #         del fixer.missingResidues[key]

    # fixer.findNonstandardResidues()
    # fixer.replaceNonstandardResidues()
    # fixer.removeHeterogens(keepWater=False)
    # fixer.findMissingAtoms()
    # fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)
    return fixer


def fromCoarsePDB(pdb_file, dna_type='B_curved', temp_name='temp'):
    """Initializes a DNA object from a pdb file containing the Coarse Grained atoms"""

    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),
                    serial=int(line[6:11]),
                    name=str(line[12:16]).strip(),
                    altLoc=str(line[16:17]),
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[21:22]),
                    resSeq=int(line[22:26]),
                    iCode=str(line[26:27]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=[0.0 if line[54:60].strip()=='' else float(line[54:60])],
                    tempFactor=[0.0 if line[60:66].strip()=='' else float(line[60:66])],
                    element=str(line[76:78]).strip(),
                    charge=str(line[78:80]).strip())

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                lines += [pdb_line(line)]
    pdb_atoms = pandas.DataFrame(lines)
    self_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                            'resname', 'chainID', 'resSeq', 'iCode',
                            'x', 'y', 'z', 'occupancy', 'tempFactor',
                            'element', 'charge']]
    # print(self.atoms.columns)
    # self.atoms.loc[:, 'chain'] = self.atoms['chainID']
    # self.atoms.loc[:, 'residue'] = self.atoms['resSeq']
    self_atoms.loc[:, 'type'] = self_atoms['name']
    # print(self.atoms.columns)
    # Initialize the system from the pdb
    return self_atoms

def pdb2table(pdb):
    """ Parses a pdb in the openmm format and outputs a table that contains all the information
    on a pdb file """
    cols = ['recname', 'serial', 'name', 'altLoc',
            'resname', 'chainID', 'resSeq', 'iCode',
            'x', 'y', 'z', 'occupancy', 'tempFactor',
            'element', 'charge']
    data = []
    for atom, pos in zip(pdb.topology.atoms(), pdb.positions):
        residue = atom.residue
        chain = residue.chain
        pos = pos.value_in_unit(simtk.unit.angstrom)
        data += [dict(zip(cols, ['ATOM', int(atom.id), atom.name, '',
                                 residue.name, chain.id, int(residue.id), '',
                                 pos[0], pos[1], pos[2], 0, 0,
                                 atom.element.symbol, '']))]
    atom_list = pandas.DataFrame(data)
    atom_list = atom_list[cols]
    atom_list.index = atom_list['serial']
    return atom_list

def CoarseGrain(pdb_table):
    """ Selects DNA atoms from a pdb table and returns a table containing only the coarse-grained atoms for 3SPN2"""
    masses = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "P": 30.973762, }
    CG = {"O5\'": 'P', "C5\'": 'S', "C4\'": 'S', "O4\'": 'S', "C3\'": 'S', "O3\'": 'P',
            "C2\'": 'S', "C1\'": 'S', "O5*": 'P', "C5*": 'S', "C4*": 'S', "O4*": 'S',
            "C3*": 'S', "O3*": 'P', "C2*": 'S', "C1*": 'S', "N1": 'B', "C2": 'B', "O2": 'B',
            "N2": 'B', "N3": 'B', "C4": 'B', "N4": 'B', "C5": 'B', "C6": 'B', "N9": 'B',
            "C8": 'B', "O6": 'B', "N7": 'B', "N6": 'B', "O4": 'B', "C7": 'B', "P": 'P',
            "OP1": 'P', "OP2": 'P', "O1P": 'P', "O2P": 'P', "OP3": 'P', "HO5'": 'P',
            "H5'": 'S', "H5''": 'S', "H4'": 'S', "H3'": 'S', "H2'": 'S', "H2''": 'S',
            "H1'": 'S', "H8": 'B', "H61": 'B', "H62": 'B', 'H2': 'B', 'H1': 'B', 'H21': 'B',
            'H22': 'B', 'H3': 'B', 'H71': 'B', 'H72': 'B', 'H73': 'B', 'H6': 'B', 'H41': 'B',
            'H42': 'B', 'H5': 'B', "HO3'": 'P'}
    cols = ['recname', 'serial', 'name', 'altLoc',
            'resname', 'chainID', 'resSeq', 'iCode',
            'x', 'y', 'z', 'occupancy', 'tempFactor',
            'element', 'charge', 'type']
    temp = pdb_table.copy()

    # Select DNA residues
    temp = temp[temp['resname'].isin(['DA', 'DT', 'DG', 'DC'])]

    # Group the atoms by sugar, phosphate or base
    temp['group'] = temp['name'].replace(CG)
    temp = temp[temp['group'].isin(['P', 'S', 'B'])]

    # Move the O3' to the next residue
    for c in temp['chainID'].unique():
        sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"]
        temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"] = list(sel)[1:] + [-1]
        sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"]
        temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"] = list(sel)[1:] + [-1]
    temp = temp[temp['resSeq'] > 0]

    # Calculate center of mass
    temp['element']=temp['element'].str.strip()
    temp['mass'] = temp.element.replace(masses).astype(float)
    temp[['x', 'y', 'z']] = (temp[['x', 'y', 'z']].T * temp['mass']).T[['x', 'y', 'z']]
    temp = temp[temp['element'] != 'H']  # Exclude hydrogens
    Coarse = temp.groupby(['chainID', 'resSeq', 'resname', 'group']).sum().reset_index()
    Coarse[['x', 'y', 'z']] = (Coarse[['x', 'y', 'z']].T / Coarse['mass']).T[['x', 'y', 'z']]

    # Set pdb columns
    Coarse['recname'] = 'ATOM'
    Coarse['name'] = Coarse['group']
    Coarse['altLoc'] = ''
    Coarse['iCode'] = ''
    Coarse['charge'] = ''
    # Change name of base to real base
    mask = (Coarse.name == 'B')
    Coarse.loc[mask, 'name'] = Coarse[mask].resname.str[-1]  # takes last letter from the residue name
    Coarse['type'] = Coarse['name']
    # Set element (depends on base)
    Coarse['element'] = Coarse['name'].replace({'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'G': 'C', 'C': 'O'})
    # Remove P from the beggining
    drop_list = []
    for chain in Coarse.chainID.unique():
        sel = Coarse[Coarse.chainID == chain]
        drop_list += list(sel[(sel.resSeq == sel.resSeq.min()) & sel['name'].isin(['P'])].index)
    Coarse = Coarse.drop(drop_list)
    # Renumber
    Coarse.index = range(len(Coarse))
    Coarse['serial'] = Coarse.index
    return Coarse[cols]

def main():
    pdb_file = sys.argv[1]
    traj_file = sys.argv[2]
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    configuration_file = f'{__location__}/open3SPN2/3SPN2.conf'
    config = configparser.ConfigParser()
    config.read(configuration_file)

    self_config = {}
    for c in config.sections():
        self_config.update({c: parseConfigTable(config[c])})

    bond_definition = self_config['Bonds']
    # print(bond_definition.DNA)
    bond_types = bond_definition[bond_definition['DNA'] == 'B_curved']
    # print(bond_types)

    # Now parse PDB files
    # pdb = fixPDB(pdb_file)
    # pdb_table = pdb2table(pdb)

    # self_atoms = CoarseGrain(pdb_table)
    self_atoms = fromCoarsePDB(pdb_file)
    traj_atoms = fromCoarsePDB(traj_file)
    # print(traj_atoms.iloc[10071])

    # print(pdb_table)

    index = {}
    cr_list = set()  # Chain residue list
    for i, atom in self_atoms.iterrows():
        index.update({(atom['chainID'], atom['resSeq'], atom['name']): i})
        cr_list.update([(atom['chainID'], atom['resSeq'])])
    cr_list = list(cr_list)
    cr_list.sort()

    data = []

    for i, ftype in bond_types.iterrows():
        # print(bond_type)
        ai = ftype['i']
        aj = ftype['j']
        s1 = ftype['s1']
        for c, r in cr_list:
            k1 = (c, r, ai)
            k2 = (c, r + s1, aj)
            # print(k1)
            if k1 in index and k2 in index:
                data += [[i, index[k1], index[k2]]]
    data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj'])
    self_bonds = data.merge(bond_types, left_on='name', right_index=True)

    x1 = self_atoms.loc[self_bonds['aai']][['x', 'y', 'z']]
    x2 = self_atoms.loc[self_bonds['aaj']][['x', 'y', 'z']]
    self_bonds['r0'] = np.diag(sdist.cdist(x1, x2))/10

    # print(self_bonds)
    total_energy = 0

    for i, b in self_bonds.iterrows():
    # Units converted from
        parameters = [b['r0'],
                        b['Kb2'] / _df ** 2 * _ef,
                        b['Kb3'] / _df ** 3 * _ef,
                        b['Kb4'] / _df ** 4 * _ef]
        #print(parameters)
        # atom1 = traj_atoms.loc[traj_atoms['serial'] == int(b['aai'])+1].iloc[0]
        atom1 = traj_atoms.iloc[int(b['aai'])]
        atom2 = traj_atoms.iloc[int(b['aaj'])]
        # print(atom1)
        xyz_atom1 = [atom1['x'], atom1['y'], atom1['z']]        
        xyz_atom2 = [atom2['x'], atom2['y'], atom2['z']]
        # print(xyz_atom2)
        r = vabs(vector(xyz_atom1, xyz_atom2)) * _df
        distance = r - b['r0']
        # print(r - b['r0'] * _df)
        # energy = r - b['r0'] * _df
        energy = parameters[1]*(r-parameters[0])**2+parameters[3]*(r-parameters[0])**4
        print(int(b['aai']), int(b['aaj']), energy, xyz_atom2)
        total_energy += energy
    print(total_energy)


if __name__ == '__main__':
    main()