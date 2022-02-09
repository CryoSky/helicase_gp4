#!/usr/bin/env python
"""
This package implements a coarse-grained ATP force field that used in helicase system.
"""

__author__ = 'Shikai Jin'
__version__ = '1.1.0'
__date__ = "Mar-28-2021"

import simtk.openmm.app
import simtk.openmm
import simtk.unit as unit
import configparser
import numpy as np
import itertools
import scipy.spatial.distance as sdist
import os
import pdbfixer
import pandas
import subprocess
import nose

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
# print("Now loading the cgligand package from %s" %str(__location__))
_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor
_df = 1 * unit.angstrom / unit.nanometer  # distance scaling factor
_af = 1 * unit.degree / unit.radian  # angle scaling factor
xml = f'{__location__}/cgligand.xml'
_ligandResidues = ['ATP', 'ADP', 'GTP', 'TTP', 'MG', 'ATM', 'ADM']

def parseConfigTable(config_section):
    """Parses a section of the configuration file as a table.
       This function is used to parse the cgligand.conf file for the interaction."""

    def readData(config_section, a):
        """Filters comments and returns values as a list."""
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


def parsePDB(pdb_file):
    """Transforms the pdb file into a pandas table for easy access and data editing."""

    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),  # record name
                    serial=int(line[6:11]),          # atom serial number
                    name=str(line[12:16]).strip(),   # atom name
                    altLoc=str(line[16:17]),         # alternate location indicator
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[21:22]),
                    resSeq=int(line[22:26]),         # residue sequence number
                    iCode=str(line[26:27]),          # code for insertion of residues
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=1.0 if line[54:60].strip() == '' else float(line[54:60]), # set to 1.0 because Plumed RMSD need 1.0
                    tempFactor=1.0 if line[60:66].strip() == '' else float(line[60:66]),
                    element=str(line[76:78]),        # element symbol, right-justified
                    charge=str(line[78:80]))         # charge on the atom, right-justified

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                lines += [pdb_line(line)]
    pdb_atoms = pandas.DataFrame(lines)
    pdb_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                           'resname', 'chainID', 'resSeq', 'iCode',
                           'x', 'y', 'z', 'occupancy', 'tempFactor',
                           'element', 'charge']]
    return pdb_atoms


def fixPDB(pdb_file):
    """Uses the pdbfixer library to fix a pdb file, replacing non standard residues, removing
    hetero-atoms and adding missing hydrogens. The input is a pdb file location,
    the output is a fixer object, which is a pdb in the openawsem format.
    Manual on https://raw.githubusercontent.com/pandegroup/pdbfixer/master/Manual.html"""
    fixer = pdbfixer.PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in list(keys):
        chain_tmp = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
            del fixer.missingResidues[key]

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms() # Only run when the SEQ section in PDB file contains info of sequence
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def writeAfterFixer(fixer):
    """Write the PDB file after running of the pdbfixer."""
    simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open('fixer_output.pdb', 'w'))


def pdb2table(pdb):
    """Parses a pdb in the openmm format and outputs a table that contains all the information
    on a pdb file. Definition based on official PDB format."""
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


class Ligand(object):
    """ A ligand object that interacts with protein and DNA/RNA."""

    def __init__(self, periodic=True):
        """Initializes an ligand object"""
        self.periodic = periodic

    def __repr__(self):
        return f'CGLIGAND ligand object ({len(self.atoms)} atoms)'
        # print the sequence and the identity of the RNA object

    def parseConfigurationFile(self, configuration_file=f'{__location__}/cgligand.conf'):
        """Reads the configuration file for the force field. The default configuration file is openna.conf
        and it contains most of the parameters used in the simulation."""
        self.configuration_file = configuration_file
        config = configparser.ConfigParser()
        config.read(configuration_file)

        # Parse all sections of the configuration file
        self.config = {}
        for c in config.sections():
            self.config.update({c: parseConfigTable(config[c])})

        # Assign main sections to variables
        self.particle_definition = self.config['Particles']
        self.bond_definition = self.config['Bonds']

    def writePDB(self, pdb_file='clean.pdb'):
        """ Writes a minimal version of the pdb file needed for openmm """
        # Compute chain field
        if type(self.atoms['chainID'].iloc[0]) is not str:
            chain_ix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
            self.atoms['chainID'] = [chain_ix[i - 1] for i in self.atoms['chainID']]

        # Compute element fields
        element_ix = {'PA': 'P', 'PB': 'P', 'PG': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'C': 'O', 'G': 'C', 'U': 'S'}  # Elements choosen to keep VMD colors
        self.atoms.loc[:, 'element'] = [element_ix[atomType] for atomType in self.atoms['name']]

        # Write pdb file
        with open(pdb_file, 'w+') as pdb:
            for i, atom in self.atoms.iterrows():
                pdb_line = f'ATOM  {i + 1:>5} {atom["name"]:^4} {atom.resname:<3} {atom.chainID}{atom.resSeq:>4}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.tempFactor:>6.2f}' + ' ' * 10 + f'{atom.element:>2}' + ' ' * 2 # modified to meet RNA requirement
                assert len(pdb_line) == 80, 'An item in the atom table is longer than expected'
                pdb.write(pdb_line + '\n')
        self.pdb_file = pdb_file
        return pdb_file
    
    def computeTopology(self, temp_name='temp'):
        """ Creates tables of bonds, angles and dihedrals with their respective parameters (bonded interactions).
        3SPN2.C requires a template structure to calculate the equilibrium bonds, angles and dihedrals.
        If template_from_structure is True, it will try to compute the equilibrium geometry using X3DNA.
        If template_from_structure is False, then the initial structure is expected to be the equilibrium geometry"""
        # Parse configuration file if not already done
        try:
            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()

        na_type = self.na_type
        #if na_type not in self.angle_definition['DNA'].unique():
        #    raise DNATypeError(self)

        # Rewrite index in case it is not ordered, important!
        self.atoms.index = range(len(self.atoms))


        self.template_atoms = self.atoms
        # Make an index to build the topology
        index = {}
        cr_list = set()  # Chain residue list
        for i, atom in self.atoms.iterrows():
            index.update({(atom['chainID'], atom['resSeq'], atom['name']): i})
            cr_list.update([(atom['chainID'], atom['resSeq'])])
        cr_list = list(cr_list)
        cr_list.sort()
        # max_chain = self.atoms['chain'].max()
        # max_residue = self.atoms['resSeq'].max()
        assert len(index) == len(self.atoms), 'Atom index was repeated'


        bond_types = self.bond_definition[self.bond_definition['ligand'] == na_type]
                
        # print(bond_types)
        # print(index)

        # Make a table with bonds
        data = []
        for i, ftype in bond_types.iterrows():
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            s1 = ftype['s1']
            for c, r in cr_list:
                k1 = (c, r, ai)
                k2 = (c, r + s1, aj)
                if k1 in index and k2 in index:
                    data += [[i, index[k1], index[k2]]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj'])
        self.bonds = data.merge(bond_types, left_on='name', right_index=True)

        if na_type == 'Ligand':
            # Make default distances the same as the initial distance
            x1 = self.template_atoms.loc[self.bonds['aai']][['x', 'y', 'z']]
            x2 = self.template_atoms.loc[self.bonds['aaj']][['x', 'y', 'z']]
            self.bonds['r0'] = np.diag(sdist.cdist(x1, x2))


# @ symbol in Python https://docs.python.org/3/reference/compound_stmts.html#index-30
# https://docs.python.org/3/library/functions.html#staticmethod
# A class method is a method that is bound to a class rather than its object. It doesn't require creation of a class instance, much like staticmethod.
#The difference between a static method and a class method is:
#Static method knows nothing about the class and just deals with the parameters
#Class method works with the class since its parameter is always the class itself.

    @classmethod
    def fromCoarsePDB(cls, pdb_file, na_type='Ligand', temp_name='temp'):
        """Initializes a RNA object from a pdb file containing the Coarse Grained atoms"""
        self = cls()

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
                        occupancy=1.0 if line[54:60].strip() == '' else float(line[54:60]),
                        tempFactor=1.0 if line[60:66].strip() == '' else float(line[60:66]),
                        element=str(line[76:78]).strip(),
                        charge=str(line[78:80]).strip())

        with open(pdb_file, 'r') as pdb:
            lines = []
            for line in pdb:
                if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                    lines += [pdb_line(line)]
        pdb_atoms = pandas.DataFrame(lines)
        self.atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                                'resname', 'chainID', 'resSeq', 'iCode',
                                'x', 'y', 'z', 'occupancy', 'tempFactor',
                                'element', 'charge']]
        # print(self.atoms.columns)
        # self.atoms.loc[:, 'chain'] = self.atoms['chainID']
        # self.atoms.loc[:, 'residue'] = self.atoms['resSeq']
        self.atoms.loc[:, 'type'] = self.atoms['name']
        # print(self.atoms.columns)
        # Initialize the system from the pdb
        self.na_type = na_type
        self.parseConfigurationFile()
        self.computeTopology(temp_name=temp_name)
        self.pdb_file = pdb_file
        return self
    
    @classmethod
    def fromPDB(cls, pdb_file, na_type='Ligand', output_pdb='clean.pdb', temp_name='temp'):
        """Creates a DNA object from a complete(atomistic) pdb file"""
        self = cls()
        pdb = fixPDB(pdb_file)
        pdb_table = pdb2table(pdb)

        self.atoms = self.CoarseGrain(pdb_table)
        self.na_type = na_type
        self.parseConfigurationFile()
        self.computeTopology(temp_name=temp_name)
        self.writePDB(output_pdb)
        #self.atomistic_model=temp
        return self
    
    # https://stackabuse.com/pythons-classmethod-and-staticmethod-explained/
    @staticmethod
    def CoarseGrain(pdb_table, ligand_name='ATP'):
        """ Selects RNA atoms from a pdb table and returns a table containing only the coarse-grained atoms for 3SPN2"""
        masses = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "P": 30.973762, "Mg": 24.305}
        CG = {"C5\'": 'S', "C5*": 'S', "C4\'": 'S', "C4*": 'S', "C3\'": 'S', "C3*": 'S', 
              "C2\'": 'S', "C2*": 'S', "C1\'": 'S', "C1*": 'S',
              "O4\'": 'S', "O4*": 'S', "O3\'": 'S', "O2\'": 'S',
              "H5'": 'S', "H5''": 'S', "H4'": 'S', "H3'": 'S', "H2'": 'S', "H2''": 'S', "H1'": 'S',
              "C6": 'B', "C5M": 'B', "C5": 'B', "C4": 'B', "C2": 'B',
              "O4": 'B', "O2": 'B', "N3": 'B', "N1": 'B',
              "PA": 'PA', "O1A": 'PA', "O2A": 'PA', "O3A": 'PA',
              "PB": 'PB', "O1B": 'PB', "O2B": 'PB', "O3B": 'PB',
              "PG": 'PG', "O1G": 'PG', "O2G": 'PG', "O3G": 'PG',
              "MG": 'MG'
              }
        cols = ['recname', 'serial', 'name', 'altLoc',
                'resname', 'chainID', 'resSeq', 'iCode',
                'x', 'y', 'z', 'occupancy', 'tempFactor',
                'element', 'charge', 'type']
        temp = pdb_table.copy()
        temp = temp[temp['resname'].isin(_ligandResidues)]
        base_name = ligand_name[0]
        

        # Group the atoms by sugar, phosphate or base
        temp['group'] = temp['name'].replace(CG) # Replace the true atom to group name and add in a new column group
        
        temp = temp[temp['group'].isin(['PA', 'PB', 'PG', 'S', 'B', 'MG'])]
        #print(temp)
        # Calculate center of mass
        temp['element'] = temp['element'].str.strip()
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
        Coarse['type'] = Coarse['name']
        # Set element (depends on base)
        Coarse['name'] = Coarse['name'].replace({'B': base_name})
        Coarse['element'] = Coarse['name'].replace({'PA': 'P', 'PB': 'P', 'PG': 'P', 'S': 'H', 'A': 'N', 'U': 'S', 'G': 'C', 'C': 'O', 'MG': 'MG'})
        Coarse['resname'] = Coarse['resname'].replace({'TTP': 'ATP', 'MG': 'MG'})
        # Remove P from the beggining
        drop_list = []
        for chain in Coarse.chainID.unique():
            sel = Coarse[Coarse.chainID == chain]
        # Renumber
        Coarse.index = range(len(Coarse))
        Coarse['serial'] = Coarse.index
        return Coarse[cols]
