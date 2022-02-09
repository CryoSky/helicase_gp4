#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2020-Nov-8, latest modified on 2021-Mar-29
# A final loading of whole helicase system: protein+DNA+ATP+Mg
# Example in Linux: python this_code.py 0.02

__author__ = 'Shikai Jin'
__version__ = '1.1.0'
__date__ = "Mar-29-2021"

import ffcgligand
import pdbfixer
import ffAWSEM
import open3SPN2
import sys

input_file = sys.argv[1]

fix=pdbfixer.PDBFixer(filename=input_file)
#print(fix.__dict__)
complex_table=ffcgligand.pdb2table(fix)
#print(complex_table)
ligand_atoms=ffcgligand.Ligand.CoarseGrain(complex_table)
#print(ligand_atoms)
dna_atoms=open3SPN2.DNA.CoarseGrain(complex_table)
protein_atoms=ffAWSEM.Protein.CoarseGrain(complex_table)
#print(ligand_atoms)
#print(dna_atoms)

#Merge the models
import pandas
Coarse=pandas.concat([protein_atoms, dna_atoms, ligand_atoms],sort=False)
Coarse.index=range(len(Coarse))
serial_new = [i + 1 for i in list(Coarse.index)]
Coarse.serial=serial_new

#Save the protein_sequence
from Bio.PDB.Polypeptide import three_to_one
_AWSEMresidues=['IPR','IGL','NGP']
protein_data=Coarse[Coarse.resname.isin(_AWSEMresidues)].copy()
resix = (protein_data.chainID + '_' + protein_data.resSeq.astype(str))
res_unique = resix.unique()
protein_data['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
protein_sequence=[r.iloc[0]['real_resname'] for i, r in protein_data.groupby('resID')]
protein_sequence_one = [three_to_one(a) for a in protein_sequence]

# Deal with protein issue
with open('protein.seq','w+') as ps:
    ps.write(''.join(protein_sequence_one))

# Create a merged PDB
def writePDB(atoms,pdb_file):
    with open(pdb_file, 'w+') as pdb:
        for i, atom in atoms.iterrows():
            pdb_line = f'{atom.recname:<6}{atom.serial:>5} {atom["name"]:^4}{atom.altLoc:1}'+\
                       f'{atom.resname:<3} {atom.chainID:1}{atom.resSeq:>4}{atom.iCode:1}   '+\
                       f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' +\
                       f'{1:>6.2f}{1:>6.2f}'+' ' * 10 +\
                       f'{atom.element:>2}{atom.charge:>2}'
            assert len(pdb_line) == 80, f'An item in the atom table is longer than expected ({len(pdb_line)})\n{pdb_line}'
            pdb.write(pdb_line + '\n')
writePDB(Coarse,'clean.pdb')
