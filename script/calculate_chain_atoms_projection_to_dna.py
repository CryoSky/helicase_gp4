#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Jul-23, latest modified on 2021-Jul-23

import MDAnalysis
import sys
input_file = sys.argv[1]
chain_id = sys.argv[2]
u = MDAnalysis.Universe(input_file, format='PDB')
import numpy as np
dna_atoms = u.select_atoms('segid G and name S')
chain_f_atoms = u.select_atoms('segid %s' %chain_id)
#new_atoms = dna_atoms[0:-1:10]
new_atoms = dna_atoms[20:41:5]
#print(new_atoms)
reference_position = dna_atoms[20].position
#print(len(new_atoms))
for ts in u.trajectory:
    r = []
    for i in range(len(new_atoms)):
        try:
            distance = new_atoms[i].position - new_atoms[i+1].position
            r.append(distance)
        except:
            pass
    z_vector = np.average(r, axis=0) # Represents the direction from nucleotide 21 to 41 
    z_vector_norm = z_vector / np.linalg.norm(z_vector)
    #print(z_vector_norm) # Only z_vector_norm useful here
    
    chain_f_zaxis = []
    for each_f_atom in chain_f_atoms:
        # print(each_f_atom.position)
        relative_position = each_f_atom.position - reference_position
        dot_sum = np.dot(relative_position, z_vector_norm)
        # new_coord = dot_sum*z_vector_norm
        # print(np.linalg.norm(dot_sum)) # np.linalg.norm calculates Frobenius Norm
        # print(np.linalg.norm(new_coord)) That is nonsense
        chain_f_zaxis.append(np.linalg.norm(dot_sum))
    print(np.mean(chain_f_zaxis))
