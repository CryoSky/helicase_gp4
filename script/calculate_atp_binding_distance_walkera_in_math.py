#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Oct-14, latest modified on 2021-Oct-14

import MDAnalysis
import sys
input = sys.argv[1]
u = MDAnalysis.Universe(input, format='PDB')
import numpy as np

# atp_atoms = u.select_atoms('segid M')


# atp_atoms_com = atp_atoms.center_of_mass()
# print(atp_atoms_com)

for ts in u.trajectory:
    walkera_atoms = u.select_atoms('segid A and resid 56-57')
    arginine_finger_atoms = u.select_atoms('segid F and resid 259')
    atp_atoms = u.select_atoms('segid M')

    walkera_atoms_com = walkera_atoms.center_of_mass()
    arginine_finger_atoms_com = arginine_finger_atoms.center_of_mass()
    atp_atoms_com = atp_atoms.center_of_mass()
    # print(atp_atoms_com)

    xyz_1 = np.array(walkera_atoms_com)
    xyz_2 = np.array(arginine_finger_atoms_com)
    xyz_3 = np.array(atp_atoms_com)

    dist = np.linalg.norm(xyz_3-xyz_1)
    #ba = xyz_1 - xyz_2
    #bc = xyz_3 - xyz_2

    #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    #angle = np.arccos(cosine_angle)
    print(dist)

#new_atoms = dna_atoms[0:-1:10]
#new_atoms = dna_atoms[20:41:5]
#print(new_atoms)
# reference_position = dna_atoms[20].position
# #print(len(new_atoms))
# for ts in u.trajectory:
#     r = []
#     for i in range(len(new_atoms)):
#         try:
#             distance = new_atoms[i].position - new_atoms[i+1].position
#             r.append(distance)
#         except:
#             pass
#     z_vector = np.average(r, axis=0) # Represents the direction from nucleotide 21 to 41 
#     z_vector_norm = z_vector / np.linalg.norm(z_vector)
#     #print(z_vector_norm) # Only z_vector_norm useful here
    
#     chain_f_zaxis = []
#     for each_f_atom in chain_f_atoms:
#         # print(each_f_atom.position)
#         relative_position = each_f_atom.position - reference_position
#         dot_sum = np.dot(relative_position, z_vector_norm)
#         # new_coord = dot_sum*z_vector_norm
#         # print(np.linalg.norm(dot_sum)) # np.linalg.norm calculates Frobenius Norm
#         # print(np.linalg.norm(new_coord)) That is nonsense
#         chain_f_zaxis.append(np.linalg.norm(dot_sum))
#     print(np.mean(chain_f_zaxis))
