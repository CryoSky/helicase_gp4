#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Oct-14, latest modified on 2021-Oct-14

import MDAnalysis
import sys
import numpy as np


input_file = sys.argv[1]
u = MDAnalysis.Universe(input_file, format='PDB')

def distance(xyz_1, xyz_2):
    # xyz_1 = np.array((x1, y1, z1))
    # xyz_2 = np.array((x2, y2, z2))

    dist = np.linalg.norm(xyz_2-xyz_1)
    return dist

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def calc_dihedral(u1, u2, u3, u4):
    """ Calculate dihedral angle method. From bioPython.PDB
    (adapted to np.array)
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    [-pi, pi].
    """

    a1 = u2 - u1
    a2 = u3 - u2
    a3 = u4 - u3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm

    return rad


for ts in u.trajectory:
    chain_e_atoms = u.select_atoms('segid E and name CA')
    chain_f_atoms = u.select_atoms('segid F and name CA')
    dna_atoms = u.select_atoms('segid G and name S')

    com1 = chain_e_atoms.center_of_mass()
    com2 = chain_f_atoms.center_of_mass()

    tip1 = u.select_atoms('segid E and resid 202 and name CA')[0].position
    #print(tip1)
    tip2 = u.select_atoms('segid F and resid 202')[0].position

    ref_point_a = com1

    minim_dist = 99

    ref_point_b = np.array([0,0,0])
    ref_point_c = dna_atoms[40].position

    for each_nucleotide in dna_atoms:
        if minim_dist > distance(com1, each_nucleotide.position):
            minim_dist = distance(com1, each_nucleotide.position)
            ref_point_b = each_nucleotide.position

    print(calc_dihedral(tip1, com1, com2, tip2))
    #print(ref_point_b)

    # ab = (ref_point_a - ref_point_b)
    # ab = ab / np.linalg.norm(ab)

    # cb = (ref_point_c - ref_point_b)
    # cb = cb / np.linalg.norm(cb)
    # print(cb)

    #r1 = loopd1_atoms[-1].position - loopd1_atoms[0].position
    #r2 = loopd2_atoms[-1].position - loopd2_atoms[0].position

    # print(atp_atoms_com)

    #xyz_1 = np.array(r1)
    #xyz_2 = np.array(r2)


    #angle = angle_between(xyz_1, xyz_2)
    #ba = xyz_1 - xyz_2
    #bc = xyz_3 - xyz_2

    #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    #angle = np.arccos(cosine_angle)
    #print(angle)

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
