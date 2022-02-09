#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-May-24, latest modified on 2021-May-24

__author__ = 'Shikai Jin'
__version__ = '0.0.4'
__date__ = "Jul-5-2021"

import argparse
import sys
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import math
import matplotlib 
import itertools
matplotlib.use('Agg') 
# Change global font 
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches

reds = cm.get_cmap('Reds', 12)


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

def check_if_native(xyz_CAi, xyz_CAj):
    if vabs(vector(xyz_CAi, xyz_CAj)) < 12.0:
        return True
    else:
        return False

    
def compute_contactmap(ca_atoms_pdb, dna_atoms_pdb):
    #print(sorted(ca_atoms_pdb))
    width1 = sorted(ca_atoms_pdb)[-1] # Don't use length because we may segment domain
    begin1 = sorted(ca_atoms_pdb)[0]
    width2 = sorted(dna_atoms_pdb)[-1] # Don't use length because we may segment domain
    begin2 = sorted(dna_atoms_pdb)[0]    
    contactmap = np.zeros((width2, width1))
    for i in range(width2):
        for j in range(width1):
            if i+1 in dna_atoms_pdb.keys() and j+1 in ca_atoms_pdb.keys():
                if check_if_native(dna_atoms_pdb[i+1], ca_atoms_pdb[j+1]): # Remember here the index start from 0
                    contactmap[i][j] = 1.0
                    # contactmap[j][i] = contactmap[i][j]
    # print  np.max(contactmap)
    return begin1, contactmap

# def merge_matrix(contactmap1):
#     # if len(contactmap1) >= len(contactmap2):
#     #     n_atoms = len(contactmap1)
#     # else:
#     #     n_atoms = len(contactmap2)
#     (n_atomsx, n_atomsy) = np.shape(contactmap1)
#     matrix = np.zeros((n_atomsx, n_atomsy))
#     #print(n_atoms)
#     for i in range(n_atomsx):
#         for j in range(n_atomsy):
#             try:
#                 matrix[i][j] = contactmap1[i][j]
#             except: #Avoid error, when two matrices have different size, error occur when loop out of smaller one
#                 continue

#     for i in range(n_atomsx, 0, -1):
#         for j in range(n_atomsy):
#             try:
#                 matrix[i][j] = contactmap1[i][j]
#             except: #Avoid error, when two matrices have different size, error occur when loop out of smaller one
#                 continue
#     return matrix

def draw_contactmap(matrix, begin, xname, yname, title):
    n_atoms = 0
    (n_atomsx, n_atomsy) = np.shape(matrix)
    #np.savetxt('test.csv', matrix)
    print(n_atomsx)
    print(n_atomsy)
    # if n_atomsx == n_atomsy:
    #     n_atoms = n_atomsx
    # else:
    #     sys.exit("matrix is not square")
        
    #plt.figure(figsize=(30,30))
    plt.figure()
    current_axis = plt.gca()
    for i in range(n_atomsx):
        for j in range(n_atomsy):
            if matrix[i][j] > 0:
                rect = patches.Rectangle((j+1,i+1), 1, 1, edgecolor='k', facecolor=reds(matrix[i][j]/50), linewidth=0.5)
                # if i >= j:
                #     rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='red', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='darkorange', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='none', facecolor='r', linewidth=0.5)
                # else:
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='mediumspringgreen', linewidth=0.5)
                #     rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='royalblue', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='none', facecolor='b', linewidth=0.5)
                current_axis.add_patch(rect)

    plt.xlim(1, 286)
    plt.ylim(1, 100)
    #plt.xlabel(xname, fontsize=2.5*(n_atoms/10))
    plt.xlabel(xname, fontsize=10)
    plt.ylabel(yname, fontsize=10)
    if title != 'None':
        plt.title(title, fontsize=15)
    #plt.xlabel(xname, fontsize=40)
    #plt.ylabel(yname, fontsize=40)
    #plt.title(title, fontsize=40) 

    #plt.xticks(fontsize=25)
    #plt.xticks(fontsize=1.8*(n_atoms/10))
    #plt.yticks(fontsize=25)
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('%s.png' %title, bbox_inches='tight', dpi=300)
    # plt.show()

def pdb_load_protein(pdb_file, chain_label, start_text, end_text):
    verbose = True
    if start_text == None:
        start_pair = [None]
    else:
        start_pair = start_text.split(',')
    if end_text == None:
        end_pair = [None]
    else:
        end_pair = end_text.split(',')
        
    if len(start_pair) != len(end_pair):
        sys.exit("The pair number of start point and end point are different")
    
    
    p = PDBParser(PERMISSIVE=1)  # Useful when find errors in PDB
    s = p.get_structure('pdb', pdb_file)
    chains = s[0].get_list()  # Compatible for multiple chains, but only for first model
    chain_id = 0
    ca_atoms_pdb = {}  # A new dictionary to record the C_alpha atom coordinate in pdb file
    pdb_chain_id = []

    for chain in chains:
        #print(dir(chain))
        if chain.id == chain_label:
            for pair_index in range(len(start_pair)):
                chain_id = chain_id + 1
                first_residue = chain.get_unpacked_list()[0]
                last_residue = chain.get_unpacked_list()[-1]
                sequence_id_flag = 0
                if end_pair[pair_index] is None:  # Default end point is last residue
                    end = int(last_residue.get_id()[1])
                    # print("End value is %s" %end)
                    # print(chain[353]['CA'].get_coord())
                else:
                    end = int(end_pair[pair_index])
                if start_pair[pair_index] is None:  # Some fxxking pdbs start from -7
                    start = int(first_residue.get_id()[1])
                else:
                    start = int(start_pair[pair_index])
                for res in chain:
                    is_regular_res = (res.has_id('CA') and res.has_id('O')) or (
                        res.get_id()[1] == last_residue and res.has_id(
                            'CA'))  # Some pdbs lack the O atom of the last one residue...interesting
                    hetero_flag = res.get_id()[
                        0]  # Get a list for each residue, include hetero flag, sequence identifier and insertion code
                    # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
                    if (hetero_flag == ' ' or hetero_flag == 'H_MSE' or hetero_flag == 'H_M3L' or hetero_flag == 'H_CAS') \
                            and is_regular_res:
                        # The residue_id indicates that there is not a hetero_residue or water ('W')
                        sequence_id = res.get_id()[1]
                        # print(sequence_id)
                        if sequence_id_flag == 0:
                            last_sequence_id = sequence_id
                        else:
                            if last_sequence_id != int(sequence_id) - 1:
                                if verbose:
                                    print ("WARNING: the residues between %s and %s are lost" % (last_sequence_id, sequence_id))
                                else:
                                    pass
                                # Some fxxking pdbs lost residues halfway
                            last_sequence_id = sequence_id
                            # print(res.get_id())
                        if sequence_id >= start and sequence_id <= end:
                            ca_atoms_pdb[sequence_id] = res['CA'].get_coord() # Biopython only gets first CA if occupancy is not 1
                            # print(ca_atoms_pdb.keys())
                            # ca_atoms_pdb.append(res['CA'].get_coord())
                            pdb_chain_id.append(chain_id)
                        sequence_id_flag = sequence_id_flag + 1
        else:
            chain_id = chain_id + 1
    #print(ca_atoms_pdb.keys())
    return ca_atoms_pdb

def pdb_load_dna(pdb_file, chain_label, start_text, end_text):
    verbose = True
    if start_text == None:
        start_pair = [None]
    else:
        start_pair = start_text.split(',')
    if end_text == None:
        end_pair = [None]
    else:
        end_pair = end_text.split(',')
        
    if len(start_pair) != len(end_pair):
        sys.exit("The pair number of start point and end point are different")
    
    
    p = PDBParser(PERMISSIVE=1)  # Useful when find errors in PDB
    s = p.get_structure('pdb', pdb_file)
    chains = s[0].get_list()  # Compatible for multiple chains, but only for first model
    chain_id = 0
    ca_atoms_pdb = {}  # A new dictionary to record the C_alpha atom coordinate in pdb file
    pdb_chain_id = []

    for chain in chains:
        #print(dir(chain))
        if chain.id == chain_label:
            for pair_index in range(len(start_pair)):
                chain_id = chain_id + 1
                first_residue = chain.get_unpacked_list()[0]
                last_residue = chain.get_unpacked_list()[-1]
                sequence_id_flag = 0
                if end_pair[pair_index] is None:  # Default end point is last residue
                    end = int(last_residue.get_id()[1])
                    # print("End value is %s" %end)
                    # print(chain[353]['CA'].get_coord())
                else:
                    end = int(end_pair[pair_index])
                if start_pair[pair_index] is None:  # Some fxxking pdbs start from -7
                    start = int(first_residue.get_id()[1])
                else:
                    start = int(start_pair[pair_index])
                for res in chain:
                    is_regular_res = True
                    if is_regular_res:
                        # The residue_id indicates that there is not a hetero_residue or water ('W')
                        sequence_id = res.get_id()[1]
                        # print(sequence_id)
                        if sequence_id_flag == 0:
                            last_sequence_id = sequence_id
                        else:
                            if last_sequence_id != int(sequence_id) - 1:
                                if verbose:
                                    print ("WARNING: the residues between %s and %s are lost" % (last_sequence_id, sequence_id))
                                else:
                                    pass
                                # Some fxxking pdbs lost residues halfway
                            last_sequence_id = sequence_id
                            # print(res.get_id())
                        if sequence_id >= start and sequence_id <= end:
                            ca_atoms_pdb[sequence_id] = res['S'].get_coord() # Biopython only gets first CA if occupancy is not 1
                            # print(ca_atoms_pdb.keys())
                            # ca_atoms_pdb.append(res['CA'].get_coord())
                            pdb_chain_id.append(chain_id)
                        sequence_id_flag = sequence_id_flag + 1
        else:
            chain_id = chain_id + 1
    #print(ca_atoms_pdb.keys())
    return ca_atoms_pdb

def main():
    parser = argparse.ArgumentParser(
        description="This script draws contactmap for protein and DNA region.")
    parser.add_argument("PDB_filename1", help="The file name of input pdb1", type=str)
    parser.add_argument("csv_filename", help="The file name of csv", type=str)
    parser.add_argument("--start1", help="The start residue of protein1", type=str)
    parser.add_argument("--end1", help="The end residue of protein1", type=str)
    parser.add_argument("--x_label", help="The label name of x axis, corresponding to PDB2", default = 'X axis', type=str)
    parser.add_argument("--y_label", help="The label name of y axis, corresponding to PDB1", default = 'Y axis', type=str)
    parser.add_argument("--title", help="The title of whole image", default = 'None', type=str)
    
    args = parser.parse_args()
    pdb1_file = args.PDB_filename1
    start1 = args.start1
    end1 = args.end1
    x_label = args.x_label
    y_label = args.y_label
    title = args.title
    csv_file = args.csv_filename

    if pdb1_file[-4:].lower() != ".pdb":
        sys.exit("It must be a pdb file.")

    full_contactmap = np.zeros((100, 286))

    n_snapshot, n_lines, n_atoms = get_pdbfile_line_atoms(pdb1_file)
    for i in range(int(n_snapshot)):
        pdb_part = get_pdbfile_i(pdb1_file, i, n_lines)
        with open('temp.pdb', 'w') as fwrite:
            fwrite.writelines(pdb_part)
        ca_atoms_pdb1 = pdb_load_protein('temp.pdb', 'F', start1, end1)
        dna_atoms_pdb1 = pdb_load_dna(pdb1_file, 'G', '1', '100')
        begin1, contactmap1 = compute_contactmap(ca_atoms_pdb1, dna_atoms_pdb1)
        full_contactmap = full_contactmap + contactmap1
        #print(full_contactmap)
    np.savetxt(csv_filename, full_contactmap)
    draw_contactmap(full_contactmap, begin1, x_label, y_label, title)
    

if __name__ == '__main__':
    main()
