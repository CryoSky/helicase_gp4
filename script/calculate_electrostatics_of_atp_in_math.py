#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-May-3, latest modified on 2021-May-3
# Output the final electrostatics energy of a given ATP in protein+DNA+ligand system
# Example in Linux: python this_code.py

import sys
import argparse
import itertools
import numpy as np
import math

__author__ = 'Shikai Jin'
__version__ = '0.0.2'
__date__ = "May-11-2021"

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

def load_protein_sequences(fasta="../crystal_structure.fasta"): # Need to modify the comment in fasta file
    seq_dic = {}
    chain = None
    with open(fasta) as f:
        for line in f:
            if line[0] == ">":
                # assert line[:19] == ">CRYSTAL_STRUCTURE:" # No need for assertion anymore
                if chain is not None:
                    seq_dic[chain] = seq
                chain = line[19]
                seq = ""
            else:
                seq += line.replace("\n", "")
        seq_dic[chain] = seq
    return seq_dic

def compute_electrostatics(pdb_part, protein_sequence, atp_atom_index):
    # First screen all CB atoms
    # Dict, key is chain ID, value is a small list with residue ID and CB atoms coord
    total_energy = 0
    total_cb_atoms = {}
    for each_ligand in atp_atom_index:
        for line in pdb_part:
            if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
                if (line[12:16].strip() == 'CB'):
                    residue_index = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom = [x, y, z]
                    try:
                        total_cb_atoms[line[21:22].strip()][residue_index]=atom
                        #print(str({residue_index: atom}))
                    except:
                        total_cb_atoms[line[21:22].strip()]={residue_index: atom}
                if int(line[6:11]) == each_ligand:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ligand_atom_position = [x, y, z]
                    print("The ligand position is:" + str(ligand_atom_position))

        for key in protein_sequence.keys():
            # try:
            sequence = protein_sequence[key]
            index = 1
            
            for _, residue in enumerate(sequence):
                if residue == 'G':
                    index += 1
                    continue
                protein_atom_position = total_cb_atoms[key][index]
                if residue == 'K' or residue == 'R':
                    r = vabs(vector(ligand_atom_position, protein_atom_position))
                    this_energy = 1 * -2 * math.exp(-r/12) / 0.2348 / r
                elif residue == 'D' or residue == 'E':
                    r = vabs(vector(ligand_atom_position, protein_atom_position))
                    this_energy = -1 * -2 * math.exp(-r/12) / 0.2348 / r
                else:
                    this_energy = 0
                total_energy += this_energy
                index += 1
        # except:
        #     sys.exit("The chain ID %s in PDB is not exist in protein sequence. Change the 20th char of fasta header line to chain ID." %key)
    return total_energy


def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the electrostatics value of all protein atoms and a ligand atom in a trajectory file.")
    parser.add_argument("trajectory", help="The file name of trajectory in pdb format.", type=str)
    parser.add_argument("protein_sequence", help="The file name of protein sequence", type=str)
    parser.add_argument("atp_atom_index", help="The index of ligand atom", type=str)
    parser.add_argument("output", help="The file name of output",
                        type=str)
    args = parser.parse_args()
    trajectory_file = args.trajectory
    protein_sequence = args.protein_sequence
    atp_atom_index = args.atp_atom_index
    #topology = args.topology
    output = args.output

    atp_atom_index = atp_atom_index.split(',')
    
    n_snapshot, n_lines, n_atoms = get_pdbfile_line_atoms(trajectory_file)
    sequence_data = load_protein_sequences(protein_sequence) # load sequences
    out = open(output, "w")
    for i in range(int(n_snapshot)):
        pdb_part = get_pdbfile_i(trajectory_file, i, n_lines)
        energy = compute_electrostatics(pdb_part, sequence_data, atp_atom_index)
        out.write(str(round(energy, 3)) + "\n")
    out.close()


if __name__ == '__main__':
    main()
