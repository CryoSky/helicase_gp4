#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Jun-14, latest modified on 2021-Jun-14
# Pick out the index for each wham file that matches the basin region in 2D map
# Example in Linux: python this_code.py

import sys
import argparse
import itertools
import numpy as np
import math

__author__ = 'Shikai Jin'
__version__ = '0.0.1'
__date__ = "Jun-14-2021"

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
        for i,line in enumerate(text_file):
            if i >= line1-1 and i <= line2:
                stdout.append(line)
#        for line in itertools.islice(text_file, line1-1, line2):
#            stdout.append(line)
    return stdout

def get_pdbfile_i(pdbfile, i, n_lines):
    line_start = 1 + n_lines * i
    line_end = n_lines * (i + 1)
    pdb_part = getline_range(pdbfile, line_start, line_end)
    return pdb_part

wham_file = sys.argv[1]
wham_basin_center_first = sys.argv[2]
wham_basin_center_first = float(wham_basin_center_first)
wham_basin_range_first = sys.argv[3]
wham_basin_range_first = float(wham_basin_range_first)
wham_basin_center_second = sys.argv[4]
wham_basin_center_second = float(wham_basin_center_second)
wham_basin_range_second = sys.argv[5]
wham_basin_range_second = float(wham_basin_range_second)
simulation_file = sys.argv[6]
simulation_prefix = simulation_file[:-4]
simulation_prefix = simulation_prefix.split('/')[-1]
#print(simulation_prefix)

frames_index = []
with open(wham_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        line = line.split()
        if float(line[1]) < wham_basin_center_first + wham_basin_range_first and float(line[1]) > wham_basin_center_first - wham_basin_range_first:
            if float(line[2]) < wham_basin_center_second + wham_basin_range_second and float(line[2]) > wham_basin_center_second - wham_basin_range_second:
                print(line[0])
                frames_index.append(int(line[0]))

n_snapshot, n_lines, n_atoms = get_pdbfile_line_atoms(simulation_file)
for i in range(int(n_snapshot)):
    pdb_part = get_pdbfile_i(simulation_file, i, n_lines)
    if i+1 in frames_index:
        #print(pdb_part[0])
        with open(simulation_prefix+f'_{i+1}.pdb', 'w') as fwrite:
            fwrite.writelines(pdb_part)
