#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Dec-27, latest modified on 2021-Dec-27

__author__ = 'Shikai Jin'
__version__ = '0.0.1'
__date__ = "Dec-27-2021"

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


def draw_contactmap(matrix, maximum_num, xname, yname, title):
    n_atoms = 0
    (n_atomsx, n_atomsy) = np.shape(matrix)
    # np.savetxt('test.csv', matrix)
    # print(n_atomsx)
    # print(n_atomsy)
    # if n_atomsx == n_atomsy:
    #     n_atoms = n_atomsx
    # else:
    #     sys.exit("matrix is not square")
        
    #plt.figure(figsize=(30,30))
    plt.figure()
    current_axis = plt.gca()
    for i in range(n_atomsx):
        for j in range(n_atomsy):
            if matrix[i][j] > maximum_num*0.4:
                print("{:d} {:d} {:.3f}".format(i+1, j+1, matrix[i][j] / maximum_num))
                rect = patches.Rectangle((j+1,i+1), 1, 1, edgecolor='red', facecolor='red', linewidth=0.5) # color can be reds(matrix[i][j]/50)
                # if i >= j:
                #     rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='red', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='darkorange', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='none', facecolor='r', linewidth=0.5)
                # else:
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='mediumspringgreen', linewidth=0.5)
                #     rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='k', facecolor='royalblue', linewidth=0.5)
                #     #rect = patches.Rectangle((i+1,j+1), 1, 1, edgecolor='none', facecolor='b', linewidth=0.5)
                current_axis.add_patch(rect)

    #plt.xlim(1, 483)
    #plt.ylim(1, 51)
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
    plt.savefig('%s' %title, bbox_inches='tight', dpi=300)
    # plt.show()


data = np.zeros((197, 483))

input_file = sys.argv[1]
output_file = sys.argv[2]
maximum_num = int(sys.argv[3])

with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        line = line.split()
        # print(line)
        try:
            data[int(line[0])][int(line[1])] += 1
        except:
            pass

xname = 'Chain F'
yname = 'DNA'
title = 'None'

draw_contactmap(data, maximum_num, xname, yname, output_file)
