import itertools
import argparse
import math
import numpy as np

def vector(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]


def vabs(a):
    return math.sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))

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


def get_ca_s_atoms(pdb_part, chain_id):
    index = 1
    chain_data = []  # Initialize the chain
    for line in pdb_part:
        try:
            if line.split()[0] != 'ENDMDL':
                if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
                    if (line[12:16].strip() == 'CA') and (line[21:22].strip() == chain_id):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        atom = [x, y, z]
                        chain_data.append(atom)
                        index += 1
        except:
            pass

    return chain_data

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

def get_rg(pdb_part):
    atoms_data = []


    for line in pdb_part:
        if (line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atom = [x, y, z]
            atoms_data.append(atom)

    num_atoms = len(atoms_data)
    rg_sum = 0
    for ia in range(0, num_atoms):
        for ja in range(ia+1, num_atoms):
            rv = vector(atoms_data[ia], atoms_data[ja])
            rg_sum += pow(rv[0], 2) + pow(rv[1], 2) + pow(rv[2], 2)
    rg_sum = math.sqrt(rg_sum/(num_atoms**2))
    return rg_sum


def compute_traj(trajectory_file, output):
    n_snapshot, n_lines, n_atoms = get_pdbfile_line_atoms(trajectory_file)
    print(n_snapshot)
    print(n_lines)
    print(n_atoms)
    out = open(output, "w")
    for i in range(int(n_snapshot)):
        pdb_part = get_pdbfile_i(trajectory_file, i, n_lines)
        com = get_rg(pdb_part)
        out.write(str(round(com, 3)) + "\n")
    out.close()

def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the Q interface value between a topology and trajectory pdb file. The trajectory must be clean. Cannot change residue index, must be same")
    parser.add_argument("trajectory", help="The file name of trajectory in dcd format, usually converted from dump.lammpstrj in AWSEM", type=str)
    #parser.add_argument("topology", help="The file name of topology for align, any frame of trajectory is fine", type=str)
    parser.add_argument("output", help="The file name of output",
                        type=str)
    args = parser.parse_args()
    trajectory = args.trajectory
    #topology = args.topology
    output = args.output

    compute_traj(trajectory, output)


if __name__ == '__main__':
    main()
