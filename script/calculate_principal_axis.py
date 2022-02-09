import numpy as np
import argparse
import MDAnalysis

__author__ = 'Shikai Jin'
__version__ = '0.0.1'
__date__ = "Aug-29-2021"


def get_axis(xyz):
    coord = np.array(xyz, float)
    center = np.mean(coord, axis=0)
    coord = coord - center
    # #This is the covariance matrix
    # inertia = np.dot(coord.transpose(), coord)
    # V = inertia / (len(coord) - 1)
    V = np.cov(coord.T)
    e_values, e_vectors = np.linalg.eig(V)
    # #order eigen values (and eigen vectors)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    return axis1


def main():
    parser = argparse.ArgumentParser(
        description="This script computes the crossing angle between the principal axes of chain F and chain E.")
    parser.add_argument("input_file", help="The file name of the input file.", type=str)
    parser.add_argument("output_file", help="The file name of output file.", type=str)
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    xyz_f = []
    xyz_e = []

    # with open(input_file, 'r') as pdb_file:
    #     for line in pdb_file:
    #         if line[12:16].strip() == "CA":
    #             if line[21] == 'E':
    #                 # extract x, y, z coordinates for carbon alpha atoms
    #                 x = float(line[30:38].strip())
    #                 y = float(line[38:46].strip())
    #                 z = float(line[46:54].strip())
    #                 xyz_e.append([x, y, z])
    #             elif line[21] == 'F':
    #                 # extract x, y, z coordinates for carbon alpha atoms
    #                 x = float(line[30:38].strip())
    #                 y = float(line[38:46].strip())
    #                 z = float(line[46:54].strip())
    #                 xyz_f.append([x, y, z])

    u = MDAnalysis.Universe(input_file, format='PDB')

    chain_e_atoms = u.select_atoms('segid E and name CA')
    chain_f_atoms = u.select_atoms('segid F and name CA')

    for ts in u.trajectory:
        xyz_e = chain_e_atoms.positions
        xyz_f = chain_f_atoms.positions
        axis1_e = get_axis(xyz_e)
        axis1_f = get_axis(xyz_f)

        xyz_1 = axis1_e
        if xyz_1[2] < 0:
            xyz_1 = -xyz_1
        xyz_2 = np.array((0,0,0))
        xyz_3 = axis1_f
        if xyz_3[2] < 0:
            xyz_3 = -xyz_3

        #dist = np.linalg.norm(xyz_2-xyz_1)
        ba = xyz_1 - xyz_2
        bc = xyz_3 - xyz_2

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        angle = angle * 180 / np.pi
        print(angle)

if __name__ == '__main__':
    main()