import itertools
import argparse
from Bio.PDB.PDBParser import PDBParser
import math

def vector(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]


def vabs(a):
    return math.sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))

def compute_qinterface_for_reference(reference_file):
    structure_interactions = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', reference_file)
    model = structure[0]
    chain_start_i = 0  # record the start residue of each chain
    count_i = 0  # count on the current index of the residue inside each chain
    ic = 0  # ic indicates the index of each chain
    chain_list = list(model.get_chains())

    chain_a_id = 'A'
    chain_b_id = 'F'

    for chain_index1 in range(0, len(chain_list)):
        chain_i = chain_list[chain_index1]
        if chain_i.full_id[-1] == chain_a_id:
            for temp1 in range(0, chain_index1):
                chain_start_i += int(math.fsum(1 for _ in chain_list[temp1].get_residues())) # Calculates the sum index from beginning to now
            count_i = 0 # Initialize after each cycle that read chain i
            ic += 1

            chain_start_j = 0
            jc = 0
            count_j = 0
            chain_i_count_flag = 0 # prevent for multiple counting

            for chain_index2 in range(0, len(chain_list)): # Don't scan repetitive chains interactions
                chain_j_count_flag = 0
                chain_j = chain_list[chain_index2]
                if chain_j.full_id[-1] == chain_b_id:
                    for temp2 in range(0, chain_index2):
                        chain_start_j += int(math.fsum(1 for _ in chain_list[temp2].get_residues()))
                    #print(chain_start_j)
                    jc += 1

                    for i, residue_i in enumerate(chain_i.get_residues()):
                        if chain_i_count_flag == 0:
                            count_i += 1
                            count_j = 0
                        for j, residue_j in enumerate(chain_j.get_residues()):
                            if chain_j_count_flag == 0:
                                count_j += 1

                            ca_i = residue_i['CA']
                            ca_j = residue_j['CA']

                            r_ijN = abs(ca_i - ca_j)
                            #print(r_ijN)
                            if r_ijN >= 12: continue
                            sigma_ij = ((abs(i - j) + 1) ** 0.15)  # 0.1 nm = 1 A
                            gamma_ij = 1.0

                            # if reference_chain_name != "ALL" and (chain.id not in reference_chain_name):
                            #	continue
                            i_index = count_i
                            j_index = count_j
                            structure_interaction = [i_index, j_index, gamma_ij, r_ijN, sigma_ij]
                            structure_interactions.append(structure_interaction)

    return structure_interactions

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

def compute_qinterface_for_traj(traj_ca_atoms_1, traj_ca_atoms_2, reference_qinterface_data):
    counts = 0
    q_sum = 0

    for each_line in reference_qinterface_data:
        try:
            #print("first is:" + str(traj_ca_atoms_1[each_line[0]-1]))
            #print("second is:" + str(traj_ca_atoms_2[each_line[1]-1]))
            rij_traj = vabs(vector(traj_ca_atoms_1[each_line[0]-1], traj_ca_atoms_2[each_line[1]-1]))
            q_sum += math.exp(-(rij_traj - each_line[3]) ** 2 / (2 * (each_line[4] ** 2)))
            counts += 1
        except:
            pass
    print(q_sum/counts)
    return q_sum/counts


def compute_traj(reference_file, trajectory_file, output):
    reference_qinterface_data = compute_qinterface_for_reference(reference_file)
    #print(reference_qinterface_data)
    n_snapshot, n_lines, n_atoms = get_pdbfile_line_atoms(trajectory_file)
    print(n_snapshot)
    print(n_lines)
    print(n_atoms)
    out = open(output, "w")
    for i in range(int(n_snapshot)):
        pdb_part = get_pdbfile_i(trajectory_file, i, n_lines)
        traj_ca_atoms_1 = get_ca_s_atoms(pdb_part, 'A')
        traj_ca_atoms_2 = get_ca_s_atoms(pdb_part, 'F')
        q = compute_qinterface_for_traj(traj_ca_atoms_1, traj_ca_atoms_2, reference_qinterface_data)
        out.write(str(round(q, 3)) + "\n")
    out.close()

def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the Q interface value between a topology and trajectory pdb file. The trajectory must be clean. Cannot change residue index, must be same")
    parser.add_argument("trajectory", help="The file name of trajectory in dcd format, usually converted from dump.lammpstrj in AWSEM", type=str)
    parser.add_argument("topology", help="The file name of topology for align, any frame of trajectory is fine", type=str)
    parser.add_argument("output", help="The file name of output",
                        type=str)
    args = parser.parse_args()
    trajectory = args.trajectory
    topology = args.topology
    output = args.output

    compute_traj(topology, trajectory, output)


if __name__ == '__main__':
    main()
