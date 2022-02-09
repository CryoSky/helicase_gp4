from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from Bio.PDB.PDBParser import PDBParser
import math


def read_reference_structure_for_qdiff(pdb_file, chain_a, chain_b, a=0.1, removeDNAchains=True): # a is converter
    structure_interactions = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    model = structure[0]
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY',
                       'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    proteinResidues += ["NGP", "IGL", "IPR"]
    rnaResidues = ['A', 'G', 'C', 'U', 'I'] # Leave for future use
    dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

    chain_start_i = 0  # record the start residue of each chain
    count_i = 0  # count on the current index of the residue inside each chain
    ic = 0  # ic indicates the index of each chain
    chain_list = list(model.get_chains()) # given the number of chains, for 6n7n is 6

    for index1 in range(0, len(chain_list)):
        chain_i = chain_list[index1]
        if chain_i.full_id[-1] == chain_a:
            for temp1 in range(0, index1):
                chain_start_i += int(math.fsum(1 for _ in chain_list[temp1].get_residues()))
            #print(chain_start_i)
            count_i = 0 # Initialize after each cycle that read chain i
            ic += 1
            if removeDNAchains and np.alltrue([a.get_resname().strip() in dnaResidues for a in chain_i.get_residues()]):
                print(f"chain {chain_i.id} is a DNA chain. will be ignored for Q evaluation")
                continue
            elif removeDNAchains and np.alltrue(
                    [a.get_resname().strip() not in proteinResidues for a in chain_i.get_residues()]):
                print(f"chain {chain_i.id} is a ligand chain. will be ignored for Q evaluation")
                continue
            current_chaini_length = int(math.fsum(1 for _ in chain_i.get_residues()))
            #chain_start_j = chain_start_i + current_chaini_length # The chain j read starts from the shift of chain i plus the current chain length

            chain_start_j = 0
            jc = 0
            count_j = 0
            chain_i_count_flag = 0
            for index2 in range(0, len(chain_list)): # Don't scan repetitive chains interactions
                chain_j_count_flag = 0
                chain_j = chain_list[index2]
                if chain_j.full_id[-1] == chain_b:
                    for temp2 in range(0, index2):
                        chain_start_j += int(math.fsum(1 for _ in chain_list[temp2].get_residues()))
                    #print(chain_start_j)
                    count_j = 0
                    jc += 1
                    if removeDNAchains and np.alltrue(
                            [a.get_resname().strip() in dnaResidues for a in chain_j.get_residues()]):
                        print(f"chain {chain_j.id} is a DNA chain. will be ignored for Q evaluation")
                        continue
                    elif removeDNAchains and np.alltrue(
                            [a.get_resname().strip() not in proteinResidues for a in chain_j.get_residues()]):
                        print(f"chain {chain_j.id} is a ligand chain. will be ignored for Q evaluation")
                        continue
                    for i, residue_i in enumerate(chain_i.get_residues()):
                        if chain_i_count_flag == 0:
                            count_i += 1
                        for j, residue_j in enumerate(chain_j.get_residues()):
                            if chain_j_count_flag == 0:
                                count_j += 1
                            ca_i = residue_i['CA']
                            ca_j = residue_j['CA']
                            #print(ca_j.get_serial_number())

                            r_ijN = abs(ca_i - ca_j) / 10.0 * nanometers  # convert to nm
                            #if r_ijN >= 1.2 * nanometers: continue
                            sigma_ij = a * ((abs(i - j) + 1) ** 0.15)  # 0.1 nm = 1 A
                            gamma_ij = 1.0
                            #i_index = oa.ca[i + chain_start_i]
                            #j_index = oa.ca[j + chain_start_j]
                            i_index = ca_i.get_serial_number() - 1
                            j_index = ca_j.get_serial_number() - 1
                            #print(ca_j.get_serial_number() - j_index=1)
                            structure_interaction = [i_index, j_index, gamma_ij, r_ijN, sigma_ij]
                            structure_interactions.append(structure_interaction)
                        chain_j_count_flag = 1
                    chain_i_count_flag = 1
    return structure_interactions

def q_diff_value(pdb_file_1, pdb_file_2, chain_a = 'A', chain_b = 'F', forceGroup=31):
    structure_interactions_1 = read_reference_structure_for_qdiff(pdb_file_1, chain_a, chain_b)
    #print("length of structure 1")
    #print(len(structure_interactions_1))
    structure_interactions_2 = read_reference_structure_for_qdiff(pdb_file_2, chain_a, chain_b)
    #print("length of structure 2")
    #print(len(structure_interactions_2))
    structure_interactions_all = []
    for index in range(len(structure_interactions_1)):
        structure_interactions_all.append([structure_interactions_1[index][0], structure_interactions_1[index][1], [structure_interactions_1[index][2], structure_interactions_1[index][3], structure_interactions_2[index][3], structure_interactions_1[index][4]]])
    #structure_interactions_all = pd.merge(structure_interactions_1, structure_interactions_2, on=['i_index','j_index'])
    #print(structure_interactions_all)

    normalization = len(structure_interactions_all)
    qdiffvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*(exp(-(r-r_ijN1)^2/(2*sigma_ij^2))-exp(-(r-r_ijN2)^2/(2*sigma_ij^2)))")
    qdiffvalue.addPerBondParameter("gamma_ij")
    qdiffvalue.addPerBondParameter("r_ijN1")
    qdiffvalue.addPerBondParameter("r_ijN2")
    qdiffvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions_all:
        qdiffvalue.addBond(*structure_interaction) #Asterisks for unpacking into function call
    qdiffvalue.setForceGroup(forceGroup)
    #print(qdiffvalue.getEnergyFunction())
    return qdiffvalue

def q_diff1_value(pdb_file_1, pdb_file_2, chain_a = 'A', chain_b = 'F', forceGroup=30):
    structure_interactions_1 = read_reference_structure_for_qdiff(pdb_file_1, chain_a, chain_b)
    #print("length of structure 1")
    #print(len(structure_interactions_1))
    structure_interactions_2 = read_reference_structure_for_qdiff(pdb_file_2, chain_a, chain_b)
    #print("length of structure 2")
    #print(len(structure_interactions_2))
    structure_interactions_all = []
    for index in range(len(structure_interactions_1)):
        structure_interactions_all.append([structure_interactions_1[index][0], structure_interactions_1[index][1], [structure_interactions_1[index][2], structure_interactions_1[index][3], structure_interactions_2[index][3], structure_interactions_1[index][4]]])
    #structure_interactions_all = pd.merge(structure_interactions_1, structure_interactions_2, on=['i_index','j_index'])
    #print(structure_interactions_all)

    normalization = len(structure_interactions_all)
    qdiffvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*(exp(-(r_ijN1-r_ijN1)^2/(2*sigma_ij^2))-exp(-(r_ijN1-r_ijN2)^2/(2*sigma_ij^2)))")
    qdiffvalue.addPerBondParameter("gamma_ij")
    qdiffvalue.addPerBondParameter("r_ijN1")
    qdiffvalue.addPerBondParameter("r_ijN2")
    qdiffvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions_all:
        qdiffvalue.addBond(*structure_interaction) #Asterisks for unpacking into function call
    qdiffvalue.setForceGroup(forceGroup)
    #print(qdiffvalue.getEnergyFunction())
    return qdiffvalue

def q_diff2_value(pdb_file_1, pdb_file_2, chain_a = 'A', chain_b = 'F', forceGroup=29):
    structure_interactions_1 = read_reference_structure_for_qdiff(pdb_file_1, chain_a, chain_b)
    #print("length of structure 1")
    #print(len(structure_interactions_1))
    structure_interactions_2 = read_reference_structure_for_qdiff(pdb_file_2, chain_a, chain_b)
    #print("length of structure 2")
    #print(len(structure_interactions_2))
    structure_interactions_all = []
    for index in range(len(structure_interactions_1)):
        structure_interactions_all.append([structure_interactions_1[index][0], structure_interactions_1[index][1], [structure_interactions_1[index][2], structure_interactions_1[index][3], structure_interactions_2[index][3], structure_interactions_1[index][4]]])
    #structure_interactions_all = pd.merge(structure_interactions_1, structure_interactions_2, on=['i_index','j_index'])
    #print(structure_interactions_all)

    normalization = len(structure_interactions_all)
    qdiffvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*(exp(-(r_ijN2-r_ijN1)^2/(2*sigma_ij^2))-exp(-(r_ijN2-r_ijN2)^2/(2*sigma_ij^2)))")
    qdiffvalue.addPerBondParameter("gamma_ij")
    qdiffvalue.addPerBondParameter("r_ijN1")
    qdiffvalue.addPerBondParameter("r_ijN2")
    qdiffvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions_all:
        qdiffvalue.addBond(*structure_interaction) #Asterisks for unpacking into function call
    qdiffvalue.setForceGroup(forceGroup)
    #print(qdiffvalue.getEnergyFunction())
    return qdiffvalue

def qdiff_interface_term(oa, q0, reference_pdb_file_1, reference_pdb_file_2, chain_a = 'A', chain_b = 'F', k_qbias=10000, qbias_min_seq_sep=3,
                         qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.8 * nanometers, forceGroup=5):
    qdiffbias = CustomCVForce(f"""energy;
                            energy=0.5*k_qbias*(qdiff-q0)^2;
                            qdiff=(qdouble-q2)/(q1-q2)""")

    q_double = q_diff_value(reference_pdb_file_1, reference_pdb_file_2, chain_a = 'A', chain_b = 'F')
    q_diff1 = q_diff1_value(reference_pdb_file_1, reference_pdb_file_2, chain_a='A', chain_b='F')
    q_diff2 = q_diff2_value(reference_pdb_file_1, reference_pdb_file_2, chain_a='A', chain_b='F')
    qdiffbias.addCollectiveVariable("qdouble", q_double)
    qdiffbias.addCollectiveVariable("q1", q_diff1)
    qdiffbias.addCollectiveVariable("q2", q_diff2)
    qdiffbias.addGlobalParameter("k_qbias", k_qbias)
    qdiffbias.addGlobalParameter("q0", q0)
    qdiffbias.setForceGroup(forceGroup)
    return qdiffbias
