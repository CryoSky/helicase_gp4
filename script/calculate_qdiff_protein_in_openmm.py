#!/usr/bin/env python3
import os
import sys
import mdtraj as md
import math

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    #print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *
import ffAWSEM

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
    #print(structure_interactions_1[0])
    structure_interactions_2 = read_reference_structure_for_qdiff(pdb_file_2, chain_a, chain_b)
    #print("length of structure 2")
    #print(len(structure_interactions_2))
    #print(structure_interactions_2[0])
    structure_interactions_all = []
    for index in range(len(structure_interactions_1)):
        structure_interactions_all.append([structure_interactions_1[index][0], structure_interactions_1[index][1], [structure_interactions_1[index][2], structure_interactions_1[index][3], structure_interactions_2[index][3], structure_interactions_1[index][4]]])
    #structure_interactions_all = pd.merge(structure_interactions_1, structure_interactions_2, on=['i_index','j_index'])
    #print(structure_interactions_all)

    normalization = len(structure_interactions_all)
    qdiffvalue = CustomBondForce(f"""energy;
                                    temp=100;
                                    energy=(1/{normalization})*gamma_ij*(exp(-(r-r_ijN1)^2/(2*sigma_ij^2))-exp(-(r-r_ijN2)^2/(2*sigma_ij^2)))""")
    qdiffvalue.addPerBondParameter("gamma_ij")
    qdiffvalue.addPerBondParameter("r_ijN1")
    qdiffvalue.addPerBondParameter("r_ijN2")
    qdiffvalue.addPerBondParameter("sigma_ij")
    #print(structure_interactions_all[0])
    for structure_interaction in structure_interactions_all:
        qdiffvalue.addBond(*structure_interaction) #Asterisks for unpacking into function call
    qdiffvalue.setForceGroup(forceGroup)
    #print(qdiffvalue.getEnergyFunction())
    return qdiffvalue


import simtk.openmm
pdb = simtk.openmm.app.PDBFile('clean.pdb')
top = pdb.topology
coord = pdb.positions
forcefield = simtk.openmm.app.ForceField('awsem.xml','3SPN2.xml', 'cgligand.xml')
s = forcefield.createSystem(top)

with open('protein.seq') as ps:
    protein_sequence_one=ps.readlines()[0]
protein = ffAWSEM.Protein.fromCoarsePDB('clean.pdb', sequence=protein_sequence_one)
protein.periodic = False

def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = []
    if computeQ:
        #print("addforce now")
        forces.append(q_diff_value('../model_001_protein_cg.pdb', "../model_050_protein_cg.pdb"))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces


forces = set_up_forces(protein, computeQ=True)
for i, force in enumerate(forces):
    s.addForce(force)

collision_rate = 5.0 / picoseconds

input = sys.argv[1]
output = sys.argv[2]
pdb_trajectory = md.load(input, stride=1)

platform = Platform.getPlatformByName('CUDA')
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(topology, s, integrator, platform)

with open(output, 'w') as fwrite:
    for step in range(len(pdb_trajectory)):
        e = []
        simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
        state = simulation.context.getState(getEnergy=True, groups={31})
        termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        e.append(termEnergy)
        print(e)
        fwrite.writelines("{0:<8.2f}\n".format(e[0]))
