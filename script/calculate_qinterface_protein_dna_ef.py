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

def read_reference_structure_for_q_interface_calculation(oa, pdb_file, reference_chain_name="ALL", min_seq_sep=3,
                                                         max_seq_sep=np.inf, contact_threshold=0.95 * nanometers,
                                                         Qflag=0, a=0.1, removeDNAchains=True):
    # Added by Shikai on Oct/3/2020.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
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
        if chain_i.full_id[-1] == 'E':
        #if chain_i.full_id[-1] != None: # This line is currently invalid
            #chain_start_i += count_i
            for temp1 in range(0, index1):
                chain_start_i += int(math.fsum(1 for _ in chain_list[temp1].get_residues()))
            print(chain_start_i)
            count_i = 0 # Initialize after each cycle that read chain i
            ic += 1
            if removeDNAchains and np.alltrue([a.get_resname().strip() in dnaResidues for a in chain_i.get_residues()]):
                print(f"chain {chain_i.id} is a DNA chain. will be ignored for Q evaluation")
                continue
            elif removeDNAchains and np.alltrue(
                    [a.get_resname().strip() not in proteinResidues for a in chain_i.get_residues()]):
                print(f"chain {chain_i.id} is a ligand chain. will be ignored for Q evaluation")
                continue
            # print(chain)

            #print(chain_start_i)
            current_chaini_length = int(math.fsum(1 for _ in chain_i.get_residues()))
            #chain_start_j = chain_start_i + current_chaini_length # The chain j read starts from the shift of chain i plus the current chain length
            chain_start_j = 0
            jc = 0
            count_j = 0
            chain_i_count_flag = 0
            for index2 in range(0, len(chain_list)): # Don't scan repetitive chains interactions
                chain_j_count_flag = 0
                chain_j = chain_list[index2]
                if chain_j.full_id[-1] == 'F':
                #if chain_j.full_id[-1] != None: # Currently useless here
                    #chain_start_j += count_j
                    for temp2 in range(0, index2):
                        chain_start_j += int(math.fsum(1 for _ in chain_list[temp2].get_residues()))
                    print(chain_start_j)
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
                    #with open('rijn.txt', 'a') as fwrite: # Leave for debug
                    #    fwrite.write("Now the first chain is %d chain \n" % index1)
                    #    fwrite.write("Now the second chain is %d chain \n" % index2)
                    for i, residue_i in enumerate(chain_i.get_residues()):
                        if chain_i_count_flag == 0:
                        #  print(i, residue_i)
                            count_i += 1
                        for j, residue_j in enumerate(chain_j.get_residues()):
                            if chain_j_count_flag == 0:
                                count_j += 1
                        # if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                        # if abs(ic - jc) == 1 and abs(
                        #         count_j - count_i) >= min_seq_sep:  # In q interface we only count the residues in the neighboring chains and the two residue number are more than 3?
                            ca_i = residue_i['CA']
                            ca_j = residue_j['CA']

                            r_ijN = abs(ca_i - ca_j) / 10.0 * nanometers  # convert to nm
                            #r_ijN = abs(ca_i - ca_j) * angstroms
                            # if Qflag == 1 and r_ijN >= contact_threshold: continue
                            if r_ijN >= 1.2 * nanometers: continue
                            #with open('rijn.txt', 'a') as fwrite: # Leave for debug
                            #    fwrite.write(str(residue_i.get_id()[1]) + ' ' + str(residue_j.get_id()[1]) + ' ' + str(r_ijN.value_in_unit(angstrom)) + '\n')
                            sigma_ij = a * ((abs(i - j) + 1) ** 0.15)  # 0.1 nm = 1 A
                            #sigma_ij_wrong = a * (abs(i - j)  ** 0.15)  # This for testing when openmm meet 0, it will convert to a small value to enable calculation (from Wei)
                            #if sigma_ij_wrong == 0:
                            #    print("True we meet 0")
                            gamma_ij = 1.0

                            # if reference_chain_name != "ALL" and (chain.id not in reference_chain_name):
                            #	continue
                            i_index = oa.ca[i + chain_start_i]
                            j_index = oa.ca[j + chain_start_j]
                            structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                            structure_interactions.append(structure_interaction)
                        chain_j_count_flag = 1
                    chain_i_count_flag = 1

    print("Done reading")
    #with open('test.txt', 'a') as fopen:
    #    for structure_interaction in structure_interactions:
    #        fopen.write(str(structure_interaction))  # Asterisks for unpacking into function call
    #        fopen.write('\n')
    print(len(structure_interactions))
    return structure_interactions
    

def q_interface_value(oa, reference_pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf,
                      contact_threshold=0.95 * nanometers, forceGroup=4):
    ### Modified by Mingchen to compute canonical QW/QO

    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_q_interface_calculation(oa, reference_pdb_file,
                                                                                  reference_chain_name=reference_chain_name,
                                                                                  min_seq_sep=min_seq_sep,
                                                                                  max_seq_sep=max_seq_sep,
                                                                                  contact_threshold=contact_threshold,
                                                                                  Qflag=0)
    # print(len(structure_interactions))
    # print(structure_interactions)
    # create bond force for q calculation
    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)  # Asterisks for unpacking into function call
    qvalue.setForceGroup(forceGroup)
    print(qvalue.getEnergyFunction())
    return qvalue


import simtk.openmm
pdb = simtk.openmm.app.PDBFile('6n7n_merged_new_new_cg_swap.pdb')
top = pdb.topology
coord = pdb.positions
forcefield = simtk.openmm.app.ForceField('awsem.xml','3SPN2.xml')
s = forcefield.createSystem(top)

#dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb')
with open('protein.seq') as ps:
    protein_sequence_one=ps.readlines()[0]
protein = ffAWSEM.Protein.fromCoarsePDB('6n7n_merged_new_new_cg_swap.pdb', sequence=protein_sequence_one)
#dna.periodic=False
protein.periodic = False

#input_pdb_filename = '6n7n_merged_new_new_cg_swap.pdb'
#protein=ffAWSEM.Protein.fromCoarsePDB('6n7n_merged_new_new_cg_swap.pdb',sequence='crystal_structure.fasta')
#protein.periodic=False

#chain = getAllChains("6n7n_merged_new_new_cg_swap.pdb")
#oa = OpenMMAWSEMSystem(input_pdb_filename, chains='ABCDEF', k_awsem=1.0, xml_filename=f"{OPENAWSEM_LOCATION}/awsem.xml", seqFromPdb=read_fasta('crystal_structure.fasta'), includeLigands=False)  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
#print(oa.ca[0])
#print(f"using force setup file from {forceSetupFile}")
#structure_interactions = read_reference_structure_for_q_interface_calculation(oa, "crystal_structure-cleaned.pdb")
#print(type(structure_interactions[-1]))
#q_value = q_interface_value(oa, "6n7t_merged_new_new_cg_swap.pdb")
#print(q_value.getNumBonds())


def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = []
    if computeQ:
        forces.append(q_interface_value(oa, "6n7n_merged_new_new_cg_swap.pdb", reference_chain_name="ABCDEF", forceGroup=1))
        # forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", residueIndexGroup=list(range(0, 15)), forceGroup=1))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces


forces = set_up_forces(protein, computeQ=True)
for i, (force) in enumerate(forces):
    s.addForce(force)

#protein.addForcesWithDefaultForceGroup(forces)

collision_rate = 5.0 / picoseconds

#pdb_trajectory = md.load('./movie_short.pdb', stride=1)
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
        state = simulation.context.getState(getEnergy=True, groups={1})
        termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        e.append(termEnergy)
        fwrite.writelines("{0:<8.2f}\n".format(e[0]))
