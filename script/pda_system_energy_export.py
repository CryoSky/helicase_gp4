#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2021-Apr-22, latest modified on 2021-Apr-22
# Output the final energy of protein+DNA+ligand system
# Version 0.0.2 Comment out the print of time step
# Version 1.0.0 Mingchen mentioned you should include ALL of your energy terms
# Example in Linux: python this_code.py 0.02

__author__ = 'Shikai Jin'
__version__ = '1.0.0'
__date__ = "Jul-26-2021"

import ffcgligand
import pdbfixer
import ffAWSEM
import pandas
import open3SPN2
import sys
import simtk.openmm
import LigandforceTerm
import numpy as np
import mdtraj as md

input_traj = sys.argv[1]
pdb_trajectory = md.load(input_traj, top=input_traj, stride=1)
input_topology = sys.argv[2]
q_bin_min = sys.argv[3]
q_bin_min = float(q_bin_min)

# Create the merged system, here pdb refers all atoms in the input
pdb = simtk.openmm.app.PDBFile(input_topology)
top = pdb.topology
coord = pdb.positions
forcefield = simtk.openmm.app.ForceField('awsem.xml', './open3SPN2/3SPN2.xml', 'cgligand.xml')
s = forcefield.createSystem(top)

# Get the three different parts of clean pdb
atp = ffcgligand.Ligand.fromCoarsePDB(input_topology)
atp.parseConfigurationFile(configuration_file=f'./cgligand.conf')
dna = open3SPN2.DNA.fromCoarsePDB(input_topology)
with open('protein.seq') as ps:
    protein_sequence_one = ps.readlines()[0]
protein = ffAWSEM.Protein.fromCoarsePDB(input_topology,sequence=protein_sequence_one)
#print(protein.atom_lists)
protein.periodic = False
atp.periodic = False
dna.periodic = False


# Set up forces for the whole system
forces = {}

# Add ligand only forces
ligand_forces = dict(Ligandbone=LigandforceTerm.Ligandbone)
for ligand_force_name in ligand_forces:
    print(ligand_force_name)
    ligand_force_current = ligand_forces[ligand_force_name](atp)
    s.addForce(ligand_force_current)
forces.update({ligand_force_name: ligand_force_current})

#Add 3SPN2 DNA forces
for dna_force_name in open3SPN2.forces:
    print(dna_force_name)
    dna_force_current = open3SPN2.forces[dna_force_name](dna)
    if dna_force_name in ['BasePair','CrossStacking']:
        dna_force_current.addForce(s)
    elif dna_force_name in ['Exclusion', 'Electrostatics']: # Keep the exclusion list consistent in CUDA/OpenCL drive
        LigandforceTerm.addNonBondedExclusions(atp, dna_force_current)
        s.addForce(dna_force_current)
    else:
        s.addForce(dna_force_current)
    forces.update({dna_force_name: dna_force_current})

#Add AWSEM protein forces
openAWSEMforces = dict(Connectivity=ffAWSEM.functionTerms.basicTerms.con_term,
                       Chain=ffAWSEM.functionTerms.basicTerms.chain_term,
                       Chi=ffAWSEM.functionTerms.basicTerms.chi_term,
                       Excl=ffAWSEM.functionTerms.basicTerms.excl_term_v2,
                       rama=ffAWSEM.functionTerms.basicTerms.rama_term,
                       rama_pro=ffAWSEM.functionTerms.basicTerms.rama_proline_term,
                       rama_ss=ffAWSEM.functionTerms.basicTerms.rama_ssweight_term,
                       contact=ffAWSEM.functionTerms.contactTerms.contact_term,
                       beta1 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_1,
                       beta2 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_2,
                       beta3 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_3,
                       pap1 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_1,
                       pap2 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_2,
                       fragment=ffAWSEM.functionTerms.templateTerms.fragment_memory_term,
                       #qinterface = ffAWSEM.functionTerms.qInterfaceTerms.qbias_interface_term,
                       qdiff = ffAWSEM.functionTerms.qDiffTerms.qdiff_interface_term,
                       )
protein.setup_virtual_sites(s)
for protein_force_name in openAWSEMforces:
    print(protein_force_name)
    if protein_force_name in ['contact']:
        protein_force_current = openAWSEMforces[protein_force_name](protein, withExclusion=False, periodic=False)
        #print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna, protein_force_current)
        LigandforceTerm.addNonBondedExclusions(atp, protein_force_current)
        #print(force.getNumExclusions())
    elif protein_force_name in ['Excl']:
        protein_force_current = openAWSEMforces[protein_force_name](protein)
        #print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna, protein_force_current)
        LigandforceTerm.addNonBondedExclusions(atp, protein_force_current)
        #print(force.getNumExclusions())
    elif protein_force_name in ['fragment']:
        protein_force_current = openAWSEMforces[protein_force_name](protein, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=True)
    elif protein_force_name in ['qinterface']:
        protein_force_current = openAWSEMforces[protein_force_name](protein, q_bin_min, reference_pdb_file="../model_050_protein_cg_correct.pdb", reference_chain_name="ALL", k_qbias=10000)
    elif protein_force_name in ['qdiff']:
        protein_force_current = openAWSEMforces[protein_force_name](protein, q_bin_min, reference_pdb_file_1="../model_001_protein_cg_correct.pdb", reference_pdb_file_2="../model_050_protein_cg_correct.pdb", k_qbias=1000)
    else:
        protein_force_current = openAWSEMforces[protein_force_name](protein)
    s.addForce(protein_force_current)
    forces.update({protein_force_name: protein_force_current})


#Add protein-ligand interaction forces
protein_ligand_forces=dict(
    ExclusionProteinLigand=LigandforceTerm.ExclusionProteinLigand,
    ElectrostaticsProteinLigand=LigandforceTerm.ElectrostaticsProteinLigand,
    )

# external_distance_bias_forces=dict(
LigandDistanceRestraint1=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint2=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint3=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint4=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint5=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint6=LigandforceTerm.LigandDistanceRestraint
)

chain_match = {'A':'M', 'B':'H', 'C':'I', 'D':'J' ,'E':'K', 'F':'L'}


for proligand_force_name, chain_pair in zip(external_distance_bias_forces, chain_match.items()):
    print(proligand_force_name)
    proligand_force = external_distance_bias_forces[proligand_force_name](atp, protein, chain_pair)
    s.addForce(proligand_force)
    forces.update({proligand_force_name: proligand_force})


for proligand_force_name in protein_ligand_forces:
    print(proligand_force_name)
    if proligand_force_name == 'LigandDistanceRestraint':
        proligand_force = protein_ligand_forces[proligand_force_name](atp, protein)
    else:
        proligand_force = protein_ligand_forces[proligand_force_name](atp, protein)
        open3SPN2.addNonBondedExclusions(dna, proligand_force)
    s.addForce(proligand_force)
    forces.update({proligand_force_name: proligand_force})

#Add protein-DNA interaction forces
for prodna_force_name in open3SPN2.protein_dna_forces:
    print(prodna_force_name)
    prodna_force = open3SPN2.protein_dna_forces[prodna_force_name](dna, protein)
    LigandforceTerm.addNonBondedExclusions(atp, prodna_force)
    s.addForce(prodna_force)
    forces.update({prodna_force_name: prodna_force})

# Now set up the simulation system
temperature=300 * simtk.openmm.unit.kelvin
platform_name='CPU'
#platform_name='CUDA'

integrator = simtk.openmm.LangevinIntegrator(temperature, 1 / simtk.openmm.unit.picosecond, 2 * simtk.openmm.unit.femtoseconds)
platform = simtk.openmm.Platform.getPlatformByName(platform_name)
simulation = simtk.openmm.app.Simulation(top, s, integrator, platform)
simulation.context.setPositions(coord)
energy_unit=simtk.openmm.unit.kilojoule_per_mole
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(energy_unit)
# The energy of the starting structure
# print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

energies = {}
# for force_name, force in forces.items():
#     #force = force(atp)
#     group=force.getForceGroup()
#     state = simulation.context.getState(getEnergy=True, groups=2**group)
#     energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

# for force_name in forces.keys():
#     print(force_name, round(energies[force_name],6),energy_unit.get_symbol())

#for step in range(1):
for step in range(len(pdb_trajectory)):
    # print(pdb_trajectory[step].time)
    simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    #print('TotalEnergy',round(energy,6),energy_unit.get_symbol())
    print(round(energy,6))
    energies = {}
#     for force_name, force in forces.items():
#         group=force.getForceGroup()
#         state = simulation.context.getState(getEnergy=True, groups=2**group)
#         energies[force_name] = state.getPotentialEnergy().value_in_unit(energy_unit)

#     for force_name in forces.keys():
#         output_energy_value = energies[force_name]
#         print(force_name, round(output_energy_value,8),energy_unit.get_symbol())
