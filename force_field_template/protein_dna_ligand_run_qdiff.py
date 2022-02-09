#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2020-Nov-8, latest modified on 2021-Jul-21
# A final loading of whole helicase system: protein+DNA+ATP+Mg
# Example in Linux: python this_code.py 0.02

# Version 1.1.1 Add LigandDistanceRestraint to deal with problem in movie 30/31
# Version 1.1.3 Fix the bug of Ligandbond (not included MG), add strength to distance bias
# Version 1.1.4 Forget to tune the strength of DH from 1 to 10.5, also change bias strength to 2000
# Version 1.1.5 Set the RANDOM_SEED
# Version 1.1.6 Fixed the bug of RNADOM_SEED, rewrite k_pdexcl part
# Version 1.1.7 After discussion with Mingchen, change # steps and k_qdiff
__author__ = 'Shikai Jin'
__version__ = '1.1.7'
__date__ = "Jul-21-2021"

import ffcgligand
import pdbfixer
import ffAWSEM
import pandas
import open3SPN2
import sys
import simtk.openmm
import LigandforceTerm
import numpy as np

# Input value for q diff term
q_bin_min = sys.argv[1]
q_bin_min = float(q_bin_min)

# Create the merged system, here pdb refers all atoms in the input
pdb = simtk.openmm.app.PDBFile('clean.pdb')
top = pdb.topology
coord = pdb.positions
forcefield = simtk.openmm.app.ForceField('awsem.xml', './open3SPN2/3SPN2.xml', 'cgligand.xml')
s = forcefield.createSystem(top)

# Get the three different parts of clean pdb
atp = ffcgligand.Ligand.fromCoarsePDB('clean.pdb')
atp.parseConfigurationFile(configuration_file=f'./cgligand.conf')
dna = open3SPN2.DNA.fromCoarsePDB('clean.pdb')
with open('protein.seq') as ps:
    protein_sequence_one = ps.readlines()[0]
protein = ffAWSEM.Protein.fromCoarsePDB('clean.pdb',sequence=protein_sequence_one)
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
    elif dna_force_name in ['Exclusion', 'Electrostatics']: # Keep the exclusion list consistent in OpenCL/OpenCL drive
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
        protein_force_current = openAWSEMforces[protein_force_name](protein, q_bin_min, reference_pdb_file_1="../model_001_protein_cg_correct.pdb", reference_pdb_file_2="../model_050_protein_cg_correct.pdb", k_qbias=2000)
    else:
        protein_force_current = openAWSEMforces[protein_force_name](protein)
    s.addForce(protein_force_current)
    forces.update({protein_force_name: protein_force_current})


#Add protein-ligand interaction forces
protein_ligand_forces=dict(
    ExclusionProteinLigand=LigandforceTerm.ExclusionProteinLigand,
    ElectrostaticsProteinLigand=LigandforceTerm.ElectrostaticsProteinLigand,
    )

external_distance_bias_forces=dict(
LigandDistanceRestraint1=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint2=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint3=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint4=LigandforceTerm.LigandDistanceRestraint,
LigandDistanceRestraint5=LigandforceTerm.LigandDistanceRestraint,
#LigandDistanceRestraint6=LigandforceTerm.LigandDistanceRestraint
)

#chain_match = {'A':'M', 'B':'H', 'C':'I', 'D':'J' ,'E':'K', 'F':'L'}
chain_match = {'A':'M', 'B':'H', 'C':'I', 'D':'J' ,'E':'K'}

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
    #print(force.getExclusionParticles(0))
    #for i in range(force.getNumExclusions()):
    #    print(force.getExclusionParticles(i))
    #print(force.getNumExclusions())
    s.addForce(proligand_force)
    forces.update({proligand_force_name: proligand_force})

#Add protein-DNA interaction forces
for prodna_force_name in open3SPN2.protein_dna_forces:
    print(prodna_force_name)
    prodna_force = open3SPN2.protein_dna_forces[prodna_force_name](dna, protein)
    LigandforceTerm.addNonBondedExclusions(atp, prodna_force)
    #print(force.getExclusionParticles(0))
    #for i in range(force.getNumExclusions()):
    #    print(force.getExclusionParticles(i))
    #print(force.getNumExclusions())
    s.addForce(prodna_force)
    forces.update({prodna_force_name: prodna_force})


# Now set up the simulation system
temperature=300 * simtk.openmm.unit.kelvin
#platform_name='Reference'
platform_name='OpenCL'
RANDOM_SEED=19951002
integrator = simtk.openmm.LangevinIntegrator(temperature, 1 / simtk.openmm.unit.picosecond, 2 * simtk.openmm.unit.femtoseconds)
integrator.setRandomNumberSeed(RANDOM_SEED)
platform = simtk.openmm.Platform.getPlatformByName(platform_name)
simulation = simtk.openmm.app.Simulation(top, s, integrator, platform)
simulation.context.setPositions(coord)
simulation.context.setVelocitiesToTemperature(temperature, RANDOM_SEED)
energy_unit=simtk.openmm.unit.kilojoule_per_mole
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(energy_unit)
# The energy of the starting structure
print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

# Check the detail energy of each term
energies = {}
for force_name, force in forces.items():
    #force = force(atp)
    group=force.getForceGroup()
    state = simulation.context.getState(getEnergy=True, groups=2**group)
    energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

for force_name in forces.keys():
    print(force_name, round(energies[force_name],6),energy_unit.get_symbol())


# Set up simulation coordinates
dcd_reporter=simtk.openmm.app.DCDReporter(f'output.dcd', 1000)
pdb_reporter=simtk.openmm.app.PDBReporter(f'output.pdb', 1000)
energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 1000, step=True,time=True,
                                                   potentialEnergy=True, temperature=True)
checkpoint_reporter=simtk.openmm.app.CheckpointReporter(f'checkpoint.chk', 1000)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(energy_reporter)
simulation.reporters.append(checkpoint_reporter)
#Run simulation
simulation.minimizeEnergy()
#simulation.context.setVelocitiesToTemperature(temperature)
#simulation.step(500000)
#for i in range(100): 
#    print(i) 
#    positions = simulation.context.getState(getPositions=True).getPositions() 
#    print(min(positions)) 
#    print(max(positions#)) 
#    simulation.step(1)

for i in range(1,21):
    #simulation.minimizeEnergy()
    simulation.step(1000)
    simulation.context.setParameter('k_pdexcl', 10**(-5*(20-i)/20))
    #update_1 = open3SPN2.protein_dna_forces['ExclusionProteinDNA'](protein, dna)
    k_output = simulation.context.getParameter('k_pdexcl')
    print("Current k_excl value is: %s" %str(k_output))
    #state = simulation.context.getState(getEnergy=True, groups=14)
    #print(state.getPotentialEnergy().value_in_unit(energy_unit))
simulation.step(780000)
