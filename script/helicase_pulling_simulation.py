import open3SPN2
import simtk.openmm
import simtk.openmm.app
import simtk.unit
import sys
import numpy as np
import matplotlib.pyplot as plt
from simtk.openmm.app.metadynamics import Metadynamics
import matplotlib
%matplotlib inline

# Initialize the DNA from a sequence.
# DNA type can be changed to 'A' or 'B'

input_structure = 'clean.pdb'

#seq='ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
dna=open3SPN2.DNA.fromCoarsePDB(input_structure, dna_type='B')
#seq='ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
#dna=open3SPN2.DNA.fromSequence(seq,dna_type='B')

# Compute the topology for the DNA structure.
# Since the dna was generated from the sequence using X3DNA,
# it is not necesary to recompute the geometry.

dna.computeTopology(template_from_X3DNA=False)

# Create the system.
# To sets.append(simtk.openmm.app.CheckpointReporter('checkpnt_1.chk', 2000))
#Run simulation
simulation.step(100000) periodic boundary conditions (periodicBox=[50,50,50]).
# The periodic box size is in nanometers.
dna.periodic=False
s=open3SPN2.System(dna, periodicBox=None)

#Add 3SPN2 forces
#s.add3SPN2forces(verbose=True)
import simtk.unit as unit

def pulling_term(oa, k_pulling=4.184, forceDirect="x", appliedToResidue=1, forceGroup=19):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    
    k_pulling *= 1.0
    pulling = simtk.openmm.openmm.CustomExternalForce(f"(-{k_pulling})*({forceDirect})")
    print(oa)
    for i in range(oa.natoms):
        if appliedToResidue == "LAST":
            appliedToResidue = oa.nres
        if appliedToResidue == "FIRST":
            appliedToResidue = 1
        if oa.resi[i] == (appliedToResidue-1):
            pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    pulling.setForceGroup(forceGroup)
    return pulling

#cv2 = simtk.openmm.openmm.CustomExternalForce(f"k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
cv2 = simtk.openmm.openmm.CustomExternalForce(f"-k*x")

cv2.addGlobalParameter("k", 20)
particle_list = [1,2,3,4,5,6,7,8,9,10]
for atom in particle_list:
    cv2.addParticle(atom) # Atom index starts from 0
s.addForce(cv2)
#psi = simtk.openmm.app.metadynamics.BiasVariable(cv2, -1, 1, 0.05, True, gridWidth=20)

forces = dict(Bond=open3SPN2.Bond,
              Angle=open3SPN2.Angle,
              Stacking=open3SPN2.Stacking,
              Dihedral=open3SPN2.Dihedral,
              BasePair=open3SPN2.BasePair,
              CrossStacking=open3SPN2.CrossStacking,
              Exclusion=open3SPN2.Exclusion,
              Electrostatics=open3SPN2.Electrostatics,
             )

verbose=True
for force_name in forces:
    if verbose:
        print(force_name)
    force = forces[force_name](s.dna)
    if force_name in ['BasePair', 'CrossStacking']:
        force.addForce(s)
    elif force_name in ['Electrostatics']:
        force = forces[force_name](s.dna,salt_concentration=150*unit.millimolar)
        s.addForce(force)
    else:
        s.addForce(force)
pulling = pulling_term(s, k_pulling=5*4.184, forceDirect="x", appliedToResidue=1)
force.addForce(pulling)
#Initialize Molecular Dynamics simulations
s.initializeMD(temperature=300 * simtk.unit.kelvin,platform_name='CUDA')
simulation=s.simulation

#Set initial positionsdef __getattr__(self, item):
simulation.context.setPositions(s.coord.getPositions())

energy_unit=simtk.openmm.unit.kilojoule_per_mole
#Total energy
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(energy_unit)
print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

#Detailed energy
energies = {}
for force_name, force in s.forces.items():
    group=force.getForceGroup()
    state = simulation.context.getState(getEnergy=True, groups=2**group)
    energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

for force_name in s.forces.keys():
    print(force_name, round(energies[force_name],6),energy_unit.get_symbol())

#Add simulation reporters
dcd_reporter=simtk.openmm.app.DCDReporter(f'output_pulling.dcd', 2000)
pdb_reporter=simtk.openmm.app.PDBReporter(f'output_pulling.pdb', 2000)
energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 2000, step=True,time=True, potentialEnergy=True, temperature=True)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(energy_reporter)
simulation.minimizeEnergy()
simulation.reporters.append(simtk.openmm.app.CheckpointReporter('checkpnt_1.chk', 2000))
#Run simulation
simulation.step(100000)