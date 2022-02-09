# 2020-Nov-5 Ver 0.1.0
# 2020-Oct-11 Ver 1.0.0
# No need to add object like class Force(object) in Python 3
# 2021-May-17 Ver 1.1.5 Change the dielectric constant based on other cg models, change back ATP hydrolysis questions
# 2021-Aug-13 Ver 1.2.0 Big problem solved! When use globalparameter remember to change name!

__author__ = 'Shikai Jin'
__version__ = '1.2.0'
__date__ = "Aug-13-2021"

import simtk.openmm.app
import simtk.openmm
import simtk.unit as unit
import numpy as np
import itertools
import pandas

_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor
_df = 1 * unit.angstrom / unit.nanometer  # distance scaling factor
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_complement = {'DA': 'DT', 'DT': 'DA', 'DG': 'DC', 'DC': 'DG', 'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_rnaResidues = ['A', 'C', 'U', 'G']
_ligandResidues = ['ATP', 'ADP', 'GTP', 'TTP', 'MG', 'ATM', 'ADM']
_proteinResidues = ['IPR', 'IGL', 'NGP']

class Force():
    """ Wrapper for the openMM force. """

    def __init__(self, ligand):
        self.periodic = ligand.periodic
        self.force = None
        self.ligand = ligand

        # Followed by Carlos's design, every forces should have reset and defineInteraction functions.
        # Define the dna force
        self.reset()

        # Define the interaction pairs
        self.defineInteraction()

    #Called when an attribute lookup has not found the attribute in the usual places 
    #(i.e. it is not an instance attribute nor is it found in the class tree for self). 
    # name is the attribute name. This method should return the (computed) attribute value or raise an AttributeError exception.
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        elif 'force' in self.__dict__:
            try:
                return getattr(self.force, attr)
            except:
                pass
        else:
            if '__repr__' in self.__dict__:
                raise AttributeError(f"type object {str(self)} has no attribute {str(attr)}")
            else:
                raise AttributeError()


class Ligandbone(Force, simtk.openmm.CustomBondForce):
    def __init__(self, ligand, force_group=29):
        self.force_group = force_group
        super().__init__(ligand)

    def getParameterNames(self):
        self.perInteractionParameters = []
        self.GlobalParameters = []
        for i in range(self.force.getNumPerBondParameters()):
            self.perInteractionParameters += [self.force.getPerBondParameterName(i)]
        for i in range(self.force.getNumGlobalParameters()):
            self.GlobalParameters += [self.force.getGlobalParameterName(i)]
        return [self.perInteractionParameters, self.GlobalParameters]

    def reset(self):
        bondForce = simtk.openmm.CustomBondForce("Kb2*(r-r0)^2+Kb3*0+Kb4*0")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        bondForce.setUsesPeriodicBoundaryConditions(self.periodic)
        bondForce.setForceGroup(self.force_group)
        self.force = bondForce

        # epsilon = 2
        # delta = 0.25
        # backbone = 0.76 * 8.4

    def defineInteraction(self):
        for i, b in self.ligand.bonds.iterrows():
            # Units converted from
            parameters = [b['r0'] * _df,
                          b['Kb2'] / _df ** 2 * _ef,
                          b['Kb3'] / _df ** 3 * _ef,
                          b['Kb4'] / _df ** 4 * _ef]
            #parameters = [2, 0.76*8.4 * _df, 0.25 * _df]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)

def addNonBondedExclusions(ligand, force, OpenCLPatch=True):
    is_ligand = ligand.atoms['resname'].isin(_ligandResidues)
    atoms = ligand.atoms.copy()
    selection = atoms[is_ligand]
    for (i, atom_a), (j, atom_b) in itertools.combinations(selection.iterrows(), r=2):
        if j < i:
            i, j = j, i
            atom_a, atom_b = atom_b, atom_a
        # Neighboring residues
        if atom_a.chainID == atom_b.chainID and (abs(atom_a.resSeq - atom_b.resSeq) <= 1):
            force.addExclusion(i, j)
            # print(i, j)
        # Base-pair residues
        elif OpenCLPatch and (atom_a['name'] in _complement.keys()) and (atom_b['name'] in _complement.keys()) and (
                atom_a['name'] == _complement[atom_b['name']]):
            force.addExclusion(i, j)
            # print(i, j)

class ProteinLigandForce(Force):
    def __init__(self, ligand, protein):
        self.protein = protein
        super().__init__(ligand)


class ExclusionProteinLigand(ProteinLigandForce):
    """ Protein-Ligand exclusion potential"""
    def __init__(self, ligand, protein, k=1, force_group=30):
        self.k = k
        self.force_group = force_group
        super().__init__(ligand, protein)

    def reset(self):
        k = self.k
        exclusionForce = simtk.openmm.CustomNonbondedForce(f"""{k}*energy;
                                                            energy=(4*epsilon*((sigma/r)^12-(sigma/r)^6)-offset)*step(cutoff-r);
                                                            offset=4*epsilon*((sigma/cutoff)^12-(sigma/cutoff)^6);
                                                            sigma=0.5*(sigma1+sigma2); 
                                                            epsilon=sqrt(epsilon1*epsilon2);
                                                            cutoff=sqrt(cutoff1*cutoff2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.addPerParticleParameter('cutoff')
        exclusionForce.setCutoffDistance(1.55)
        # exclusionForce.setUseLongRangeCorrection(True)
        exclusionForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        self.force = exclusionForce

    def defineInteraction(self):
        particle_definition = self.ligand.config['Protein-Ligand particles']
        ligand_particle_definition=particle_definition[(particle_definition['molecule'] == 'Ligand')]
        protein_particle_definition = particle_definition[(particle_definition['molecule'] == 'Protein')]

        # Merge Ligand and protein particle definitions
        particle_definition = pandas.concat([ligand_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name
        self.particle_definition = particle_definition

        is_ligand = self.ligand.atoms['resname'].isin(_ligandResidues)
        is_protein = self.ligand.atoms['resname'].isin(_proteinResidues)
        atoms = self.ligand.atoms.copy()
        atoms['is_ligand'] = is_ligand
        atoms['is_protein'] = is_protein
        atoms['epsilon']=np.nan
        atoms['radius']=np.nan
        atoms['cutoff'] = np.nan
        ligand_list = []
        protein_list = []
        for i, atom in atoms.iterrows():
            if atom.is_ligand:
                param = particle_definition.loc['Ligand' + atom['name']]
                parameters = [param.epsilon * _ef,
                              param.radius * _df,
                              param.cutoff * _df]
                ligand_list += [i]
            elif atom.is_protein:
                param = particle_definition.loc['Protein' + atom['name']]
                parameters = [param.epsilon * _ef,
                              param.radius * _df,
                              param.cutoff * _df]
                protein_list += [i]
            else:
                #print(f'Residue {i} not included in protein-ligand exclusions')
                parameters = [0, .1,.1]
            atoms.loc[i, ['epsilon', 'radius', 'cutoff']] = parameters
            self.atoms = atoms
            self.force.addParticle(parameters) # The number of particles for which you set 
            # parameters must be exactly equal to the number of particles in the System
        self.force.addInteractionGroup(ligand_list, protein_list)
        # addExclusions

        addNonBondedExclusions(self.ligand, self.force)

class ElectrostaticsProteinLigand(ProteinLigandForce):
    """Ligand-protein and protein-protein electrostatics."""
    def __init__(self, ligand, protein, k=1, force_group=31): # Igonred: 10.5 comes from the previous calculation that difference between 1 ATP and 1 ADP contributes to 5.3 kJ. 5.3*10.5=57 to the ATP hydrolysis energy.
        self.k = k
        self.force_group = force_group
        super().__init__(ligand, protein)

    def reset(self):
        #dielectric = 78 # e * a use dielectric permittivity of water
        dielectric = 8 # Based on other papers see 20210512 slide
        #print(dielectric)
        # Debye length
        Na = simtk.unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = 1.2 * unit.nanometer # np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        denominator_pa = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator_pa = denominator_pa.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator_pa)
        k = self.k
        electrostaticForce = simtk.openmm.CustomNonbondedForce(f"""k_electro_pa*energy;
                                                                energy=q1*q2*exp(-r/inter_dh_length_pa)/inter_denominator_pa/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('k_electro_pa', k) # For global parameters remember to rename them!
        electrostaticForce.addGlobalParameter('inter_dh_length_pa', ldby)
        electrostaticForce.addGlobalParameter('inter_denominator_pa', denominator_pa)

        electrostaticForce.setCutoffDistance(4)
        if self.periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        electrostaticForce.setForceGroup(self.force_group)
        self.force = electrostaticForce
        #print(self.protein.__dict__.keys())
        #print(self.ligand.__dict__.keys())

    def defineInteraction(self):
        particle_definition = self.ligand.config['Protein-Ligand particles']
        ligand_particle_definition=particle_definition[(particle_definition['molecule'] == 'Ligand')]
        protein_particle_definition = particle_definition[(particle_definition['molecule'] == 'Protein')]

        # Merge Ligand and protein particle definitions
        particle_definition = pandas.concat([ligand_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name
        self.particle_definition = particle_definition

        # Open Sequence dependent electrostatics
        sequence_electrostatics = self.ligand.config['Sequence dependent electrostatics']
        sequence_electrostatics.index = sequence_electrostatics.resname


        # Select only dna and protein atoms
        is_ligand = self.protein.atoms['resname'].isin(_ligandResidues)
        is_protein = self.protein.atoms['resname'].isin(_proteinResidues)
        atoms = self.protein.atoms.copy()
        atoms['is_ligand'] = is_ligand
        atoms['is_protein'] = is_protein
        ligand_list = []
        protein_list = []

        for i, atom in atoms.iterrows():
            if atom.is_ligand:
                param = particle_definition.loc['Ligand' + atom['name']]
                charge = param.charge
                parameters = [charge]
                if charge != 0:
                    ligand_list += [i]
                    #print(atom.chainID, atom.resSeq, atom.resname, atom['name'], charge)
            elif atom.is_protein:
                atom_param = particle_definition.loc['Protein' + atom['name']]
                seq_param = sequence_electrostatics.loc[atom.real_resname]
                charge = atom_param.charge * seq_param.charge
                parameters = [charge]
                if charge != 0:
                    protein_list += [i]
                    #print(atom.chainID, atom.resSeq, atom.resname, atom['name'], charge)
            else:
                #print(f'Residue {i} not included in protein-ligand electrostatics')
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(ligand_list, protein_list)
        # self.force.addInteractionGroup(protein_list, protein_list) #protein-protein electrostatics should be included using debye Huckel Terms

        # addExclusions
        addNonBondedExclusions(self.ligand, self.force)


class LigandDistanceRestraint(ProteinLigandForce):
    """Ligand-protein distance bias."""
    def __init__(self, ligand, protein, chain_pair=('A', 'M'), k=100, d0=1.2, force_group=25):
        self.k = k
        self.d0 = d0
        self.force_group = force_group
        self.chain_pair = chain_pair
        super().__init__(ligand, protein)

    def reset(self):
        self.d0 = self.d0 *_df  # convert to nm
        #print(self.d0)
        distanceForce = simtk.openmm.CustomCentroidBondForce(2, f"""k_constraint*energy;
                                                                    energy=0.5*(distance(g1,g2)-d0)^2;
                                                            """)
        distanceForce.addPerBondParameter('k_constraint') # Use global parameter maybe a good idea?
        distanceForce.addPerBondParameter('d0')
        distanceForce.setUsesPeriodicBoundaryConditions(self.periodic)
        self.force = distanceForce

    def defineInteraction(self):
        particle_definition = self.ligand.config['Protein-Ligand particles']
        ligand_particle_definition = particle_definition[(particle_definition['molecule'] == 'Ligand')]
        protein_particle_definition = particle_definition[(particle_definition['molecule'] == 'Protein')]

        # Merge ligand and protein particle definitions
        particle_definition = pandas.concat([ligand_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name
        self.particle_definition = particle_definition

        is_ligand = self.ligand.atoms['resname'].isin(_ligandResidues)
        is_protein = self.ligand.atoms['resname'].isin(_proteinResidues)
        atoms = self.ligand.atoms.copy()
        atoms['is_ligand'] = is_ligand
        atoms['is_protein'] = is_protein

        atoms1 = atoms[atoms['resSeq'].isin([241, 272])] # Corresponding to the residue 504 and 535 in Yang's paper
        atoms2 = atoms[atoms.is_ligand]

        # chain_match = {'A':'M', 'B':'H', 'C':'I', 'D':'J' ,'E':'K', 'F':'L'}

        protein_chain = self.chain_pair[0] # Must choose protein chain as keys then ligand chain as values
        ligand_chain = self.chain_pair[1]

        #print(matched_index)
        atoms_protein_thischain = atoms1[atoms1['chainID'] == protein_chain]
        atoms_ligand_thischain = atoms2[atoms2['chainID'] == ligand_chain]

        #print(atoms_ligand_thischain)

        protein_atoms_index = atoms_protein_thischain.index.values.tolist()
        protein_atoms_mass = []
        ligand_atoms_index = atoms_ligand_thischain.index.values.tolist()
        ligand_atoms_mass = []

        for index, each_atom in atoms_protein_thischain.iterrows():
            #print(type(each_atom))
            atom_name = each_atom.at['name']
            atom_mass = particle_definition.loc[particle_definition['name'] == atom_name, 'mass'].iloc[0]
            protein_atoms_mass.append(atom_mass)

        for index, each_atom in atoms_ligand_thischain.iterrows():
            atom_name = each_atom.at['name']
            atom_mass = particle_definition.loc[particle_definition['name'] == atom_name, 'mass'].iloc[0]
            ligand_atoms_mass.append(atom_mass)

        parameters = [5000, 1.0]
        self.force.addGroup(protein_atoms_index, protein_atoms_mass)
        self.force.addGroup(ligand_atoms_index, ligand_atoms_mass)
        self.force.addBond([0, 1], parameters)
