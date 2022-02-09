import simtk.openmm.app
import simtk.openmm
import sys
import pdbfixer
from simtk.openmm.app import PDBFile

def fixPDB(pdb_file):
    """Uses the pdbfixer library to fix a pdb file, replacing non standard residues, removing
    hetero-atoms and adding missing hydrogens. The input is a pdb file location,
    the output is a fixer object, which is a pdb in the openawsem format."""
    fixer = pdbfixer.PDBFixer(filename=pdb_file, )
    # fixer.findMissingResidues()
    # chains = list(fixer.topology.chains())
    # keys = fixer.missingResidues.keys()
    # for key in list(keys):
    #     chain_tmp = chains[key[0]]
    #     if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
    #         del fixer.missingResidues[key]

    # fixer.findNonstandardResidues()
    # fixer.replaceNonstandardResidues()
    # fixer.removeHeterogens(keepWater=False)
    # fixer.findMissingAtoms()
    # fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)
    return fixer
def main():
    input_file = sys.argv[1]
    fixer = fixPDB(input_file)
    PDBFile.writeFile(fixer.topology, fixer.positions, open('test.pdb', 'w'))

if __name__ == '__main__':
    main()