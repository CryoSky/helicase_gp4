import sys
import pandas
import numpy as np

input_pdb = sys.argv[1]
index1 = sys.argv[2]
index2 = sys.argv[3]

def parsePDB(pdb_file):
    """Transforms the pdb file into a pandas table for easy access and data editing."""

    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),  # record name
                    serial=int(line[6:11]),          # atom serial number
                    name=str(line[12:16]).strip(),   # atom name
                    altLoc=str(line[16:17]),         # alternate location indicator
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[21:22]),
                    resSeq=int(line[22:26]),         # residue sequence number
                    iCode=str(line[26:27]),          # code for insertion of residues
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=1.0 if line[54:60].strip() == '' else float(line[54:60]), # set to 1.0 because Plumed RMSD need 1.0
                    tempFactor=1.0 if line[60:66].strip() == '' else float(line[60:66]),
                    element=str(line[76:78]),        # element symbol, right-justified
                    charge=str(line[78:80]))         # charge on the atom, right-justified

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                lines += [pdb_line(line)]
    pdb_atoms = pandas.DataFrame(lines)
    pdb_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                           'resname', 'chainID', 'resSeq', 'iCode',
                           'x', 'y', 'z', 'occupancy', 'tempFactor',
                           'element', 'charge']]
    return pdb_atoms

pdb_lines = parsePDB(input_pdb)

index1 = int(index1)
index2 = int(index2)

atom1 = pdb_lines[pdb_lines['serial']==index1]
atom2 = pdb_lines[pdb_lines['serial']==index2]

x1 = atom1.x
y1 = atom1.y
z1 = atom1.z

x2 = atom2.x
y2 = atom2.y
z2 = atom2.z

xyz_1 = np.array((x1, y1, z1))
xyz_2 = np.array((x2, y2, z2))

dist = np.linalg.norm(xyz_2-xyz_1)
#print(dist)

r = dist/10
epsilon = 0.08333
sigma = 0.62
cutoff = 1.55
offset = -0.00135968
if r <= 1.55:
    energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)-offset
else:
    energy = 0
#print(energy)
print("%8d %8.3f %8.3f" %(index2, dist, energy))
