import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

total_energy = []
ligand_energy = []
dna_bond_energy = []
with open(input_file, 'r') as fopen:
    lines = fopen.readlines()[78:]
    for line in lines:
        if line[0] == '-':
            total_energy.append(float(line))
        elif line[:24] == 'LigandDistanceRestraint1':
            line = line.split()
            ligand_energy.append(float(line[-2]))
        elif line[:4] == 'Bond':
            line = line.split()
            dna_bond_energy.append(float(line[-2]))
print(len(total_energy))
print(len(ligand_energy))
print(len(dna_bond_energy))

difference = []
zip_object = zip(total_energy, ligand_energy, dna_bond_energy)
for list1_i, list2_i, list3_i in zip_object:
    difference.append(list1_i)
    #difference.append(list1_i-list2_i)



with open(output_file, 'w') as fwrite:
    fwrite.writelines('%f\n' % l for l in difference)
