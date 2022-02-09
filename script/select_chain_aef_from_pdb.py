import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

new_data = []

with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        if line[0:5] == 'MODEL' or line[0:6] == 'ENDMDL':
            new_data.append(line)
        elif line[0:6] == 'HETATM':
            if line[21] in ['A', 'E', 'F']:
                new_data.append(line)

with open(output_file, 'w') as fwrite:
    fwrite.writelines(new_data) 
            
