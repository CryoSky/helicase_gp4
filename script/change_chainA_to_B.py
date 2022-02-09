import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
data = []

with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        if line[21] == 'M':
            #print(line[21])
            data.append(line[:21] + 'L' + line[22:])
        else:
            data.append(line)

with open(output_file, 'w') as fwrite:
    fwrite.writelines(data)
