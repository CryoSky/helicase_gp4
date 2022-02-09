import sys

input_file = sys.argv[1]
shift = sys.argv[2]
output_file = sys.argv[3]

data = []
with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        line = line.split()
        line[1] = float(line[1]) + float(shift)
        data.append([line[0], line[1], line[2]])
#print(data)

with open(output_file, 'w') as fwrite:
    for each_data in data:
        fwrite.writelines("%s %4.2f %s\n" %(each_data[0], each_data[1], each_data[2]))
