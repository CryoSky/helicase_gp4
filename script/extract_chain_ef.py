import sys
input_file = sys.argv[1]
output_file = sys.argv[2]


data = []
with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        try:
            if line[21]=='E' or line[21]=='F':
                data.append(line)
        except:
            pass

with open(output_file, 'w') as fwrite:
    fwrite.writelines(data)
