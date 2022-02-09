import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
data = []


with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    header = lines[0:5]
    data.extend(header)
    for line in lines[5:]:
        check = line.split()
        try:
            if check[8] == 'A,' and check[16] == 'B,':
                #print(line)
                data.append(line)
        except:
            data.append(line)
print(data)
with open(output_file, 'w') as fwrite:
    fwrite.writelines(data)
