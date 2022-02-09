import sys
input = sys.argv[1]
output = sys.argv[2]


out_lines = []

with open(input, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        #print(line[13:15])
        if line[13:15] == 'CB':
            new_line = line[0:77]+'B'+line[78:-1]+'\n'
            #print(new_line)
        else:
            new_line = line
        out_lines.append(new_line)

with open(output, 'w') as fwrite:
    fwrite.writelines(out_lines)
