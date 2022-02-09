import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
data = []
with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    # for line in lines:
    #     data.append(line)


for each_line in lines:
    if each_line[21] == 'L':
        # if each_line[12:16].strip() in ['PA', 'PB']:
        #     data.append(each_line[0:6] + "{:5d}".format(int(each_line[6:11])) + each_line[11:17] + 'ADP' + each_line[20:])
        # elif each_line[12:16].strip() == 'PG':
        #     pass
        # elif each_line[12:16].strip() in ['A', 'S']:
        #     data.append(each_line[0:6] + "{:5d}".format(int(each_line[6:11]) - 1) + each_line[11:17] + 'ADP' + each_line[20:])
        # else:
        #     data.append(each_line[0:6] + "{:5d}".format(int(each_line[6:11]) - 1) + each_line[11:])
        pass
    elif each_line[21] == 'M':
        data.append(each_line[0:6] + "{:5d}".format(int(each_line[6:11]) - 1) + each_line[11:])
    else:
        data.append(each_line)

with open(output_file, 'w') as fwrite:
    fwrite.writelines(data)
