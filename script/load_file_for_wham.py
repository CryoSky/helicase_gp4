import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]

data1 = []
data2 = []



with open(input1, 'r') as fopen1:
    lines = fopen1.readlines()
    for line in lines:
        data1.append(float(line))

with open(input2, 'r') as fopen2:
    lines = fopen2.readlines()
    for line in lines:
        data2.append(float(line))

#new_data1 = data1[-250:]
#new_data2 = data2[-250:]

file_len = len(data1)
with open(output, 'w') as fwrite:
    for i in range(1, file_len+1):
        fwrite.writelines("%d %.2f %.2f\n" %(i, data1[i-1], data2[i-1]))
