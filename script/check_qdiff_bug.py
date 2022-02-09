import sys
input_1 = sys.argv[1]
input_2 = sys.argv[2]
counts = 1
list_1 = []
list_2 = []
with open (input_1, 'r') as fread1:
    for line in fread1.readlines():
        list_1.append(line)

with open (input_2, 'r') as fread2:
    for line in fread2.readlines():
        list_2.append(line)

print(len(list_1))
print(len(list_2))

print("Now start debug!")

for i in range(len(list_1)):
    test = float(list_1[i]) - float(list_2[i])
    if test < 0:
        print(counts)
        #print("%f %d" %(test, counts))
    counts += 1