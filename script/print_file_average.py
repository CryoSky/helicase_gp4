import sys
import numpy as np

input_file = sys.argv[1]
data = []

with open(input_file, 'r') as fopen:
    lines = fopen.readlines()
    for line in lines:
        data.append(float(line))

print(np.average(data))
