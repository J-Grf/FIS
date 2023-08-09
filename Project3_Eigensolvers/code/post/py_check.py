#!/usr/local/bin/python3.9
import numpy as np
from scipy.sparse import csc_matrix
import re

# read msr matrix
path = "../../input/test3_sym.txt"
JM = 0
VM = 0
i = 0
arraysize = 0
dim = 0
flag =''

with open(path) as f:
    data = f.readlines()
    flag = re.findall(r"([s|n])", data[0])[0]
    #print(flag)
    S = re.findall(r"([-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?\d+)?)", data[1])
    dim = int(S[0])
    arraysize = int(S[1])
    JM = np.zeros(arraysize)
    VM = np.zeros(arraysize)
    for line in data[2:]:
        Snumbers = re.findall(r"([-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?\d+)?)", line)
        JM[i] = int(Snumbers[0])
        VM[i] = float(Snumbers[1])
        i+=1

print(JM)
print(VM)

Matrix = np.zeros(shape = (dim, dim))
for i,v in enumerate(VM[:dim]):
    Matrix[i][i] = v
print(Matrix)

tmplen = 0
print(flag)
for i in range(dim):
    i1 = int(JM[i])
    i2 = int(JM[i+1])
    tmplen = i2 - i1
    for j in range(int(tmplen)):
        if abs(VM[i1-1]) > 1E-20:
            Matrix[i][int(JM[i1-1]-1)] = VM[i1-1]p
            if flag == 's':
                Matrix[int(JM[i1-1]-1)][i] = VM[i1-1]
            i1 += 1

print(Matrix)
MatrixCSC = csc_matrix(Matrix, dtype = float)
    