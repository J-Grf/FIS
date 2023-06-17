#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import math

path = '../main/' 
plotdir = '../../plots/'
pgfdir = '../../report/pgf/'
files = ['maxError.txt', 'grid.txt', 'infNorm.txt', 'u_est.txt']
data = {}
print('hello')
fs = 15

for f in files:
    data[f] = []
    with open(path + f) as inp:
        content = inp.readlines()
        for line in content:
            data[f].append(float(line))
print(data)

for f in ['maxError.txt', 'infNorm.txt']:
    idx = [i for i in range(len(data[f]))]
    name = f.split('.txt')[0]
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()
    ax.semilogy(idx, data[f])
    plt.savefig(plotdir + name + ".pdf", dpi = 100)

f = 'u_est.txt'
X, Y = np.meshgrid(data['grid.txt'], data['grid.txt'])

""" XEX = np.arange(0, 1, 1/1000)
YEX = np.arange(0, 1, 1/1000)
XEX, YEX = np.meshgrid(XEX, YEX) """
size = len(data['grid.txt'])
U = np.zeros(shape = (size, size))
Error = U.copy()
for i in range(size):
    for j in range(size):
        U[i][j] = data[f][i * size + j]
        Error[i][j] = abs(np.sin(2 * X[i][j] * np.pi ) * np.sin(2 * Y[i][j] * np.pi) - U[i][j])

#U_ex = np.sin(2 * XEX * np.pi ) * np.sin(2 * YEX * np.pi) 

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.set_title("U")
surf = ax.plot_surface(X, Y, U, cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)

ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.set_title("Error")
surf2 = ax.plot_surface(X, Y, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(surf2, shrink=0.5, aspect=5, pad = 0.15)

plt.savefig(plotdir + "u.pdf", dpi=100)