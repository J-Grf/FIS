#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import sys

path = '../main/' 
plotdir = '../../plots/'
pgfdir = '../../report/pgf/'
U_data = ['maxError.txt', 'infNorm.txt', 'u.txt','u_est.txt', 'u_ex.txt', 'u_ex_fine.txt', 'Restriction.txt', 'Prolongation.txt', 'grid_coarse.txt', 'grid_fine.txt', 'grid.txt']
data = {}
fs = 15
N = 4

for f in U_data:
    data[f] = []
    try:
        with open(path + f) as inp:
            content = inp.readlines()
            for line in content:
                data[f].append(float(line))
    except FileNotFoundError:
        continue
print(data)

if(int(sys.argv[1]) == 1):
    for f in ['maxError.txt', 'infNorm.txt']:
        idx = [i for i in range(len(data[f]))]
        name = f.split('.txt')[0]
        fig = plt.figure(name)
        fig.clear()
        ax = fig.gca()
        ax.semilogy(idx, data[f])
        plt.savefig(plotdir + name + ".pdf", dpi = 100)
    
    X, Y = np.meshgrid(data['grid.txt'], data['grid.txt'])
    size = len(data['grid.txt'])
    U_ex = np.zeros(shape = (size, size))
    U_est = U_ex.copy()
    Error = U_ex.copy()
    for i in range(size):
        for j in range(size):
            U_est[i][j] = data['u_est.txt'][i * size + j]
            U_ex[i][j] = np.sin(2 * X[i][j] * np.pi ) * np.sin(2 * Y[i][j] * np.pi)
            Error[i][j] = abs(U_ex[i][j] - U_est[i][j])
    
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.set_title("U_ex")
    surf = ax.plot_surface(X, Y, U_ex, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)
    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.set_title("U_est")
    surf = ax.plot_surface(X, Y, U_est, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)
    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.set_title("Error")
    surf = ax.plot_surface(X, Y, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)

    plt.savefig(plotdir + "u_sg.pdf", dpi=100)

if(int(sys.argv[1]) == 2):
    X, Y = np.meshgrid(data['grid.txt'], data['grid.txt'])
    X_c, Y_c = np.meshgrid(data['grid_coarse.txt'], data['grid_coarse.txt'])
    X_f, Y_f = np.meshgrid(data['grid_fine.txt'], data['grid_fine.txt'])

    """ XEX = np.arange(0, 1, 1/1000)
    YEX = np.arange(0, 1, 1/1000)
    XEX, YEX = np.meshgrid(XEX, YEX) """

    #Check Restriction
    size = len(data['grid_coarse.txt'])
    size_f = 2 * size - 1

    print("Restriction: size: " + str(size) + " size_f: " + str(size_f))
    U_ex = np.zeros(shape = (size_f, size_f))
    U_res = np.zeros(shape = (size, size))
    Error = np.zeros(shape = (size, size))

    for i in range(size_f):
        for j in range(size_f):
            U_ex[i][j]  = data['u_ex.txt'][i * size_f + j]

    for i in range(size):
        for j in range(size):
            U_res[i][j] = data['Restriction.txt'][i * size + j]
            Error[i][j] = U_ex[2*i][2*j] - U_res[i][j]

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.set_title("U_ex")
    surf = ax.plot_surface(X, Y, U_ex, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)

    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.set_title("U_res")
    surf2 = ax.plot_surface(X_c, Y_c, U_res, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf2, shrink=0.5, aspect=5, pad = 0.15)

    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.set_title("Error")
    surf2 = ax.plot_surface(X_c, Y_c, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf2, shrink=0.5, aspect=5, pad = 0.15)

    plt.savefig(plotdir + "restriction.pdf", dpi=100)

    #Check Prolongation
    size = len(data['grid.txt'])
    size_f = 2 * size - 1
    print("Prolongation: size: " + str(size) + " size_f: " + str(size_f))
    U_ex = np.zeros(shape = (size_f, size_f))
    U_res = np.zeros(shape = (size_f, size_f))
    Error = np.zeros(shape = (size_f, size_f))

    for i in range(size_f):
        for j in range(size_f):
            U_ex[i][j]  = data['u_ex_fine.txt'][i * size_f + j]
            U_res[i][j] = data['Prolongation.txt'][i * size_f + j]
            Error[i][j] = U_ex[i][j] - U_res[i][j]

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.set_title("U_ex")
    surf = ax.plot_surface(X_f, Y_f, U_ex, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)

    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.set_title("U_res")
    surf2 = ax.plot_surface(X_f, Y_f, U_res, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf2, shrink=0.5, aspect=5, pad = 0.15)

    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.set_title("Error")
    surf2 = ax.plot_surface(X_f, Y_f, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf2, shrink=0.5, aspect=5, pad = 0.15)

    plt.savefig(plotdir + "Prolongation.pdf", dpi=100)

if(int(sys.argv[1]) == 3):
    X, Y = np.meshgrid(data['grid.txt'], data['grid.txt'])

    size = len(data['grid.txt'])
    U_ex = np.zeros(shape = (size, size))
    U_est = U_ex.copy()
    Error = U_ex.copy()
    for i in range(size):
        for j in range(size):
            U_est[i][j] = data['u.txt'][i * size + j]
            U_ex[i][j]  = data['u_ex.txt'][i * size + j]
            Error[i][j] = abs(U_ex[i][j] - U_est[i][j])
    
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.set_title("U_ex")
    surf = ax.plot_surface(X, Y, U_ex, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)
    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.set_title("U_est")
    surf = ax.plot_surface(X, Y, U_est, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)
    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.set_title("Error")
    surf = ax.plot_surface(X, Y, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, pad = 0.15)

    plt.savefig(plotdir + "u_mg.pdf", dpi=100)
