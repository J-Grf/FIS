#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import sys
import pgf

path = '../main/' 
plotdir = '../../plots/'
pgfdir = '../../report/plots/'
U_data = ['maxError.txt', 'infNorm.txt', 'u.txt','u_est.txt', 'u_ex.txt', 
          'u_ex_fine.txt', 'Restriction.txt', 'Prolongation.txt', 'grid_coarse.txt', 
          'grid_fine.txt', 'grid.txt']
data = {}
fs = 15
N = 4
n = sys.argv[2]

for f in U_data:
    data[f] = []
    try:
        with open(path + f) as inp:
            content = inp.readlines()
            for line in content:
                data[f].append(float(line))
    except FileNotFoundError:
        continue
#print(data)

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
    ax = fig.add_subplot(1, 3, 1, projection='3d', aspect='auto')
    ax.set_title(r"$\mathbf{U}_{ex}$")
    surf = ax.plot_surface(X, Y, U_ex, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
    
    ax = fig.add_subplot(1, 3, 2, projection='3d', aspect='auto')
    ax.set_title(r"$\mathbf{U}_{est}$")
    surf = ax.plot_surface(X, Y, U_est, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #fig.colorbar(surf, shrink=0.5, aspect=3, pad = 0.15)
    
    ax = fig.add_subplot(1, 3, 3, projection='3d', aspect='auto')
    ax.set_title("Error")
    surf = ax.plot_surface(X, Y, Error, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.update_layout(height=800, width=800)
    #fig.colorbar(surf, shrink=0.5, aspect=3, pad = 0.17)
    #fig.subplots_adjust(wspace=1, hspace=1)
    

if(int(sys.argv[1]) == 4):
    X, Y = np.meshgrid(data['grid.txt'], data['grid.txt'])

    size = len(data['grid.txt'])
    U_ex = np.zeros(shape = (size, size))
    U_est = np.zeros(shape = (size, size))
    Error = np.zeros(shape = (size, size))
    maxError = 0
    for i in range(size):
        for j in range(size):
            U_est[i][j] = data['u.txt'][i * size + j]
            U_ex[i][j]  = data['u_ex.txt'][i * size + j]
            Error[i][j] = abs(U_ex[i][j] - U_est[i][j])
            if(Error[i][j] > maxError):
                maxError = Error[i][j]
    print("maxError: " + str(maxError))

    s = ["u_ex_mg_" + n, "u_est_mg_" + n, "error_" + n]

    for i, obj in enumerate([U_ex, U_est, Error]):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        surf = ax.plot_surface(X, Y, obj, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        if i < 2: 
            ax.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        else:
            if n == '7':
                ax.set_zticks([0.0, 0.00005, 0.0001, 0.00015, 0.0002])
            else:
                ax.set_zticks([0.0, 0.003, 0.006, 0.009, 0.012])
        fig.tight_layout()
        plt.savefig(plotdir + s[i] + ".pdf", dpi=100)
        pgf.savePgf(pgfdir + s[i] + ".pgf", factor=0.8)

    