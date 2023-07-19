#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import os
import pgf
import sys
import subprocess

timingsDir = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/rawData/"
plotDir = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/plots/"
pgfDir  = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/report/plots/"
maxNu = 50
maxGamma = 10
maxIt = 0
timingsFile = ""
Type = ""
n = 0
fs = 15

def removeFileIfExists(path):
    if os.path.isfile(path):
        os.remove(path)
    with open(path,'w') as f:
        f.write("n  time in seconds  nu1  nu2  gamma \n")

def runMG(maxIt, Type, n):
    if Type == "NU":
        for i in range(1, maxIt + 1):
            print("h")
            subprocess.run(["./main.exe", "mg", "n", str(n), str(i), "1", "2", Type], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/code/main")
    elif Type == "GAMMA":
        for i in range(1, maxIt + 1):
            subprocess.run(["./main.exe", "mg", "n", str(n), "1", "1", str(i), Type], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/code/main")

def createPlots(Type, timingsFile, n, maxIt):
    #read timings-file
    times = []
    with open(timingsFile, 'r') as f:
        data = f.readlines()
        for line in data:
            match = re.findall(r'(\d+[.]\d+)', line)
            if len(match) >= 1:
                #print(times)
                times.append(float(match[0]))
    
    #plot
    name = "Timings"
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()
    ax.plot([i for i in range(maxIt)], times, label = r"$n = $" + str(n))

    ax.set_ylim(ymin=0, ymax=1.1 * max(times))
    ax.set_xlim(xmin=0, xmax=maxIt+1)
    xlab = ""
    if Type == "NU":
        xlab = r"$\nu_1$"
    elif Type == "GAMMA":
        xlab = r"$\gamma$"

    ax.set_xlabel(xlab, fontsize=fs)
    ax.set_ylabel("runtime [s]", rotation = 0, fontsize = fs, labelpad = 35)
    
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + "Timings_" + Type + "_" + str(n) + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + "Timings_" + Type + "_" + str(n) + ".pgf", factor = 0.8)

n = int(sys.argv[2])
if sys.argv[1] == "NU":
    timingsFile = timingsDir + "timings_nu_" + str(n) + ".txt"
    Type = "NU"
    maxIt = maxNu
elif sys.argv[1] == "GAMMA":
    timingsFile = timingsDir + "timings_gamma_" + str(n) + ".txt"
    Type = "GAMMA"
    maxIt = maxGamma
else:
    print("not known!")

removeFileIfExists(timingsFile)
runMG(maxIt, Type, n)
createPlots(Type, timingsFile, n, maxIt)
