#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import pgf
import os
import subprocess

path = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/rawData/timings_Lanczos.txt"
inp = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/input/cg_test_msr.txt"
plotDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/plots/"
pgfDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/report/plots/"
m = [30, 50, 75, 100]
fs = 15
regexFloat = re.compile(r'([-+]?[0-9]*[.][0-9]+(?:[eE][-+]?\d+)?)')

def removeFileIfExists(path):
    if os.path.isfile(path):
        os.remove(path)
    """ with open(path,'w') as f:
        f.write("n  time in seconds  nu1  nu2  gamma \n") """

removeFileIfExists(path)

for i in m:
    subprocess.run(["./main.exe", "LANC", inp, str(i)], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/code/main")

#read timings-file
times = []
lambdas = []
errors = []
with open(path, 'r') as f:
    data = f.readlines()
    for line in data:
        match = re.findall(regexFloat, line)
        if len(match) >= 1:
            print(times)
            times.append(float(match[0]))
            lambdas.append(float(match[1]))
            errors.append(float(match[2]))

#print(lambdas)

def PlotSemiLog(name, x, y, xlab, ylab, x_min = None):
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()

    ax.semilogy(x, y)
    #ax.set_ylim(ymax = 100, ymin=1E-12)
    if x_min != None:
        ax.set_xlim(xmin = x_min)
    
    ax.set_xlabel(xlab, fontsize=fs)
    ax.set_ylabel(ylab, fontsize=fs, labelpad=40, rotation = 0)
    
    if "Lambdas" in name:
        plt.axhline(y = 9.5986080894852857E3, color = 'r', linestyle='-')
        ax.set_yticks([9.584E3, 9.586E3, 9.588E3, 9.59E3, 9.592E3, 9.594E3, 9.596E3, 9.598E3])
    elif "Timings" in name:
        ax.set_yticks([4E-2, 6E-2, 8E-2, 1E-1, 2E-1])
        ax.set_yticklabels([4E-2, 6E-2, 8E-2, 1E-1, 2E-1])
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + name + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + name + ".pgf", factor = 0.8)

PlotSemiLog("Lanczos_Timings", m, times, r"$m$", r"$t \ [s]$")
PlotSemiLog("Lanczos_Lambdas", m, lambdas, r"$m$", r"$\lambda_{max}$")
PlotSemiLog("Lanczos_Errors", m, errors, r"$m$", r"$|\lambda_{cg}^{ex} - \lambda^{(k)}|$")

    