#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import pgf
import os
import subprocess
import sys

path = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/rawData/timings_Lanczos.txt"
inp = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/input/cg_test_msr.txt"
plotDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/plots/"
pgfDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/report/plots/"
m = [30, 50, 75, 100]
eps = [1E0, 1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11, 1E-12, 1E-13, 1E-14]
fs = 15
regexFloat = re.compile(r'([-+]?[0-9]*[.][0-9]+(?:[eE][-+]?\d+)?)')

def removeFileIfExists(path):
    if os.path.isfile(path):
        os.remove(path)
    """ with open(path,'w') as f:
        f.write("n  time in seconds  nu1  nu2  gamma \n") """

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

def PlotLogLog(name, x, y, xlab, ylab):
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()

    ax.loglog(x, y)
    ax.set_xticks(eps)
    #ax.set_ylim(ymax = 100, ymin=1E-12)
    
    ax.set_xlabel(xlab, fontsize=fs)
    ax.set_ylabel(ylab, fontsize=fs, labelpad=40, rotation = 0)
    
    if "Errors" in name:
        ax.set_yticks([1E-10, 1E-8, 1E-6, 1E-4, 1E-2, 1E0, 1E2])
        ax.set_yticklabels([1E-10, 1E-8, 1E-6, 1E-4, 1E-2, 1E0, 1E2])
        ax.set_xticks([1E0, 1E-2, 1E-4, 1E-6, 1E-8, 1E-10, 1E-12, 1E-14])
        ax.set_xticklabels([1E0, 1E-2, 1E-4, 1E-6, 1E-8, 1E-10, 1E-12, 1E-14])
            
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + name + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + name + ".pgf", factor = 0.8)

if sys.argv[1] == "1":

    removeFileIfExists(path)

    for i in m:
        subprocess.run(["./main.exe", "LANC", inp, str(i), "-1"], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/code/main")

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

    PlotSemiLog("Lanczos_Timings", m, times, r"$m$", r"$t \ [s]$")
    PlotSemiLog("Lanczos_Lambdas", m, lambdas, r"$m$", r"$\lambda_{max}$")
    PlotSemiLog("Lanczos_Errors", m, errors, r"$m$", r"$|\lambda_{cg}^{ex} - \lambda^{(k)}|$")

if sys.argv[1] == "2":

    m = sys.argv[2]
    pathM = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/rawData/timings_Lanczos_" + m + ".txt"
    removeFileIfExists(pathM)
    for i in eps:
        subprocess.run(["./main.exe", "LANC", inp, m, str(i)], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/code/main")

    times = []
    eps = []
    lambdas = []
    errors = []
    with open(pathM, 'r') as f:
        data = f.readlines()
        for line in data:
            match = re.findall(regexFloat, line)
            if len(match) >= 1:
                eps.append(float(match[0]))
                times.append(float(match[1]))
                lambdas.append(float(match[2]))
                errors.append(float(match[3]))
    #print(errors)
    
    PlotLogLog("Lanczos_Timings" + str(m), eps, times, r"$\varepsilon$", r"$t \ [s]$")
    PlotLogLog("Lanczos_Lambdas" + str(m), eps, lambdas, r"$\varepsilon$", r"$\lambda_{max}$")
    PlotLogLog("Lanczos_Errors" + str(m), eps, errors, r"$\varepsilon$", r"$|\lambda_{cg}^{ex} - \lambda^{(k)}|$")
