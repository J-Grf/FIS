#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import pgf

plotDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/plots/"
pgfDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/report/plots/"

regexFloat = re.compile(r'([-+]?[0-9]*[.][0-9]+(?:[eE][-+]?\d+)?)')
fs = 15

def PlotSemiLog(name, x, y):
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()

    ax.semilogy(x, y)
    ax.set_ylim(ymax = 1, ymin=1E-11)
    ax.set_xlim(xmin = 1)
    ax.set_xlabel(r"$k$", fontsize=fs)
    ax.set_ylabel(r"|\lambda^{(k)} - \lambda^{(k-1)}|", fontsize=fs, labelpad=28, rotation = 0)
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + name + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + name + ".pgf", factor = 0.8)


it = []
time = []
error = []
with open("../../rawData/powerIterationConvergence.txt") as f:
    content = f.readlines()
    for line in content:
        #print(line)
        match = regexFloat.findall(line)
        if(len(match) > 0):
            it.append(float(match[0]))
            time.append(float(match[1]))
            error.append(float(match[2]))


PlotSemiLog("PowerIt_1", it, error)
PlotSemiLog("PowerIt_2", time, error)


