#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import pgf

plotDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/plots/"
pgfDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/report/plots/"

regexFloat = re.compile(r'([-+]?[0-9]*[.][0-9]+(?:[eE][-+]?\d+)?)')
fs = 12

def PlotSemiLog(name, x, y, xlab, x_min, ylab):
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()

    ax.semilogy(x, y)
    ax.set_ylim(ymax = 100, ymin=1E-12)
    ax.set_xlim(xmin = x_min)
    ax.set_xlabel(xlab, fontsize=fs)
    ax.set_ylabel(ylab, fontsize=fs, labelpad=35, rotation = 0)
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + name + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + name + ".pgf", factor = 0.8)


it = []
time = []
diff = []
error = []
i = 0
with open("../../rawData/powerIterationConvergence.txt") as f:
    content = f.readlines()
    for line in content:
        print(line)
        match = regexFloat.findall(line)
        if(len(match) > 0):
            i+=1
            it.append(i)
            time.append(float(match[0]))
            diff.append(float(match[1]))
            error.append(float(match[2]))

ylab1 = r"$|\lambda^{(k)} - \lambda^{(k-1)}|$"
PlotSemiLog("PowerIt_1_diff", it, diff, r"$k$", 1, ylab1)
PlotSemiLog("PowerIt_2_diff", time, diff, r"$t \ [s]$", 0, ylab1)

ylab2 = r"$|\lambda_{cg}^{ex} - \lambda^{(k)}|$"
PlotSemiLog("PowerIt_1_error", it, error, r"$k$", 1, ylab2)
PlotSemiLog("PowerIt_2_error", time, error, r"$t \ [s]$", 0, ylab2)


