#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import pgf
import sys

plotDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/plots/"
pgfDir = "/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/report/plots/"

regexFloat = re.compile(r'([-+]?[0-9]*[.]?[0-9]+(?:[eE][-+]?\d+)?)')
fs = 12

def PlotSemiLog(name, x, y, xlab, x_min, ylab, x_max):
    fig = plt.figure(name)
    fig.clear()
    ax = fig.gca()

    ax.semilogy(x, y)
    ax.set_ylim(ymax = 1E8, ymin=1E-10)
    #ax.set_ylim(ymax = 1E4, ymin = 1E-8)
    ax.set_xlim(xmin = x_min, xmax = x_max)
    ax.set_xlabel(xlab, fontsize=fs)
    ax.set_ylabel(ylab, fontsize=fs, labelpad=40, rotation = 0)
    ax.grid()
    plt.tight_layout()
    plt.savefig(plotDir + name + ".pdf", dpi = 100)
    pgf.savePgf(pgfDir + name + ".pgf", factor = 0.8)


it = []
time = []
diff = []
error = []
i = 0
with open(sys.argv[1]) as f:
    content = f.readlines()
    for line in content:
        print(line)
        match = regexFloat.findall(line)
        if(len(match) > 0):
            i+=1
            it.append(i)
            time.append(float(match[1]))
            diff.append(float(match[2]))
            error.append(float(match[3]))

ylab1 = r"$|\lambda^{(k)} - \lambda^{(k-1)}|$"
print(diff)
PlotSemiLog("PowerIt_test_convergence", it, diff, r"$k$", 1, ylab1, 800)
#PlotSemiLog("PowerIt_2_diff", time, diff, r"$t \ [s]$", 0, ylab1, 3.5)

#ylab2 = r"$|\lambda_{cg}^{ex} - \lambda^{(k)}|$"
#PlotSemiLog("PowerIt_1_error", it, error, r"$k$", 1, ylab2, 2200)
#PlotSemiLog("PowerIt_2_error", time, error, r"$t \ [s]$", 0, ylab2, 3.5)


