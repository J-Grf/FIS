#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import re
import os
import pgf
import sys

hfont = {'fontname':'Helvetica'}
rootdir = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/code/main/"
plotDir = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/plots/"
pgfDir  = "/Users/johannes/Desktop/FIS/projects/Project2_Multigrid/report/plots/"
fs = 15

nu = int(sys.argv[1])
regex = re.compile(r'residualRatio_(\d+)_' + str(nu) + r'_(\d+).txt')
print(str(regex))

fileMap = {}
n = nu1 = nu2 = 0 
for root, dirs, files in os.walk(rootdir):
  for file in files:
    match = regex.match(file)

    if match:
       fileMap[file] = []
       fileMap[file].append(match.group(1))
       fileMap[file].append(str(nu))
       fileMap[file].append(match.group(2))

print(fileMap)

data = {}
regexFloat = re.compile(r'([-+]?[0-9]*[.][0-9]+(?:[eE][-+]?\d+)?)')
for file in fileMap.keys():
    data[file] = []
    with open(rootdir + file) as inp:
        content = inp.readlines()
        data[file].append(1.0)
        for line in content:
            #print(line)
            match = regexFloat.findall(line)
            if(len(match) > 0):
                data[file].append(float(match[0]))
    #print(data)
print(data)
name = "ResidualPlot"
fig = plt.figure(name)
fig.clear()
ax = fig.gca()

legendmap = {}
for k, v in data.items():
    color = ""
    if fileMap[k][0] == str(4):
        color = "C0"
    else:
        color = "C1"
    idx = [i for i in range(len(v))]
    print(idx)
    print(v)
    legendmap[fileMap[k][0]], = ax.semilogy(idx, v, label = r"$n = $" + fileMap[k][0], color = color)

ax.set_ylim(ymax = 1, ymin=1E-11)
ticks = [1.0, 1E-2, 1E-4, 1E-6, 1E-8, 1E-10, 1E-12]
ax.yaxis.set_ticks(ticks)
#ax.yaxis.set_ticklabels(ticks)
ax.set_xlim(xmin=0, xmax=12)
ax.set_xlabel(r"$m$", fontsize=fs)
ax.set_ylabel(r"$\frac{||r^{(m)}||_{\infty}}{||r^{(0)}||_{\infty}}$", rotation = 0, fontsize = fs, labelpad = 28)

#handles, labels = ax.get_legend_handles_labels()
#ax.legend([legendmap['4'], legendmap['7']], [r"$n = 4$", r"$n = 7$"])

leg = ax.legend([legendmap['4'], legendmap['7']], [r"$n = 4$", r"$n = 7$"], loc='upper right', framealpha=1.0, fontsize=12)
leg.get_frame().set_edgecolor('k')

fileId = ""
if(nu == 1):
    fileId = "_nu1_1"
else:
    fileId = "_nu1_2"
#ax.set_title(r"Residuals for $n=4$ and $n=7$, $\nu_1 = \nu_2 = 1$")
ax.grid()
plt.tight_layout()
plt.savefig(plotDir + name + fileId + ".pdf", dpi = 100)
pgf.savePgf(pgfDir + name + fileId + ".pgf", factor = 0.8)
