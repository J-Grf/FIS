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

def removeFileIfExists(path):
    if os.path.isfile(path):
        os.remove(path)
    """ with open(path,'w') as f:
        f.write("n  time in seconds  nu1  nu2  gamma \n") """

removeFileIfExists(path)

for i in m:
    subprocess.run(["./main.exe", "LANC", inp, str(m)], check=True, cwd="/Users/johannes/Desktop/FIS/projects/Project3_Eigensolvers/code/main")

#read timings-file
times = []
with open(path, 'r') as f:
    data = f.readlines()
    for line in data:
        match = re.findall(r'(\d+[.]\d+)', line)
        if len(match) >= 1:
            print(times)
            times.append(float(match[0]))

fig = plt.figure("Lanczos_Timings")
fig.clear()
ax = fig.gca()
ax.semilogy(m, times)
#ax.set_xlim()
#ax.set_ylim()
ax.set_xlabel(r"t [s]", fontsize=fs)
ax.set_ylabel(r"|\lambda^{(k)} - \lambda^{(k-1)}|", fontsize=fs, labelpad=28, rotation = 0)

ax.grid()
plt.tight_layout()
plt.savefig(plotDir + "Timings_Lanczos.pdf", dpi = 100)
pgf.savePgf(pgfDir + "Timings_Lanczos.pgf", factor = 0.8)

    