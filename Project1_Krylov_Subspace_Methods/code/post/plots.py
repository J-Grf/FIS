#!/usr/local/bin/python3.9

import matplotlib.pyplot as plt
import pgf
import re

path = '../../rawData/rawData2/' 
plotdir = '../../plots/plots2/'
pgfdir = '../../report/pgf/'
files = ['DotPKrylov.txt', 'relResiduals_NONE.txt', 'relResiduals_JACOBI.txt', 'relResiduals_GAUSSSEIDEL.txt', 'relResiduals_ILU.txt', 'eANorm.txt', 'rNorm.txt']
data = {}
print('hello')

for f in files:
    data[f] = []
    with open(path + f) as inp:
        content = inp.readlines()
        for line in content:
            data[f].append(float(line))
print(data)


#plot relative residuals against iteration index GMRES
maxIdx = max([len(data['relResiduals_NONE.txt']), len(data['relResiduals_JACOBI.txt']), len(data['relResiduals_GAUSSSEIDEL.txt']), len(data['relResiduals_ILU.txt'])])
fig = plt.figure('Residuals')
fig.clear()
ax = fig.gca()

for R in ['relResiduals_NONE.txt', 'relResiduals_JACOBI.txt', 'relResiduals_GAUSSSEIDEL.txt', 'relResiduals_ILU.txt']:
    print("len " + str(len(data[R])))
    idx = [i for i in range(len(data[R]))]
    l = R.split('_')[1].split('.txt')[0]
    l += ", " + r"$\tilde{m} = $ " + str(len(data[R]))
    print(l)
    ax.semilogy(idx, data[R], label = l)
    ax.set_xlabel(r'$k$')
    ax.set_ylabel(r'$\frac{||\mathbf{r}_k||_2}{||\mathbf{r}_0||_2}$', rotation = 0, labelpad = 18)
    fs = 12
    leg = ax.legend(loc='upper right', framealpha=1.0, fontsize=fs)
    leg.get_frame().set_edgecolor('k')

ax.set_xlim(xmin=0, xmax=500)
ax.set_ylim(ymin=10E-9)
ax.grid()
plt.savefig(plotdir + "relResiduals.pdf", dpi = 100)
pgf.savePgf(pgfdir + "relResiduals.pgf", factor = 0.8)

#plot othogonality of Krylov vectors
file = 'DotPKrylov.txt'
fig = plt.figure('DotP')
fig.clear()
ax = fig.gca()
ax.semilogy([i for i in range(len(data[file][1:]))], data[file][1:])
ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$(\mathbf{v}_1, \mathbf{v}_k)$', rotation = 0, labelpad = 18)
ax.grid()
ax.set_xlim(xmin=0, xmax=500)
plt.savefig(plotdir + "DotPKrylovVectors.pdf", dpi = 100)
pgf.savePgf(pgfdir + "DotPKrylovVectors.pgf", factor = 0.8)

#plot norms
fig = plt.figure('Norms')
fig.clear()
ax = fig.gca()
for R in ['rNorm.txt', 'eANorm.txt']:
    idx = [i for i in range(len(data[R]))]
    l = ""
    if R == "eANorm.txt":
        l = r'$||\mathbf{e}||_A$'
    else:
        l = r'$||\mathbf{r}||_2$'

    ax.semilogy(idx, data[R], label = l)
    ax.set_xlabel(r'$k$')
    ax.set_ylabel(r'$||\mathbf{e}||_A, ||\mathbf{r}||_2$', rotation = 0, labelpad = 18)
    fs = 12
    leg = ax.legend(loc='upper right', framealpha=1.0, fontsize=fs)
    leg.get_frame().set_edgecolor('k')

ax.grid()
ax.set_xlim(xmin=0, xmax=25000)
ax.set_ylim(ymax = 10E3, ymin=10E-4)
plt.savefig(plotdir + "CG_Norms.pdf", dpi = 100)
pgf.savePgf(pgfdir + "CG_Norms.pgf", factor = 0.8)

file = 'timings.txt'
fig = plt.figure('Timings')
fig.clear()
ax = fig.gca()
with open("/Users/johannes/Desktop/FIS/projects/Project1_Krylov_Subspace_Methods/code/main/" + file) as f:
    data = f.readline()
    timingsStr = re.findall(r"(\d+\.\d+)",data)

timings = []
for i in range(len(timingsStr)):
    timings.append(float(timingsStr[i]))

ax.plot([i for i in range(20, len(timingsStr) + 20)], timings)
ax.set_xlabel(r'$m$')
ax.set_ylabel(r'CPU-time [s]', rotation = 0, labelpad = 35)
ax.grid()
ax.set_xlim(xmin=0, xmax=480)
ax.set_xticks([0,20,100,200,300,400])
plt.axvline(20, color="k", linestyle="--", linewidth=0.5)
plt.savefig(plotdir + "Timings.pdf", dpi = 100)
pgf.savePgf(pgfdir + "Timings.pgf", factor = 0.8)