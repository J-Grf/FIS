import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        # Adjust to your LaTex-Engine
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "axes.unicode_minus": False,
        "figure.autolayout": True,
    }
)

TW = 5.90666 #textwidth
TH = 8.20633 #textheight

def savePgf(path, factor=0.55):
    plt.gcf().set_size_inches(TH * factor, TW * factor)
    # Fixes cropped labels
    plt.tight_layout()
    # Save as pgf
    plt.savefig(path) #bbox_inches='tight')

#To optionally set font
""" font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 10}
matplotlib.rc('font', **font) """