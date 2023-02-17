import matplotlib.pyplot as mpl
import numpy as np
from pathlib import Path
from scipy.io import loadmat
from lmfit.models import LinearModel
from scipy.stats import linregress

import pylibrary.plotting.plothelpers as PH
# sns.set_style(rc={"pdf.fonttype": 42})
mpl.style.use("~/.matplotlib/figures.mplstyle")

diskpath = Path('/Volumes/Pegasus_002/VCN-SBEM-Data/matfiles')

mode = "no_outlier"

def fitline(ax:object, x:np.array, y:np.array):
    # lm_A = LinearModel()
    # pars_A = lm_A.guess(largest_inputs, x=soma_SA)
    # out_A = lm_A.fit(largest_inputs, pars_A, x=soma_SA)
    # print(out_A.fit_report())
    lr_A = linregress(x, y)
    print(lr_A)
    ax.plot(x, lr_A.intercept + lr_A.slope*x, 'r', label='fitted line')
    ax.text(1.0, 0.01, f"{r2:s}={lr_A.rvalue**2:5.3f} p={lr_A.pvalue:5.3f}", 
            fontdict = {"fontsize": 8, "ha":"right", "va":'bottom'}, transform=ax.transAxes)


soma_sa_input = Path('somaSA_and_LargestInput.mat')
x = loadmat(Path(diskpath, soma_sa_input))
print(x.keys())
largest_inputs = x['Largest_Inputs'].squeeze()
large_input_count = x['LargeInput_Count'].squeeze()
soma_SA_bylarge = x['SomaSA_Covered_ByLarge'].squeeze()
soma_SA = x['soma_SAs'].squeeze()

if mode == 'no_outlier':
    no_outliers = np.where(soma_SA <= 1900.)
    print(no_outliers)
    largest_inputs = largest_inputs[no_outliers]
    large_input_count = large_input_count[no_outliers]
    soma_SA_bylarge = soma_SA_bylarge[no_outliers]
    soma_SA = soma_SA[no_outliers]

CD_x_pct = 100.*soma_SA_bylarge/soma_SA

P = PH.regular_grid(2, 2, figsize=(7, 5), labelposition=(-0.15, 1.05),
    order="rowsfirst",
    verticalspacing=0.15, horizontalspacing=0.15, 
    margins={"leftmargin": 0.1, "rightmargin": 0.05,
            "topmargin": 0.1, "bottommargin": 0.1,},
    panel_labels=["A", "B", "C", "D"],
            )
mum2 = r"($\mu m^2$)"
r2 = r"$r^2$"

P.axdict["A"].scatter(soma_SA, largest_inputs, s=9, c='k', marker='o', edgecolors = 'w', linewidths=0.5)
fitline(P.axdict["A"], x=soma_SA, y=largest_inputs)
P.axdict["A"].set_xlim(1000, 1600)
P.axdict["A"].set_ylim(50, 300)
P.axdict["A"].set_xlabel(f"Cell Body SA {mum2:s}")
P.axdict["A"].set_ylabel(f"Largest Input ASA {mum2:s}")

P.axdict["B"].scatter(CD_x_pct, largest_inputs, s=9, c='k', marker='o', edgecolors = 'w', linewidths=0.5)
fitline(P.axdict["B"], x=CD_x_pct, y=largest_inputs)
P.axdict["B"].set_xlim(20, 80)
P.axdict["B"].set_ylim(50, 300)
P.axdict["B"].set_xlabel(f"Percent of Soma Covered by Large Inputs")
P.axdict["B"].set_ylabel(f"Largest Input ASA {mum2:s}")

P.axdict["C"].scatter(soma_SA, large_input_count, s=9, c='k', marker='o', edgecolors = 'w', linewidths=0.5)
fitline(P.axdict["C"], x=soma_SA, y=large_input_count)
P.axdict["C"].set_xlim(1000, 1600)
P.axdict["C"].set_ylim(0, 15)
PH.set_axes_ticks(P.axdict["C"], yticks=[0, 5, 10, 15], yticks_str=["0", "5", "10", "15"], 
    y_minor=[1,2,3,4,6,7,8,9,11,12,13,14])
P.axdict["C"].set_xlabel(f"Cell Body SA {mum2:s}")
P.axdict["C"].set_ylabel(f"Number of Large Inputs")

P.axdict["D"].scatter(CD_x_pct, soma_SA, s=9, c='k', marker='o', edgecolors = 'w', linewidths=0.5)
fitline(P.axdict["D"], x=CD_x_pct, y=soma_SA)
P.axdict["D"].set_xlim(20, 80)
P.axdict["D"].set_ylim(1000, 1600)
P.axdict["D"].set_xlabel(f"Percent of Soma Covered by Large Inputs")
P.axdict["D"].set_ylabel(f"Cell Body SA {mum2:s}")

for ax in P.axdict:
    PH.nice_plot(P.axdict[ax], direction="outward", position=-0.03, ticklength=3.0)

mpl.savefig(f"/Volumes/Pegasus_002/VCN-SBEM-Data/SBEM-paper Figures/Figure2/Figure2_supp/Figure2_Supplemental1_{mode:s}.pdf")
mpl.savefig(f"/Volumes/Pegasus_002/VCN-SBEM-Data/SBEM-paper Figures/Figure2/Figure2_supp/Figure2_Supplemental1_{mode:s}.png")
mpl.show()