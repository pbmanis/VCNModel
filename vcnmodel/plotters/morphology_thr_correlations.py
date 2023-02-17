"""
Plot two panels showing how threshold changes with either AIS length
or dendrite surface area.
This replaces two plots previously made in Prism (Hillock-AIS-Threshold.pfzx)
aisdata and aisthrdata are taken from the data tables in that file

The threshold data as printed out by model_run2.py into thrrun.txt has been accumulated
by hand into the folder "VCN_ThresholdTest_Results" in the VCN-SBEM-Data folder.
The data in aisthrdata is a restructured vrsion of the 10.05.2021 and 08.26.2022 runs listed
in that file.

"""
from pylibrary.plotting import plothelpers as PH
import pandas as pd
import io
import matplotlib.pyplot as mpl
import numpy as np
import scipy as sp
import scipy.stats as SPS
import scipy.interpolate as SCINT

import vcnmodel.group_defs as GRPDEF


aisdata="""cell,DendAreas,OrigThr,OrigStd,MeshInfl,NewStandardized
2,4674.18,0.655,0.530,nan,0.617
5,3755.39,1.410,0.560,nan,0.583
6,4130.73,0.570,0.515,0.589,nan
9,3380.13,0.405,0.465,0.439,nan
10,4060.61,0.585,0.565,0.594,nan
11,3137.8,0.370,0.405,0.398,nan
13,3263.11,0.435,0.415,0.445,nan
17,3709.79,0.550,0.495,0.483,nan
18,3893.267,0.350,0.440,0.581,nan
30,3989.83,0.545,0.580,0.563,nan
"""


aisthrdata="""
AIS_Length,2,5,6,9,10,11,13,17,18,30
10,0.850,0.802,0.778,0.722,0.841,0.616,0.609,0.628,0.691,0.873
12,0.763,0.719,0.7,0.648,0.753,0.556,0.552,0.569,0.623,0.781
14,0.694,0.655,0.639,0.591,0.686,0.509,0.506,0.522,0.569,0.709
16,0.639,0.605,0.589,0.544,0.631,0.472,0.469,0.484,0.527,0.652
18,0.598,0.566,0.552,0.509,0.589,0.442,0.439,0.455,0.492,0.608
20,0.561,0.530,0.517,0.478,0.552,0.416,0.413,0.427,0.463,0.569
22,0.528,0.498,0.488,0.45,0.519,0.392,0.389,0.403,0.436,0.534
24,0.498,0.472,0.461,0.425,0.491,0.37,0.369,0.381,0.413,0.505
"""

somavsdendritearea= """cellname,soma_area,dend_area
BC02, 1446.05, 4674.18
BC05, 1343.55, 3755.39
BC06, 1464.27, 3755.39
BC09, 1340.36, 3380.13
BC10, 1388.7, 4060.61
BC11, 1288.45, 3137.8
BC13, 1305.31, 3263.11
BC17, 1357.62, 3709.79
BC18, 1292.47, 3893.267
BC30, 1508.36, 3989.83
"""

aislengths = {6: 14.16,
              9: 24.16,
              10: 18.66,
              11: 22.19,
              13: 17.58,
              17: 15.07,
              18: 12.58,
              30: 21.37,
            }

def prepare_data(datas):
    sio = io.StringIO(datas)
    df = pd.read_table(sio, sep=",")
    return df

def plot_scatter(ax, df, cells, symbol, data=["DendAreas", "MeshInfl"]):
    colorl = []
    dfs = df[df["cell"].isin(cells)]
    for i, cell in enumerate(cells):
        bc_index = GRPDEF.gradeACells.index(int(cell))
        colorl.append(GRPDEF.sns_colors[bc_index])
    # note the table uses MeshInfl as the Threshold value (referncing the cell area correction)
    # and NewStandardized for cells 2 and 5 for the axon.
    ax.scatter(dfs[data[0]], dfs[data[1]], s=20, marker=symbol, c=colorl)

class SD():
    def __init__(self):
        self.cells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
    
    def SomavsDend(self, ax=None):
        df = prepare_data(somavsdendritearea)
        
        def make_celln(row):
            row.cell = int(row.cellname[2:])
            row.color = GRPDEF.sns_colors[self.cells.index(row.cell)]
            return row
        df["cell"] = ''
        df["color"] = 'k'
        df = df.apply(make_celln, axis=1)
        print(df.head())
        colorl = []
        dfs = df[df["cell"].isin(self.cells)]
        for i, cell in enumerate(self.cells):
            bc_index = GRPDEF.gradeACells.index(int(cell))
            colorl.append(GRPDEF.sns_colors[bc_index])
        # note the table uses MeshInfl as the Threshold value (referncing the cell area correction)
        # and NewStandardized for cells 2 and 5 for the axon.
        symbol = 'o'
        if ax is None:
            f, ax = mpl.subplots(1,1)
        cells0 = GRPDEF.get_cells_in_group(1)
        plot_scatter(ax, df, cells0, symbol=GRPDEF.group_symbols[0], data=['soma_area', 'dend_area'])
        cells1 = GRPDEF.get_cells_in_group(2)
        plot_scatter(ax, df, cells1, symbol=GRPDEF.group_symbols[1], data=['soma_area', 'dend_area'])

        # ax.scatter(dfs["soma_area"], dfs["dend_area"], s=20, marker=symbol, c=colorl)
        res = SPS.linregress(dfs["soma_area"], dfs["dend_area"])
        xl = np.linspace(np.min(dfs["soma_area"]), np.max(dfs["soma_area"]))
        yl = res.slope * xl + res.intercept
        ax.plot(xl, yl, 'k', linewidth=0.5)
        ax.set_ylim(3000, 5000)
        ax.set_ylabel(r"Dendrite Area (${\mu m^2}$)")
        ax.set_xlim(1250, 1550)
        ax.set_xlabel(r"Soma Area (${\mu m^2}$)")
        r2 = r"${r^2}$"
        ax.text(1550, 3100, f"p = {res.pvalue:.3f}, {r2:s} = {res.rvalue**2:.3f} ", fontsize=8, ha="right")
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)


class AIS():
    def __init__(self):
        pass
    
    def DendvsThr(self, ax=None):
        df = prepare_data(aisdata)

        show = False
        if ax is None:
            f, ax = mpl.subplots(1,1)
            show = True
        colorl = []
        # map colors
        cells0 = GRPDEF.get_cells_in_group(1)
        plot_scatter(ax, df, cells0, symbol=GRPDEF.group_symbols[0])
        cells1 = GRPDEF.get_cells_in_group(2)
        plot_scatter(ax, df, cells1, symbol=GRPDEF.group_symbols[1])
        # note the table uses MeshInfl as the Threshold value (referncing the cell area correction)
        # and NewStandardized for cells 2 and 5 for the axon.
#        ax.scatter(df["DendAreas"], df["MeshInfl"], s=12, marker=marks, c=colorl)
        # 2 and 5 are plotted separately here
        dfs = df[df["MeshInfl"].isnull()]
        colorl = []
        for cell in [2, 5]:
            mark = GRPDEF.get_group_symbol(cell)
            bc_index = GRPDEF.gradeACells.index(int(cell))
            colorl.append(GRPDEF.sns_colors[bc_index])
            ax.scatter(dfs["DendAreas"], dfs["NewStandardized"], 
                s=20, marker=mark, edgecolors=colorl, facecolors=["w", "w"], linewidths=1)
        dfn = df[df["MeshInfl"].notnull()]
        res = sp.stats.linregress(dfn["DendAreas"], dfn["MeshInfl"])
        xp = np.linspace(3000, 4500, 100)
        # print(res.slope)
        # print(res.intercept)
        # print(f"R-squared: {res.rvalue**2:.6f}")
        # print(f"P value: {res.pvalue:.5f}")
        ax.plot(xp, res.intercept + res.slope*xp, 'k-', linewidth=0.5)
        ax.text(x=1.0, y=0.05, s=f"p = {res.pvalue:.5f}, r$^2$ = {res.rvalue**2:5.3f}",
            horizontalalignment="right", fontsize=8,
            transform=ax.transAxes)
        ax.set_xlim(3000, 5000)
        ax.set_ylim(0.35, 0.65)
        ax.set_xlabel("Dendrite Area ($\mu m^2$)")
        ax.set_ylabel("Threshold (nA)")
        xtpos = np.arange(3000, 5001, 500)
        ytpos = np.arange(0.35, 0.66, 0.05)
        PH.set_axes_ticks(ax=ax,
            xticks = xtpos,
            xticks_str = [f"{int(x):d}" for x in xtpos],
           # xticks_pad:Union[List, None]=None,
            #x_minor = np.arange(0, 25, 2),
            major_length = 3.0,
            minor_length =1.5,
            yticks = ytpos, 
            yticks_str=[f"{x:4.2f}" for x in ytpos], 
            yticks_pad=[1]*7,
           # y_minor=[0.3, 0.5, 0.7, 0.9],
            fontsize=8,
        )
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)
        if show:
            mpl.show()

    def AISLengthThr(self, ax=None, show=False):
        df = prepare_data(aisthrdata)
        df_true = prepare_data(aisdata) # for the true length

        if ax is None:
            f, ax = mpl.subplots(1,1)
            show = True
        ax.set_clip_on(False)
        colnames = df.columns
        # print("colnames: ", colnames)
        # print(df["AIS_Length"])
        syms = ['o', 'o', '^', 's', 'v', 'd', 'o', 'o', '^', 'v']
        facecolors = GRPDEF.sns_colors
        # edgecolors = [None]*len(colors)

        for i, coln in enumerate(colnames):
            if i == 0:
                continue
            celln = int(coln)
            bc_index = GRPDEF.gradeACells.index(celln)
            color = GRPDEF.sns_colors[bc_index]
            sym = syms[bc_index]
            if celln in [13, 18, 30]:
                mk_face = "none"
                mk_edge = facecolors[bc_index]
            else:
                mk_face = facecolors[bc_index]
                mk_edge = mk_face
            
            print(coln, bc_index, color, sym, mk_face, mk_edge)
            bc_num = f"BC{int(coln):02d}"
            ax.plot(df["AIS_Length"].values, df[coln].values, 
            # remove symbols for revision - only plot symbols for the actual AIS length on the same data
                # sym, markerfacecolor=mk_face, markeredgecolor = mk_edge, markeredgewidth=1,  
                linestyle='-', color=color, linewidth=1,
                label=bc_num)
            if celln not in [2, 5]:
                la_y = SCINT.lagrange(df["AIS_Length"].values, df[coln].values)
                ax.plot(aislengths[celln], la_y(aislengths[celln]), 'o', color=GRPDEF.sns_colors[bc_index], markersize=4)

        ax.set_xlim(0, 25)
        ax.set_xlabel(f"AIS Length ($\mu$m)")
        ax.set_ylabel(f"Threshold (nA)")
        ax.set_ylim(0.2, 1.0)
        xtpos = np.arange(0, 26, 8)
        ytpos = [0.2, 0.4, 0.6, 0.8, 1.0]
        PH.set_axes_ticks(ax=ax,
            xticks = xtpos,
            xticks_str = [f"{int(x):d}" for x in xtpos],
           # xticks_pad:Union[List, None]=None,
            x_minor = np.arange(0, 25, 2),
            major_length = 3.0,
            minor_length =1.5,
            yticks = ytpos, 
            yticks_str=[f"{x:3.1f}" for x in ytpos], 
            yticks_pad=[1]*7,
            y_minor=[0.3, 0.5, 0.7, 0.9],
            fontsize=8,
        )
        # ax.legend(fontsize=6, loc="upper left", ncol=1, frameon=False)
        # P.axdict[pan[5]].legend(handles=custom_legend, handlelength=1, loc="upper right", fontsize=7, labelspacing=0.33, markerscale=0.5)
        ax.legend(loc="upper left", fontsize=6, ncol=1, frameon=False, labelspacing=0.33, markerscale=0.5)
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)
        if show:
            mpl.show()


if __name__ == "__main__":
    # A = AIS()
    # df = prepare_data(aisdata)
    # print(df.head())
    # # A.DendvsThr()
    # A.AISLengthThr(show=True)
    S = SD()
    S.SomavsDend()
    mpl.show()