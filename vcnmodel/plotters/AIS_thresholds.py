"""
Plot two panels showing how threshold changes with either AIS length
or dendrite surface area.
This replaces two plots previously made in Prism (Hillock-AIS-Threshold.pfzx)
aisdata and aisthrdata are taken from the data tables in that file
"""
from pylibrary.plotting import plothelpers as PH
import pandas as pd
import io
import matplotlib.pyplot as mpl
import numpy as np
import scipy as sp


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
AIS_Length,6,9,10,11,13,17,18,30
10,0.778,0.722,0.841,0.616,0.609,0.628,0.691,0.873
12,0.7,0.648,0.753,0.556,0.552,0.569,0.623,0.781
14,0.639,0.591,0.686,0.509,0.506,0.522,0.569,0.709
16,0.589,0.544,0.631,0.472,0.469,0.484,0.527,0.652
18,0.552,0.509,0.589,0.442,0.439,0.455,0.492,0.608
20,0.517,0.478,0.552,0.416,0.413,0.427,0.463,0.569
22,0.488,0.45,0.519,0.392,0.389,0.403,0.436,0.534
24,0.461,0.425,0.491,0.37,0.369,0.381,0.413,0.505
"""

# seaborn default palette, first 10 colors
colors = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
]
gradeACells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]

class AIS():
    def __init__(self):
        pass

    def prepare_data(self, datas):
        sio = io.StringIO(datas)
        df = pd.read_table(sio, sep=",")
        return df

    def DendvsThr(self, ax=None):
        df = self.prepare_data(aisdata)
        show = False
        if ax is None:
            f, ax = mpl.subplots(1,1)
            show = True
        colorl = []
        # map colors
        cells = df['cell'].values
        for i, cell in enumerate(cells):
            bc_index = gradeACells.index(int(cell))
            colorl.append(colors[bc_index])
        ax.scatter(df["DendAreas"], df["MeshInfl"], s=12, marker='o', c=colorl)
        # 2 and 5 are plotted separately here
        dfs = df[df["MeshInfl"].isnull()]
        colorl = []
        for cell in [2, 5]:
            bc_index = gradeACells.index(int(cell))
            colorl.append(colors[bc_index])
            ax.scatter(dfs["DendAreas"], dfs["NewStandardized"], s=20, marker='o', edgecolors=colorl, facecolors=["w", "w"], linewidths=1)
        dfn = df[df["MeshInfl"].notnull()]
        res = sp.stats.linregress(dfn["DendAreas"], dfn["MeshInfl"])
        xp = np.linspace(3000, 4500, 100)
        print(res.slope)
        print(res.intercept)
        print(f"R-squared: {res.rvalue**2:.6f}")
        print(f"P value: {res.pvalue:.5f}")
        ax.plot(xp, res.intercept + res.slope*xp, 'k-', linewidth=0.5)
        ax.text(x=1.0, y=0.1, s=f"p = {res.pvalue:.5f}, r$^2$={res.rvalue**2:5.3f}",
            horizontalalignment="right", fontsize=9,
            transform=ax.transAxes)
        ax.set_xlim(3000, 5000)
        ax.set_ylim(0.35, 0.65)
        ax.set_xlabel("Dendritic Area ($\mu m^2$)")
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

    def AISLengthThr(self, ax=None,):
        df = self.prepare_data(aisthrdata)
        show = False
        if ax is None:
            f, ax = mpl.subplots(1,1)
            show = True
        ax.set_clip_on(False)
        colnames = df.columns
        # print("colnames: ", colnames)
        # print(df["AIS_Length"])
        syms = ['o', 'o', '^', 's', 'v', 'd', 'o', 'o', '^', 'v']
        facecolors = colors
        # edgecolors = [None]*len(colors)

        for i, coln in enumerate(colnames):
            if i == 0:
                continue
            bc_index = gradeACells.index(int(coln))
            color = colors[bc_index]
            sym = syms[bc_index]
            if int(coln) in [13, 18, 30]:
                mk_face = "none"
                mk_edge = facecolors[bc_index]
            else:
                mk_face = facecolors[bc_index]
                mk_edge = mk_face
            print(coln, bc_index, color, sym, mk_face, mk_edge)
            bc_num = f"BC{int(coln):02d}"
            ax.plot(df["AIS_Length"].values, df[coln].values, 
                sym, markerfacecolor=mk_face, markeredgecolor = mk_edge, markeredgewidth=1,  
                linestyle='-', color=color, linewidth=1,
                label=bc_num)
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
        ax.legend(fontsize=7, loc="upper left", ncol=1, frameon=False)
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)
        if show:
            mpl.show()


if __name__ == "__main__":
    A = AIS()
    df = A.prepare_data(aisdata)
    print(df.head())
    #A.DendvsThr()
    A.AISLengthThr()