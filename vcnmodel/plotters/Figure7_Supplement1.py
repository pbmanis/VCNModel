"""Figure 7 Supplement 1
Replot of VS-SAM-AN-BUcnmodel-test.pzfx

Just plot the ANF data, 100% modulation, all on one plot vs frequency
"""

import io
from pathlib import Path

import matplotlib.pyplot as mpl
import matplotlib.patches
import numpy as np
import pandas as pd
import seaborn as sns
import pylibrary.plotting.plothelpers as PH
from scipy.io import loadmat
# from lmfit.models import LinearModel
from scipy.stats import linregress

# sns.set_style(rc={"pdf.fonttype": 42})
mpl.style.use("~/.matplotlib/figures.mplstyle")

ANFVSModData="""Freq,dBSPL,VSAN,
50, 0, 0.26
50, 5, 0.27
50, 10,	0.76
50, 15,	0.68
50, 20,	0.54
50, 25,	0.41
50, 30,	0.25
50, 35,	0.18
50, 40,	0.08
50, 45,	0.06
50, 50,	0.05
50, 55,	0.06
50, 60,	0.06
50, 65,	0.02
50, 70,	0.00
100, 0, 0.11
100, 5, 0.38
100, 10, 0.76
100, 15, 0.77
100, 20, 0.71
100, 25, 0.58
100, 30, 0.39
100, 35, 0.28
100, 40, 0.18
100, 45, 0.13
100, 50, 0.06
100, 55, 0.08
100, 60, 0.07
100, 65, 0.06
100, 70, 0.06
300, 0, 0.2
300, 5, 0.16
300, 10, 0.75
300, 15, 0.82
300, 20, 0.79
300, 25, 0.74
300, 30, 0.58
300, 35, 0.41
300, 40, 0.22
300, 45, 0.15
300, 50, 0.16
300, 55, 0.14
300, 60, 0.15
300, 65, 0.15
300, 70, 0.15
750, 0, 0.11
750, 5, 0.13
750, 10, 0.33
750, 15, 0.49
750, 20, 0.58
750, 25, 0.49
750, 30, 0.3
750, 35, 0.22
750, 40, 0.09
750, 45, 0.02
750, 50, 0.03
750, 55, 0.09
750, 60, 0.07
750, 65, 0.12
750, 70, 0.13
"""


class ANF_VS():
    def __init__(self):
        pass

    def prepare_data(self, datas):
        sio = io.StringIO(datas)
        df = pd.read_table(sio, sep=",")
        return df


    def VS_Mod(self, ax=None,):
        df = self.prepare_data(ANFVSModData)
        show = False
        if ax is None:
            f, ax = mpl.subplots(1,1)
            show = True
        ax.set_clip_on(False)
        colnames = df.columns

        freqs = list(set(df.Freq.values))

        ax.set_xlim(0, 70)
        ax.set_xlabel(f"Stimulus (dB SPL)")
        ax.set_ylabel(f"Vector Strength")
        ax.set_ylim(0.0, 1.0)

        rect1 = matplotlib.patches.Rectangle((14, 0), 2, 1.0, color="goldenrod", linewidth=None, alpha=0.3)
        ax.add_patch(rect1)
        rect2 = matplotlib.patches.Rectangle((29, 0), 2, 1.0, color="#cccccc", linewidth=None, alpha=0.5)
        ax.add_patch(rect2)
        palette=sns.color_palette("deep")

        for i, fr in enumerate(freqs):
            dff = df[df['Freq'] == fr]
            ax.plot(dff["dBSPL"], dff["VSAN"], color=palette[i], 
                linewidth=1.0, marker='o', markersize=6.0,
                label=f"{int(fr):d}",
                clip_on=False
            )
        xtpos = np.arange(0, 71, 10)
        ytpos = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        PH.set_axes_ticks(ax=ax,
            xticks = xtpos,
            xticks_str = [f"{int(x):d}" for x in xtpos],
           # xticks_pad:Union[List, None]=None,
            x_minor = np.arange(5, 70, 10),
            major_length = 3.0,
            minor_length =1.5,
            yticks = ytpos, 
            yticks_str=[f"{x:3.1f}" for x in ytpos], 
            yticks_pad=[1]*7,
            y_minor=[0.1, 0.3, 0.5, 0.7, 0.9],
            fontsize=8,
        )
        ax.legend(fontsize=9, loc="upper right", ncol=1, frameon=False, title="100% SAM (Hz)")
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)

        mpl.savefig(f"/Volumes/Pegasus_002/VCN-SBEM-Data/SBEM-paper Figures/Figure7/Figure7_supp/Figure7_Supplemental1_V3.pdf")
        mpl.savefig(f"/Volumes/Pegasus_002/VCN-SBEM-Data/SBEM-paper Figures/Figure7/Figure7_supp/Figure7_Supplemental1_V3.png")
        mpl.show()


def main():
    ANF = ANF_VS()
    df = ANF.VS_Mod()


if __name__ == '__main__':
    main()