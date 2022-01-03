"""
SAM-VS plot generated graphs that summarize the Sinusoidal Amplitude Modulation
simulation results.


"""
from dataclasses import dataclass, field
from pathlib import Path
import string
import typing
from typing import Union
import io
import re
import lmfit
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit import Model
from matplotlib import pyplot as mpl
import toml

from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as ST
from pylibrary.tools import utility as PU
from vcnmodel import cell_config as cell_config

import VS_data_15dB   # VS data from 12-3 to 12-8 2021
import VS_data_30dB   # VS_data from 8_16_2021

config = toml.load(open("wheres_my_data.toml", "r"))

"""
Plot the vector strength measures as a function of frequency for VCNModel
Reads the data from VS_data_xxdB.py

9 Oct 2020 pbm
12 2021 pbm - plot VS for both 15 and 30 dB on same figure.

"""

# datas = {'100': data100Hz, '400': data400Hz}  # just keep adding...
panels = {50: ["A", "B"], 100: ["C","D"], 200: ['E', 'F'], 
        300: ['G', 'H'], 400: ['I', "J"], 500: ['K', 'L'], 750: ["M", "N"], 1000:["O", "P"]}


# plotwhat = "phasesd"
# plotwhat = "phasesdtime"

def reset_style():
    """
    Reset the figure style to use an AI-compatible font
    """
    sns.set_style(rc={"pdf.fonttype": 42})
    mpl.style.use("~/.matplotlib/figures.mplstyle")
    

def setup_figure():
        xfig = 8.0
        yfig = 8.0
        yn = 2
        xn = 4
        ymar = 0.65
        xmar = 0.65
        xfigw = xfig - xmar*2
        yfigh = yfig - ymar*2
        rowh = (yfigh-ymar - (yn-1)*ymar)/yn
        rows = np.linspace(yfigh - rowh, ymar, yn)
        cols = np.linspace(xmar, (xfigw- xmar), xn)
        colw = (xfigw-(xn-1)*xmar)/xn

        sizer = {}
        for i, l in enumerate(string.ascii_uppercase[:xn*yn]):
            icol = i % 4
            irow = int(i/4)
            sizer[l] = {
                "pos": [cols[icol], colw, rows[irow], rowh],
                "labelpos": (-0.15, 1.02),
                "noaxes": False,
            }
        # P = PH.arbitrary_grid(
        #     sizer,
        #     order="columnsfirst",
        #     units="in",
        #     figsize=(xfig, yfig),
        #     label=True,
        #     # showgrid=True,
        # )
        P = PH.regular_grid(2, 4, order="rowsfirst", units="in",
            horizontalspacing=0.5, verticalspacing=0.5,
            margins ={'bottommargin': 0.5, 'leftmargin': 0.5, 'rightmargin': 0.5, 'topmargin': 0.5},
            figsize=(8, 5), label=True,
            showgrid=False,
            )
        return P

def vs_an_freq(df):
    """
    get average an vector strength for these runs, and plot a horizontal line.
    """
    freqs = set(df['frequency'])
 
    vs_by_freq = {key: [] for key in freqs}
    
    for c in df["Configuration"]:
        dx1 = df.loc[df['Configuration']==c]
        for f in dx1["frequency"]:
        #df. loc[((df['a'] > 1) & (df['b'] > 0)) | ((df['a'] < 1) & (df['c'] == 100))]
            dfc = dx1.loc[(df["frequency"]==f)]['AN_VS'].values
            vs_by_freq[f].extend(dfc)

    for f in vs_by_freq.keys():
        vs_by_freq[f] = np.mean(vs_by_freq[f])
    return vs_by_freq

def plot_summary(measure, df=None, P=None, axname=None, hue='Cell', legend=False):
    """
    Make a summary plot 
    """
    if measure == 'phasesdtime':
        dataframe['phasesdtime'] = 1e3*(df['phasesd']/np.pi)/df['frequency']

    sns.boxplot(
        x="Configuration",
        y=measure,
        data=df,
        showcaps=False,
        boxprops={"facecolor": "None"},
        showfliers=False,
        whiskerprops={"linewidth": 0},
        ax=P.axdict[axname],
    )
    scp = sns.swarmplot(x="Configuration", y=measure, data=df, hue=hue, ax=P.axdict[axname],
        clip_on=False)
 
    if measure == 'phase':
        scp.set_ylim(0.0, 2.0*np.pi)
        scp.set_ylabel('Phase (radians)')
    else:
        scp.set_ylim(0.0, 2.0)
    if legend:
        scp.legend(
        loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=2, 
        fontsize=figstyle.Legend['fontsize']*0.7, markerscale=0.5,
        fancybox=False, shadow=False, facecolor='w',
        labelspacing=0.25
    )
    else:
        scp.legend().remove()

    PH.talbotTicks(P.axdict[axname], axes='y', pointSize=figstyle.Ticks['labelsize'],
        tickPlacesAdd={'x': 0, 'y': 2}, floatAdd={'x': 0, 'y': 2})
    mpl.setp(P.axdict[axname].get_xticklabels(), ha="right", rotation=30)


def plot_summary2(measure, df=None, P=None, axname=None, hue='Configuration', legend=False):
    if measure == 'phasesdtime':
        df['phasesdtime'] = 1e3*(df['phasesd']/np.pi)/df['frequency']
    sns.boxplot(
        x="frequency",
        y=measure,
        data=df,
        showcaps=False,
        boxprops={"facecolor": "None"},
        showfliers=False,
        whiskerprops={"linewidth": 0},
        ax=P.axdict[axname],
    )
    # print(dataframe[hue])
    scp = sns.swarmplot(x="frequency", y=measure, data=df, hue=hue, ax=P.axdict[axname],
        clip_on=False)
        # size=[40, 400], sizes="maxArea")
    # mpl.yscale("log")
    if measure == 'phase':
        scp.set_ylim(0.0, 2.0)
        scp.set_ylabel('Phase (radians)')
    else:
        scp.set_ylim(0.0, 1.0)
    if legend:
        scp.legend(
        loc="lower right", bbox_to_anchor=(1.0, 0.0), ncol=2, 
        fontsize=figstyle.Legend['fontsize']*0.7, markerscale=0.5,
        fancybox=False, shadow=False, facecolor='w',
        labelspacing=0.25
    )
    else:
        scp.legend().remove()

    PH.talbotTicks(P.axdict[axname], axes='y', pointSize=figstyle.Ticks['labelsize'],
        tickPlacesAdd={'x': 0, 'y': 2}, floatAdd={'x': 0, 'y': 2})
    mpl.setp(P.axdict[axname].get_xticklabels(), ha="right", rotation=30)

"""
Plot selected cells for a poster, 1 column, two rows
"""

class Plot_Strength_Groups():
    def __init__(self, plotwhat, df=None, P=None, panels=None, legend=None):
        self.df = df
        self.ax1 = P.axdict[panels[0]]
        self.ax2 = P.axdict[panels[1]]
        self.legend = legend
        print("what: ", plotwhat)
        self.make_plots()

    def _scalefun(self, x):
        return [f"{y:.2f}" for y in x]

    def _getCell(self, cell):
        return(f"VCN_c{int(cell):02d}")

    # df = df.sort_values('Cell')
    def make_plots(self):
        labels = set(self.df['frequency'])
        labels = [str(int(l)) for l in sorted(labels)]
        sizes = (2, 6)
        self.df = self.df[self.df["Configuration"].isin(["all"])]
        # self.ax1.set_yscale('squareroot')
        scp1 = sns.lineplot(data=self.df, 
                x="frequency",
                y=plotwhat,
                hue="strength",
                marker='o',
                markersize=9,
                hue_order=['subthreshold', 'one', 'two'],
                sizes=sizes,
                size="strength",
                ax=self.ax1,
                palette="colorblind",
                # estimator=None,
        )
        if self.legend != 1:
            scp1.legend().remove()
            
        self.ax1.plot(self.df['frequency'][0:7], self.df['AN_VS'][0:7], 'r-', linewidth=1.5)  # only 7 frequencies

        # sns.lineplot(data=self.df,
        #             x="frequency", y="AN_VS",
        #             hue_order=['subthreshold', 'one', 'two'],
        #                  sizes=sizes,
        # #             size="strength",
        #              ax=self.ax1, palette="husl",
        #             # estimator=None,
        # )


        self.ax1.set_xlim(40, 1100)
        self.ax1.set_xscale("log")
        self.ax1.set_ylim((0, 1))

        # ax2.set_yscale('squareroot')
        scp2 = sns.lineplot(data=self.df, 
            x="frequency", y='phasesd',
            hue="strength", marker='o', markersize=9,
                        hue_order=['subthreshold', 'one', 'two'],
        #                 sizes=sizes,
                        size="strength",
                         ax=self.ax2, palette="colorblind",
                    # estimator=None,
        )
        if self.legend != 2:
            scp2.legend().remove()


        self.ax2.set_xlim(40, 1100)
        self.ax2.set_xscale("log")
        self.ax2.set_ylim((0, 2*np.pi))


reset_style()
maxtc = 0
maxvs = 1.0
figstyle = ST.styler(journal='CerebralCortex', figuresize='single', font='Arial')
plabels = [[panels[p][0], panels[p][1]] for p in panels]
plabels = PU.flatten(plabels)

def individual_plots(plotwhat, df=None, P=None, panels=None):

    sns.set()
    hue = "Cell"
    # # df['phase'] = df['phase']/(np.pi)  # convert to radians

    spc = re.compile("[ ;,\t\f\v]+") # format replacing all spaces or commas or semicolons or form feeds or vertical feeds with tabs
    for i, k in enumerate(panels.keys()):
        # data = datas[k] # clean up the table
        # datas[k] = re.sub(spc, ',', datas[k])  # just do it, make comma sep
        dff = df.loc[df['frequency'] == int(k)]
        pa = vs_an_freq(dff)
        if i == 0:
            legend = True
        else:
            legend = False
        plot_summary(dff, plotwhat, P, panels[int(k)][0], hue=hue, legend=False)
        PH.referenceline(P.axdict[panels[k][0]], pa[k], [-0.5,3], 'k')
        plot_summary(dff, 'phase', P, panels[int(k)][1], hue=hue, legend=legend)

        if panels[k][0] not in ['I']:
            P.axdict[panels[k][0]].xaxis.label.set_visible(False)
            P.axdict[panels[k][0]].tick_params(labelbottom=False)
        if panels[k][1] not in ['J']:
            P.axdict[panels[k][1]].xaxis.label.set_visible(False)
            P.axdict[panels[k][1]].tick_params(labelbottom=False)

        tfont = {'fontsize': 11,
         'fontweight' : 'normal',
         'verticalalignment': 'baseline',
         'horizontalalignment': 'center'}

        P.axdict[panels[k][0]].set_title(f"{int(k):d}Hz", fontdict=tfont)
    return P

def plot_freq_manipulation(plotwhat, df=None, P=None, panels=None, legend=None):
    """
    Swarm plots of VS and phase/phasetime of each cell, with the
    different manipulations of the inputs (all, largestonly, nolargest)
    """

    sns.set()
    hue = "Cell"
    # df['phase'] = df['phase']/(np.pi)  # convert to radians
    # df['phasesd'] = df['phasesd']/np.pi  # also the phase SD
    if plotwhat == 'phasesdtime':
        df['phasesdtime'] = 1e3*(df['phasesd']/np.pi)/df['frequency']

    labels=['50', '100','200','400','750', '1000']

    hue='Configuration'
    scp1 = sns.swarmplot(x="frequency", y=plotwhat, data=df, hue=hue,
            ax=P.axdict[panels[0]],
            clip_on=False, s=3)
    if legend == 1:
        scp1.legend(
            loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1,
            fontsize=figstyle.Legend['fontsize']*0.7, markerscale=0.5,
            fancybox=False, shadow=False, facecolor='w',
            labelspacing=0.25
        )
    frx = vs_an_freq(df)
    frx = {k: frx[k] for k in sorted(frx)}
    w = 0.5
    for i, f in enumerate(frx.keys()):
        xv = [i-w, i+w]
        PH.referenceline(P.axdict[panels[0]],
            frx[f], xv, linewidth=1.25, color='r')

    scp2 = sns.swarmplot(x="frequency", y="phase", data=df, hue=hue, 
            ax=P.axdict[panels[1]],
            clip_on=False, s=3)
    if legend != 2:
        scp2.legend().remove()
    # P.axdict[panels[0]].set_xlim(40, 1100)
    # P.axdict[panels[1]].set_xlim(40, 1100)


def convert_data(data):
    sio = io.StringIO(data.data)
    df = pd.read_table(sio, sep=",")

    def _label(d):
        if d in [5, 10, 2, 6, 30]:
            return('subthreshold')
        elif d in [9, 13, 18]:
            return('one')
        elif d in [11,17]:
            return 'two'
        else:
            raise ValueError(f"Failed to map cell number {str(d):s} {str(type(d)):s} to strength")
    print(df['Cell'].values)
    df['strength'] = [_label(int(x)) for i, x in enumerate(df['Cell']) if x != "Cell"]
    df = df.sort_values(["Cell", "frequency", "Configuration"], ascending = (True, True, True))
    df = df.reset_index()
    return df

P = setup_figure()  # make overall plot
# P = individual_plots()
datas_15 = convert_data(VS_data_15dB)
datas_30 = convert_data(VS_data_30dB)
plotwhat = "VectorStrength"

plot_freq_manipulation(plotwhat, df=datas_15, P=P, panels=['A', 'B'])
PSG1 = Plot_Strength_Groups(plotwhat, df=datas_15, P=P, panels=['C', 'D'], legend=1)
plot_freq_manipulation(plotwhat, df=datas_30, P=P, panels=['E', 'F'])
PSG2 = Plot_Strength_Groups(plotwhat, df=datas_30, P=P, panels=['G', 'H'], legend=None)

# save_file = f"Fig_M6_VS_SAM_Supplmental.pdf"
# P.figure_handle.text(
#     0.99,
#     0.01,
#     f"All Cells",
#     fontsize=10,
#     horizontalalignment="right",
#     verticalalignment="bottom",
# )
# ftxt = r"${savefile:s}$"
# P.figure_handle.text(
#     0.99,
#     0.99,
#     ftxt,  # .replace('_', '\_'),
#     transform=P.figure_handle.transFigure,
#     horizontalalignment="right",
#     verticalalignment="top",
# )
# mpl.savefig(
#     Path(config["baseDataDirectory"], "Figures", save_file),
#     metadata={
#         "Creator": "Paul Manis",
#         "Author": "Paul Manis",
#         "Title": "SBEM Project Figure 6 Modeling Supplemental : VS to SAM, All cells",
#     },
# )

mpl.show()
