import io
import re
import lmfit
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit import Model
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as ST
from pylibrary.tools import utility as PU
# import plotnine as PN
from vcnmodel import cell_config as cell_config

import VS_data
"""
Plot the vector strength measures as a function of frequency for VCNModel
Reads the data from VS_data.py

9 Oct 2020 pbm

"""

# datas = {'100': data100Hz, '400': data400Hz}  # just keep adding...
panels = {100: ["A", "B"], 200: ["C","D"], 400: ['E', 'F'], 750: ['G', 'H'], 1000: ["I", "J"]}

def vs_an_freq(df):
    # get average an vector strength for these runs, and plot a horizontal line.
    freqs = set(df['frequency'])
 
    vs_by_freq = {key: [] for key in freqs}
    
    for c in df["Configuration"]:
        dx1 = df.loc[df['Configuration']==c]
        for f in dx1["frequency"]:
        #df. loc[((df['a'] > 1) & (df['b'] > 0)) | ((df['a'] < 1) & (df['c'] == 100))]
            dfc = dx1.loc[(df["frequency"]==f)]['ANVS'].values
            vs_by_freq[f].extend(dfc)

    for f in vs_by_freq.keys():
        vs_by_freq[f] = np.mean(vs_by_freq[f])
    return vs_by_freq

def plot_summary(dataframe, measure, P, axname, hue='Cell', legend=False):
    sns.boxplot(
        x="Configuration",
        y=measure,
        data=dataframe,
        showcaps=False,
        boxprops={"facecolor": "None"},
        showfliers=False,
        whiskerprops={"linewidth": 0},
        ax=P.axdict[axname],
    )
    scp = sns.swarmplot(x="Configuration", y=measure, data=dataframe, hue=hue, ax=P.axdict[axname],
        clip_on=False)
 
    if measure == 'phase':
        scp.set_ylim(0.0, 2.0)
        scp.set_ylabel('Phase (radians)')
    else:
        scp.set_ylim(0.0, 1.0)
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


def plot_summary2(dataframe, measure, P, axname, hue='Configuration', legend=False):
    sns.boxplot(
        x="frequency",
        y=measure,
        data=dataframe,
        showcaps=False,
        boxprops={"facecolor": "None"},
        showfliers=False,
        whiskerprops={"linewidth": 0},
        ax=P.axdict[axname],
    )
    # print(dataframe[hue])
    scp = sns.swarmplot(x="frequency", y=measure, data=dataframe, hue=hue, ax=P.axdict[axname],
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

maxtc = 0        
maxvs = 1.0
figstyle = ST.styler(journal='CerebralCortex', figuresize='single', font='Arial')
plabels = [[panels[p][0], panels[p][1]] for p in panels]
plabels = PU.flatten(plabels)
P = PH.regular_grid(
    rows=len(panels),
    cols=2,
    order="rowsfirst",
    margins={
        "bottommargin": 0.25,
        "leftmargin": 0.1,
        "rightmargin": 0.05,
        "topmargin": 0.05,
    },
    horizontalspacing=0.15,
    verticalspacing=0.04,
    figsize=(6, 7),
    labelposition=(-0.1, 1.05),
    panel_labels=plabels,
    fontweight=figstyle.get_fontweights(),
    fontsize=figstyle.get_fontsizes()
)
sns.set()
hue = "Cell"
datas = VS_data.data
sio = io.StringIO(datas)
df = pd.read_table(sio, sep=",")
df['phase'] = df['phase']/(np.pi)  # convert to radians

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
    plot_summary(dff, 'VectorStrength', P, panels[int(k)][0], hue=hue, legend=False)
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


P2 = PH.regular_grid(
    rows=1,
    cols=2,
    order="rowsfirst",
    margins={
        "bottommargin": 0.25,
        "leftmargin": 0.1,
        "rightmargin": 0.05,
        "topmargin": 0.15,
    },
    horizontalspacing=0.15,
    figsize=(8, 4),
    labelposition=(-0.1, 1.05),
    panel_labels=['A', 'B'],
    fontweight=figstyle.get_fontweights(),
    fontsize=figstyle.get_fontsizes()
)
sns.set()
hue = "Cell"
datas = VS_data.data
sio = io.StringIO(datas)
df = pd.read_table(sio, sep=",")
df['phase'] = df['phase']/(np.pi)  # convert to radians
# df['frequency'] = int(df['frequency'])
bins = [-np.inf, 150, 250, 450, 800, 1200]

labels=['100','200','400','750', '1000']

# df['frequency'] = pd.cut(df.frequency, bins=bins, labels=labels)
# print(df['frequency'])
hue='Configuration'
# plot_summary2(dff, 'VectorStrength', P2, 'A', hue=hue, legend=False)
# plot_summary2(dff, 'phase', P2, 'B', hue=hue, legend=True)
# P2.axdict['A'].scatter(df['frequency'], df['VectorStrength'])
scp1 = sns.swarmplot(x="frequency", y="VectorStrength", data=df, hue=hue, ax=P2.axdict['A'],
        clip_on=False, s=3)
scp1.legend(
        loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1, 
        fontsize=figstyle.Legend['fontsize']*0.7, markerscale=0.5,
        fancybox=False, shadow=False, facecolor='w',
        labelspacing=0.25
    )
frx = vs_an_freq(df)
print(frx)
frx = {k: frx[k] for k in sorted(frx)}
w = 0.5
for i, f in enumerate(frx.keys()):
    xv = [i-w, i+w]
    PH.referenceline(P2.axdict['A'], frx[f], xv, linewidth=1.25, color='r')

scp2 = sns.swarmplot(x="frequency", y="phase", data=df, hue=hue, ax=P2.axdict['B'],
        clip_on=False, s=3)
scp2.legend().remove()
mpl.show()
