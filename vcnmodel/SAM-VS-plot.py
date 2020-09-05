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
from vcnmodel import cell_config as cell_config

import VS_data

# all runs 20 dB SPL 8/28, 8/29 2020 pbm 100ms pip
# '5','largestonly',0.,0,0.,0.,1.0,141.0200,7,
# data100Hz = """Cell,Configuration,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,maxArea,ninputs,
# 2,all,0.8776,1596,0.5797,0.0813,2458.1971,0.0000e+00,174.0100,6,
# 2,largestonly,0.8970,423,0.6739,0.0742,680.7433,1.5082e-148,174.0100,6,
# 2,removelargest,0.9433,1491,0.5986,0.0544,2653.7081,0.0000e+00,174.0100,6,
# 5,all,0.9875,1104,0.6053,0.0252,2153.1209,0.0000e+00,141.0200,7,
# 5,removelargest,0.9881,590,0.6168,0.0246,1152.1698,6.4493e-251,141.0200,7,
# 6,all,0.7635,1728,0.5785,0.1169,2014.3841,0.0000e+00,159.0100,6,
# 6,largestonly,0.8906,674,0.6712,0.0766,1069.1356,6.9209e-233,159.0100,6,
# 6,removelargest,0.8400,1633,0.5860,0.0940,2304.7505,0.0000e+00,159.0100,6,
# 9,all,0.5737,2079,0.6044,0.1678,1368.5723,6.5811e-298,245.7600,8,
# 9,largestonly,0.7399,1632,0.6358,0.1235,1787.0277,0.0000e+00,245.7600,8,
# 9,removelargest,0.7381,1753,0.5868,0.1240,1910.1438,0.0000e+00,245.7600,8,
# 10,all,0.8937,1586,0.5771,0.0754,2533.5692,0.0000e+00,121.5700,10,
# 10,largestonly,0.9329,44,0.7292,0.0593,76.5788,2.3503e-17,121.5700,10,
# 10,removelargest,0.9477,1533,0.5870,0.0521,2753.8836,0.0000e+00,121.5700,10,
# 11,all,0.5131,2202,0.6094,0.1839,1159.5771,1.5887e-252,213.0900,7,
# 11,largestonly,0.7409,1604,0.6279,0.1232,1761.1645,0.0000e+00,213.0900,7,
# 11,removelargest,0.6391,1918,0.5952,0.1506,1566.9475,0.0000e+00,213.0900,7,
# 13,all,0.7132,1789,0.5762,0.1309,1819.7771,0.0000e+00,146.5700,7,
# 13,largestonly,0.8666,1165,0.6499,0.0852,1749.6419,0.0000e+00,146.5700,7,
# 13,removelargest,0.8504,1623,0.5800,0.0906,2347.6915,0.0000e+00,146.5700,7,
# 17,all,0.6341,1945,0.5908,0.1519,1564.1857,0.0000e+00,278.3200,7,
# 17,largestonly,0.7475,1595,0.6252,0.1214,1782.6515,0.0000e+00,278.3200,7,
# 17,removelargest,0.7427,1737,0.5988,0.1228,1916.2568,0.0000e+00,278.3200,7,
# 30,all,0.7967,1684,0.5832,0.1073,2137.9377,0.0000e+00,172.1400,9,
# 30,largestonly,0.8789,918,0.6747,0.0809,1418.2080,1.0966e-308,172.1400,9,
# 30,removelargest,0.9473,1531,0.5921,0.0524,2747.5264,0.0000e+00,172.1400,9,
# """
#
# data400Hz="""Cell,Configuration,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,maxArea,ninputs,
# 2,all,0.7545,2288,0.7486,0.11952604.8635,0.0000e+00,174.0100,6,
# 2,largestonly,0.6511,165,0.9012,0.1474139.8866,4.2073e-31,174.0100,6,
# 2,removelargest,0.7573,1482,0.7845,0.11871699.9993,0.0000e+00,174.0100,6,
# 5,all,0.9111,168,0.7135,0.0687278.9283,2.7008e-61,141.0200,7,
# 5,removelargest,0.9171,46,0.7346,0.066277.3833,1.5719e-17,141.0200,7,
# 6,all,0.7836,2762,0.7247,0.11123391.6776,0.0000e+00,159.0100,6,
# 6,largestonly,0.6297,365,0.8996,0.1531289.4815,1.3801e-63,159.0100,6,
# 6,removelargest,0.7536,2357,0.7574,0.11972677.0857,0.0000e+00,159.0100,6,
# 9,all,0.7128,3600,0.7300,0.13103658.2071,0.0000e+00,245.7600,8,
# 9,largestonly,0.6976,2072,0.7494,0.13512016.9487,0.0000e+00,245.7600,8,
# 9,removelargest,0.7263,2836,0.7748,0.12732992.0559,0.0000e+00,245.7600,8,
# 10,all,0.7713,2548,0.7913,0.11473031.8213,0.0000e+00,121.5700,10,
# 10,largestonly,0.9519,3,0.0706,0.05005.4367,6.5982e-02,121.5700,10,
# 10,removelargest,0.7826,2143,0.8200,0.11142625.2039,0.0000e+00,121.5700,10,
# 11,all,0.7138,3823,0.7062,0.13073896.0910,0.0000e+00,213.0900,7,
# 11,largestonly,0.6926,2042,0.7377,0.13641958.8187,0.0000e+00,213.0900,7,
# 11,removelargest,0.7107,3110,0.7344,0.13153141.2778,0.0000e+00,213.0900,7,
# 13,all,0.7470,3083,0.7313,0.12163440.7632,0.0000e+00,146.5700,7,
# 13,largestonly,0.6475,1005,0.8557,0.1484842.6051,1.0731e-183,146.5700,7,
# 13,removelargest,0.7513,2504,0.7582,0.12042826.6539,0.0000e+00,146.5700,7,
# 17,all,0.7506,3092,0.6804,0.12063483.7467,0.0000e+00,278.3200,7,
# 17,largestonly,0.7057,2036,0.7292,0.13292028.1244,0.0000e+00,278.3200,7,
# 17,removelargest,0.7202,2440,0.7110,0.12902530.9454,0.0000e+00,278.3200,7,
# 30,all,0.7444,2679,0.7874,0.12232968.9036,0.0000e+00,172.1400,9,
# 30,removelargest,0.7485,2089,0.8388,0.12112340.8211,0.0000e+00,172.1400,9,
# 30,largestonly,0.6601,584,0.9440,0.1451508.9383,3.0581e-111,172.1400,9,
# """

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
