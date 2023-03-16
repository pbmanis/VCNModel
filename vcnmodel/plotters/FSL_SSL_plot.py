import io

import lmfit
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit import Model
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as ST
# all runs 20 dB SPL 8/28, 8/29 2020 pbm 100ms pip

# Plot the FSL and SSl for each of tehse conditions, with different sets of input configuration.

data = """
Cell,Configuration,FSL,FSLSD,SSL,SSLSD,maxArea
2,all,4.103,1.265,10.451,3.488,174.0100
2,largestonly,14.972,21.533,43.593,41.588,174.0100
2,removelargest,4.639,2.7861,4.312,7.685,174.0100
5,all,7.347,19.499,27.533,20.728,141.0200
5,removelargest,5.269,5.995,27.658,18.482,141.0200
6,all,4.926,4.610,11.876,6.353,159.0100
6,largestonly,8.101,6.505,23.253,15.966,159.0100
6,removelargest,4.694,2.780,11.035,4.193,159.0100
9,all,4.247,1.421,9.345,3.489,245.7600
9,largestonly,4.341,0.865,9.912,2.270,245.7600
9,removelargest,4.268,1.436,8.950,1.789,245.7600
10,all,3.985,0.720,11.007,3.544,121.5700
10,largestonly,60.230,67.625,108.258,48.031,121.5700
10,removelargest,4.050,0.337,11.581,4.523,121.5700
11,largestonly,4.392,0.919,9.902,2.043,213.0900
11,all,4.263,1.491,9.271,3.267,213.0900
11,removelargest,4.097,1.070,8.645,2.173,213.0900
13,all,4.320,1.694,9.597,2.829,146.5700
13,largestonly,4.991,1.860,12.665,3.006,146.5700
13,removelargest,4.082,1.135,9.694,2.538,146.5700
17,all,4.216,1.278,10.499,4.829,278.3200
17,largestonly,4.457,1.432,10.245,2.440,278.3200
17,removelargest,4.299,2.424,9.800,3.297,278.3200
30,all,4.183,1.478,9.734,2.905,172.1400
30,largestonly,7.881,6.792,18.763,9.587,172.1400
30,removelargest,4.324,1.513,10.957,3.393,172.1400
"""


def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


sio = io.StringIO(data)
df = pd.read_table(sio, sep=",")
maxv = 20.0
figstyle = ST.styler(journal='Generic', figuresize='single', font='Arial')

P = PH.regular_grid(
    rows=1,
    cols=2,
    margins={
        "bottommargin": 0.25,
        "leftmargin": 0.1,
        "rightmargin": 0.05,
        "topmargin": 0.15,
    },
    horizontalspacing=0.15,
    figsize=(6, 4),
    labelposition=(-0.1, 1.05),
    panel_labels=["A", "B"],
    fontweight=figstyle.get_fontweights(),
    fontsize=figstyle.get_fontsizes()
)
sns.set()
sns.boxplot(
    x="Configuration",
    y="FSLSD",
    data=df,
    showcaps=False,
    boxprops={"facecolor": "None"},
    showfliers=False,
    whiskerprops={"linewidth": 0},
    ax=P.axdict["B"],
)
hue = "maxArea"
scp = sns.swarmplot(x="Configuration", y="FSLSD", data=df, hue=hue, ax=P.axdict["B"])
# mpl.yscale("log")
scp.set_ylim(0.0, maxv)
scp.set_ylabel("SD FSL (ms)")
PH.talbotTicks(P.axdict["B"], axes='y', pointSize=figstyle.Ticks['labelsize'])
mpl.setp(P.axdict["B"].get_xticklabels(), ha="right", rotation=30)
# scp.legend(loc='upper center',
#     bbox_to_anchor=(0.1, 1.0), ncol=2,
#     fontsize=8, markerscale=0.5)

sns.boxplot(
    x="Configuration",
    y="FSL",
    data=df,
    showcaps=False,
    boxprops={"facecolor": "None"},
    showfliers=False,
    whiskerprops={"linewidth": 0},
    ax=P.axdict["A"],
)
scpf = sns.swarmplot(x="Configuration", y="FSL", data=df, hue=hue, ax=P.axdict["A"])
scpf.set_ylim(0.0, maxv)
scpf.set_ylabel("FSL (ms)")
# mpl.yscale("log")
scpf.legend_.remove()
scp.legend(
    loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=2, 
    fontsize=figstyle.Legend['fontsize']*0.8, markerscale=0.5,
    fancybox=False, shadow=False, facecolor='w',
    labelspacing=0.2
)
PH.talbotTicks(P.axdict["A"], axes='y', pointSize=figstyle.Ticks['labelsize'])
mpl.setp(P.axdict["A"].get_xticklabels(), ha="right", rotation=30)


mpl.show()
