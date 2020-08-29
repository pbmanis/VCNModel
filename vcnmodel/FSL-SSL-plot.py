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

data = """
Cell,Configuration,FSL,FSLSD,SSL,SSLSD
'2','all',4.103,1.265,10.451,3.488
'2','removelargest',4.639,2.7861,4.312,7.685
'2','largestonly',14.972,21.533,43.593,41.588
'5','all',7.347,19.499,27.533,20.728
'5','removelargest',5.269,5.995,27.658,18.482
'5','largestonly',nan,nan,nan
'6','all',4.926,4.610,11.876,6.353
'6','largestonly',8.101,6.505,23.253,15.966
'6','removelargest',4.694,2.780,11.035,4.193
'9','all',4.247,1.421,9.345,3.489
'9','largestonly',4.341,0.865,9.912,2.270
'9','removelargest',4.268,1.436,8.950,1.789
'10','all',3.985,0.720,11.007,3.544
'10','largestonly',60.230,67.625,108.258,48.031
'10','removelargest',4.050,0.337,11.581,4.523
'11','largestonly',4.392,0.919,9.902,2.043
'11','all',4.263,1.491,9.271,3.267
'11','removelargest',4.097,1.070,8.645,2.173
'13','all',4.320,1.694,9.597,2.829
'13','largestonly',4.991,1.860,12.665,3.006
'13','removelargest',4.082,1.135,9.694,2.538
'17','all',4.216,1.278,10.499,4.829
'17','largestonly',4.457,1.432,10.245,2.440
'17','removelargest',4.299,2.424,9.800,3.297
'30','all',4.183,1.478,9.734,2.905
'30','largestonly',7.881,6.792,18.763,9.587
'30','removelargest',4.324,1.513,10.957,3.393
"""


def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


sio = io.StringIO(data)
df = pd.read_table(sio, sep=",")
maxv = 20.0
figstyle = ST.styler(journal='CerebralCortex', figuresize='single', font='Arial')

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
scp = sns.swarmplot(x="Configuration", y="FSLSD", data=df, hue="Cell", ax=P.axdict["B"])
# mpl.yscale("log")
scp.set_ylim(0.0, maxv)

PH.talbotTicks(P.axdict["B"], axes='y', pointSize=figstyle.Ticks['labelsize'])
mpl.setp(P.axdict["B"].get_xticklabels(), ha="right", rotation=30)
# scp.legend(loc='upper center',
#     bbox_to_anchor=(0.1, 1.0), ncol=2,
#     fontsize=8, markerscale=0.5)
scp.legend_.remove()
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
scpf = sns.swarmplot(x="Configuration", y="FSL", data=df, hue="Cell", ax=P.axdict["A"])
scpf.set_ylim(0.0, maxv)
# mpl.yscale("log")

scpf.legend(
    loc="upper center", bbox_to_anchor=(0.4, 1.0), ncol=2, 
    fontsize=figstyle.Legend['fontsize']*0.8, markerscale=0.5,
    fancybox=False, shadow=False, facecolor='w'
)
PH.talbotTicks(P.axdict["A"], axes='y', pointSize=figstyle.Ticks['labelsize'])
mpl.setp(P.axdict["A"].get_xticklabels(), ha="right", rotation=30)

# gmodel = Model(boltz)
#
# gmodel.set_param_hint("A", value=1, min=0.0, max=1.0, vary=True)
# gmodel.set_param_hint("vhalf", value=100., min=10., max=300.)
# gmodel.set_param_hint(
#     "k", value=5, min=0.01, max=200.0, vary=True
# )
# gparams = gmodel.make_params()
#
# result = gmodel.fit(
#     df.Eff, method="nedler", params=gparams, x=df.ASA,
# )
#
# mpl.text(200., result.params['A'].value-0.05, F"Max Eff: {result.params['A'].value:.2f} (1$\sigma$ = {result.params['A'].stderr:.2f})", fontsize=8,
#     verticalalignment='top')
# uni = r"$\mu m^2$"
# print(result.fit_report())
# mpl.text(150., result.params['A'].value/2.,
#         f"ASA at half max: {result.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result.params['vhalf'].stderr:.1f})\n" +
#         f"Slope: {result.params['k'].value:.1f} (1$\sigma$ = {result.params['k'].stderr:.1f})",
#     fontsize=8)
#
# xfit = np.linspace(0., 300., 300)
# bfit = gmodel.eval(params=result.params, x=xfit)
# ax.plot(xfit, bfit, 'b-')
# ax.set_xlabel(f"ASA ({uni:s})")
# ax.set_ylabel('Efficacy (Bushy spikes/input spikes)')

mpl.show()
