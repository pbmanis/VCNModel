
import io
import re
from typing import Union

import lmfit
import matplotlib
import numpy as np
import scipy.stats as SPS
import pandas as pd
import seaborn as sns
from lmfit import Model
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as STY

matplotlib.rcParams["mathtext.fontset"] = "stixsans"
# VCNc02*    1446.05    4674.18    0 14.41
# VCNc05*    1343.55    3755.39    0 3.06
data = """
Cell	SomaArea	DendriteArea	Group AIS
VCNc06	1464.27	3755.39	0 14.16
VCNc09	1340.36	3380.13	1 24.16
VCNc10	1388.7	4060.61	0 18.66
VCNc11	1288.45	3137.8	1 22.19
VCNc13	1305.31	3263.11	1 17.58
VCNc17	1357.62	3709.79	0 15.07
VCNc18	1292.47	3893.267	0 12.58
VCNc30	1508.36	3989.83	0 21.37
"""



spc = re.compile("[ ;,\t\f\v]+")  # format replacing all spaces with tabs
dataiter = re.finditer(spc, data)
data = re.sub(spc, ",", data)

sio = io.StringIO(data)
df = pd.read_table(sio, sep=",")
print(df.head())

df_0 = df[df.Group == 0]
df_1 = df[df.Group == 1]

u, p = SPS.mannwhitneyu(df_0.DendriteArea, df_1.DendriteArea)
print(u, p)
u, p = SPS.mannwhitneyu(df_0.AIS, df_1.AIS)
print(u, p)
print(df.sort_values(by=["Group"]))



f, ax = mpl.subplots(1, 1)
ax.margins(x=1, y=1)
ax = sns.barplot(
    x="Group",
    y="AIS",
    data=df,
    hue=[f"{c:s}" for c in df.Cell],
#    sizes=(200, 200),
    clip_on=False,
)
mpl.show()