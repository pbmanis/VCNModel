import numpy as np
import io
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH
import seaborn as sns
import pandas as pd
import lmfit
from lmfit import Model

sns.set()

data = """
ASA	Eff	Cell
174.01	0.3984375	'2'
124.89	0	'2'
108.53	0	'2'
101.16	0	'2'
80.73	0	'2'
45.89	0	'2'
141.02	0	'5'
106.09	0	'5'
99.89	0	'5'
93.86	0	'5'
93.52	0	'5'
89.85	0	'5'
53.78	0	'5'
159.01	0.515625	'6'
140.81	0.272	'6'
132.17	0.194915254237288	'6'
109.34	0	'6'
97.97	0.0078125	'6'
94.22	0	'6'
245.76	0.817460317460317	'9'
132.27	0.548148148148148	'9'
107.92	0.2578125	'9'
74.43	0	'9'
71.5	0	'9'
62.98	0	'9'
54.81	0	'9'
35.03	0	'9'
121.57	0	'10'
112.88	0.00719424460431655	'10'
71.71	0	'10'
71.3	0	'10'
69.63	0	'10'
66.74	0	'10'
51.8	0	'10'
46.39	0	'10'
42.08	0	'10'
40.09	0	'10'
213.09	0.748091603053435	'11'
156.383	0.725352112676056	'11'
134.46	0.60431654676259	'11'
89.76	0.136690647482014	'11'
77.04	0.0230769230769231	'11'
36.76	0	'11'
35.67	0	'11'
146.57	0.641221374045801	'13'
94.11	0.0704225352112676	'13'
91.38	0.0503597122302158	'13'
90.88	0.0287769784172662	'13'
77.02	0	'13'
57.42	0	'13'
56.98	0	'13'
278.32	0.755725190839695	'17'
261.49	0.788732394366197	'17'
105.06	0.00719424460431655	'17'
62.36	0	'17'
41.04	0	'17'
38.19	0	'17'
36.75	0	'17'
172.14	0.553030303030303	'30'
88.27	0	'30'
79.11	0	'30'
78.5	0	'30'
66.28	0	'30'
63.15	0	'30'
59.21	0	'30'
42.12	0	'30'
41.92	0	'30'
"""

def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))

sio = io.StringIO(data)
df = pd.read_table(sio, sep='\t')
ax = sns.scatterplot(x='ASA', y="Eff", data=df, hue='Cell')

gmodel = Model(boltz)

gmodel.set_param_hint("A", value=1, min=0.0, max=1.0, vary=True)
gmodel.set_param_hint("vhalf", value=100., min=10., max=300.)
gmodel.set_param_hint(
    "k", value=5, min=0.01, max=200.0, vary=True
)
gparams = gmodel.make_params()

result = gmodel.fit(
    df.Eff, method="nedler", params=gparams, x=df.ASA,
)

mpl.text(200., result.params['A'].value-0.05, F"Max Eff: {result.params['A'].value:.2f} (1$\sigma$ = {result.params['A'].stderr:.2f})", fontsize=8, 
    verticalalignment='top')
uni = r"$\mu m^2$"
print(result.fit_report())
mpl.text(150., result.params['A'].value/2., 
        f"ASA at half max: {result.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result.params['vhalf'].stderr:.1f})\n" +
        f"Slope: {result.params['k'].value:.1f} (1$\sigma$ = {result.params['k'].stderr:.1f})",
    fontsize=8)

xfit = np.linspace(0., 300., 300)
bfit = gmodel.eval(params=result.params, x=xfit)
ax.plot(xfit, bfit, 'b-')
ax.set_xlabel(f"ASA ({uni:s})")
ax.set_ylabel('Efficacy (Bushy spikes/input spikes)')

mpl.show()

