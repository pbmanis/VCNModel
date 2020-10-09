import numpy as np
import io
import re

import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH
import seaborn as sns
import pandas as pd
import lmfit
from lmfit import Model

sns.set()

data = """
ASA	Eff	Cell    SDRatio
174.01	0.3984375	2 3.2324
124.89	0	2        3.2324
108.53	0	2        3.2324
101.16	0	2        3.2324
80.73	0	2        3.2324
45.89	0	2        3.2324
141.02	0	5   2.7951
106.09	0	5   2.7951
99.89	0	5   2.7951
93.86	0	5   2.7951
93.52	0	5   2.7951
89.85	0	5   2.7951
53.78	0	5   2.7951
159.01	0.5156250000000000	6   2.8210
140.81	0.27200000000000	6   2.8210
132.17	0.194915254237288	6   2.8210
109.34	0	6   2.8210
97.97	0.0078125	6   2.8210
94.22	0	6   2.8210
245.76	0.817460317460317	9   2.5218
132.27	0.548148148148148	9   2.5218
107.92	0.2578125	9   2.5218
74.43	0	9   2.5218
71.5	0	9   2.5218
62.98	0	9   2.5218
54.81	0	9   2.5218
35.03	0	9   2.5218
121.57	0	10
112.88	0.00719424460431655	10  2.9240
71.71	0	10  2.9240
71.3	0	10  2.9240
69.63	0	10  2.9240
66.74	0	10  2.9240
51.8	0	10  2.9240
46.39	0	10  2.9240
42.08	0	10  2.9240
40.09	0	10  2.9240
213.09	0.748091603053435	11  2.4353
156.383	0.725352112676056	11  2.4353
134.46	0.60431654676259	11  2.4353
89.76	0.136690647482014	11  2.4353
77.04	0.0230769230769231	11  2.4353
36.76	0	11  2.4353
35.67	0	11  2.4353
146.57	0.641221374045801	13  2.4999
94.11	0.0704225352112676	13  2.4999
91.38	0.0503597122302158	13  2.4999
90.88	0.0287769784172662	13  2.4999
77.02	0	13  2.4999
57.42	0	13  2.4999
56.98	0	13  2.4999
278.32	0.755725190839695	17  3.0123
261.49	0.788732394366197	17  3.0123
105.06	0.00719424460431655	17  3.0123
62.36	0	17  3.0123
41.04	0	17  3.0123
38.19	0	17  3.0123
36.75	0	17  3.0123
172.14	0.553030303030303	30  2.6451
88.27	0	30  2.6451
79.11	0	30  2.6451
78.5	0	30  2.6451
66.28	0	30  2.6451
63.15	0	30  2.6451
59.21	0	30  2.6451
42.12	0	30  2.6451
41.92	0	30  2.6451
"""

def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))

def boltz_resid(p, x, data):
    return(p['A'] * (1.0 / (1.0 + np.exp(-(x - p['vhalf']) / p['k']))) - data)


spc = re.compile("[ ;,\t\f\v]+") # format replacing all spaces with tabs
dataiter = re.finditer(spc, data)
data = re.sub(spc, ',', data)

sio = io.StringIO(data)
df = pd.read_table(sio, sep=',')
print(df.head())
f, ax = mpl.subplots(1,1)
ax.margins(x=1, y=1)
ax = sns.scatterplot(x="ASA", y="SDRatio", data=df, hue=[f"VCN\_c{c:02d}" for c in df.Cell], clip_on=False)
ax.set_xlim(0, 300)
ax.set_ylim(2, 3.05)

ax.legend(
        loc="lower right", bbox_to_anchor=(1.0, .0), ncol=2, 
        fontsize=9, markerscale=0.5,
        fancybox=False, shadow=False, facecolor='w',
        labelspacing=0.2
    )
mpl.show()
exit()

ax = sns.scatterplot(x='ASA', y="Eff", data=df, hue=[f"VCN\_c{c:02d}" for c in df.Cell])

ax.legend(
        loc="lower right", bbox_to_anchor=(1.0, 0.0), ncol=2, 
        fontsize=9, markerscale=0.5,
        fancybox=False, shadow=False, facecolor='w',
        labelspacing=0.2
    )


gmodel = Model(boltz)

gmodel.set_param_hint("A", value=1, min=0.0, max=1.0, vary=True)
gmodel.set_param_hint("vhalf", value=100., min=10., max=300.)
gmodel.set_param_hint(
    "k", value=5, min=0.01, max=200.0, vary=True
)
gparams = gmodel.make_params()
print('gparams: ', gparams.pretty_print())

weights = np.ones(len(df.ASA))
# weights[df.Eff<0.1] = 0.3

result_brute = gmodel.fit(
    df.Eff, method="brute", params=gparams, 
    x=df.ASA,
    weights=weights,
)
for p in result_brute.params:
    result_brute.params[p].stderr = 0. # abs(res2.params[p].value * 0)

print('\nBrute fit: ')
print(result_brute.params.pretty_print())
print('-'*80)
xfit = np.linspace(0., 300., 300)
bfit = gmodel.eval(params=result_brute.params, x=xfit)

methods = ['leastsq', 'least_squares', 'differential_evolution', 'basin_hopping', 'ampgo', 'nelder', 'lbfgsb', 'powell', 'cg',
     'cobyla', 'bfgs', #'trust-exact', 'trust-krylov', 'trust-constr', 
     #'dogleg', 
     'slsqp', 'shgo', 'dual_annealing']
# need jacobian: 'newton','tnc''trust-ncg'
resdict = {key: None for key in methods} 
ix = np.argsort(df.Eff)
y = df.ASA[ix]
x = df.Eff[ix]
weights = np.ones(len(df.Eff))# * np.random. *(df.Eff/np.max(df.Eff))**2
for meth in resdict.keys():
    result_LM = gmodel.fit(
        x, method=meth, params=result_brute.params, # gparams, 
        x=y,
        weights=weights,
        max_nfev = 20000, # fit_kws={'maxfev': 5000}
        jac = '3-point',
    
    )
    print(f'\n{meth:s} fit: ')
    print(result_LM.params.pretty_print())
    print('-'*80)
    xfit = np.linspace(0., 300., 300)
    lev_fit = gmodel.eval(params=result_LM.params, x=xfit)
    # lmfit.report_fit(result_brute.params, min_correl=0.5)

    for p in result_LM.params:
        if result_LM.params[p].stderr == None:
            result_LM.params[p].stderr = 0. # abs(res2.params[p].value * 0)
    resdict[meth] = result_LM
# ci, trace = lmfit.conf_interval(mini, res2, sigmas=[2, 2, 2], trace=True)
# lmfit.printfuncs.report_ci(ci)

mpl.text(200., result_LM.params['A'].value-0.05, F"LM: Max Eff: {result_LM.params['A'].value:.2f} (1$\sigma$ = {result_LM.params['A'].stderr:.2f})", fontsize=8,
    verticalalignment='top')

mpl.text(0., 0.8, F"Brute: Max Eff: {result_brute.params['A'].value:.2f} (1$\sigma$ = {result_brute.params['A'].stderr:.2f})", fontsize=8,
    verticalalignment='top')

uni = r"$\mu m^2$"
print('='*80)
mpl.text(150., 0.5*result_LM.params['A'],
        f"LM ASA at half max: {result_LM.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result_LM.params['vhalf'].stderr:.1f})\n" +
        f"Slope: {result_LM.params['k'].value:.1f} (1$\sigma$ = {result_LM.params['k'].stderr:.1f})",
    fontsize=8)


mpl.text(0., 0.7,
        f"Br ASA at half max: {result_brute.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result_brute.params['vhalf'].stderr:.1f})\n" +
        f"Slope: {result_brute.params['k'].value:.1f} (1$\sigma$ = {result_brute.params['k'].stderr:.1f})",
    fontsize=8)

# xd = df.ASA
# breff = gmodel.eval(params=result_brute.params, x=xd)
# brlms = np.sum((breff - df.Eff)**2)
# lmeff = gmodel.eval(params=result_LM.params, x=xd)
# lmlms = np.sum((lmeff - df.Eff)**2)
mpl.text(0., 0.6, f"Br $\chi^2$: {result_brute.chisqr:.3f}\nLM $\chi^2$:{result_LM.chisqr:.3f}", fontsize=8)

ax.plot(xfit, bfit, 'b-', label="Brute")
for m in resdict.keys():
    fit = gmodel.eval(params=resdict[m].params, x=xfit)
    ax.plot(xfit, fit, label=meth)
ax.set_xlabel(f"ASA ({uni:s})")
ax.set_ylabel('Efficacy (Bushy spikes/input spikes)')

mpl.show()

