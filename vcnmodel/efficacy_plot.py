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
159.01	0.51562500000000	6   2.8210
140.81	0.27200000000000	6   2.8210
132.17	0.19491525423729	6   2.8210
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



data_Full = """
N	Cell	syn#	nout	nin	Eff	ASA	SDRatio	nsites
0	2	0	98	1128	0.08688	174.01	3.23	139
1	2	1	1	1087	0.00092	124.89	3.23	100
2	2	2	0	1131	0.00000	108.53	3.23	87
3	2	3	0	1123	0.00000	101.16	3.23	81
4	2	4	0	1061	0.00000	80.73	3.23	65
5	2	5	0	1130	0.00000	45.89	3.23	37
0	5	0	0	1104	0.00000	141.02	2.80	113
1	5	1	0	1112	0.00000	106.09	2.80	85
2	5	2	0	1120	0.00000	99.89	2.80	80
3	5	3	0	1094	0.00000	93.86	2.80	75
4	5	4	0	1108	0.00000	93.52	2.80	75
5	5	5	0	1119	0.00000	89.85	2.80	72
6	5	6	0	1096	0.00000	53.78	2.80	43
0	6	0	212	1128	0.18794	159.01	2.82	127
1	6	1	75	1087	0.06900	140.81	2.82	113
2	6	2	45	1131	0.03979	132.17	2.82	106
3	6	3	3	1123	0.00267	109.34	2.82	87
4	6	4	2	1061	0.00189	97.97	2.82	78
5	6	5	0	1130	0.00000	94.22	2.82	75
0	9	0	821	1113	0.73765	245.76	2.52	196
1	9	1	299	1083	0.27608	132.27	2.52	106
2	9	2	100	1142	0.08757	107.92	2.52	86
3	9	3	1	1145	0.00087	74.43	2.52	59
4	9	4	1	1100	0.00091	71.50	2.52	57
5	9	5	1	1080	0.00093	62.98	2.52	50
6	9	6	0	1100	0.00000	54.81	2.52	44
7	9	7	0	1090	0.00000	35.03	2.52	28
0	10	0	5	1110	0.00450	121.57	2.92	97
1	10	1	1	1111	0.00090	112.88	2.92	90
2	10	2	0	1095	0.00000	71.71	2.92	57
3	10	3	0	1083	0.00000	71.30	2.92	57
4	10	4	0	1101	0.00000	69.63	2.92	56
5	10	5	0	1086	0.00000	66.74	2.92	53
6	10	6	0	1110	0.00000	51.80	2.92	41
7	10	7	0	1075	0.00000	46.39	2.92	37
8	10	8	0	1144	0.00000	42.08	2.92	34
9	10	9	0	1125	0.00000	40.09	2.92	32
0	11	0	788	1104	0.71377	213.09	2.44	170
1	11	1	612	1112	0.55036	156.38	2.44	125
2	11	2	439	1120	0.39196	134.46	2.44	107
3	11	3	55	1094	0.05027	89.76	2.44	72
4	11	4	15	1108	0.01354	77.04	2.44	62
5	11	5	0	1119	0.00000	36.76	2.44	29
6	11	6	0	1096	0.00000	35.67	2.44	29
0	13	0	426	1104	0.38587	146.57	2.50	117
1	13	1	24	1112	0.02158	94.11	2.50	75
2	13	2	15	1120	0.01339	91.38	2.50	73
3	13	3	12	1094	0.01097	90.88	2.50	73
4	13	4	4	1108	0.00361	77.02	2.50	62
5	13	5	0	1119	0.00000	57.42	2.50	46
6	13	6	0	1096	0.00000	56.98	2.50	46
0	17	0	801	1104	0.72554	278.32	2.73	222
1	17	1	755	1112	0.67896	261.49	2.73	209
2	17	2	5	1120	0.00446	105.06	2.73	84
3	17	3	0	1094	0.00000	62.36	2.73	50
4	17	4	0	1108	0.00000	41.04	2.73	33
5	17	5	0	1119	0.00000	38.19	2.73	31
6	17	6	0	1096	0.00000	36.75	2.73	29
0	30	0	291	1146	0.25393	172.14	2.65	138
1	30	1	0	1074	0.00000	88.27	2.65	71
2	30	2	0	1159	0.00000	79.11	2.65	63
3	30	3	0	1123	0.00000	78.50	2.65	63
4	30	4	0	1060	0.00000	66.28	2.65	53
5	30	5	0	1075	0.00000	63.15	2.65	50
6	30	6	0	1099	0.00000	59.21	2.65	47
7	30	7	0	1080	0.00000	42.12	2.65	34
8	30	8	0	1144	0.00000	41.92	2.65	33
"""


data_NoDend = """
N	Cell	syn#	nout	nin	Eff	ASA	SDRatio	nsites
0	2	0	981	1128	0.86968	174.01	3.23	139
1	2	1	843	1087	0.77553	124.89	3.23	100
2	2	2	811	1131	0.71706	108.53	3.23	87
3	2	3	761	1123	0.67765	101.16	3.23	81
4	2	4	612	1061	0.57681	80.73	3.23	65
5	2	5	45	1130	0.03982	45.89	3.23	37
0	5	0	740	1104	0.67029	141.02	2.80	113
1	5	1	329	1112	0.29586	106.09	2.80	85
2	5	2	235	1120	0.20982	99.89	2.80	80
3	5	3	133	1094	0.12157	93.86	2.80	75
4	5	4	133	1108	0.12004	93.52	2.80	75
5	5	5	103	1119	0.09205	89.85	2.80	72
6	5	6	0	1096	0.00000	53.78	2.80	43
0	6	0	873	1128	0.77394	159.01	2.82	127
1	6	1	786	1087	0.72309	140.81	2.82	113
2	6	2	790	1131	0.69850	132.17	2.82	106
3	6	3	630	1123	0.56100	109.34	2.82	87
4	6	4	509	1061	0.47974	97.97	2.82	78
5	6	5	459	1130	0.40619	94.22	2.82	75
0	9	0	1023	1113	0.91914	245.76	2.52	196
1	9	1	838	1083	0.77378	132.27	2.52	106
2	9	2	810	1142	0.70928	107.92	2.52	86
3	9	3	605	1145	0.52838	74.43	2.52	59
4	9	4	578	1100	0.52545	71.50	2.52	57
5	9	5	437	1080	0.40463	62.98	2.52	50
6	9	6	294	1100	0.26727	54.81	2.52	44
7	9	7	18	1090	0.01651	35.03	2.52	28
0	10	0	747	1110	0.67297	121.57	2.92	97
1	10	1	692	1111	0.62286	112.88	2.92	90
2	10	2	164	1095	0.14977	71.71	2.92	57
3	10	3	188	1083	0.17359	71.30	2.92	57
4	10	4	158	1101	0.14351	69.63	2.92	56
5	10	5	116	1086	0.10681	66.74	2.92	53
6	10	6	16	1110	0.01441	51.80	2.92	41
7	10	7	3	1075	0.00279	46.39	2.92	37
8	10	8	0	1144	0.00000	42.08	2.92	34
9	10	9	0	1125	0.00000	40.09	2.92	32
0	11	0	1002	1104	0.90761	213.09	2.44	170
1	11	1	956	1112	0.85971	156.38	2.44	125
2	11	2	922	1120	0.82321	134.46	2.44	107
3	11	3	754	1094	0.68921	89.76	2.44	72
4	11	4	662	1108	0.59747	77.04	2.44	62
5	11	5	44	1119	0.03932	36.76	2.44	29
6	11	6	52	1096	0.04745	35.67	2.44	29
0	13	0	943	1104	0.85417	146.57	2.50	117
1	13	1	750	1112	0.67446	94.11	2.50	75
2	13	2	742	1120	0.66250	91.38	2.50	73
3	13	3	742	1094	0.67824	90.88	2.50	73
4	13	4	623	1108	0.56227	77.02	2.50	62
5	13	5	306	1119	0.27346	57.42	2.50	46
6	13	6	330	1096	0.30109	56.98	2.50	46
0	17	0	992	1104	0.89855	278.32	2.73	222
1	17	1	998	1112	0.89748	261.49	2.73	209
2	17	2	822	1120	0.73393	105.06	2.73	84
3	17	3	376	1094	0.34369	62.36	2.73	50
4	17	4	32	1108	0.02888	41.04	2.73	33
5	17	5	23	1119	0.02055	38.19	2.73	31
6	17	6	10	1096	0.00912	36.75	2.73	29
0	30	0	981	1146	0.85602	172.14	2.65	138
1	30	1	778	1074	0.72439	88.27	2.65	71
2	30	2	748	1159	0.64538	79.11	2.65	63
3	30	3	743	1123	0.66162	78.50	2.65	63
4	30	4	612	1060	0.57736	66.28	2.65	53
5	30	5	573	1075	0.53302	63.15	2.65	50
6	30	6	503	1099	0.45769	59.21	2.65	47
7	30	7	180	1080	0.16667	42.12	2.65	34
8	30	8	151	1144	0.13199	41.92	2.65	33
"""
def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))

def boltz_resid(p, x, data):
    return(p['A'] * (1.0 / (1.0 + np.exp(-(x - p['vhalf']) / p['k']))) - data)

data = data_NoDend

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
# exit()

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

