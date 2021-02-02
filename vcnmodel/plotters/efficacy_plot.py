import io
import re
from typing import Union

import lmfit
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit import Model
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as STY

matplotlib.rcParams["mathtext.fontset"] = "stixsans"
# matplotlib.rcParams['font.family'] = 'sans-serif'
# matplotlib.rcParams['font.sans-serif'] = ['stixsans'] #, 'Tahoma', 'DejaVu Sans',
#                                #'Lucida Grande', 'Verdana']
# matplotlib.rcParams["pdf.fonttype"] = 42
# matplotlib.rcParams["text.usetex"] = False
# sns.set_style(rc={"pdf.fonttype": 42})


"""
Plot the efficacy data into an axis
"""
# sns.set()

data_original = """
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
0	2	0	76	1128	0.06738	174.01	3.23	134
1	2	1	1	1087	0.00092	124.89	3.23	96
2	2	2	0	1131	0.00000	108.53	3.23	83
3	2	3	0	1123	0.00000	101.16	3.23	78
4	2	4	0	1061	0.00000	80.73	3.23	62
5	2	5	0	1130	0.00000	45.89	3.23	35
0	5	0	0	1104	0.00000	141.02	2.80	108
1	5	1	0	1112	0.00000	106.09	2.80	82
2	5	2	0	1120	0.00000	99.89	2.80	77
3	5	3	0	1094	0.00000	93.86	2.80	72
4	5	4	0	1108	0.00000	93.52	2.80	72
5	5	5	0	1119	0.00000	89.85	2.80	69
6	5	6	0	1096	0.00000	53.78	2.80	41
0	6	0	152	1128	0.13475	159.01	2.82	122
1	6	1	56	1087	0.05152	140.81	2.82	108
2	6	2	40	1131	0.03537	132.17	2.82	102
3	6	3	3	1123	0.00267	109.34	2.82	84
4	6	4	0	1061	0.00000	97.97	2.82	75
5	6	5	0	1130	0.00000	94.22	2.82	72
0	9	0	794	1113	0.71339	245.76	2.52	189
1	9	1	265	1083	0.24469	132.27	2.52	102
2	9	2	56	1142	0.04904	107.92	2.52	83
3	9	3	0	1145	0.00000	74.43	2.52	57
4	9	4	1	1100	0.00091	71.50	2.52	55
5	9	5	0	1080	0.00000	62.98	2.52	48
6	9	6	0	1100	0.00000	54.81	2.52	42
7	9	7	0	1090	0.00000	35.03	2.52	27
0	10	0	3	1110	0.00270	121.57	2.92	93
1	10	1	2	1111	0.00180	112.88	2.92	87
2	10	2	0	1095	0.00000	71.71	2.92	55
3	10	3	0	1083	0.00000	71.30	2.92	55
4	10	4	0	1101	0.00000	69.63	2.92	54
5	10	5	0	1086	0.00000	66.74	2.92	51
6	10	6	0	1110	0.00000	51.80	2.92	40
7	10	7	0	1075	0.00000	46.39	2.92	36
8	10	8	0	1144	0.00000	42.08	2.92	32
9	10	9	0	1125	0.00000	40.09	2.92	31
0	11	0	779	1104	0.70562	213.09	2.44	164
1	11	1	567	1112	0.50989	156.38	2.44	120
2	11	2	404	1120	0.36071	134.46	2.44	103
3	11	3	36	1094	0.03291	89.76	2.44	69
4	11	4	5	1108	0.00451	77.04	2.44	59
5	11	5	0	1119	0.00000	36.76	2.44	28
6	11	6	0	1096	0.00000	35.67	2.44	27
0	13	0	387	1104	0.35054	146.57	2.50	113
1	13	1	17	1112	0.01529	94.11	2.50	72
2	13	2	10	1120	0.00893	91.38	2.50	70
3	13	3	13	1094	0.01188	90.88	2.50	70
4	13	4	0	1108	0.00000	77.02	2.50	59
5	13	5	0	1119	0.00000	57.42	2.50	44
6	13	6	0	1096	0.00000	56.98	2.50	44
0	17	0	780	1104	0.70652	278.32	2.73	214
1	17	1	753	1112	0.67716	261.49	2.73	201
2	17	2	7	1120	0.00625	105.06	2.73	81
3	17	3	0	1094	0.00000	62.36	2.73	48
4	17	4	0	1108	0.00000	41.04	2.73	32
5	17	5	0	1119	0.00000	38.19	2.73	29
6	17	6	0	1096	0.00000	36.75	2.73	28
0	18	0	723	1113	0.64960	273.22	3.01	210
1	18	1	89	1083	0.08218	164.19	3.01	126
2	18	2	1	1142	0.00088	103.85	3.01	80
3	18	3	0	1145	0.00000	70.54	3.01	54
4	18	4	0	1100	0.00000	60.30	3.01	46
5	18	5	0	1080	0.00000	55.66	3.01	43
6	18	6	0	1100	0.00000	45.18	3.01	35
7	18	7	0	1090	0.00000	37.12	3.01	29
0	30	0	232	1146	0.20244	172.14	2.65	132
1	30	1	0	1074	0.00000	88.27	2.65	68
2	30	2	0	1159	0.00000	79.11	2.65	61
3	30	3	0	1123	0.00000	78.50	2.65	60
4	30	4	0	1060	0.00000	66.28	2.65	51
5	30	5	0	1075	0.00000	63.15	2.65	49
6	30	6	0	1099	0.00000	59.21	2.65	46
7	30	7	0	1080	0.00000	42.12	2.65	32
8	30	8	0	1144	0.00000	41.92	2.65	32
"""


data_NoDend = """
N	Cell	syn#	nout	nin	Eff	ASA	SDRatio	nsites
0	2	0	962	1128	0.85284	174.01	3.23	134
1	2	1	811	1087	0.74609	124.89	3.23	96
2	2	2	789	1131	0.69761	108.53	3.23	83
3	2	3	742	1123	0.66073	101.16	3.23	78
4	2	4	546	1061	0.51461	80.73	3.23	62
5	2	5	20	1130	0.01770	45.89	3.23	35
0	5	0	693	1104	0.62772	141.02	2.80	108
1	5	1	274	1112	0.24640	106.09	2.80	82
2	5	2	189	1120	0.16875	99.89	2.80	77
3	5	3	89	1094	0.08135	93.86	2.80	72
4	5	4	98	1108	0.08845	93.52	2.80	72
5	5	5	47	1119	0.04200	89.85	2.80	69
6	5	6	0	1096	0.00000	53.78	2.80	41
0	6	0	852	1128	0.75532	159.01	2.82	122
1	6	1	763	1087	0.70193	140.81	2.82	108
2	6	2	765	1131	0.67639	132.17	2.82	102
3	6	3	586	1123	0.52182	109.34	2.82	84
4	6	4	457	1061	0.43073	97.97	2.82	75
5	6	5	387	1130	0.34248	94.22	2.82	72
0	9	0	1013	1113	0.91015	245.76	2.52	189
1	9	1	947	1083	0.87442	132.27	2.52	102
2	9	2	922	1142	0.80736	107.92	2.52	83
3	9	3	809	1145	0.70655	74.43	2.52	57
4	9	4	783	1100	0.71182	71.50	2.52	55
5	9	5	701	1080	0.64907	62.98	2.52	48
6	9	6	606	1100	0.55091	54.81	2.52	42
7	9	7	222	1090	0.20367	35.03	2.52	27
0	10	0	717	1110	0.64595	121.57	2.92	93
1	10	1	675	1111	0.60756	112.88	2.92	87
2	10	2	133	1095	0.12146	71.71	2.92	55
3	10	3	147	1083	0.13573	71.30	2.92	55
4	10	4	141	1101	0.12807	69.63	2.92	54
5	10	5	91	1086	0.08379	66.74	2.92	51
6	10	6	2	1110	0.00180	51.80	2.92	40
7	10	7	2	1075	0.00186	46.39	2.92	36
8	10	8	0	1144	0.00000	42.08	2.92	32
9	10	9	0	1125	0.00000	40.09	2.92	31
0	11	0	1007	1104	0.91214	213.09	2.44	164
1	11	1	1005	1112	0.90378	156.38	2.44	120
2	11	2	980	1120	0.87500	134.46	2.44	103
3	11	3	860	1094	0.78611	89.76	2.44	69
4	11	4	800	1108	0.72202	77.04	2.44	59
5	11	5	267	1119	0.23861	36.76	2.44	28
6	11	6	229	1096	0.20894	35.67	2.44	27
0	13	0	848	1104	0.76812	146.57	2.50	113
1	13	1	543	1112	0.48831	94.11	2.50	72
2	13	2	502	1120	0.44821	91.38	2.50	70
3	13	3	492	1094	0.44973	90.88	2.50	70
4	13	4	259	1108	0.23375	77.02	2.50	59
5	13	5	32	1119	0.02860	57.42	2.50	44
6	13	6	52	1096	0.04745	56.98	2.50	44
0	17	0	996	1104	0.90217	278.32	2.73	214
1	17	1	997	1112	0.89658	261.49	2.73	201
2	17	2	810	1120	0.72321	105.06	2.73	81
3	17	3	341	1094	0.31170	62.36	2.73	48
4	17	4	19	1108	0.01715	41.04	2.73	32
5	17	5	10	1119	0.00894	38.19	2.73	29
6	17	6	7	1096	0.00639	36.75	2.73	28
0	18	0	1000	1113	0.89847	273.22	3.01	210
1	18	1	967	1083	0.89289	164.19	3.01	126
2	18	2	885	1142	0.77496	103.85	3.01	80
3	18	3	652	1145	0.56943	70.54	3.01	54
4	18	4	499	1100	0.45364	60.30	3.01	46
5	18	5	396	1080	0.36667	55.66	3.01	43
6	18	6	175	1100	0.15909	45.18	3.01	35
7	18	7	55	1090	0.05046	37.12	3.01	29
0	30	0	1001	1146	0.87347	172.14	2.65	132
1	30	1	790	1074	0.73557	88.27	2.65	68
2	30	2	773	1159	0.66695	79.11	2.65	61
3	30	3	754	1123	0.67142	78.50	2.65	60
4	30	4	620	1060	0.58491	66.28	2.65	51
5	30	5	607	1075	0.56465	63.15	2.65	49
6	30	6	532	1099	0.48408	59.21	2.65	46
7	30	7	142	1080	0.13148	42.12	2.65	32
8	30	8	158	1144	0.13811	41.92	2.65	32

"""


def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def boltz_resid(p, x, data):
    return p["A"] * (1.0 / (1.0 + np.exp(-(x - p["vhalf"]) / p["k"]))) - data


class EfficacyPlots(object):
    def __init__(self, parent_figure:object=None, draft=False):
        self.parent_figure = parent_figure
        self.draft = draft

    def plot_efficacy(self, datasetname:str, ax:object, loc:tuple=(0., 0., 0., 0.)):
        self.titles = ["Intact", "Intact", "No Dendrites", "No Dendrites"]
        dmap = {'Full': data_Full, "NoDend": data_NoDend}
        dataset = dmap[datasetname]
        if loc is None:
            x0 = 0
            y0 = 0
        else:
            x0 = loc[0]
            y0 = loc[2]
        if self.parent_figure is None:
            sizer = {
                "A": {
                    "pos": [1+x0, 3, 1+y0, 3],
                    "labelpos": (-0.15, 1.02),
                    "noaxes": False,
                    },
                }
            self.P = PH.arbitrary_grid(
                sizer,
                order="columnsfirst",
                units="in",
                figsize=(5, 5),
                label=True,
                # verticalspacing=0.12,
                # horizontalspacing=0.12,
                # margins={
                #     "bottommargin": 0.1,
                #     "leftmargin": 0.1,
                #     "rightmargin": 0.1,
                #     "topmargin": 0.1,
                # },
                # parent=self.parent_plot,
            )
        else:
            self.P = self.parent_figure
        self.plot_dataset(dataset,
                ax=ax) # , title=self.titles[0]) # , legend=legend_on)
        return
        # axn = [1, 3]
        # for i, data in enumerate(datasets):
        #     if i == 1:
        #         legend_on = True
        #     else:
        #         legend_on = False
        #         print(i, axn)
        #     self.plot_dataset(data, plotno=i,
        #     ax=self.P.axdict[self.axes_labels[axn[i]]], title=self.titles[axn[i]], legend=legend_on)
        # return self.P

    def plot_ASA_SD(self):
        spc = re.compile("[ ;,\t\f\v]+")  # format replacing all spaces with tabs
        dataiter = re.finditer(spc, data)
        data = re.sub(spc, ",", data)

        sio = io.StringIO(data)
        df = pd.read_table(sio, sep=",")
        print(df.head())
        f, ax = mpl.subplots(1, 1)
        ax.margins(x=1, y=1)
        ax = sns.scatterplot(
            x="ASA",
            y="SDRatio",
            data=df,
            hue=[f"VCN_c{c:02d}" for c in df.Cell],
            sizes=(200, 200),
            clip_on=False,
        )
        ax.set_xlim(0, 300)
        ax.set_ylim(2, 3.05)

        ax.legend(
            loc="lower right",
            bbox_to_anchor=(1.0, 0.0),
            ncol=2,
            fontsize=9,
            markerscale=0.5,
            fancybox=False,
            shadow=False,
            facecolor="w",
            labelspacing=0.2,
        )
        mpl.show()

    def plot_dataset(self, data, ax: object = None, title: str = None, legend:bool=True):
        spc = re.compile("[ ;,\t\f\v]+")  # format replacing all spaces with tabs
        dataiter = re.finditer(spc, data)
        data = re.sub(spc, ",", data)

        sio = io.StringIO(data)
        df = pd.read_table(sio, sep=",")
        cell_names=[f"VCN_c{c:02d}" for c in df.Cell]
        sns.scatterplot(
            x="ASA",
            y="Eff",
            data=df,
            ax=ax,
            hue=cell_names,
            size=cell_names,
            sizes=(40,40),
            legend="full",
        )
        label = [f"VCN\_c{c:02d}" for c in df.Cell]

        if legend:
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
                ncol=1,
                fontsize=8,
                markerscale=0.5,
                fancybox=False,
                shadow=False,
                facecolor="w",
                labelspacing=0.2,
            )
        ax.set_xlim(0, 350.)
        gmodel = Model(boltz)
        gmodel.set_param_hint("A", value=1, min=0.0, max=1.0, vary=True)
        gmodel.set_param_hint("vhalf", value=140.0, min=10.0, max=300.0)
        gmodel.set_param_hint("k", value=20., min=0.01, max=200.0, vary=True)
        gparams = gmodel.make_params()
        # print("gparams: ", gparams.pretty_print())

        weights = np.ones(len(df.ASA))

        if self.draft:  # include brute force fit.
            result_brute = gmodel.fit(
            df.Eff, method="brute", params=gparams, x=df.ASA, weights=weights,
            )
            for p in result_brute.params:
                result_brute.params[p].stderr = 0.0  # abs(res2.params[p].value * 0)

            print("\nBrute fit: ")
            print(result_brute.params.pretty_print())
            print("-" * 80)
            xfit = np.linspace(0.0, 300.0, 300)
            bfit = gmodel.eval(params=result_brute.params, x=xfit)
            ax.plot(xfit, bfit, "b-", label="Brute")
            
        methods = ["leastsq"]
        # methods = ['leastsq', 'least_squares', 'differential_evolution', 'basin_hopping', 'ampgo', 'nelder', 'lbfgsb', 'powell', 'cg',
        #      'cobyla', 'bfgs', #'trust-exact', 'trust-krylov', 'trust-constr',
        #      #'dogleg',
        #      'slsqp', 'shgo', 'dual_annealing']
        # need jacobian: 'newton','tnc''trust-ncg'
        resdict = {key: None for key in methods}
        ix = np.argsort(df.Eff)
        y = df.ASA[ix]
        x = df.Eff[ix]
        weights = np.ones(len(df.Eff))  # * np.random. *(df.Eff/np.max(df.Eff))**2
        for meth in resdict.keys():
            result_LM = gmodel.fit(
                x,
                method=meth,
                params=gparams,  # gparams,
                x=y,
                weights=weights,
                max_nfev=20000,  # fit_kws={'maxfev': 5000}
                jac="3-point",
            )
            print(f"\n{meth:s} fit: ")
            print(result_LM.params.pretty_print())
            print("-" * 80)
            xfit = np.linspace(0.0, 300.0, 300)
            lev_fit = gmodel.eval(params=result_LM.params, x=xfit)
            # lmfit.report_fit(result_brute.params, min_correl=0.5)

            for p in result_LM.params:
                if result_LM.params[p].stderr == None:
                    result_LM.params[p].stderr = 0.0  # abs(res2.params[p].value * 0)
            resdict[meth] = result_LM
        # ci, trace = lmfit.conf_interval(mini, res2, sigmas=[2, 2, 2], trace=True)
        # lmfit.printfuncs.report_ci(ci)
        ax.text(
            1.1, 0.8,
            f"Max Efficacy: {result_LM.params['A'].value:.2f} (1$\sigma$ = {result_LM.params['A'].stderr:.2f})",
            fontsize=7,
            verticalalignment="top", horizontalalignment="right",
            transform=ax.transAxes,
        )
        uni = r"$\mu m^2$"
        asalab = r"$ASA_{0.5}$"
        ax.text(
            1.1,
            0.4,
            f"{asalab:s}: {result_LM.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result_LM.params['vhalf'].stderr:.1f})\n",
            fontsize=7,
            verticalalignment="center", horizontalalignment="right",
            transform=ax.transAxes,
        )
        ax.text(
            1.1,
            0.35,
            f"k: {result_LM.params['k'].value:.1f} (1$\sigma$ = {result_LM.params['k'].stderr:.1f})",
            fontsize=7,
            verticalalignment="center", horizontalalignment="right",
            transform=ax.transAxes,
        )
        # return
        #
        # if self.draft:
        #     ax.text(
        #         0.0,
        #         1.0,
        #         f"Brute: Max Eff: {result_brute.params['A'].value:.2f} (1$\sigma$ = {result_brute.params['A'].stderr:.2f})",
        #         fontsize=6,
        #         verticalalignment="top",
        #     )
        #     ax.text(
        #         0.0,
        #         0.9,
        #         f"Br ASA at half max: {result_brute.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result_brute.params['vhalf'].stderr:.1f})\n"
        #         + f"Slope: {result_brute.params['k'].value:.1f} (1$\sigma$ = {result_brute.params['k'].stderr:.1f})",
        #         fontsize=6,
        #     )
        #     mpl.text(
        #         0.0,
        #         0.8,
        #         f"Br $\chi^2$: {result_brute.chisqr:.3f}\nLM $\chi^2$:{result_LM.chisqr:.3f}",
        #         fontsize=6,
            # )

        for m in resdict.keys():
            fit = gmodel.eval(params=resdict[m].params, x=xfit)
            ax.plot(xfit, fit, "k-", label=meth)
        ax.set_xlabel(f"ASA ({uni:s})")
        ax.set_ylabel("Efficacy (Bushy spikes/input spikes)")
        ax.set_ylim(0, 1.0)
        ax.set_title(title, fontsize=12, fontweight="bold")
        PH.nice_plot(ax, direction='outward', ticklength=2.)
        PH.talbotTicks(ax, density=(1.0, 1.0), insideMargin=0, tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1},
                   axrange={'x': (0, 300.), 'y':(0, 1)}, pointSize=None)
        

if __name__ == "__main__":
    sizer = {
        "A": {
            "pos": [0.5, 3, 0.5, 3],
            "labelpos": (-0.15, 1.02),
            "noaxes": True,
        },
        "B": {"pos": [4.0, 3.0, 0.5, 3.0], "labelpos": (-0.15, 1.02)},
    }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) f
    P = PH.arbitrary_grid(
        sizer,
        order="columnsfirst",
        units='in',
        figsize=(8, 4.5),
        label=True,
    )
    EFP = EfficacyPlots(parent_plot=P)
    # EFP.plot_data("B")
    EFP.plot_efficacy(dataset="Full", ax=P.axdict["B"])
    mpl.show()