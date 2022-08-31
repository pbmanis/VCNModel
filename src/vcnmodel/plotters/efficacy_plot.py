import io
import re
from typing import Union

import lmfit
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit import Model
import matplotlib.pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as STY

matplotlib.rcParams["mathtext.fontset"] = "stixsans"

"""
Plot the efficacy data into an axis

data_original is an old efficacy measurement data set
plots appeared in some posters. In the original table, cells 2 and 5 have their original reconstructed axons. From runs
done in 2020.

This is superseeded by the dataFull and dataNoDend tables below, from later sets of runs
where cells 2 and 5 have their "standardized" axons (mean axons/ais/ah) as the axon reconstructions are incomplete.
All data are from 30 dB SPL stimuli

The second set of runs were done on 8/2/2021, runANSingles, Full/NoDend, 1 sec, 5 reps, no depression 
The table for the efficacy is printed out when the dataset is selected in DataTablesVCN,
and the 'Singles' analysis is run. THe analysis displays the spikes, as well as printing the table in a format
suitable for pasting into this file
"""

# data_original = """
# ASA    Eff    Cell    SDRatio
# 174.01    0.3984375    2 3.2324
# 124.89    0    2        3.2324
# 108.53    0    2        3.2324
# 101.16    0    2        3.2324
# 80.73    0    2        3.2324
# 45.89    0    2        3.2324
# 141.02    0    5   2.7951
# 106.09    0    5   2.7951
# 99.89    0    5   2.7951
# 93.86    0    5   2.7951
# 93.52    0    5   2.7951
# 89.85    0    5   2.7951
# 53.78    0    5   2.7951
# 159.01    0.51562500000000    6   2.8210
# 140.81    0.27200000000000    6   2.8210
# 132.17    0.19491525423729    6   2.8210
# 109.34    0    6   2.8210
# 97.97    0.0078125    6   2.8210
# 94.22    0    6   2.8210
# 245.76    0.817460317460317    9   2.5218
# 132.27    0.548148148148148    9   2.5218
# 107.92    0.2578125    9   2.5218
# 74.43    0    9   2.5218
# 71.5    0    9   2.5218
# 62.98    0    9   2.5218
# 54.81    0    9   2.5218
# 35.03    0    9   2.5218
# 121.57    0    10
# 112.88    0.00719424460431655    10  2.9240
# 71.71    0    10  2.9240
# 71.3    0    10  2.9240
# 69.63    0    10  2.9240
# 66.74    0    10  2.9240
# 51.8    0    10  2.9240
# 46.39    0    10  2.9240
# 42.08    0    10  2.9240
# 40.09    0    10  2.9240
# 213.09    0.748091603053435    11  2.4353
# 156.383    0.725352112676056    11  2.4353
# 134.46    0.60431654676259    11  2.4353
# 89.76    0.136690647482014    11  2.4353
# 77.04    0.0230769230769231    11  2.4353
# 36.76    0    11  2.4353
# 35.67    0    11  2.4353
# 146.57    0.641221374045801    13  2.4999
# 94.11    0.0704225352112676    13  2.4999
# 91.38    0.0503597122302158    13  2.4999
# 90.88    0.0287769784172662    13  2.4999
# 77.02    0    13  2.4999
# 57.42    0    13  2.4999
# 56.98    0    13  2.4999
# 278.32    0.755725190839695    17  3.0123
# 261.49    0.788732394366197    17  3.0123
# 105.06    0.00719424460431655    17  3.0123
# 62.36    0    17  3.0123
# 41.04    0    17  3.0123
# 38.19    0    17  3.0123
# 36.75    0    17  3.0123
# 172.14    0.553030303030303    30  2.6451
# 88.27    0    30  2.6451
# 79.11    0    30  2.6451
# 78.5    0    30  2.6451
# 66.28    0    30  2.6451
# 63.15    0    30  2.6451
# 59.21    0    30  2.6451
# 42.12    0    30  2.6451
# 41.92    0    30  2.6451
# """


data_Full = """
N	Cell	syn#	nout	nin	Eff	ASA	SDRatio	nsites
0	2	0	85	1074	0.07914	174.01	3.23	134
1	2	1	1	1027	0.00097	124.89	3.23	96
2	2	2	0	1083	0.00000	108.53	3.23	83
3	2	3	0	1070	0.00000	101.16	3.23	78
4	2	4	0	999	    0.00000	80.73	3.23	62
5	2	5	0	1071	0.00000	45.89	3.23	35
0	5	0	15	1053	0.01425	141.02	2.80	108
1	5	1	0	1048	0.00000	106.09	2.80	82
2	5	2	0	1071	0.00000	99.89	2.80	77
3	5	3	0	1041	0.00000	93.86	2.80	72
4	5	4	0	1052	0.00000	93.52	2.80	72
5	5	5	0	1064	0.00000	89.85	2.80	69
6	5	6	0	1045	0.00000	53.78	2.80	41
0	6	0	80	1074	0.07449	159.01	2.82	122
1	6	1	30	1027	0.02921	140.81	2.82	108
2	6	2	14	1083	0.01293	132.17	2.82	102
3	6	3	2	1070	0.00187	109.34	2.82	84
4	6	4	1	999	0.00100	97.97	2.82	75
5	6	5	0	1071	0.00000	94.22	2.82	72
0	9	0	730	1059	0.68933	245.76	2.52	189
1	9	1	141	1034	0.13636	132.27	2.52	102
2	9	2	23	1090	0.02110	107.92	2.52	83
3	9	3	1	1083	0.00092	74.43	2.52	57
4	9	4	0	1038	0.00000	71.50	2.52	55
5	9	5	0	1022	0.00000	62.98	2.52	48
6	9	6	0	1048	0.00000	54.81	2.52	42
7	9	7	0	1039	0.00000	35.03	2.52	27
0	10	0	1	1062	0.00094	121.57	2.92	93
1	10	1	2	1062	0.00188	112.88	2.92	87
2	10	2	0	1036	0.00000	71.71	2.92	55
3	10	3	0	1031	0.00000	71.30	2.92	55
4	10	4	0	1056	0.00000	69.63	2.92	54
5	10	5	0	1028	0.00000	66.74	2.92	51
6	10	6	0	1049	0.00000	51.80	2.92	40
7	10	7	0	1025	0.00000	46.39	2.92	36
8	10	8	0	1088	0.00000	42.08	2.92	32
9	10	9	0	1069	0.00000	40.09	2.92	31
0	11	0	712	1053	0.67616	213.09	2.44	164
1	11	1	472	1048	0.45038	156.38	2.44	120
2	11	2	277	1071	0.25864	134.46	2.44	103
3	11	3	14	1041	0.01345	89.76	2.44	69
4	11	4	2	1052	0.00190	77.04	2.44	59
5	11	5	0	1064	0.00000	36.76	2.44	28
6	11	6	0	1045	0.00000	35.67	2.44	27
0	13	0	338	1053	0.32099	146.57	2.50	113
1	13	1	8	1048	0.00763	94.11	2.50	72
2	13	2	4	1071	0.00373	91.38	2.50	70
3	13	3	5	1041	0.00480	90.88	2.50	70
4	13	4	0	1052	0.00000	77.02	2.50	59
5	13	5	0	1064	0.00000	57.42	2.50	44
6	13	6	0	1045	0.00000	56.98	2.50	44
0	17	0	796	1053	0.75594	278.32	2.73	214
1	17	1	771	1048	0.73569	261.49	2.73	201
2	17	2	20	1071	0.01867	105.06	2.73	81
3	17	3	0	1041	0.00000	62.36	2.73	48
4	17	4	0	1052	0.00000	41.04	2.73	32
5	17	5	0	1064	0.00000	38.19	2.73	29
6	17	6	0	1045	0.00000	36.75	2.73	28
0	18	0	717	1059	0.67705	273.22	3.01	210
1	18	1	125	1034	0.12089	164.19	3.01	126
2	18	2	0	1090	0.00000	103.85	3.01	80
3	18	3	0	1083	0.00000	70.54	3.01	54
4	18	4	0	1038	0.00000	60.30	3.01	46
5	18	5	0	1022	0.00000	55.66	3.01	43
6	18	6	0	1048	0.00000	45.18	3.01	35
7	18	7	0	1039	0.00000	37.12	3.01	29
0	30	0	148	1091	0.13566	172.14	2.65	132
1	30	1	0	1012	0.00000	88.27	2.65	68
2	30	2	0	1102	0.00000	79.11	2.65	61
3	30	3	0	1071	0.00000	78.50	2.65	60
4	30	4	0	1006	0.00000	66.28	2.65	51
5	30	5	0	1026	0.00000	63.15	2.65	49
6	30	6	0	1048	0.00000	59.21	2.65	46
7	30	7	0	1031	0.00000	42.12	2.65	32
8	30	8	0	1089	0.00000	41.92	2.65	32
"""


data_NoDend = """
N	Cell	syn#	nout	nin	Eff	ASA	SDRatio	nsites
0	2	0	918	1074	0.85475	174.01	3.23	134
1	2	1	794	1027	0.77313	124.89	3.23	96
2	2	2	773	1083	0.71376	108.53	3.23	83
3	2	3	730	1070	0.68224	101.16	3.23	78
4	2	4	559	999	0.55956	80.73	3.23	62
5	2	5	55	1071	0.05135	45.89	3.23	35
0	5	0	922	1053	0.87559	141.02	2.80	108
1	5	1	844	1048	0.80534	106.09	2.80	82
2	5	2	839	1071	0.78338	99.89	2.80	77
3	5	3	801	1041	0.76945	93.86	2.80	72
4	5	4	795	1052	0.75570	93.52	2.80	72
5	5	5	798	1064	0.75000	89.85	2.80	69
6	5	6	425	1045	0.40670	53.78	2.80	41
0	6	0	807	1074	0.75140	159.01	2.82	122
1	6	1	721	1027	0.70204	140.81	2.82	108
2	6	2	720	1083	0.66482	132.17	2.82	102
3	6	3	537	1070	0.50187	109.34	2.82	84
4	6	4	401	999	0.40140	97.97	2.82	75
5	6	5	354	1071	0.33053	94.22	2.82	72
0	9	0	963	1059	0.90935	245.76	2.52	189
1	9	1	890	1034	0.86074	132.27	2.52	102
2	9	2	863	1090	0.79174	107.92	2.52	83
3	9	3	756	1083	0.69806	74.43	2.52	57
4	9	4	725	1038	0.69846	71.50	2.52	55
5	9	5	652	1022	0.63796	62.98	2.52	48
6	9	6	566	1048	0.54008	54.81	2.52	42
7	9	7	175	1039	0.16843	35.03	2.52	27
0	10	0	790	1062	0.74388	121.57	2.92	93
1	10	1	761	1062	0.71657	112.88	2.92	87
2	10	2	358	1036	0.34556	71.71	2.92	55
3	10	3	358	1031	0.34724	71.30	2.92	55
4	10	4	345	1056	0.32670	69.63	2.92	54
5	10	5	284	1028	0.27626	66.74	2.92	51
6	10	6	55	1049	0.05243	51.80	2.92	40
7	10	7	34	1025	0.03317	46.39	2.92	36
8	10	8	9	1088	0.00827	42.08	2.92	32
9	10	9	8	1069	0.00748	40.09	2.92	31
0	11	0	950	1053	0.90218	213.09	2.44	164
1	11	1	945	1048	0.90172	156.38	2.44	120
2	11	2	929	1071	0.86741	134.46	2.44	103
3	11	3	818	1041	0.78578	89.76	2.44	69
4	11	4	761	1052	0.72338	77.04	2.44	59
5	11	5	217	1064	0.20395	36.76	2.44	28
6	11	6	179	1045	0.17129	35.67	2.44	27
0	13	0	806	1053	0.76543	146.57	2.50	113
1	13	1	469	1048	0.44752	94.11	2.50	72
2	13	2	461	1071	0.43044	91.38	2.50	70
3	13	3	442	1041	0.42459	90.88	2.50	70
4	13	4	234	1052	0.22243	77.02	2.50	59
5	13	5	18	1064	0.01692	57.42	2.50	44
6	13	6	27	1045	0.02584	56.98	2.50	44
0	17	0	940	1053	0.89269	278.32	2.73	214
1	17	1	950	1048	0.90649	261.49	2.73	201
2	17	2	783	1071	0.73109	105.06	2.73	81
3	17	3	287	1041	0.27570	62.36	2.73	48
4	17	4	18	1052	0.01711	41.04	2.73	32
5	17	5	11	1064	0.01034	38.19	2.73	29
6	17	6	3	1045	0.00287	36.75	2.73	28
0	18	0	953	1059	0.89991	273.22	3.01	210
1	18	1	923	1034	0.89265	164.19	3.01	126
2	18	2	834	1090	0.76514	103.85	3.01	80
3	18	3	637	1083	0.58818	70.54	3.01	54
4	18	4	481	1038	0.46339	60.30	3.01	46
5	18	5	418	1022	0.40900	55.66	3.01	43
6	18	6	206	1048	0.19656	45.18	3.01	35
7	18	7	58	1039	0.05582	37.12	3.01	29
0	30	0	953	1091	0.87351	172.14	2.65	132
1	30	1	743	1012	0.73419	88.27	2.65	68
2	30	2	709	1102	0.64338	79.11	2.65	61
3	30	3	692	1071	0.64613	78.50	2.65	60
4	30	4	567	1006	0.56362	66.28	2.65	51
5	30	5	530	1026	0.51657	63.15	2.65	49
6	30	6	474	1048	0.45229	59.21	2.65	46
7	30	7	135	1031	0.13094	42.12	2.65	32
8	30	8	123	1089	0.11295	41.92	2.65	32

"""


def boltz(x, A, vhalf, k):
    return A * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def boltz_resid(p, x, data):
    return p["A"] * (1.0 / (1.0 + np.exp(-(x - p["vhalf"]) / p["k"]))) - data


class EfficacyPlots(object):
    def __init__(self, parent_figure:object=None, draft=False):
        self.parent_figure = parent_figure
        self.draft = draft

    def plot_efficacy(self, datasetname:str, ax:object, loc:tuple=(0., 0., 0., 0.), figuremode="full"):
        self.figuremode = figuremode
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
            )
        else:
            self.P = self.parent_figure
        self.plot_dataset(dataset, ax=ax)
        return

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

    def fit_dataset(self, df_full, sel_cells = None, ax=None):
        gmodel = Model(boltz)
        if ic is None:
            A_ic = 1
            Vh_ic = 140.0
            k_ic = 20
        else:
            A_ic = ic["A"]
            Vh_ic = ic["Vh"]
            k_ic = ic["k"]
        gmodel.set_param_hint("A", value=A_ic, min=0.0, max=1.0, vary=True)
        gmodel.set_param_hint("vhalf", value=Vh_ic, min=10.0, max=300.0)
        gmodel.set_param_hint("k", value=k_ic, min=0.01, max=200.0, vary=True)
        gparams = gmodel.make_params()
        # print("gparams: ", gparams.pretty_print())
        if sel_cells is not None:
            df = df_full[df_full["Cell"].isin(sel_cells)]

        weights = np.ones(len(df.ASA))

        # if self.draft:  # include brute force fit.
        #     result_brute = gmodel.fit(
        #     df.Eff, method="brute", params=gparams, x=df.ASA, weights=weights,
        #     )
        #     for p in result_brute.params:
        #         result_brute.params[p].stderr = 0.0  # abs(res2.params[p].value * 0)
        #
        #     print("\nBrute fit: ")
        #     print(result_brute.params.pretty_print())
        #     print("-" * 80)
        #     xfit = np.linspace(0.0, 300.0, 300)
        #     bfit = gmodel.eval(params=result_brute.params, x=xfit)
        #     ax.plot(xfit, bfit, "b-", label="Brute")
        # otherwise we can use another method...
        # methods = ['leastsq', 'least_squares', 'differential_evolution', 'basin_hopping', 'ampgo', 'nelder', 'lbfgsb', 'powell', 'cg',
        #      'cobyla', 'bfgs', #'trust-exact', 'trust-krylov', 'trust-constr',
        #      #'dogleg',
        #      'slsqp', 'shgo', 'dual_annealing']
        # need jacobian: 'newton','tnc''trust-ncg'
        # ix = np.argsort(df.Eff)
        y = df.ASA # [ix]
        x = df.Eff #[ix]
        weights = np.ones(len(df.Eff))  # * np.random. *(df.Eff/np.max(df.Eff))**2
        for meth in resdict.keys():
            result_LM = gmodel.fit(
                x,
                method=meth,
                params=gparams,  # gparams,
                x=y,
                weights=weights,
                # maxfev=20000,  # fit_kws={'maxfev': 5000}
                # jac="3-point",
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
        return result_LM, resdict, gmodel, xfit

    def plot_dataset(self, data, ax: object = None, title: str = None, legend:bool=True):
        spc = re.compile("[ ;,\t\f\v]+")  # format replacing all spaces with tabs
        dataiter = re.finditer(spc, data)
        data = re.sub(spc, ",", data)

        sio = io.StringIO(data)
        df = pd.read_table(sio, sep=",")
        cell_names=[f"VCN\_c{c:02d}" for c in df.Cell]
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
        result_LM1, resdict1, gmodel1, xfit2 = self.fit_dataset(df, sel_cells = [5, 30, 17])
        result_LM2, resdict2, gmodel2, xfit1 = self.fit_dataset(df, sel_cells = [9,11, 13])
        result_LM, resdict, gmodel, xfit = self.fit_dataset(df, sel_cells  = None)
        
        # ci, trace = lmfit.conf_interval(mini, res2, sigmas=[2, 2, 2], trace=True)
        # lmfit.printfuncs.report_ci(ci)
        return result_LM, xfit, yfit

    def plot_fit(self, df, ax, y0, method='leastsq', cells = [],
                 color='k-', ic=None):
        result_LM, xfit, yfit = self.fit_dataset(df, method, ic=ic)
        if len(cells) > 0:
            tc = str(cells)
            lab = f"{method:s}  {tc:s}"
        else:
            lab = method
        ax.plot(xfit, yfit, color, label=lab)  # all cells

        if self.figuremode == "full":
            ax.text(
                1.1, 0.8+y0,
                f"Max Efficacy: {result_LM.params['A'].value:.2f} (1$\sigma$ = {result_LM.params['A'].stderr:.2f})",
                fontsize=7,
                verticalalignment="top", horizontalalignment="right",
                transform=ax.transAxes,
            )
            uni = r"$\mu m^2$"
            asalab = r"$ASA_{0.5}$"
            ax.text(
                1.1,
                0.4+y0,
                f"{asalab:s}: {result_LM.params['vhalf'].value:.1f} {uni:s} (1$\sigma$ = {result_LM.params['vhalf'].stderr:.1f})\n",
                fontsize=7,
                color=color[0],
                verticalalignment="center", horizontalalignment="right",
                transform=ax.transAxes,
            )
            ax.text(
                1.1,
                0.2+y0,
                f"k: {result_LM.params['k'].value:.1f} (1$\sigma$ = {result_LM.params['k'].stderr:.1f})",
                fontsize=7,
                color=color[0],
                verticalalignment="center", horizontalalignment="right",
                transform=ax.transAxes,
            )

    def plot_dataset(self, data, ax: object = None, title: str = None, legend:bool=True):
        spc = re.compile("[ ;,\t\f\v]+")  # format replacing all spaces with tabs
        dataiter = re.finditer(spc, data)
        data = re.sub(spc, ",", data)

        sio = io.StringIO(data)
        df = pd.read_table(sio, sep=",")
        cell_names=[f"VCN\_c{c:02d}" for c in df.Cell]
        uni = r"$\mu m^2$"

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
            ax.plot(xfit, fit, "k-", label=m)
        ax.set_xlabel(f"ASA ({uni:s})")
        label = [f"VCN\_c{c:02d}" for c in df.Cell]
        method = "leastsq"

        cells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
        self.plot_fit(df, ax, y0=0, method=method, cells = cells, color='k-')

        cells1 = [9,11,13]
        df1 = df[df.Cell.isin(cells1)]
        self.plot_fit(df1, ax, y0=0.05, method=method, cells = cells1, color='r-')

        cells2 = [2,5,6,10,30]
        df2 = df[df.Cell.isin(cells2)]
        self.plot_fit(df2, ax, y0=0.1, method=method, cells = cells2, color='b-', ic={"A": 1.0, "Vh": 200., "k": 10})

        ax.set_xlim(0, 350.)
        ax.set_xlabel(f"Input ASA ({uni:s})")

        ax.set_ylabel("Efficacy (Bushy spikes/input spikes)")
        ax.set_ylim(0, 1.0)
        ax.set_title(title, fontsize=12, fontweight="bold")
        PH.nice_plot(ax, direction='outward', ticklength=2.)
        PH.talbotTicks(ax, density=(1.0, 1.0), insideMargin=0, tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1},
                   axrange={'x': (0, 300.), 'y':(0, 1)}, pointSize=None)
        if legend:
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
                ncol=1,
                fontsize=6,
                markerscale=0.5,
                fancybox=False,
                shadow=False,
                facecolor="w",
                labelspacing=0.2,
            )

if __name__ == "__main__":
    sizer = {
        "A": {
            "pos": [0.5, 3, 0.5, 3],
            "labelpos": (-0.15, 1.02),
            "noaxes": True,
        },
        "B": {"pos": [4.0, 3.0, 0.5, 3.0], "labelpos": (-0.15, 1.02)},
    }
    P = PH.arbitrary_grid(
        sizer,
        order="columnsfirst",
        units='in',
        figsize=(8, 4.5),
        label=True,
    )
    EFP = EfficacyPlots(parent_figure=P)
    # EFP.plot_data("B")
    EFP.plot_efficacy("Full", ax=P.axdict["B"])
    mpl.show()