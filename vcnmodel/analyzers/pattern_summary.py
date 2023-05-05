""" pattern_summary:
Compute the input patterns, and generate a summary figure
breaking the cells down by 'drivers' and 'coincidence'
for different input patterns.

This also has the routines for plotting panel 4F.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""

import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import matplotlib.pyplot as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import discriminant_analysis
from vcnmodel.util.get_data_paths import get_data_paths
from vcnmodel.util.get_data_paths import update_disk
import vcnmodel.cell_config as cell_config
import vcnmodel.plotters.figure_data as FD
import vcnmodel.plotters.plot_sims as PS
import vcnmodel.util.readmodel as readmodel
from matplotlib.lines import Line2D
from vcnmodel.util.set_figure_path import set_figure_path
from pylibrary.tools import cprint as CP
from pylibrary.plotting import plothelpers as PH
from pyrsistent import PSet
from vcnmodel.analyzers import reverse_correlation as REVCORR
import vcnmodel.group_defs as GRPDEF
from vcnmodel.analyzers.analyzer_data_classes import PData
from vcnmodel.plotters import \
    figure_data as FD  # table of simulation runs used for plotting figures

cprint = CP.cprint

PSC = PS.PlotSims(parent=None)

# @dataclass
# class PData:
#     """
#     data class for some parameters that control what we read
#     """

#     gradeA: list = field(default_factory=GRPDEF.grAList)
#     default_modelName: str = "XM13_nacncoop"
#     soma_inflate: bool = True
#     dend_inflate: bool = True
#     basepath: str = ""  # config["baseDataDirectory"]
#     renderpath: str = ""  # " str(Path(self.config["codeDirectory"], "Renderings"))
#     revcorrpath: str = ""  # config["revcorrDataDirectory"]
#     thiscell: str = ""


def get_synaptic_info(gbc: str, add_inputs="none") -> tuple:
    SC = cell_config.CellConfig(add_inputs=add_inputs)
    # print(SC.VCN_Inputs)
    if isinstance(gbc, int):
        gbc_string = f"VCN_c{gbc:02d}"
    elif isinstance(gbc, str) and len(gbc) <= 2:
        gbc_string = f"VCN_c{int(gbc):02d}"
    elif isinstance(gbc, str) and gbc.startswith("VCN_"):
        gbc_string = gbc
    elif isinstance(gbc, str) and gbc.startswith("BC"):
        gbc_string = gbc.replace("BC", "VCN_c")
    syninfo = SC.VCN_Inputs[gbc_string]
    return (SC, syninfo)


def read_pkl_file(datafile) -> Union[tuple, None]:
    """
    Reads pickle file entry to get the overall info about the
    run .
    """
    try:
        with open(datafile, "rb") as fh:
            d = pickle.load(fh, fix_imports=True, encoding="latin1")
    except:
        import vcnmodel.util.fixpicklemodule as FPM

        try:
            with open(datafile, "rb") as fh:
                d = FPM.pickle_load(fh)  # , fix_imports=True) # , encoding="latin1")
        except (ModuleNotFoundError, IOError, pickle.UnpicklingError):
            return None
    return d  # just the info we need


def revcorr(gbc, filename, PD, protocol="runANPSTH", revcorrtype="RevcorrSPKS"):
    RM = readmodel.ReadModel()
    MD = RM.get_data(filename, PD=PD, protocol=protocol)
    if not MD.success:
        return None
    # (d, AR, SP, RM, RCP, RCD) = res
    PD.thiscell = gbc
    RCP = MD.RCP
    RCD = MD.RCD
    RCP.algorithm = revcorrtype
    SC, syninfo = get_synaptic_info(gbc)

    # print(f"    compute_revcorr: Preparing for computation for: {str(gbc):s}")
    RCP.ninputs = len(syninfo[1])
    # print("ninputs from syninfo: ", RCP.ninputs)
    RCP.ninputs  # save for trace viewer.
    RCD.sites = np.zeros(RCP.ninputs)
    for isite in range(RCP.ninputs):  # precompute areas
        area = syninfo[1][isite][0]
        if area > RCP.amax:
            RCP.amax = area
        RCD.sites[isite] = int(np.around(area * SC.synperum2))
    cprint("g", f"    compute_revcorr # of inputs: {RCP.ninputs:d}")

    """
        2. set up parameters
    """

    if isinstance(MD.SI, dict):
        min_time = (MD.SI["pip_start"] + 0.025) * 1000.0  # push onset out of the way
        max_time = (MD.SI["pip_start"] + MD.SI["pip_duration"]) * 1000.0
        F0 = MD.SI["F0"]
        dB = (MD.SI["dB"],)
    else:
        if isinstance(MD.RI.pip_start, list):
            start = MD.RI.pip_start[0]
        else:
            start = MD.RI.pip_start
        min_time = (start + 0.025) * 1000.0
        max_time = (start + MD.RI.pip_duration) * 1000.0
    expt = MD.RI.Spirou
    # print("    compute_revcorr experiment type: ", expt)

    RCD.sv_trials = []
    RCD.C = [None] * RCP.ninputs
    RCD.CB = [None] * RCP.ninputs
    RCD.STTC = [None] * RCP.ninputs
    RCD.max_coin_rate = 0.0
    RCD.npost_spikes = 0  # count spikes used in analysis
    RCD.mean_post_interval = 0.0
    RCD.mean_pre_intervals = np.zeros(RCP.ninputs)
    RCD.nsp_avg = 0
    RCP.min_time = min_time  # driven window without onset
    RCP.max_time = max_time
    RCD.pre_w = [-2.7, -0.5]  # sets the values for the pre window

    """
    3. Prepare storage arrays
    """
    nbins = int((max_time - min_time) / RCP.binw)

    RCD.CBT = np.arange(RCP.minwin, RCP.maxwin, RCP.binw)  # refcorr time base
    for isite in range(RCP.ninputs):
        RCD.CB[isite] = np.zeros_like(RCD.CBT)
        ncpts = len(np.arange(RCP.minwin, -RCP.minwin, RCP.binw))
        RCD.C[isite] = np.zeros(ncpts)
        RCD.STTC[isite] = np.zeros(ncpts)

    """
    4. Do the calculations.
    """
    syn_ASA = np.array([syninfo[1][isite][0] for isite in range(RCP.ninputs)])
    revcorr_params = {
        "deselect": False,
        "threshold": 400.0,
        "window": [-2.7, -0.5],
        "ASAs": syn_ASA,
    }
    MD2, filtered_stx_by_trial = REVCORR.revcorr(
        MD, nbins, revcorrtype, ASAs=syn_ASA, revcorr_params=revcorr_params
    )
    RCP, RCD, PAT = REVCORR.spike_pattern_analysis(
        MD
    )  # perfrom an analysis of input spike patters
    return RCP, RCD, PAT


def get_pattern_data(dataset:str="Spont"):

    datapaths = get_data_paths()
    
    fig_data = FD.figure_revcorr
    db = dataset
    group = {"MixedMode": GRPDEF.MixedMode, "Coincidence": GRPDEF.Coincidence}

    compares = pd.DataFrame(columns=["Cell", "PatternName", "Percent", "#Spikes", "Group"])
    for gbc in GRPDEF.grAList():
        BC = f"VCN_c{gbc:02d}"
        cprint("g", f"\nBushy cell: {BC:s}")
        fn = fig_data[gbc][db]
        filename = Path(
            datapaths["cellDataDirectory"],
            f"VCN_c{gbc:02d}",
            "Simulations",
            "AN",
            fn + ".pkl",
        )
        print("Getting compares from: ", filename)
        pklf = read_pkl_file(filename)
        # PSC.compute_revcorr(
        #     P=None,
        #     gbc=str(gbc),
        #     fn=pklf.files[0],
        #     PD = PData(),
        #     protocol="runANPSTH",
        #     revcorrtype="RevcorrSPKS",
        # )
        if gbc in group["MixedMode"]:
            gname = "Mixed Mode"
        else:
            gname = "Coincidence"
        # maybe disk has changed, so adjust the filename to read the data from

        pfs = update_disk(pklf.files[0])  # update disk if needed
        RCP, RCD, PAT = revcorr(gbc, pfs, PD=PData(gradeA=GRPDEF.gradeACells))
        for t in PAT.filttable.keys():
            pat = PAT.filttable[t]
            # pat.print()
            pct = 100.0 * pat.sumtrue / (pat.sumtrue + pat.sumfalse)
            Nspikes = pat.sumtrue + pat.sumfalse
            df = pd.DataFrame(
                [{"Cell": gbc, "PatternName": t, "Percent": pct, "#Spikes":Nspikes, "Group": gname}]
            )
            compares = pd.concat([compares, df], ignore_index=True)
    return compares




def Figure5_Supplemental2_Patterns(reanalyze=False, dataset:str="Spont"):
    assert dataset in ["Spont", "30dB"]
    compare_file = f"pattern_compares_{dataset:s}.pkl"
    if reanalyze:
        compares = get_pattern_data(dataset=dataset)
        compares.to_pickle(compare_file)

    compares = pd.read_pickle(compare_file)
    print("compares", compares.head())
    # print(compares.head(10))
    nrows = 4
    ncols = 5
    f, ax = mpl.subplots(nrows, ncols)
    f.set_size_inches(9, 7)
    ax = ax.ravel()
    pats = [
        "1st_largest_alone",
        "2nd_largest_alone",
        "3rd_largest_alone",
        "4th_largest_alone",
        "5th_largest_alone",
        # end of row 1
        "1st_largest+others",
        "2nd_largest+others",
        "3rd_largest+others",
        "4th_largest+others",
        "5th_largest+others",
        # end of row 2
        "1st+2nd_alone",
        "1st+2nd+others",
        "1st+2nd+3rd+others",
        "",
        "",
        # end of row 3
        "not_largest",
        "not_2_largest",
        "not_3_largest",
        "",
        "",
        # end of row 4
    ]
    pnames = sorted(set(pats))
    pdatanames = sorted(set(compares["PatternName"].values))
    # print(pnames)
    # print(pdatanames)

    for p in pats:
        if p == "":
            continue
        if p not in pdatanames:
            raise ValueError(f"Pattern {p:s} is missing from pdatanames")

    # create your own color array
    bar_colors = [(0.8, 1.0, 0.8, 0.05), (0.8, 0.8, 1.0, 0.05)]

    # add color array to set_palette
    # function of seaborn
    sns.set_palette(bar_colors, desat=0.2)
    ylims = [25, 75, 25, 75]
    yticks = [[0, 5, 10, 15, 20, 25],
               [0, 25, 50, 75],
               [0, 5, 10, 15, 20, 25],
               [0, 25, 50, 75],
               ]
    yticks_str = [[f"{y:2d}" for y in yticks[yi]] for yi in range(len(yticks))]
    compares.loc[:, "Group"] = compares.apply(lambda row: GRPDEF.remap_Group(row.Cell) , axis=1)
    for i, axl in enumerate(pats):
        irow = i // ncols
        icol = i % ncols
        # print(compares[compares["PatternName"] == pats[i]].head(25))
        if len(axl) > 0:
            axt = sns.boxplot(
                data=compares[compares["PatternName"] == pats[i]],
                x="Group",
                y="Percent",
                width=0.5,
                ax=ax[i],
            )
            # # adding transparency to colors
            # for patch in axt.artists:
            #     r, g, b, a = patch.get_facecolor()
            #     # patch.set_facecolor((r, g, b, 0.3))
            #     patch.set(alpha=None, facecolor=((r,g,b,0.3)))

            swp = sns.swarmplot(
                data=compares[compares["PatternName"] == pats[i]],
                x="Group",
                y="Percent",
                hue="Cell",
                palette=GRPDEF.sns_colors,
                ax=ax[i],
                size=4.5,
            )

            if i < len(pats) - 1:
                swp.legend_.remove()
            # else:
            #     f.legend(labelspacing=0.18, fontsize=6.5, loc="lower right") # (1.5, 0.))
            #     for lh in swp.legend_.legendHandles:
            #         lh.set_alpha(1)
            #         lh._sizes = [10]
            ax[i].set_title(pats[i].replace("_", " "), fontsize=9)

            ax[i].set_ylim(0, ylims[irow])
            PH.set_axes_ticks(ax=ax[i], 
                yticks = yticks[irow],
                yticks_str = yticks_str[irow], 
            fontsize=8)
        else:
            PH.noaxes(ax[i])

        # print(irow, icol)
        if irow == 2 and icol == 3:
            # print("**********************************")
            custom_legend = []
            for ic, cell in enumerate([2, 5, 6, 9, 10, 11, 13, 17, 18, 30]):
                custom_legend.append(
                Line2D([0], [0], marker='o', color='w', markerfacecolor=GRPDEF.sns_colors[ic], markersize=5, label=f"{FD.BC_name:s}{cell:02d}"))
            
            ax[i].legend(handles=custom_legend, handlelength=1, loc="center left", fontsize=7, title="Cell", labelspacing=0.33)
        # if irow == 2 and icol == 3:
        #     ax[i].text(0.7, 0.95, dataset, fontsize=8)
    dsnames = {"Spont": "Spontaneous Activity", "30dB": "30dB SPL Tone"}
    mpl.suptitle(dsnames[dataset], fontsize=10, fontweight="bold")
    mpl.tight_layout()

    fn = set_figure_path(fignum=5, filedescriptor="Patterns_V3", suppnum=2)
    mpl.savefig(fn)
    fn2 = fn.with_suffix(".png")
    mpl.savefig(fn2)
    mpl.show()





def Figure5E_pattern_plot(axin=None, dataset="Spont", mode:str='mmcd', cell_legend:bool=True):
    """Plot the fraction of inputs for specific subsets of
    input patterns for cells 5, 30, 9 and 17

    Parameters
    ----------
    axin : matplotlib axes instance (default: None)
        If none, make a stand-alone plot
        else, place plot into the specified axes object.
    
    dataset: Which dataset to compute from (default: "Spont")

    mode : which plotting mode to use. (default: mmcd)
        Plot either multiple sets of patterns ("multi"), or plot
        just the ones that include a specific input as the 
        largest ("mmcd" mode)
    
    cell_legend: Enable or disable putting the cell markers up in 
        the legend. If False, the ASA size legend is moved up to the
        top of the legend space


    Raises
    ------
    ValueError
        _description_
    """
    assert mode in ['mmcd', 'multi']
    assert dataset in ["Spont", "15dB", "30dB"]
    compare_file = f"pattern_compares_{dataset:s}.pkl"
    analyzed = Path(compare_file)
    if not analyzed.is_file():
        compares = get_pattern_data(dataset)
        compares.to_pickle(compare_file)

    compares = pd.read_pickle(compare_file)
    print("pattern comparisons:\n", compares.head(10))


    compares.loc[:, "Group"] = compares.apply(lambda row: GRPDEF.get_BC_Group(row.Cell) , axis=1)

    if axin is None:
        f, ax = mpl.subplots(1, 1)
        f.set_size_inches(5, 5)
        asa_scale = 1.0
    else:
        ax = axin
        asa_scale = 0.6
    
    if mode == 'multi': # this was the mode for the early paper development
        pats = [
            "1st_largest_alone",
            "1st_largest+others",
            "2nd_largest_alone",
            "2nd_largest+others",
            "3rd_largest_alone",
            "3rd_largest+others",
            "4th_largest+others",
            "5th_largest+others",
        ]

    elif mode == 'mmcd':
        pats = [  # plot for just the largest
    #        "1st_largest_alone",
            "1st_largest+others",
    #        "2nd_largest_alone",
            "2nd_largest+others",
    #        "3rd_largest_alone",
            "3rd_largest+others",
            "4th_largest+others",
            "5th_largest+others",
        ]

    if mode == "multi":  # a panel with mixed up groups.
        c_cells = [5, 30, 9, 17]
        
        # print(compares.values)
        df = compares[compares["PatternName"].isin(pats) & compares["Cell"].isin(c_cells)]
        markers = []
        for c in c_cells:
            markers.append(GRPDEF.get_group_symbol(c))
        df.loc[:, "Marker"] = df.apply(lambda row: GRPDEF.get_group_symbol(row.Group), axis=1)

        # remap colors for the selected cells
        colors = GRPDEF.sns_colors
        subset_palette = [colors[1], colors[3], colors[7], colors[9]]
        sns.set_palette(subset_palette)
        sns.pointplot(
            x="PatternName",
            y="Percent",
            hue="Cell",
            style="Markers", # "Group",
            data=df,
            # markers=markers,
            scale=0.75,
            kind="swarm",
            ax=ax,
        )
    
    elif mode == "mmcd":
        c_cells = GRPDEF.gradeACells
        i_cells = range(len(c_cells))
        cumulative = {cell: [] for cell in c_cells}
        df = compares[compares["PatternName"].isin(pats) & compares["Cell"].isin(c_cells)].copy()
        # mask = df["Group"] == "Coincidence"
        # df["marker"] = None
        df.loc[:, "Marker"] = df.apply(lambda row: GRPDEF.get_group_symbol(row.Group), axis=1)
        df.loc[:, "Input"] = df.apply(lambda row: int(row.PatternName[0]), axis=1)
        # function for the lambda in apply:
        def get_ASA(cell, input, scale=1.0):
            SC, syninfo = get_synaptic_info(cell)
            return syninfo[1][input-1][0]/scale
        
        df.loc[:, "ASA"] = df.apply(lambda row: get_ASA(row.Cell, row.Input) , axis=1)
        maxASA = np.max(df["ASA"].values)*asa_scale
        df.loc[:, "ScaledASA"] = df.apply(lambda row: get_ASA(row.Cell, row.Input, scale=maxASA), axis=1)
        def get_pos(input, group):
            if group == "Coincidence":
                pos = -0.25+np.random.uniform(low=-0.125, high=0.125)
            else:
                pos = 0.25+np.random.uniform(low=-0.125, high=0.125)
            return(input+pos)
        
        df.loc[:, "pn"] = df.apply(lambda row: get_pos(row.Input, row.Group), axis=1)  # x axis position/offset
        print(f"\n{'GBC##':<6s} {'1st':^6s}  {'2nd':^6s}  {'3rd':^6s}  {'4th':^6s}  {'5th':^6s}  {'sum':^6s}")
        # marker_list = []

        for ic, c in enumerate(c_cells):
            for cpat in pats:
                cumulative[c].append(float(df.loc[(df['Cell'] == c) & (df['PatternName'] == cpat)]['Percent'].values[0]))
            cout = " ".join([f" {x:6.3f}" for x in cumulative[c]])
          #  print(f"{FD.BC_name:s}{c:02d}: {str(cout):s}  {np.sum(cumulative[c]):6.3f}")
            cumulative[c] = np.array(cumulative[c])/np.sum(cumulative[c])
        sns.set_palette(GRPDEF.sns_colors)
        asa_sizes = tuple(np.array([50, 200])*asa_scale)
        fg = sns.scatterplot(x="pn", 
            y="Percent", hue="Cell", size="ASA",
            sizes=asa_sizes, alpha=0.6, palette=GRPDEF.sns_colors,
            style="Group",
            markers={1: GRPDEF.get_group_symbol('Coincidence'),
                     2: GRPDEF.get_group_symbol('MixedMode'),
                    # 2: GRPDEF.get_group_symbol('First-in')
                    },
            # x_jitter=True,
            data=df,
            ax=ax,
            clip_on = False)

    # custom legend
    legend2_y = 1.0
    legend1_y = 1.2
    if cell_legend:
        legend2_y = 0.7
        custom_legend = []
        legorder = [0, 1, 2, 4, 6, 9, 3, 5, 7, 8, None]
        colors = GRPDEF.sns_colors
        for i in range(11):
            legn = legorder[i]
            if legn is None:
                l = Line2D([], [], linewidth=0, label="")
                custom_legend.append(l)
            else:
                cell = c_cells[legn]
                if cell in GRPDEF.MixedMode:
                    marker = 'D'
                    markersize = 4
                else:
                    marker = 'o'
                    markersize=5
                # print("Cell: ", cell, "Marker: ", df[df["Cell"]==cell]["Marker"].values)
                l = Line2D([], [], marker=marker, # df[df["Cell"]==cell]["Marker"].values[0], 
                    color=colors[legn], markerfacecolor=colors[legn], linewidth = 0, markersize=markersize, label=f"{FD.BC_name:s}{cell:02d}")
            custom_legend.append(l)

        legend1 = ax.legend(handles=custom_legend, handlelength=1, 
                loc="upper center", bbox_to_anchor=(0.8, legend1_y),
                fontsize=6, labelspacing=0.35, ncol=2, title="     Cell   \n CD      MM "
                )
        mpl.setp(legend1.get_title(), fontsize='7')

        ax.add_artist(legend1)
    
    legend2 = []
    diams = [50, 100, 200, 250]
    for d in diams:
        l = Line2D([], [], marker='o', color="grey", linewidth=0, markersize=np.sqrt(np.array(d))*asa_scale, label=f"{d:>3d}")
        legend2.append(l)
    legend2 = ax.legend(handles=legend2, labels=diams,
            loc="upper center", bbox_to_anchor=(0.8, legend2_y),
            labelspacing=0.9,
            title="ASA", fontsize=6)

    ax.add_artist(legend2)
    mpl.setp(ax.get_legend().get_title(), fontsize='7')

    ax.set_xlabel("Input Patterns") # 
    ax.set_ylabel("% of spikes evoked by input pattern")
    ax.set_clip_on(False)
    ax.set_ylim(0, 100)

    if axin is not None:
        ax.legend().remove()
    # remap the database names to a more user-friendly name
    namemap = {"1st_largest_alone": "$1^{st}$",
              "1st_largest+others": "$1^{st}+$",
              "2nd_largest_alone": "$2^{nd}$",
              "2nd_largest+others": "$2^{nd}+$",
              "3rd_largest_alone": "$3^{rd}$",
              "3rd_largest+others": "$3^{rd}+$",
              "4th_largest+others": "$4^{th}+$",
              "5th_largest+others": "$5^{th}+$",
    }
    tl = []
    for i, name in enumerate(namemap):
        if i not in [1, 3, 5, 6, 7]:
            continue
        t = namemap[name].replace("_", "\n")
        tl.append(t)
    PH.set_axes_ticks(ax=ax, xticks=[1,2,3,4,5], xticks_str=tl, yticks=[0, 25, 50, 75, 100],
        yticks_str = ['0', '25', '50', '75', '100'], 
        y_minor = [5, 10, 15, 20, 30, 35, 40, 45, 55, 60, 65, 70, 80, 85, 90, 95],
       fontsize=8)
        # ax.set_xticklabels(tl, rotation=0, fontdict={"fontsize": 6.5, "fontweight": "normal"})
    PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)
    if axin is None:
        mpl.tight_layout()
        mpl.show()


if __name__ == "__main__":
    Figure5_Supplemental3_Patterns(reanalyze=False, dataset="Spont")  # supplemental plot for Figure 5
    Figure5E_pattern_plot(dataset="Spont", mode='mmcd')
