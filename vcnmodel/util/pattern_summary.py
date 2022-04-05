""" pattern_summary:
Compute the input patterns, and generate a summary figure
breaking the cells down by 'drivers' and 'coincidence'
for different input patterns

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

import numpy as np
import pandas as pd
import seaborn as sns
import toml
import matplotlib.pyplot as mpl
import vcnmodel.cell_config as cell_config
import vcnmodel.plotters.figure_data as FD
import vcnmodel.plotters.plot_sims as PS
import vcnmodel.util.readmodel as readmodel
from pylibrary.tools import cprint as CP
from pyrsistent import PSet
from vcnmodel.analyzers import reverse_correlation as REVCORR

cprint = CP.cprint

PSC = PS.PlotSims(parent=None)


def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """
    return [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]


@dataclass
class PData:
    """
    data class for some parameters that control what we read
    """

    gradeA: list = field(default_factory=grAList)
    default_modelName: str = "XM13_nacncoop"
    soma_inflate: bool = True
    dend_inflate: bool = True
    basepath: str = ""  # config["baseDataDirectory"]
    renderpath: str = ""  # " str(Path(self.config["codeDirectory"], "Renderings"))
    thiscell: str = ""


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


def get_pattern_data():

    where_is_data = Path("wheres_my_data.toml")
    if where_is_data.is_file():
        datapaths = toml.load("wheres_my_data.toml")
    else:
        raise ValueError()
    basepath = datapaths["cellDataDirectory"]
    fig_data = FD.figure_revcorr
    db = "Spont"
    #    group = {"First-in": [9, 11, 13, 17], "Coincidence": [2, 5, 6, 10, 18, 30]}
    group = {"First-in": [9, 11, 13, 17, 18], "Coincidence": [2, 5, 6, 10, 30]}

    compares = pd.DataFrame(columns=["Cell", "PatternName", "Percent", "Group"])
    for gbc in grAList():
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
        pklf = read_pkl_file(filename)
        # PSC.compute_revcorr(
        #     P=None,
        #     gbc=str(gbc),
        #     fn=pklf.files[0],
        #     PD = PData(),
        #     protocol="runANPSTH",
        #     revcorrtype="RevcorrSPKS",
        # )
        if gbc in group["First-in"]:
            gname = "First-in"
        else:
            gname = "Coincidence"
        RCP, RCD, PAT = revcorr(gbc, pklf.files[0], PD=PData())
        for t in PAT.filttable.keys():
            pat = PAT.filttable[t]
            REVCORR.print_pattern(t, pat)
            pct = 100.0 * pat.sumtrue / (pat.sumtrue + pat.sumfalse)
            df = pd.DataFrame(
                [{"Cell": gbc, "PatternName": t, "Percent": pct, "Group": gname}]
            )
            compares = pd.concat([compares, df], ignore_index=True)
    return compares


compares = get_pattern_data()
compares.to_pickle('compares.pkl')

compares = pd.read_pickle("compares.pkl")
# print(compares.head(10))
f, ax = mpl.subplots(4, 3)
f.set_size_inches(5, 8)
ax = ax.ravel()
pats = [
    "1st_largest_alone",
    "1st_largest+others",
    "1st+2nd_alone",
    "2nd_largest_alone",
    "2nd_largest+others",
    "1st+2nd+others",
    "3rd_largest_alone",
    "3rd_largest+others",
    "1st+2nd+3rd+others",
    "not_largest",
    "not_2_largest",
    "not_3_largest",
]
pnames = sorted(set(pats))
pdatanames = sorted(set(compares["PatternName"].values))
assert pnames == pdatanames

# create your own color array
bar_colors = [(0.8, 0.8, 0.8, 0.2), (0.0, 0.0, 1.0, 0.5)]

# add color array to set_palette
# function of seaborn
sns.set_palette(bar_colors, desat=0.2)

for i, axl in enumerate(pats):
    # with sns.set_palette(bar_colors):
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
        ax=ax[i],
        size=3,
        # color="black",
    )
    # print(dir(swp))
    # print(dir(swp.legend_))
    if i < len(pats) - 1:
        swp.legend_.remove()
    else:
        # print(dir(swp.legend_))
        swp.legend(labelspacing=0.18, fontsize=6.5)
        for lh in swp.legend_.legendHandles:
            lh.set_alpha(1)
            lh._sizes = [10]
            # print(dir(lh_label))
    ax[i].set_title(pats[i], fontsize=9)
    ax[i].set_ylim(0, 100)
mpl.tight_layout()
where_is_data = Path("wheres_my_data.toml")
if where_is_data.is_file():
    datapaths = toml.load("wheres_my_data.toml")
else:
    raise ValueError()
mpl.savefig(
    Path(
        datapaths["baseDataDirectory"],
        datapaths["figureDirectory"],
        "Figure4",
        "Figure4_supp",
        "Figure4-Supplemental3_Patterns.pdf",
    )
)
mpl.show()
