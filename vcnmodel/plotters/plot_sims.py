import argparse
import dataclasses
import datetime
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import numpy as np
import seaborn
from ephys.ephysanalysis import MakeClamps, RmTauAnalysis, SpikeAnalysis
import matplotlib.colors
import matplotlib.colorbar
from matplotlib import pyplot as mpl
from matplotlib import rc
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as PLS
from pylibrary.tools import cprint as CP
import vcnmodel.model_params
from vcnmodel import cell_config as cell_config
from vcnmodel import spikestatistics as SPKS
from cnmodel.util import vector_strength
cprint = CP.cprint

"""
Functions to compute some results and plot the
simulation results from model_run2
Reads the new format filenames

This is written as a set of functions that can be called
from elsewhere.
Included is a parser so that the plots
can be called from the command line as well.

Some of the parameters must be instantiated by creating an instance
of PData that is passed into the routnines.

Wrapper for various analysis functions, handles multiple cells.


Command line usage:

usage: plot_sims.py [-h] [-p {AN,an,IO,IV,iv,gifnoise}]
                    [-a {traces,PSTH,revcorr,SAC,tuning,singles}]
                    [-M MODELTYPE] [-s]
                    [-e {None,delays,largestonly,removelargest,mean,allmean,twolargest}]
                    [--dendritemode {normal,passive,active}] [-d DBSPL]
                    [-r NREPS] [-c]
                    cell [cell ...]

Plot GBC results

positional arguments:
  cell                  Select the cell(s) or 'A' for all(no default)

optional arguments:
  -h, --help            show this help message and exit
  -p {AN,an,IO,IV,iv,gifnoise}, --protocol {AN,an,IO,IV,iv,gifnoise}
                        Select the protocol (default: IV) from: ['AN', 'an',
                        'IO', 'IV', 'iv', 'gifnoise']
  -a {traces,PSTH,revcorr,SAC,tuning,singles}, --analysis {traces,PSTH,revcorr,SAC,tuning,singles}
                        Select the analysis type (default: traces) from:
                        ['traces', 'PSTH', 'revcorr', 'SAC', 'tuning',
                        'singles']
  -M MODELTYPE, --modeltype MODELTYPE
                        Select the model type (default XM13_nacncoop) from:
                        ['mGBC', 'XM13', 'RM03', 'XM13_nacncoop']
  -s, --scaled          use scaled data or not
  -e {None,delays,largestonly,removelargest,mean,allmean,twolargest},
                --experiment {None,delays,largestonly,removelargest,mean,allmean,twolargest}
                        Select the experiment type from: [None, 'delays',
                        'largestonly', 'removelargest', 'mean', 'allmean',
                        'twolargest']
  --dendritemode {normal,passive,active}
                        Choose dendrite table (normal, active, passive)
  -d DBSPL, --dB DBSPL  Select the models at specific intensity
  -r NREPS, --nreps NREPS
                        Select the models with # reps
  -c, --check           Just check selection criteria and return

"""


# import pylibrary.tools.utility as PU
# make a shortcut for each of the clases
AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()
rc("text", usetex=False)

modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
experimenttypes = [
    None,
    "delays",
    "largestonly",
    "removelargest",
    "mean",
    "allmean",
    "twolargest",
]
modetypes = ["find", "singles", "IO", "multi"]
analysistypes = ["traces", "PSTH", "revcorr", "SAC", "tuning", "singles"]
dendriteChoices = [
    "normal",
    "passive",
    "active",
]

orient_cells = {
    2: [140.0, 0.0, -144.0],
    6: [140.0, -59.0, -12.0],
    5: [140.0, -46.0, 121.0],
    9: [140.0, -74.0, 18.0],
    11: [140.0, -2.0, -181.0],
    10: [140.0, 5.0, -35.0],
    13: [140.0, -22.0, 344.0],
    17: [140.0, -158.0, 39.0],
    30: [140.0, -134.0, -181.0],
}


def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """

    return [2, 5, 6, 9, 10, 11, 13, 17, 30]


@dataclass
class PData:
    """
    data class for some parameters that control what we read
    """

    gradeA: list = field(default_factory=grAList)
    default_modelName: str = "XM13_nacncoop"
    soma_inflate: bool = True
    dend_inflate: bool = True
    basepath: str = "/Users/pbmanis/Desktop/Python/VCN-SBEM-Data"
    renderpath: str = "/Users/pbmanis/Desktop/Python/vcnmodel/Renderings"
    thiscell: str = ""


def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


@dataclass
class RevCorrPars:
    ntrials: int = 1
    ninputs: int = 1

    # clip trace to avoid end effects
    min_time: float = 10.0  # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time: float = 250.0  # this window needs to be at least as long as maxwin
    binw: float = 0.1
    minwin: float = -5
    maxwin: float = 2.5
    amax: float = 0.0


@dataclass
class RevCorrData:
    C: list = field(default_factory=def_empty_list)
    TC: list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)
    tx: np.array = field(default_factory=def_empty_np)
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(default_factory=def_empty_np)
    sv_all: np.array = field(default_factory=def_empty_np)
    sv_avg: np.array = field(default_factory=def_empty_np)
    sites: np.array = field(default_factory=def_empty_np)
    nsp: int = 0
    max_coin_rate: float = 0


def norm(p: Union[list, np.ndarray], n: int) -> np.ndarray:
    """
    Simple function to normalize the n'th point of p
    by the min and max
    """
    pmin = np.min(p)
    pmax = np.max(p)
    return (p[n] - pmin) / float(pmax - pmin)


def twinax(fig: object, ax1: object, pos: float = 0.0) -> object:
    """
    Create a 'twin' axis on the right side of a plot
    Note: pyqtgraph.plotting.styles can also do an inset
    which may be used instead
    """
    ax2 = fig.add_axes(ax1.get_position(True), sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.set_offset_position("right")
    ax2.tick_params(direction="in", length=5.0, width=1.0, labelsize=6)
    ax2.spines["right"].set_position(("data", pos))
    ax2.spines["left"].set_color("none")
    ax2.spines["top"].set_color("none")
    # ax2.set_autoscalex_on(ax1.get_autoscalex_on())
    # ax1.yaxis.tick_left()
    ax2.xaxis.set_visible(False)
    ax2.patch.set_visible(False)
    #    PH.adjust_spines(ax2, distance=0.)
    return ax2


def get_changetimestamp():
    # trip filemode based on date of simulatoin
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate, "%Y-%m-%d-%H:%M")
    changetimestamp = datetime.datetime.timestamp(dts)
    return changetimestamp


def clean_spiketimes(spikeTimes, mindT=0.7):
    """
    Clean up spike time array, removing all less than mindT
    spikeTimes is a 1-D list or array
    mindT is difference in time, same units as spikeTimes
    If 1 or 0 spikes in array, just return the array
    """
    if len(spikeTimes) > 1:
        dst = np.diff(spikeTimes)
        st = np.array(spikeTimes[0])  # get first spike
        sok = np.where(dst > mindT)
        st = np.append(st, [spikeTimes[s + 1] for s in sok])
        # print st
        spikeTimes = st[~np.isnan(st)]
    return spikeTimes


def get_data_file(
    fn: Union[str, Path], changetimestamp: object, PD: dataclass
) -> Union[None, tuple]:
    """
    Get a data file, and also parse information from the file
    for display
    """
    fnp = Path(fn)
    fns = str(fn)
    ivdatafile = None
    if not fnp.is_file():
        cprint('r', f"   File: {str(fnp):s} NOT FOUND")
        return None
    else:
        print(f"   File {str(fnp):s} found.")
    mtime = fnp.stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d-%H:%M")
    print(f"pgbcivr2: Checking file: {fnp.name:s} [{timestamp_str:s}]")
    # print(mtime, changetimestamp)
    if mtime > changetimestamp:
        filemode = "vcnmodel.v1"
    else:
        filemode = "vcnmodel.v0"
    print("pgbcivr2: file mode: ", filemode)
    with (open(fnp, "rb")) as fh:
        d = pickle.load(fh)
    # print(d.keys())
    # if "Params" not in list(d.keys()):
    #               print("File missing Params; need to re-run")
    #               continue
    # print(d['Results'][0]['inputSpikeTimes'])
    #

    if filemode in ["vcnmodel.v0"]:
        # print(d['runInfo'].keys())
        par = d["runInfo"]
        par["soma_inflation"] = False
        par["dendrite_inflation"] = False
        if fns.find("soma=") > -1:
            par["soma_inflation"] = True
        if fns.find("dend=") > -1:
            par["dendrite_inflation"] = True
        par["soma_autoinflate"] = False
        par["dendrite_autoinflate"] = False
    elif filemode in ["vcnmodel.v1"]:
        try:
            par = d["Params"]
        except ValueError:
            try:
                par = d["self.Params"]
            except ValueError:
                raise ValueError("File missing Params; need to re-run")
        if isinstance(par, vcnmodel.model_params.Params):
            par = dataclasses.asdict(par)
    # print('pgbcivr2: Params: ', par)
    print(PD)
    print(par)
    if PD.soma_inflate and PD.dend_inflate:
        if par["soma_inflation"] > 0.0 and par["dendrite_inflation"] > 0.0:
            ivdatafile = Path(fn)
            stitle = "Soma and Dend scaled"
            print(stitle)
    elif PD.soma_inflate and not PD.dend_inflate:
        if par["soma_inflation"] > 0.0 and par["dendrite_inflation"] < 0.0:
            ivdatafile = Path(fn)
            stitle = "Soma only scaled"
            print(stitle)
    elif PD.dend_inflate and not PD.soma_inflate:
        if par["soma_inflation"] < 0.0 and par["dendrite_inflation"] > 0.0:
            ivdatafile = Path(fn)
            stitle = "Dend only scaled"
            print(stitle)
    elif not par["soma_autoinflate"] and not par["dendrite_autoinflate"]:
        print("\nConditions x: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)
        ivdatafile = Path(fn)
        stitle = "No scaling (S, D)"
        print(stitle)
    else:
        return None

    if ivdatafile is None:
        print("no file matching conditions : ", str(ivdatafile))
        return None

    if not ivdatafile.is_file():
        print("no file? : ", str(ivdatafile))
        return None

    print("\npgbcivr2: datafile to read: ", str(ivdatafile))
    return par, stitle, ivdatafile, filemode, d


def analyze_data(ivdatafile: Union[Path, str], filemode, protocol: str) -> tuple:
    """
    Provide basic spike detection, shape analysis, and
    IV analysis if appropriate
    We use ephys.acq4read.read_pfile to read the pickled data
    file into acq4 format, then we can use the ephys
    analysis tools to analyze the data
    """
    AR.read_pfile(ivdatafile, filemode=filemode)
    bridge_offset = 0.0
    threshold = -32.0  # mV
    tgap = 0.0  # gap before fittoign taum
    RM.setup(AR, SP, bridge_offset=bridge_offset)
    SP.setup(
        clamps=AR,
        threshold=threshold,
        refractory=0.0001,
        peakwidth=0.002,
        interpolate=True,
        verify=True,
        mode="peak",
    )

    SP.set_detector("Kalluri")  # spike detector

    # AR.tstart = 1.0
    # AR.tend = 0.4*1e3
    # AR.tdur = 0.399*1e3
    # print(AR.tstart, AR.tend, AR.tdur, AR.sample_rate)
    SP.analyzeSpikes()
    SP.analyzeSpikeShape()

    RMA = None
    if protocol == "IV":
        # for trial in range(len(AR.traces)):
        #     mpl.plot(AR.time_base, AR.traces[trial]*1e3, linewidth=0.5)
        # mpl.show()
        # exit()
        SP.fitOne(function="fitOneOriginal")
        RM.analyze(
            rmpregion=[0.0, AR.tstart - 0.001],
            tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
            to_peak=True,
            tgap=tgap,
        )

        RMA = RM.analysis_summary
    # print("Analysis: RMA: ", RMA)

    return AR, SP, RMA


def plot_traces(
    ax: object, fn: Union[Path, str], PD: dataclass, protocol: str,
) -> tuple:
    changetimestamp = get_changetimestamp()
    x = get_data_file(fn, changetimestamp, PD)
    mtime = Path(fn).stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d-%H:%M")
    if x is None:
        print("No simulation found that matched conditions")
        print(fn)
        return
    # unpack x
    inx = str(fn).find('_Syn')
    synno = None
    if inx > 0:
        synno = int(fn[inx+4:inx+7])

    par, stitle, ivdatafile, filemode, d = x
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
    ntr = len(AR.traces)  # number of trials
    v0 = -160.0
    deadtime = 50.
    trstep = 25.0 / ntr
    inpstep = 5.0 / ntr
    sz = 50.0 / ntr
    noutspikes = 0
    ninspikes = 0
    for trial in range(len(AR.traces)):
        ax.plot(AR.time_base, AR.traces[trial] * 1e3, linewidth=0.5)
        ax.plot(
            AR.time_base[SP.spikeIndices[trial]],
            AR.traces[trial][SP.spikeIndices[trial]] * 1e3,
            "ro",
            markersize=2.5,
        )
 
        sinds = np.array(SP.spikeIndices[trial])*AR.sample_rate[trial]
        noutspikes += len(np.argwhere(sinds > deadtime))
        if protocol in ["AN", "runANSingles"]:
            if (trial in list(d['Results'].keys())
                and "inputSpikeTimes" in list(d["Results"][trial].keys())):
                    spkt = d["Results"][trial]["inputSpikeTimes"]
            elif "inputSpikeTimes" in list(d['Results'].keys()):
                 spkt = d["Results"]["inputSpikeTimes"][trial]
            # print('input spike trains: ', len(spkt))
            # print('spkt: ', spkt)
            tr_y = trial * (trstep + len(spkt) * inpstep)
            if synno is None:
                for ian in range(len(spkt)):
                    vy = v0 + tr_y * np.ones(len(spkt[ian])) + inpstep * ian
                    ax.scatter(spkt[ian], vy, s=sz, marker="|", linewidths=0.35)
            else:
                ian = synno
                vy = v0 + tr_y * np.ones(len(spkt[ian])) + inpstep * ian
                ax.scatter(spkt[ian], vy, s=sz, marker="|", linewidths=0.35)
                ninspikes += len(spkt[ian] > deadtime)

                # print(len(vy), vy)
            #                 print(spkt[ian])
            ax.set_ylim(-140.0, 40.0)
        else:
            ax.set_ylim(-200.0, 50.0)
    # print('Nout/Nin: ', float(noutspikes)/ninspikes)
    ax.set_xlim(0.080, np.max(AR.time_base))

    ftname = str(ivdatafile.name)
    ip = ftname.find("_II_") + 4
    ftname = ftname[:ip] + "...\n" + ftname[ip:]
    toptitle = f"{ftname:s}"
    if protocol == "IV":
        toptitle += f"\nRin={RMA['Rin']:.1f} Mmhm1  Taum={RMA['taum']:.2f} ms"

        secax = PLS.create_inset_axes([0.4, 0, 0.4, 0.4], ax)
        secax.plot(
            RM.ivss_cmd * 1e12,
            RM.ivss_v * 1e3,
            "ks-",
            markersize=3,
            markerfacecolor="k",
            zorder=10,
            clip_on=False,
        )
        ltz = np.where(RM.ivss_cmd <= 0)[0]
        secax.plot(
            RM.ivss_cmd[ltz] * 1e12,
            RM.ivpk_v[ltz] * 1e3,
            "ko-",
            markersize=3,
            markerfacecolor="w",
            zorder=10,
            clip_on=False,
        )
        PH.crossAxes(
            secax, xyzero=[0.0, -60.0], limits=[-1.0 * 1e3, -150, 1.0 * 1e3, -25.0]
        )
        PH.talbotTicks(
            secax,
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.02,
            pointSize=6,
            tickPlacesAdd={"x": 0, "y": 0},
            floatAdd={"x": 0, "y": 0},
        )
    PH.calbar(
        ax,
        calbar=[20.0, -160.0, 10.0, 20.0],
        unitNames={"x": "ms", "y": "mV"},
        fontsize=9,
    )
    if RMA is not None:
        PH.referenceline(ax, RMA["RMP"])
        ax.text(
            -1.0,
            RMA["RMP"],
            f"{RMA['RMP']:.1f}",
            verticalalignment="center",
            horizontalalignment="right",
            fontsize=9,
        )
    toptitle += f"\n{timestamp_str:s}"
    ax.set_title(toptitle, fontsize=5)
    return(synno, noutspikes, ninspikes)


def plot_revcorr_map(
    P,
    pgbc,
    inputlist,
    ntrials,
    C,
    TC,
    st,
    tx,
    ti_avg,
    sv_all,
    sv_avg,
    sites,
    nsp,
    max_coin_rate,
    window,
):
    pass


def plot_revcorr2(ax: object, PD: dataclass, RCP: dataclass, RCD: dataclass):
    seaborn.set_style("ticks")
    # secax = twinax(P.figure_handle, ax, pos=maxwin)
    secax = PLS.create_inset_axes([0, 0.5, 1, 0.5], ax)
    PH.noaxes(secax, "xy")

    secax.set_facecolor((1, 1, 1, 0))
    secax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # ax.set_facecolor((0.7, 0.7, 0.7))
    summarySiteTC = {}
    for isite in range(
        RCP.ninputs
    ):  # range(ninputs):  # for each ANF input (get from those on first trial)
        # print('svall: ', RCD.sv_all, type(RCD.sv_all), RCD.sv_all.shape)
        if RCD.sv_all.shape == ():
            continue
        stepsize = int(RCD.sv_all.shape[0] / 20)
        if stepsize > 0:
            sel = list(range(0, RCD.sv_all.shape[0], stepsize))
        else:
            sel = list(range(0, RCD.sv_all.shape[0], 1))
        sel = list(range(0, RCD.sv_all.shape[0], 1))
        # if stelen(sel) < 10:
        #     if sv_all.shape[0] > 10:
        #         sel = list(range(10))
        #     else:
        #         sel = list(range(sv_all.shape[0]))

        if RCD.C[isite] is not None:
            # print('plotting C')
            nc = int(len(RCD.C[isite]) / 2)
            RCD.TC = RCD.TC / len(RCD.st)
            summarySiteTC[isite] = RCD.TC
            # print(RCD.sites, isite)
            color = mpl.cm.viridis(norm(RCD.sites, isite))
            # print('color', color)
            ax.plot(
                RCD.tx,
                RCD.C[isite][:nc],
                color=color,
                label=("Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))),
                linewidth=1.5,
                zorder=5,
            )
            RCD.max_coin_rate = np.max((RCD.max_coin_rate, np.max(RCD.C[:nc])))

        if (
            isite == 0
        ):  # only plot the first time through - the APs are the same no matter the trial
            # print('Plotting V')
            for t in sel:
                print(len(RCD.sv_all))
                secax.plot(
                    RCD.ti_avg, RCD.sv_all[t], color="#666666", linewidth=0.2, zorder=1
                )

            secax.plot(RCD.ti_avg, RCD.sv_avg, color="k", linewidth=0.75, zorder=2)
            PH.calbar(
                secax,
                calbar=[1.0, -10, 1.0, 20.0],
                axesoff=True,
                orient="right",
                unitNames={"x": "ms", "y": "mV"},
            )
            PH.referenceline(secax, -60.0)
            seaborn.despine(ax=secax)
    print(f"Total spikes plotted: {RCD.nsp:d}")
    # for trial in range(RCP.ntrials):
    #     stx = RCD.st[trial]
    secax.plot([0.0, 0.0], [-120.0, 10.0], "r", linewidth=0.5)
    # print('finished inputs')
    seaborn.despine(ax=ax)
    ax.set_ylabel("Rate of coincidences/bin (Hz)", fontsize=10)
    ax.set_xlabel("T (ms)", fontsize=10)
    ax.set_xlim((RCP.minwin, RCP.maxwin))
    if 0.2 < RCD.max_coin_rate < 1.0:
        ax.set_ylim(0, 1)
    elif 0.2 >= RCD.max_coin_rate:
        ax.set_ylim(0, 0.25)
    elif RCD.max_coin_rate > 1.0:
        ax.set_ylim(0, 2.0)
    secax.set_ylim([-70.0, 10.0])
    secax.set_xlim((RCP.minwin, RCP.maxwin))
    # secax.set_ylabel('Vm', rotation=-90., fontsize=10)
    secax.tick_params(direction="in", length=5.0, width=1.0, labelsize=9)
    ax.tick_params(direction="in", length=5.0, width=1.0, labelsize=9)
    # ax2.set_ylabel('Total Correlation W=%.1f-%0.1f'% (tcwidth[0], tcwidth[1]), fontsize=12)
    # ax2.set_ylim(0, 1.0) # maxtc*1.2)
    # PH.talbotTicks(ax2, axes='xy',
    #                density=(1.0, 1.0), insideMargin=0.05, pointSize=10,
    #                tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1})

    # a = re_self.search(fn[p])
    # b = re_c10.search(fn[p])
    # if a is not None:
    #     inp = 'VCN_c09'
    # elif b is not None:
    #     inp = 'VCN_c10'
    # else:
    #     inp = "Undefined"
    # ax.set_title(f'VCN_{p:s} input={inp:s} \
    # [{int(np.min(sites)):d}-{int(np.max(sites)):d}]\nAmax={amax:.1f}', y=0.9, x=0.02,
    #     horizontalalignment='left', fontsize=6)
    return summarySiteTC


def get_synaptic_info(gbc: str) -> tuple:
    SC = cell_config.CellConfig()
    syninfo = SC.VCN_Inputs[gbc]
    return (SC, syninfo)


def get_data(
    fn: Union[Path, str], PD: dataclass, changetimestamp, protocol
) -> Union[None, tuple]:

    X = get_data_file(fn, changetimestamp, PD)
    # mtime = Path(fn).stat().st_mtime
    # timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d-%H:%M")
    if X is None:
        print("No simulation found that matched conditions")
        print("Looking for file: ", fn)
        return None
    # unpack x
    par, stitle, ivdatafile, filemode, d = X

    # 2. find spikes
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
    # set up analysis parameters and result storage
    RCP = RevCorrPars()
    RCD = RevCorrData()
    RCD.st = SP.spikeIndices
    trials = list(d["Results"].keys())
    RCP.ntrials = len(trials)
    RCD.npost = 0  # number of postsynaptic spikes
    RCD.npre = 0  # number of presynaptic spikes
    for i, tr in enumerate(trials):
        trd = d["Results"][tr]  # trial data
        ti = trd["time"]
        # sv = trd["somaVoltage"]
        # dt = ti[1] - ti[0]
        # st[tr] = PU.findspikes(ti, sv, thr, dt=dt, mode='threshold')
        # st[tr] = clean_spiketimes(st[tr])
        RCD.npost += len(RCD.st[tr])
        for n in range(len(trd["inputSpikeTimes"])):  # for each sgc
            RCD.npre += len(trd["inputSpikeTimes"][n])
    RCD.ti = ti
    print(f"Detected {RCD.npost:d} Post spikes")
    print(f"Detected {RCD.npre:d} Presynaptic spikes")
    print("# trials: ", RCP.ntrials)

    # clip trace to avoid end effects
    RCP.max_time = (
        np.max(ti) - RCP.min_time
    )  # this window needs to be at least as long as maxwin
    RCD.tx = np.arange(RCP.minwin, 0, RCP.binw)
    return (d, AR, SP, RMA, RCP, RCD)


def compute_revcorr(
    ax: object,
    gbc: str,
    fn: Union[str, Path],
    PD: object,
    protocol: str,
    thr: float = -20.0,
    width: float = 4.0,
) -> Union[None, tuple]:

    changetimestamp = get_changetimestamp()
    #
    # 1. Gather data
    #
    SC, syninfo = get_synaptic_info(gbc)

    res = get_data(fn, PD, changetimestamp, protocol)
    if res is None:
        return None
    (d, AR, SP, RMA, RCP, RCD) = res
    RCP.ninputs = len(syninfo[1])
    RCD.sites = np.zeros(RCP.ninputs)
    for isite in range(RCP.ninputs):  # precompute areas
        area = syninfo[1][isite][0]
        if area > RCP.amax:
            RCP.amax = area
        RCD.sites[isite] = int(np.around(area * SC.synperum2))

    print("ninputs: ", RCP.ninputs)
    maxtc = 0
    si = d["Params"]
    ri = d["runInfo"]
    # print(si)
    # print(ri)
    # print(isinstance(si, dict))
    # print(isinstance(ri, dict))
    if isinstance(si, dict):
        min_time = (si["pip_start"] + 0.025 )*1000. # push onset out of the way
        max_time = (si["pip_start"]+si["pip_duration"])*1000.
        F0 = si['F0']
        dB = si["dB"],
    else:
        if isinstance(ri.pip_start, list):
            start = ri.pip_start[0]
        else:
            start = ri.pip_start
        min_time = (start+0.025)*1000.
        max_time = (start+ri.pip_duration)*1000.
    RCD.sv_sites = []
    RCD.C = [None] * RCP.ninputs
    RCD.max_coin_rate = 0.0
    nspk_plot = 0
    # spksplotted = False
    RCP.min_time = min_time  # driven window without onset
    RCP.max_time = max_time

    for isite in range(RCP.ninputs): # for each ANF input
        # print('isite: ', isite)
        firstwithspikes = False
        # RCD.sv_avg = np.zeros(1)  # these are allocated on first trial for each input run
        # RCD.sv_all = np.zeros(1)
        RCD.nsp = 0
        # n_avg = 0

        for trial in range(RCP.ntrials):  # sum across trials
            stx = AR.time_base[SP.spikeIndices[trial]]
            stx = stx[
                (stx > RCP.min_time) & (stx < RCP.max_time)
            ]  # more than minimum delay
            anx = d["Results"][trial]["inputSpikeTimes"][isite]
            anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
            if len(stx) == 0 or len(anx) == 0:
                # print('empty array for stx or anx')
                continue
            # andirac = np.zeros(int(200.0 / RCP.binw) + 1)
            if not firstwithspikes:
                firstwithspikes = True
                RCD.C[isite] = SPKS.correlogram(
                    stx, anx, width=-RCP.minwin, bin=RCP.binw, T=None
                )
                RCD.TC = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
                if np.isnan(RCD.TC):
                    RCD.TC = 0.0
            else:
                RCD.C[isite] += SPKS.correlogram(
                    stx, anx, width=-RCP.minwin, bin=RCP.binw, T=None
                )
                tct = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
                if ~np.isnan(tct):
                    RCD.TC = RCD.TC + tct

            # accumulate postsynaptic spike waveforms
            if RCD.nsp == 0:  # first spike in trace
                reltime = np.around(RCD.ti, 5) - np.around(stx[0], 5)
                areltime = np.argwhere(
                    (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
                ).squeeze()
                RCD.sv_avg = d["Results"][trial]["somaVoltage"][areltime]
                RCD.ti_avg = RCD.ti[0 : len(areltime)] + RCP.minwin
                RCD.sv_all = np.zeros((RCP.ntrials, RCD.ti_avg.shape[0]))
                RCD.nsp += 1

            else:  # rest of spikes
                for n in range(1, len(stx)):
                    reltime = np.around(RCD.ti, 5) - np.around(stx[n], 5)
                    areltime = np.argwhere(
                        (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
                    ).squeeze()
                    if len(areltime) > len(RCD.sv_avg):
                        areltime = areltime[0 : len(RCD.sv_avg)]
                    if len(areltime) < len(RCD.sv_avg):
                        nextend = len(RCD.sv_avg) - len(areltime)
                        areltime = np.append(
                            areltime,
                            np.arange(areltime[-1] + 1, areltime[-1] + nextend + 1),
                        )
                    if trial == 0:
                        RCD.sv_avg = d["Results"][trial]["somaVoltage"][areltime]
                        RCD.sv_all = np.zeros((RCP.ntrials, RCD.ti_avg.shape[0]))  # initialize the array
                        RCD.ti_avg = RCD.ti[0 : len(areltime)] + RCP.minwin
                    else:
                        RCD.sv_avg += d["Results"][trial]["somaVoltage"][areltime]
                    RCD.sv_all[trial] =  d["Results"][trial]["somaVoltage"][areltime]
                    RCD.nsp += 1

            nspk_plot += RCD.nsp
            RCD.sv_sites.append(RCD.sv_all)
        # RCD.C[isite] /= RCP.ntrials
        if RCD.nsp > 0:
            RCD.sv_avg /= RCD.nsp
            # print('trial: ', i, sv_avg)
            # C = C + SPKS.spike_triggered_average(st[trial]*ms, andirac*ms, max_interval=width*ms, dt=binw*ms)
        # print('RCD: ', RCD.TC)
        if isinstance(RCD.TC, float) and RCD.TC > maxtc:
            maxtc = RCD.TC
        else:
            maxtc = 1.0

    pre_w = [-2.7, -0.7]

    summarySiteTC = plot_revcorr2(ax, PD, RCP, RCD)
    # return summarySiteTC, RCD.sites

    ################
    # now some pariwise, etc. stats on events prior to a spike
    ################

    # tind = np.where((RCD.tx > pre_w[0]) & (RCD.tx < pre_w[1]))[0]
    pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
    participation = np.zeros(RCP.ninputs)
    nperspike = []
    nspikes = 0
    for trial in range(RCP.ntrials):  # accumulate across all trials
        spks = AR.time_base[SP.spikeIndices[trial]]
        for s in spks:  # for each postsynaptic spike
            if (
                s < RCP.min_time or s > RCP.max_time
            ):  # trim to those only in a response window
                continue
            # print('pre: ', s)
            nspikes += 1
            nps = 0
            for isite in range(RCP.ninputs):  # test inputs in a window prior
                anxi = d["Results"][trial]["inputSpikeTimes"][
                    isite
                ]  # input spike times for one input
                anxi = anxi[
                    (anxi > RCP.min_time) & (anxi < RCP.max_time)
                ]  # trim to those only in a response window
                ani = anxi - s
                # print('ani: ', ani)
                nevi = len(
                    np.where((ani >= pre_w[0]) & (ani <= pre_w[1]))[0]
                )  # count spikes in ith window
                if nevi > 0:
                    participation[isite] += 1
                    nps += 1
                for jsite in range(isite + 1, RCP.ninputs):

                    anxj = d["Results"][trial]["inputSpikeTimes"][jsite]
                    anxj = anxj[(anxj > RCP.min_time) & (anxj < RCP.max_time)]
                    anj = anxj - s
                    nevj = len(np.where((anj >= pre_w[0]) & (anj <= pre_w[1]))[0])
                    if isite != jsite:
                        if nevj > 0 and nevi > 0:
                            # print(nevi, nevi)
                            pairwise[isite, jsite] += 1
                    else:
                        if nevj > 0:
                            pairwise[isite, jsite] += 1
            nperspike.append(nps)

    s_pair = np.sum(pairwise)
    if s_pair > 0.0:
        pairwise /= s_pair
    psh = pairwise.shape
    pos = np.zeros((psh[0], psh[1], 2))

    for i in range(RCP.ninputs):
        for j in range(RCP.ninputs):
            # print(f"{pairwise[i,j]:.3f}", end='  ')
            pos[i, j, 0] = i + 1
            pos[i, j, 1] = j + 1
        # print()
    # print(nperspike)
    import scipy.stats

    # print(np.unique(nperspike, return_counts=True))
    nperspike = [n for n in nperspike if n != 0]
    nperspike = scipy.stats.itemfreq(nperspike).T
    # print('nperspike counts: ', nperspike)
    # nperspike = np.array(np.unique(nperspike, return_counts=True))/nspikes
    # properly fill out output
    # xnspike = np.arange(RCP.ninputs)
    ynspike = np.zeros(RCP.ninputs)
    for j, i in enumerate(nperspike[0]):
        # print(i, j, nperspike[1,j])
        ynspike[i - 1] = nperspike[1, j]

    ynspike = np.cumsum(ynspike / nspikes)

    # print(RCD.sites)
    # print(pos)
    maxp = np.max(pairwise)
    print(pairwise)
    PSum = PH.regular_grid(
        rows=2,
        cols=2,
        order="rowsfirst",
        figsize=(6, 6),
        showgrid=False,
        verticalspacing=0.08,
        horizontalspacing=0.08,
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.1,
            "rightmargin": 0.1,
            "topmargin": 0.1,
        },
        label=["A", "B", "C", "D"],
        labelposition=(-0.05, 1.05),
    )
    sax = PSum.axdict
    # f, sax = mpl.subplots(3,1)
    # f.set_size_inches( w=3.5, h=9)
    sax["A"].plot(np.arange(RCP.ninputs) + 1, RCD.sites, "bo")
    print('pairwise: ', pairwise)
    colormap = "plasma"
    if s_pair > 0.0:
        pclip = np.clip(pairwise, np.min(np.where(pairwise > 0)), np.max(pairwise))
        pclip[np.where(pclip == 0)] = np.nan
        pclipped = pclip - np.nanmin(pclip)
        sax["C"].scatter(
            pos[:, :, 0], pos[:, :, 1], s=200 * pairwise / maxp, c=pclipped, cmap=colormap
        )
        vmax = np.nanmax(pclip) * 100
        vmin = np.nanmin(pclip) * 100
        print("vmin, vmax: ", vmin, vmax)
    else:
        vmin = 0
        vmax = 1
    # sax['B'].plot(np.arange(RCP.ninputs)+1, participation/nspikes, 'gx')
    sax["B"].plot(RCD.sites, participation / nspikes, "gx")
 


    axcbar = PLS.create_inset_axes([0.8, 0.05, 0.05, 0.5], sax["C"])
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    ticks = np.linspace(vmin, vmax, num=4, endpoint=True)
    cm_sns = mpl.cm.get_cmap(colormap)
    c2 = matplotlib.colorbar.ColorbarBase(axcbar, cmap=cm_sns, ticks=ticks, norm=norm)

    # PH.nice_plot(sax['C'], position=-0.2)
    sax["D"].plot(np.arange(RCP.ninputs) + 1, ynspike, "m^-")

    sax["A"].set_ylim(bottom=0)
    sax["B"].set_ylim((0, 1.0))
    sax["B"].set_xlim(left=0)
    sax["D"].set_ylim(0, 1.05)
    sax["A"].set_ylabel("# Release Sites")
    sax["A"].set_xlabel("Input #")
    sax["B"].set_xlabel("# Release Sites")
    sax["B"].set_ylabel("Participation")
    sax["C"].set_ylabel("Input #")
    sax["C"].set_xlabel("Input #")
    sax["D"].set_xlabel(f"# Inputs in [{pre_w[0]:.1f} to {pre_w[1]:.1f}] before spike")

    PH.cleanAxes(PSum.axarr.ravel())
    # PH.talbotTicks(sax["C"])
    PH.talbotTicks(sax["A"])

    PH.talbotTicks(sax["B"], tickPlacesAdd={"x": 0, "y": 1}, floatAdd={"x": 0, "y": 2})
    PH.talbotTicks(sax["D"], tickPlacesAdd={"x": 0, "y": 1}, floatAdd={"x": 0, "y": 2})
    PH.talbotTicks(
        axcbar, tickPlacesAdd={"x": 0, "y": 2}, floatAdd={"x": 0, "y": 2}, pointSize=7
    )
    PSum.figure_handle.suptitle(
        f"Cell {gbc:s}", fontsize=12
    )
    mpl.show()

    return (summarySiteTC, RCD.sites)


def plot_tuning(args, filename=None, filenames=None):
    PD = PData()
    changetimestamp = get_changetimestamp()
    channel = "nacncoop"
    cols = 3
    rows = 2
    sizex = cols * 3
    sizey = rows * 2.5
    P = PH.regular_grid(
        rows,
        cols,
        order="rowsfirst",
        figsize=(sizex, sizey),
        # panel_labels=plabels,
        labelposition=(0.05, 0.95),
        margins={
            "leftmargin": 0.1,
            "rightmargin": 0.01,
            "topmargin": 0.15,
            "bottommargin": 0.15,
        },
    )

    print("...")

    truefile = {"Passive": "passive", "Canonical": "normal", "Active": "active"}
    for ic, chdist in enumerate(["Passive", "Canonical", "Active"]):
        args.dendritemode = truefile[chdist]
        if filename is None:
            fna = select_filenames(filenames, args)
            print("\n plot data from: ".join([str(f) for f in fna]))
        else:
            fna = [Path(filename)]  # just one file

        for k, fn in enumerate(fna):
            plot_traces(P.axarr[0, k], fn, PD, args.protocol)

        P.axarr[0, ic].text(50, 1.0, chdist, horizontalalignment="center")
        basen = fn.parts[-4]

        pngfile = Path(PD.renderpath, f"{basen:s}_{channel:s}_{truefile[chdist]:s}.png")
        print(pngfile)
        imaged = mpl.imread(pngfile)
        P.axarr[1, ic].imshow(imaged)


def plot_AN_response(P, fn, PD, protocol):
    changetimestamp = get_changetimestamp()
    x = get_data_file(fn, changetimestamp, PD)
    mtime = Path(fn).stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d-%H:%M")
    if x is None:
        print("No simulation found that matched conditions")
        print(fn)
        return
    # unpack x
    par, stitle, ivdatafile, filemode, d = x
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
    ntr = len(AR.traces)  # number of trials
    v0 = -160.0
    trstep = 25.0 / ntr
    inpstep = 5.0 / ntr
    sz = 50.0 / ntr
    # print(dir(AR))
    si = d["Params"]
    ri = d["runInfo"]
    # print(si)
    # print(ri)
    # print(isinstance(si, dict))
    # print(isinstance(ri, dict))
    if isinstance(si, dict):
        totaldur = (
            si["pip_start"]
            + np.max(si["pip_start"])
            + si["pip_duration"]
            + si["pip_offduration"]
        )
        soundtype = si['soundtype']
        pip_start = si["pip_start"]
        pip_duration = si["pip_duration"]
        F0 = si['F0']
        dB = si["dB"],
        fmod = si["fmod"],
        dmod = si["dmod"],
    else:
        totaldur = (ri.pip_start + np.max(ri.pip_start) + ri.pip_duration + ri.pip_offduration)
        pip_start = ri.pip_start
        pip_duration = ri.pip_duration
        soundtype = ri.soundtype
        F0 = ri.F0
        dB = ri.dB
        fmod = ri.fmod
        dmod = ri.dmod
    # print(d['Results'][0].keys())
    # print(dir(AR))

    ntr = len(AR.traces)
    allst = []
    for i in range(ntr):  # for all trails in the measure.
        trial = AR.traces[i]
        trd = d['Results'][i]
        w = trd["stimWaveform"]
        stb = trd["stimTimebase"]
        P.axdict["A"].plot(
            AR.time / 1000.0, AR.traces[i], linewidth=0.5
        )
        P.axdict["B"].plot(stb, w, linewidth=0.5)  # stimulus underneath
        P.axdict["C"].plot(
            trd["spikeTimes"] / 1000.0,
            i * np.ones(len(trd["spikeTimes"])),
            "o",
            markersize=2.5,
            color="b",
        )
        inputs = len(trd["inputSpikeTimes"])
        for k in range(inputs):
            tk = trd["inputSpikeTimes"][k] / 1000.0
            y = (i + 0.1 + k * 0.05) * np.ones(len(tk))
            P.axdict["E"].plot(tk, y, "|", markersize=2.5, color="r", linewidth=0.5)

        allst.extend(trd["spikeTimes"] / 1000.0)
    allst = np.array(allst)
    allst = np.sort(allst)
    # the histogram of the data
    ax = P.axdict["F"]
    if len(allst) > 0:
        P.axdict["F"].hist(allst, 100, density=True, facecolor="blue", alpha=0.75)
    else:
        P.axdict["F"].text(0.5, 0.5, "No Spikes", fontsize=14, 
        color='r', transform=ax.transAxes, horizontalalignment="center")
    
    for a in ["A", "B", "C", "E", "F"]:  # set some common layout scaling
        P.axdict[a].set_xlim((0.0, totaldur))
        P.axdict[a].set_xlabel("T (s)")

    if soundtype == "SAM":  # calculate vs and plot histogram
        if len(allst) == 0:
            P.axdict["D"].text(0.5, 0.5, "No Spikes", fontsize=14, color='r', 
            transform=ax.transAxes, horizontalalignment="center")
        else:
                  
            # combine all spikes into one array, plot PSTH
        # allst = []
       #  for trial in range(ntr):
       #      trd = d['Results'][trial]
       #      allst.extend(trd["spikeTimes"] / 1000.0)
       #  allst = np.array(allst)
       #  allst = np.sort(allst)
       #  # the histogram of the data
       #  if len(allst) > 0:
       #      P.axdict["F"].hist(allst, 100, density=True, facecolor="blue", alpha=0.75)
       #  else:
       #      P.axdict["F"].text(0.5, 0.5, "No Spikes", fontsize=14, color='r')
        # print('allst: ', allst)
            phasewin = [
                pip_start[0] + 0.25 * pip_duration,
                pip_start[0] + pip_duration,
            ]
            # print (phasewin)
            spkin = allst[np.where(allst > phasewin[0])]
            spikesinwin = spkin[np.where(spkin <= phasewin[1])]

            # set freq for VS calculation
            if soundtype == "tone":
                print(
                    "Tone: F0=%.3f at %3.1f dbSPL, cell CF=%.3f"
                    % (F0, dB, F0)
                )
            if soundtype == "SAM":
                tstring = (
                    "SAM Tone: F0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f"
                    % (
                        F0,
                        dB,
                        fmod,
                        dmod,
                        F0,
                    )
                )
                # print(tstring)
                P.figure_handle.suptitle(tstring, fontsize=10)
            # print('spikes: ', spikesinwin)
            vs = vector_strength(spikesinwin * 1000.0, fmod)  # vs expects spikes in msec
            print(" Sound type: ", soundtype)
            print(
                "AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
                % (F0, vs["r"], vs["d"] * 1e6, vs["R"], vs["p"], vs["n"])
            )
            # print(vs['ph'])
            P.axdict["D"].hist(vs["ph"], bins=2 * np.pi * np.arange(30) / 30.0)
            P.axdict["D"].set_xlim((0.0, 2 * np.pi))
            P.axdict["D"].set_title(
                "Phase (VS = {0:.3f})".format(vs["r"]),
                fontsize=9,
                horizontalalignment="center",
            )

    # P.figure_handle.savefig("AN_res.pdf")
#     plt.show()


def select_filenames(fng, args) -> Union[list, None]:
    # print(fng)
    print("Selct filename for experiment: ", args.experiment)
    match_index = []
    for ix, fn in enumerate(fng):
        ematch = True
        if args.experiment is not None:
            if str(fn).find(args.experiment) < 0:
                ematch = False

        dmatch = True
        if args.dbspl is not None:
            srch = f"_{int(args.dbspl):03d}dB_"
            if str(fn).find(srch) < 0:
                dmatch = False

        rmatch = True
        if args.nreps is not None:
            srch = f"_{int(args.nreps):03d}_"
            if str(fn).find(srch) < 0:
                rmatch = False

        mmatch = True
        if args.dendriteMode != "normal":
            srch = f"_{args.dendriteMode:s}_"
            if str(fn).find(srch) < 0:
                mmatch = False

        if dmatch and ematch and rmatch and mmatch:
            match_index.append(ix)

    fng = [fng[i] for i in match_index]
    if len(fng) == 0:
        print("No Files matching conditions found")
        return None
    # print('found (filtered):\n', '\n   '.join([str(f) for f in fng]))
    return fng


def build_parser():
    """
    create the command line parser and process the line
    """

    parser = argparse.ArgumentParser(description="Plot GBC results")
    parser.add_argument(
        dest="cell",
        action="store",
        nargs="+",
        type=str,
        default=None,
        help="Select the cell(s) or 'A' for all(no default)",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        dest="protocol",
        action="store",
        default="IV",
        choices=runtypes,
        help=("Select the protocol (default: IV) from: %s" % runtypes),
    )
    parser.add_argument(
        "-a",
        "--analysis",
        dest="analysis",
        action="store",
        default=analysistypes[0],
        choices=analysistypes,
        help=(
            f"Select the analysis type (default: {analysistypes[0]:s}) from: {str(analysistypes):s}"
        ),
    )
    parser.add_argument(
        "-M",
        "--modeltype",
        dest="modeltype",
        action="store",
        default="XM13_nacncoop",
        help=("Select the model type (default XM13_nacncoop) from: %s " % modeltypes),
    )

    parser.add_argument(
        "-s",
        "--scaled",
        dest="scaled",
        action="store_true",
        help=("use scaled data or not"),
    )

    parser.add_argument(
        "-e",
        "--experiment",
        dest="experiment",
        action="store",
        default="delays",
        choices=experimenttypes,
        help=("Select the experiment type from: %s " % experimenttypes),
    )

    parser.add_argument(
        "--dendritemode",
        dest="dendriteMode",
        default="normal",
        choices=dendriteChoices,
        help="Choose dendrite table (normal, active, passive)",
    )

    parser.add_argument(
        "-d",
        "--dB",
        dest="dbspl",
        type=float,
        action="store",
        default=None,
        help=("Select the models at specific intensity"),
    )

    parser.add_argument(
        "-r",
        "--nreps",
        dest="nreps",
        type=int,
        action="store",
        default=None,
        help=("Select the models with # reps"),
    )

    parser.add_argument(
        "-c",
        "--check",
        dest="check",
        action="store_true",
        help=("Just check selection criteria and return"),
    )

    return parser


def cmdline_display(args, PD):
    """
    Display analysis and traces
    when called from the commandline
    """
    args.protocol = args.protocol.upper()
    changetimestamp = get_changetimestamp()
    # PD.gradeA = [cn for cn in args.cell]
    print("args.cell: ", args.cell)
    if args.cell[0] in ["A", "a"]:
        pass  # (all grade a cells defined in pd dataclass)
    else:
        PD.gradeA = [int(c) for c in args.cell]
    rows, cols = PH.getLayoutDimensions(len(PD.gradeA))
    plabels = [f"VCN_c{g:02d}" for g in PD.gradeA]
    for i in range(rows * cols - len(PD.gradeA)):
        plabels.append(f"empty{i:d}")

    sizex = cols * 3
    sizey = rows * 2.5
    if rows * cols < 4:
        sizex *= 2
        sizey *= 2

    chan = "_"  # for certainity in selection, include trailing underscore in this designation

    modelName = args.modeltype
    print("Scaling: ", args.scaled)
    stitle = "unknown scaling"
    if args.scaled:
        PD.soma_inflate = True
        PD.dend_inflate = True
        stitle = "scaled"
    else:
        PD.soma_inflate = False
        PD.dend_inflate = False
        stitle = "notScaled"

    for ig, gbc in enumerate(PD.gradeA):
        basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
        pgbc = f"VCN_c{gbc:02d}"
        if args.protocol == "IV":
            name_start = f"IV_Result_VCN_c{gbc:02d}_*.p"
            args.experiment = None
        elif args.protocol == "AN":
            name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"

        print(f"Searching for:  {str(Path(basefn, name_start)):s}")
        fng = list(Path(basefn).glob(name_start))
        print(f"Found: {len(fng):d} files.")
        """ cull list by experiment and db """

        if args.check:
            return
        print("\nConditions: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)

        if args.analysis in ["traces", "revcorr"]:
            fng = select_filenames(fng, args)
            times = np.zeros(len(fng))
            for i, f in enumerate(fng):
                times[i] = f.stat().st_mtime
            # pick most recent file = this should be better managed (master file information)

            ix = np.argmax(times)
            fng = [fng[ix]]

            if ig == 0:
                P = PH.regular_grid(
                    rows,
                    cols,
                    order="rowsfirst",
                    figsize=(sizex, sizey),
                    panel_labels=plabels,
                    labelposition=(0.05, 0.95),
                    margins={
                        "leftmargin": 0.1,
                        "rightmargin": 0.01,
                        "topmargin": 0.15,
                        "bottommargin": 0.15,
                    },
                )

            ax = P.axdict[pgbc]
            for fn in fng:
                if args.analysis == "traces":
                    plot_traces(ax, fn, PD, args.protocol)
                elif args.analysis == "revcorr":
                    res = compute_revcorr(
                        ax, pgbc, fn, PD, args.protocol
                    )
                    if res is None:
                        return
        elif args.analysis == "singles":
            fna = select_filenames(fng, args)
            print(fna)

        elif args.analysis == "tuning":
            plot_tuning(args, filename=None, filenames=fng)
    if chan == "_":
        chan = "nav11"
    else:
        chan = chan[:-1]
    P.figure_handle.suptitle(
        f"Model: {modelName:s}  Na Ch: {chan:s} Scaling: {stitle:s} ", fontsize=7
    )
    mpl.show()


def getCommands():
    parser = build_parser()
    args = parser.parse_args()

    # if args.configfile is not None:
    #     config = None
    #     if args.configfile is not None:
    #         if ".json" in args.configfile:
    #             # The escaping of "\t" in the config file is necesarry as
    #             # otherwise Python will try to treat is as the string escape
    #             # sequence for ASCII Horizontal Tab when it encounters it
    #             # during json.load
    #             config = json.load(open(args.configfile))
    #             print(f"Reading JSON configuration file: {args.configfile:s}")
    #         elif ".toml" in args.configfile:
    #             print(f"Reading TOML configuration file: {args.configfile:s}")
    #             config = toml.load(open(args.configfile))
    #
    #     vargs = vars(args)  # reach into the dict to change values in namespace
    #
    #     for c in config:
    #         if c in args:
    #             # print("Getting parser variable: ", c)
    #             vargs[c] = config[c]
    #         else:
    #             raise ValueError(
    #                 f"config variable {c:s} does not match with comand parser variables"
    #             )
    #
    #     print("   ... All configuration file variables read OK")
    # now copy into the dataclass
    PD = PData()
    # PD_names = dir(PD)
    # for key, value in vargs.items():
    #      if key in parnames:
    #          # print('key: ', key)
    #          # print(str(value))
    #          exec(f"PD.{key:s} = {value!r}")
    #      # elif key in runnames:
    #      #     exec(f"runinfo.{key:s} = {value!r}")
    cmdline_display(args, PD)
    return (args, PD)


if __name__ == "__main__":
    getCommands()
