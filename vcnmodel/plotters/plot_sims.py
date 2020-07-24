import argparse
import dataclasses
import datetime
import pickle
from dataclasses import dataclass, field
from collections import OrderedDict
from pathlib import Path
from typing import Union, List, Tuple

import numpy as np
import scipy.stats
from numba import jit
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
from vcnmodel import sttc as STTC
from lmfit import Model
# try using a standard ccf
# import neo
# from quantities import s, Hz, ms
# from elephant.spike_train_generation import homogeneous_poisson_process
# import elephant.spike_train_correlation as ESTC
# from elephant.conversion import BinnedSpikeTrain

from vcnmodel import correlation_calcs as CXC
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

modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop", "XM13A_nacncoop"]
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
    algorithm: str = "RevcorrSPKS"  # save the algorithm name
    # clip trace to avoid end effects
    min_time: float = 10.0  # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time: float = 250.0  # this window needs to be at least as long as maxwin
    binw: float = 0.1
    minwin: float = -5
    maxwin: float = 2.5
    amax: float = 0.0


@dataclass
class RevCorrData:
    C: list = field(default_factory=def_empty_list)  # using Brian 1.4 correlation
    CB: list = field(default_factory=def_empty_list)  # using elephant/neo, binned correlation
    CBT: list = field(default_factory=def_empty_list)  # using elephant/neo, binned correlation
    TC: float = 0. # list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)
    tx: np.array = field(default_factory=def_empty_np)
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(default_factory=def_empty_np)
    sv_all: np.array = field(default_factory=def_empty_np)
    sv_avg: np.array = field(default_factory=def_empty_np)
    sites: np.array = field(default_factory=def_empty_np)
    nsp_avg: int = 0
    npost_spikes: int = 0
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
    # print(PD)
    # print(par)
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
    if protocol in ["IV", "runIV"]:
        protocol = "IV"
    elif protocol in ["VC", "runVC"]:
        protocol = "VC"
    print('Protocol: ', protocol)
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
        AR.traces[trial][0] = AR.traces[trial][1]
        if protocol in ["VC", "vc", "vclamp"]:
            AR.traces[trial] = AR.traces[trial].asarray() * 1e6
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
        elif protocol in ["VC", "vc", "vclamp"]:
            ax.set_ylim((-100.0, 100.0))
        else:
            ax.set_ylim(-200.0, 50.0)
            ax.set_xlim(0.080, np.max(AR.time_base))
    # print('Nout/Nin: ', float(noutspikes)/ninspikes)

    ftname = str(ivdatafile.name)
    ip = ftname.find("_II_") + 4
    ftname = ftname[:ip] + "...\n" + ftname[ip:]
    toptitle = f"{ftname:s}"
    if protocol in ["IV"]:
        toptitle += f"\nRin={RMA['Rin']:.1f} M$\Omega$  $\\tau_m$={RMA['taum']:.2f} ms"

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
    elif protocol in ["VC", "vc", "vclamp"]:
        PH.calbar(
            ax,
            calbar=[20.0, -10.0, 10.0, 10.0],
            unitNames={"x": "ms", "y": "nA"},
            fontsize=9,
        )    
    else:
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
    ax.set_title(toptitle, y=1.0, fontsize=8, verticalalignment="top")
    return(synno, noutspikes, ninspikes)


def boltzI(x, gmax, vhalf, k, E):
    return(gmax*(x-E)*(1./(1.0+np.exp(-(x-vhalf)/k))))
    # return(gmax*(1./(1.0+np.exp(-(x-vhalf)/k))))


def analyzeVC(
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
    if protocol in["VC", "runVC"]:
        protocol = "VC"
    print('Protocol: ', protocol)
    par, stitle, ivdatafile, filemode, d = x
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
        # for j in range(0, ntraces):
        # if verbose:
        #     print('    analyzing trace: %d' % (j))
        # vss[j] = np.mean(V[j,tss[0]:tss[1]])  # steady-state voltage
        # ic[j] = np.mean(I[j,tss[0]:tss[1]])  # corresponding currents
        # vm[j] = np.mean(V[j, 0:int((ts-1.0)/dt)])  # resting potential - for 1 msec prior to step
        # ax.plot(t, I[j])
        # if verbose:
        #     print('   >>> completed analyzing trace %d' % j)
    # print(dir(AR))
    # print(AR.tstart)
    # print(AR.tstart_tdur)
    # print(AR.tdur)
    # print(AR.tend)
    # print(AR.sample_rate)
    # print(AR.commandLevels)
    tss = [0,0]
    tss[0] = int(AR.tstart+(AR.tdur/4.0)/AR.sample_rate[0])  # wrong duration stored in traces - need to fix.
    tss[1] = int((AR.tstart + AR.tdur*0.5)/AR.sample_rate[0])
    I = AR.traces.asarray()
    V = AR.cmd_wave
    ntraces = np.shape(V)[0]
    # initialize all result arrays, lists and dicts
    vss = np.empty(ntraces)
    vmin = np.zeros(ntraces)
    vrmss = np.zeros(ntraces)
    vm = np.zeros(ntraces)
    ic = np.zeros(ntraces)
    # print(np.array(tss)*AR.sample_rate[0])
    for j in range(0, ntraces):
        vss[j] = np.mean(V[j,tss[0]:tss[1]])  # steady-state voltage
        ic[j] = np.mean(I[j,tss[0]:tss[1]])  # corresponding currents
        vm[j] = np.mean(V[j, 0:int((AR.tstart-1.0)/AR.sample_rate[0])])  # resting potential - for 1 msec prior to step
        ax[0].plot(AR.time, V[j,:], 'k-')

    print(vss)
    # now fit traces to variation of g = gmax * I/(V-Vr)
    gmodel = Model(boltzI)
    # def boltzI(self, x, gmax, vhalf, k, E):
#             return(gmax*(1./(np.exp((x-E)/k))))
#
    vss = np.array(vss)*1e3
    # transform to g
    Ek = -84.
    gss = ic/(vss-Ek)
    ax[1].plot(vss, ic, 'ko-', markersize=2)
    print(ic)
    gmodel.set_param_hint('gmax', value=2.e-10, min=0.,vary=True)
    gmodel.set_param_hint('vhalf', value=-38., min=-90., max=90.)
    gmodel.set_param_hint('k', value=7, min=0.05, max=200.0, vary=True)
    gmodel.set_param_hint('E', value=Ek, vary=False)
    gparams = gmodel.make_params()
    weights = np.ones_like(vss)
    for i, w in enumerate(weights):
        if vss[i] <= -100:
            weights[i] = 0.
    result = gmodel.fit(ic, method='nedler', params=gparams, x=vss, weights = weights, max_nfev=10000)
    print(result.fit_report())
    print(vss)
    print(result.best_fit)
    ax[1].plot(vss, result.best_fit, 'r-')
    mpl.show()
    
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
    
def make_patch_spines_invisible(axn):
    axn.set_frame_on(True)
    axn.patch.set_visible(False)
    for sp in axn.spines.values():
        sp.set_visible(False)

def plot_revcorr2(ax: object, PD: dataclass, RCP: dataclass, RCD: dataclass):
    seaborn.set_style("ticks")
    # secax = twinax(P.figure_handle, ax, pos=maxwin)
    secax = PLS.create_inset_axes([0, 0.5, 1, 0.5], ax)
    PH.noaxes(secax, "xy")

    secax.set_facecolor((1, 1, 1, 0))
    secax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # ax.set_facecolor((0.7, 0.7, 0.7))
    ax2 = ax.twinx()
    ax2.spines['right'].set_position(("axes", 0.9))
    summarySiteTC = {}
    maxrevcorr = 0
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
        refzero = int(RCP.minwin/RCP.binw)
        if RCD.C[isite] is not None:
            # print('plotting C')
            nc = int(len(RCD.C[isite]) / 2)
            RCD.TC = RCD.TC / len(RCD.st)
            summarySiteTC[isite] = RCD.TC
            # print(RCD.sites, isite)
            color = mpl.cm.viridis(norm(RCD.sites, isite))
            # print('color', color)
            maxrevcorr = np.max((maxrevcorr, np.max(RCD.CB[isite])))
            totalrevcorr = np.sum(RCD.CB[isite])


            make_patch_spines_invisible(ax2)
            if isite in range(RCP.ninputs):
                if RCP.algorithm == "RevcorrSPKS":
                    ax.plot(
                        RCD.tx,
                        RCD.C[isite][:nc],
                        color=color,
                        label=("Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))),
                        linewidth=1.5,
                        zorder=5,
                    )

                # print(RCD.npost_spikes)
     #            print(RCD.CB[isite])
                elif RCP.algorithm == "RevcorrSimple":
                    ax.plot(
                        RCD.CBT - RCP.maxwin, 
                        RCD.CB[isite]/float(RCD.npost_spikes),
                        color=color,
                        label=("Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))),
                        linewidth=1,
                        zorder=5,
                        alpha=1.0
                    )
                elif RCP.algorithm  == "RevcorrSTTC": # use spike time tiling method
                   ax.plot(
                       RCD.tx,
                       RCD.STTC[isite][:nc],
                       color=color,
                       label=("Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))),
                       linewidth=1.5,
                       zorder=5,
                   )


            RCD.max_coin_rate = np.max((RCD.max_coin_rate, np.max(RCD.C[:nc])))

        if (
            isite == 0
        ):  # only plot the first time through - the APs are the same no matter the input
            for t in sel:
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
    print(f"Total spikes plotted: {RCD.nsp_avg:d}")

    secax.plot([0.0, 0.0], [-120.0, 10.0], "r", linewidth=0.5)

    seaborn.despine(ax=ax)
    secax.set_ylabel("Rate of coincidences/bin (Hz)", fontsize=10)
    ax.set_xlabel("T (ms)", fontsize=10)
    ax.set_xlim((RCP.minwin, RCP.maxwin))
    if 0.2 < RCD.max_coin_rate < 1.0:
        ax.set_ylim(0, 1)
    elif 0.2 >= RCD.max_coin_rate:
        ax.set_ylim(0, 0.25)
    elif RCD.max_coin_rate > 1.0:
        ax.set_ylim(0, 2.0)
    yls = ax.get_ylim()
    # ax2.set_ylim(yls)
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

from numba import jit

@jit(parallel=True, cache=True,)
def nb_revcorr(st1, st2, binwidth, corrwindow):
    xds  = np.zeros(int((corrwindow[1]-corrwindow[0])/binwidth))
    print(xds.shape)
    # [None]*len(st1)
    for i, sp in enumerate(st1):
        diff = st2 - sp
        v = np.where((corrwindow[0] <= diff) & (diff <= corrwindow[1]))
        # print('v: ', v)
       #  if len(v) > 0:
       #      iv = [int(vx/binwidth) for vx in v]
       #      xds[iv] = 1
       #      print('xds: ', xds)# for d in diff:
        #     if corrwindow[0] <= d <= corrwindow[1]:
        #         xds[i] = d
    # print(np.array(xds))
    # print(np.array(xds).ravel().shape)
    return xds, len(st1)  # return the n postsynaptic spikes 



def revcorr(
    st1: Union[np.ndarray, List]= None,
    st2: Union[np.ndarray, List]= None,
    binwidth:float=0.1,
    datawindow:Union[List, Tuple]= [0., 100.], # time window for selecteing spikes
    corrwindow:Union[List, Tuple]= [-5., 1.], # time window to examine correlation relative to st1
    ) -> np.ndarray:
    if st1 is None or st2 is None:
        raise ValueError("coincident_spikes_correlation: reference and comparator must be defined")
    refa = st1 # [a for a in st1 if (datawindow[0] <= a <= datawindow[1])]
    refb = st2 # [b for b in st2 if (datawindow[0] <= b <= datawindow[1])]    
    # xds  = [None]*len(refa)
    xds  = np.zeros(int((corrwindow[1]-corrwindow[0])/binwidth))
    for i, sp in enumerate(refa):
        diff = refb - sp
        v = diff[np.where((corrwindow[0] <= diff) & (diff <= corrwindow[1]))]
        # print('v: ', v)
        if len(v) > 0:
            indxs = [int(vx/binwidth) for vx in v]
            xds[indxs] = xds[indxs] + 1
        # for d in diff:
        #     if corrwindow[0] <= d <= corrwindow[1]:
        #         xds.append(d)
    # print(np.array(xds))
    # print(np.array(xds).ravel().shape)
    return xds, len(st1)  # return the n postsynaptic spikes 


def compute_revcorr(
    ax: object,
    gbc: str,
    fn: Union[str, Path],
    PD: object,
    protocol: str,
    revcorrtype: str = 'RevcorrSPKS',
    thr: float = -20.0,
    width: float = 4.0,
    ) -> Union[None, tuple]:

    changetimestamp = get_changetimestamp()

    #
    # 1. Gather data
    #
    print('Getting data')
    SC, syninfo = get_synaptic_info(gbc)

    res = get_data(fn, PD, changetimestamp, protocol)
    if res is None:
        return None
    (d, AR, SP, RMA, RCP, RCD) = res
    RCP.algorithm = revcorrtype
    if RCP.algorithm == 'RevcorrSTTC':
        sttccorr = STTC.STTC()
        
    print('Preparing for computation')
    RCP.ninputs = len(syninfo[1])
    RCD.sites = np.zeros(RCP.ninputs)
    for isite in range(RCP.ninputs):  # precompute areas
        area = syninfo[1][isite][0]
        if area > RCP.amax:
            RCP.amax = area
        RCD.sites[isite] = int(np.around(area * SC.synperum2))

    print("ninputs: ", RCP.ninputs)
    maxtc = 0

    #
    # 2. set up parameters
    #
    si = d["Params"]
    ri = d["runInfo"]
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
    RCD.CB = [None] * RCP.ninputs
    RCD.STTC = [None] * RCP.ninputs
    RCD.max_coin_rate = 0.0
    RCD.npost_spikes = 0  # count spikes used in analysis
    RCD.nsp_avg = 0
    nspk_plot = 0
    # spksplotted = False
    RCP.min_time = min_time  # driven window without onset
    RCP.max_time = max_time

    #
    # 3. Prepare storage arrays
    #
    nbins = int((max_time-min_time)/RCP.binw)
    
    RCD.CBT = np.arange(RCP.minwin, RCP.maxwin, RCP.binw) # refcorr time base
    for isite in range(RCP.ninputs):
        RCD.CB[isite] = np.zeros_like(RCD.CBT)
        ncpts = len(np.arange(RCP.minwin, -RCP.minwin, RCP.binw))
        RCD.C[isite] = np.zeros(ncpts)
        RCD.STTC[isite] = np.zeros(ncpts)

    #
    # sum across trials, and sort by inputs
    #
    print('starting loop')
    start_time = datetime.datetime.now()
    for trial in range(RCP.ntrials):  # sum across trials
        stx = AR.time_base[SP.spikeIndices[trial]]
        stx = stx[  # get postsynaptic spikes and trim to analysis window
            (stx > RCP.min_time) & (stx < RCP.max_time)
        ]  
        if len(stx) == 0:
            continue
        RCD.npost_spikes += len(stx)

        # accumulate spikes and calculate average spike
        for n in range(len(stx)):
            reltime = np.around(RCD.ti, 5) - np.around(stx[n], 5)
            areltime = np.argwhere(
                (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
            ).squeeze()
            
            if RCD.nsp_avg == 0: # init arrays
                RCD.sv_avg = d["Results"][trial]["somaVoltage"][areltime]
                RCD.nsp_avg = 1
                RCD.ti_avg = RCD.ti[0 : len(areltime)] + RCP.minwin
                RCD.sv_all = np.zeros((RCP.ntrials, RCD.ti_avg.shape[0]))  # initialize the array
            else:
                if len(areltime) > len(RCD.sv_avg):
                    areltime = areltime[0 : len(RCD.sv_avg)]
                if len(areltime) < len(RCD.sv_avg):
                    nextend = len(RCD.sv_avg) - len(areltime)
                    areltime = np.append(
                        areltime,
                        np.arange(areltime[-1] + 1, areltime[-1] + nextend + 1),
                    )
                RCD.sv_avg += d["Results"][trial]["somaVoltage"][areltime]
                RCD.nsp_avg += 1
            RCD.sv_all[trial] =  d["Results"][trial]["somaVoltage"][areltime]


            nspk_plot += RCD.nsp_avg
            RCD.sv_sites.append(RCD.sv_all)

        # Now get reverse  correlation for each input        
        for isite in range(RCP.ninputs): # for each ANF input
            anx = d["Results"][trial]["inputSpikeTimes"][isite] # get input an spikes and trim
            anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
            if len(anx) == 0:
                continue
            
            if revcorrtype == "RevcorrSPKS":
                RCD.C[isite] += SPKS.correlogram(
                    stx, anx, width=-RCP.minwin, binwidth=RCP.binw, T=None
                    )
           
            elif revcorrtype == "RevcorrSimple":
                refcorr, npost = revcorr(stx, anx,
                    binwidth=  RCP.binw,
                    datawindow = [RCP.min_time, RCP.max_time],
                    corrwindow=[RCP.minwin, RCP.maxwin]
                    )
                RCD.CB[isite] = RCD.CB[isite] + refcorr

            elif revcorrtype == "RevcorrSTTC":
                sttccorr.set_spikes(RCP.binw, stx, anx, RCP.binw)
                refcorr = sttccorr.calc_ccf_sttc(corrwindow=[RCP.minwin, RCP.maxwin], binwidth=RCP.binw)
                RCD.STTC[isite] = RCD.STTC[isite] + refcorr
                 
            tct = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
            if ~np.isnan(tct):
                RCD.TC += tct

        if isinstance(RCD.TC, float) and RCD.TC > maxtc:
            maxtc = RCD.TC
        else:
            maxtc = 1.0
    if RCD.nsp_avg > 0:
        RCD.sv_avg /= RCD.nsp_avg
    
    elapsed_time = datetime.datetime.now() - start_time
    print('Time for calculation: ', elapsed_time)
    pre_w = [-2.7, -0.7]

    summarySiteTC = plot_revcorr2(ax, PD, RCP, RCD)

    ax.set_title(
        f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}",
        fontsize=11
    )
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
    PSum = PH.regular_grid(
        rows=2,
        cols=2,
        order="rowsfirst",
        figsize=(6, 6),
        showgrid=False,
        verticalspacing=0.1,
        horizontalspacing=0.1,
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.1,
            "rightmargin": 0.1,
            "topmargin": 0.15,
        },
        label=["A", "B", "C", "D"],
        labelposition=(-0.05, 1.05),
    )
    sax = PSum.axdict
    # f, sax = mpl.subplots(3,1)
    # f.set_size_inches( w=3.5, h=9)
    sax["A"].plot(np.arange(RCP.ninputs) + 1, RCD.sites, "bo")
    # print('pairwise: ', pairwise)
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
        # print("vmin, vmax: ", vmin, vmax)
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
        f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}",
        fontsize=11
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
            fna = filename  # just one file
        fna = [Path(fn) for fn in fna if str(fn).find("ASA=") < 0]
        print('fns: ', fna)
        for k, fn in enumerate(fna):
            plot_traces(P.axarr[0, k], fn, PD, args.protocol)

        P.axarr[0, ic].text(50, 1.0, chdist, horizontalalignment="center")
        basen = fn.parts[-4]

        pngfile = Path(PD.renderpath, f"{basen:s}_{channel:s}_{truefile[chdist]:s}.png")
        # print(pngfile)
        imaged = mpl.imread(pngfile)
        P.axarr[1, ic].imshow(imaged)
    return P

def setup_PSTH():
    sizer = OrderedDict(
        [
            ("A", {"pos": [0.08, 0.4, 0.81, 0.15]}),
            ("B", {"pos": [0.08, 0.4, 0.50, 0.27]}),
            ("C", {"pos": [0.08, 0.4, 0.36, 0.10]}),
            ("D", {"pos": [0.08, 0.4, 0.06, 0.24]}),
            # rhs
            ("E", {"pos": [0.55, 0.4, 0.75, 0.18]}),
            ("F", {"pos": [0.55, 0.4, 0.53, 0.18]}),
            ("G", {"pos": [0.55, 0.4, 0.30, 0.18]}),
            ("H", {"pos": [0.55, 0.4, 0.07, 0.18]}),
        ]
    )  # dict elements are [left, width, bottom, height] for the axes in the plot.
    n_panels = len(sizer.keys())
    gr = [
        (a, a + 1, 0, 1) for a in range(0, n_panels)
    ]  # just generate subplots - shape does not matter
    axmap = OrderedDict(zip(sizer.keys(), gr))
    P = PH.Plotter((n_panels, 1), order="columnsfirst",
        axmap=axmap, label=True, figsize=(8.0, 6.0))
    P.resize(sizer)  # perform positioning magic
    P.axdict["A"].set_ylabel("mV", fontsize=8)

    P.axdict["D"].set_title("Bushy Spike Raster", fontsize=9)
    P.axdict["D"].set_ylabel("Trial")

    P.axdict["B"].set_title("PSTH")
    P.axdict["B"].set_ylabel("Sp/sec")

    P.axdict["C"].set_title("Stimulus", fontsize=9)
    P.axdict["C"].set_ylabel("Amplitude (Pa)", fontsize=8)
    P.axdict["C"].set_xlabel("T (s)")

    P.axdict["E"].set_title("Phase", fontsize=8)
    P.axdict["F"].set_title("?", fontsize=8)
    P.axdict["G"].set_title("ANF Spike Raster", fontsize=9)
    P.axdict["H"].set_title("ANF PSTH", fontsize=9)

    return P


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
    # print('total, pipstart, pipdur: ', totaldur, pip_start, pip_duration)
    ntr = len(AR.traces)
    all_an_st = []
    all_bu_st = []
    for i in range(ntr):  # for all trials in the measure.
        vtrial = AR.traces[i]*1e3
        trd = d['Results'][i]
        w = trd["stimWaveform"]
        stb = trd["stimTimebase"]
        all_bu_st.extend(trd["spikeTimes"])
        # print('max bu spktimes: ', np.max(trd["spikeTimes"]))# Panel A: traces
        if i < 1:
            P.axdict["A"].plot(
                AR.time_base / 1000.0,
                vtrial,
                'k-',
                linewidth=0.5
                )
            P.axdict["A"].plot(
                AR.time_base[SP.spikeIndices[i]]/1000.,
                vtrial[SP.spikeIndices[i]],
                "ro",
                markersize=1.5,
            )
        if i == 0:
            P.axdict["C"].plot(stb, w, 'k-', linewidth=0.5)  # stimulus underneath
        P.axdict["D"].plot(
            trd["spikeTimes"] / 1000.0,
            i * np.ones(len(trd["spikeTimes"])),
            "|",
            markersize=1.5,
            color="b",
        )
        inputs = len(trd["inputSpikeTimes"])
        for k in range(inputs):
            tk =trd["inputSpikeTimes"][k]
            all_an_st.extend(tk)
            y = (i + 0.1 + k * 0.05) * np.ones(len(tk))
            if i % 10 == 0:
                P.axdict["G"].plot(tk/1000., y, "|", markersize=2.5, color="k", linewidth=0.5)
            # print('max an spktimes: ', np.max(tk))# Panel A: traces

    all_bu_st = np.sort(np.array(all_bu_st))/1000.
    all_an_st = np.sort(np.array(all_an_st))/1000.
    # the histogram of the data
    hbins = np.arange(0., np.max(AR.time/1000.), 1e-3)  # 0.5 msec bins
    print(all_bu_st.shape, np.max(all_bu_st), all_an_st.shape, np.max(all_an_st), np.max(hbins))
    if len(all_bu_st) > 0:
        P.axdict["B"].hist(all_bu_st, bins=hbins, facecolor="k", alpha=1)
        # n, bins, patch = mpl.hist(all_bu_st, hbins, facecolor="blue", alpha=1)
   #      print('Spikes in PSTH: ', np.sum(n))
   #      print(n)
   #      print(bins)
    else:
        P.axdict["B"].text(0.5, 0.5, "No Spikes", fontsize=14, 
            color='r', transform=P.axdict["C"].transAxes, horizontalalignment="center")
    if len(all_an_st) > 0:
        P.axdict["H"].hist(all_an_st, bins=hbins, facecolor="k", alpha=1)
    else:
        P.axdict["H"].text(0.5, 0.5, "No Spikes", fontsize=14, 
            color='r', transform=P.axdict["H"].transAxes, horizontalalignment="center")
     
    for a in ["A", "B", "C", "D", "F", "G", "H"]:  # set some common layout scaling
        P.axdict[a].set_xlim((0.0, totaldur))
        # P.axdict[a].set_xlabel("T (s)")

    if soundtype == "SAM":  # calculate vs and plot histogram
        if len(all_bu_st) == 0:
            P.axdict["H"].text(0.5, 0.5, "No Spikes", fontsize=14, color='r', 
            transform=P.axdict["H"].transAxes, horizontalalignment="center")
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
            spkin = all_bu_st[np.where(all_bu_st > phasewin[0])]
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
            P.axdict["E"].hist(vs["ph"], bins=2 * np.pi * np.arange(30) / 30.0)
            P.axdict["E"].set_xlim((0.0, 2 * np.pi))
            P.axdict["E"].set_title(
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
            name_start = f"IV_Result_VCN_c{gbc:02d}_inp=self_{modelName:s}*.p"
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
            P = plot_tuning(args, filename=None, filenames=fng)
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
