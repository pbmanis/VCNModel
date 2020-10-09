import argparse
import dataclasses
import datetime
import functools
import operator
import pickle
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Union

import lmfit
import matplotlib.colorbar  # type: ignore
import matplotlib.colors  # type: ignore
import numpy as np  # type: ignore
import pyperclip
import pyqtgraph as pg  # type: ignore
import scipy.stats  # type: ignore
import seaborn
# from cnmodel.util import vector_strength
from ephys.ephysanalysis import MakeClamps, RmTauAnalysis, SpikeAnalysis
from lmfit import Model  # type: ignore
from matplotlib import pyplot as mpl  # type: ignore
from matplotlib import rc  # type: ignore
from numba import jit  # type: ignore
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as PLS
from pylibrary.tools import cprint as CP
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

import vcnmodel.model_params
from vcnmodel import analysis as SPKANA
from vcnmodel import cell_config as cell_config
# from vcnmodel import correlation_calcs as CXC
from vcnmodel import spikestatistics as SPKS
from vcnmodel import sttc as STTC

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


SpirouChoices = [
    "all",
    "max=mean",
    "all=mean",
    "removelargest",
    "largestonly",
    "twolargest",
    "threelargest",
    "fourlargest",
]


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


@dataclass
class SpikeData:
    """
    Data class to hold information about each spike
    """

    trial: int = -1  # the trial the spike came from
    time_index: int = 0  # index into the time array for this spike
    dt: float = 0.025  # sample interval, msec
    start_win: float = -5.0
    end_win: float = 5.0
    waveform: np.array = None  # the waveform of this postspike, clipped
    prespikes: np.array = None  # time indices to pre spikes


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
    CB: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    CBT: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    TC: float = 0.0  # list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)
    tx: np.array = field(default_factory=def_empty_np)
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(default_factory=def_empty_np)
    sv_all: np.array = field(default_factory=def_empty_np)
    sv_avg: np.array = field(default_factory=def_empty_np)
    sites: np.array = field(default_factory=def_empty_np)
    nsp_avg: int = 0
    npost_spikes: int = 0
    npre_spikes: int = 0
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


def boltzI(x, gmax, vhalf, k, E):
    return gmax * (x - E) * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def boltzG(x, gmax, vhalf, k, E):
    return gmax * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def expdecay(x, decay, amplitude, offset):
    return offset + amplitude * np.exp(-x / decay)


def exp2decay(x, a0, a1, tau0, tau1, offset):
    return offset + a0 * np.exp(-x / tau0) + a1 * np.exp(-x / tau1)


from numba import jit


@jit(
    parallel=True, cache=True,
)
def nb_revcorr(st1, st2, binwidth, corrwindow):
    xds = np.zeros(int((corrwindow[1] - corrwindow[0]) / binwidth))

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


def winprint(func):
    """
    Wrapper decorator for functions that print to the text area
    Clears the print area first,
    and puts a line of '*' when the function returns
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        self.textclear()
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        self.textappend("*" * 80)
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print

def winprint_continuous(func):
    """
    Wrapper decorator for functions that print to the text area
    Clears the print area first,
    and puts a line of '*' when the function returns
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print
    

def time_func(func):
    """
    Decorator to ime functions.
    Place inside (after) winprint if using
    Output is to terminal.
    """

    @functools.wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        print(f"Starting : {func.__name__!r}")
        start_time = time.perf_counter()  # 1
        value = func(self, *args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_timer


class PlotSims:
    def __init__(self, parent):
        self.parent = (
            parent  # mostly to provide access to datatables elements (display)
        )
        self.firstline = True

    def textclear(self):
        if self.parent is None:
            print("parent is None")
            raise ValueError()
        else:
            self.parent.textbox.clear()

    def textappend(self, text, color="white"):
        if self.parent is None:
            cprint(color, text)  # just go straight to the terminal
        else:
            self.parent.textbox.setTextColor(self.parent.QColor(color))
            self.parent.textbox.append(text)
            self.parent.textbox.setTextColor(self.parent.QColor("white"))

    def get_data_file(
        self, fn: Union[str, Path], changetimestamp: object, PD: dataclass
    ) -> Union[None, tuple]:
        """
        Get a data file, and also parse information from the file
        for display
        """
        fnp = Path(fn)
        fns = str(fn)
        ivdatafile = None
        if self.firstline:
            if not fnp.is_file():
            # cprint('r', f"   File: {str(fnp):s} NOT FOUND")
                self.textappend(f"   File: {str(fnp):s} NOT FOUND", color="red")
                return None
            else:
                self.textappend(f"   File: {str(fnp):s} OK")
            # print(f"   File {str(fnp):s} found.")
        mtime = fnp.stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if self.firstline:
            self.textappend(f"pgbcivr2: Checking file: {fnp.name:s} [{timestamp_str:s}]")
        # print(mtime, changetimestamp)
        if mtime > changetimestamp:
            filemode = "vcnmodel.v1"
        else:
            filemode = "vcnmodel.v0"
        if self.firstline:
            self.textappend(f"pgbcivr2: file mode: {filemode:s}")
        with (open(fnp, "rb")) as fh:
            d = pickle.load(fh)

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
            self.textappend(
                f"Unable to match soma/dendrite inflation conditions", "red"
            )
            return None

        if ivdatafile is None or not ivdatafile.is_file():
            if self.flirstline:
                self.textappend(f"no file matching conditions : {str(ivdatafile):s}")
            return None

        if self.firstline:
            self.textappend(f"\npgbcivr2: datafile to read: {str(ivdatafile):s}")
        if "time" in list(d["Results"].keys()):
            d["Results"] = self._data_flip(d["Results"])
        return par, stitle, ivdatafile, filemode, d

    def _data_flip(self, data_res):
        """
        Convert from old data file format to new
        """
        # flip order to put trials first (this was an old format)
        trials = range(len(data_res["somaVoltage"]))
        for tr in trials:
            sv = data_res["somaVoltage"][tr]
            dv = data_res["dendriteVoltage"][tr]
            st = data_res["spikeTimes"][tr]
            isp = data_res["inputSpikeTimes"][tr]
            if len(data_res["stimWaveform"].shape) > 0:
                swv = data_res["stimWaveform"][tr]
            else:
                swv = None
            stb = data_res["stimTimebase"][tr]

            ti = data_res["time"][tr]

            data_res[tr] = {
                "somaVoltage": sv,
                "dendriteVoltage": dv,
                "spikeTimes": st,
                "inputSpikeTimes": isp,
                "stimWaveform": swv,
                "stimTimebase": stb,
                "time": ti,
            }
        delete = [key for key in data_res[0].keys()]  # get from saved keys
        for key in delete:
            del data_res[key]  # but delete top level
        print("_data_flip: ", data_res.keys())
        return data_res

    @time_func
    def analyze_data(
        self, ivdatafile: Union[Path, str], filemode, protocol: str
    ) -> tuple:
        """
        Provide basic spike detection, shape analysis, and
        IV analysis if appropriate
        We use ephys.acq4read.read_pfile to read the pickled data
        file into acq4 format, then we can use the ephys
        analysis tools to analyze the data
        """
        AR.read_pfile(ivdatafile, filemode=filemode)

        # AR.traces = AR.traces*1000.
        # print(AR.traces.shape)
        # f, ax = mpl.subplots(1,1)
        # for i in range(AR.traces.shape[0]):
        #     mpl.plot(AR.traces[i])
        # mpl.show()
        bridge_offset = 0.0
        threshold = -35.0  # mV
        tgap = 0.0  # gap before fittoign taum
        RM.setup(AR, SP, bridge_offset=bridge_offset)
        SP.setup(
            clamps=AR,
            threshold=threshold * 1e-3,
            refractory=0.001,
            peakwidth=0.001,
            interpolate=True,
            verify=True,
            mode="peak",
        )

        SP.set_detector("Kalluri")  # spike detector
        SP.analyzeSpikes()
        SP.analyzeSpikeShape()
        # print('AnalyzeData: ', SP.spikes)
        RMA = None
        if protocol == "IV":
            SP.fitOne(function="fitOneOriginal")
            RM.analyze(
                rmpregion=[0.0, AR.tstart - 0.001],
                tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
                to_peak=True,
                tgap=tgap,
            )

            RMA = RM.analysis_summary
        return AR, SP, RMA

    @winprint
    def plot_traces(
        self,
        ax: object,
        fn: Union[Path, str],
        PD: dataclass,
        protocol: str,
        ymin=-80.0,
        ymax=20.0,
        iax=None,
        nax=0,
        figure=None,
    ) -> tuple:

        changetimestamp = get_changetimestamp()
        x = self.get_data_file(fn, changetimestamp, PD)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if x is None:
            print("No simulation found that matched conditions")
            print(fn)
            return
        # unpack x
        inx = str(fn).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(str(fn)[inx + 4 : inx + 7])
        if protocol in ["IV", "runIV"]:
            protocol = "IV"
        elif protocol in ["VC", "runVC"]:
            protocol = "VC"
        print("Protocol: ", protocol)
        par, stitle, ivdatafile, filemode, d = x
        si = d["Params"]
        ri = d["runInfo"]
        AR, SP, RMA = self.analyze_data(ivdatafile, filemode, protocol)
        ntr = len(AR.traces)  # number of trials
        v0 = -160.0
        deadtime = 50.0
        trstep = 25.0 / ntr
        inpstep = 2.0 / ntr
        sz = 50.0 / ntr
        noutspikes = 0
        ninspikes = 0
        for trial in range(len(AR.traces)):
            AR.traces[trial][0] = AR.traces[trial][1]
            if protocol in ["VC", "vc", "vclamp"]:
                AR.traces[trial] = AR.traces[trial].asarray() * 1e6
            ax.plot(AR.time_base, AR.traces[trial] * 1e3, "k-", linewidth=0.5)
            print('trial: ', trial)
            print(d["Results"].keys())
            if "spikeTimes" in d["Results"][trial].keys():
                spikeindex = [int(t / si.dtIC) for t in d["Results"][trial]["spikeTimes"]]
            else:
                spikeindex = SP.spikeIndices[trial]
            ax.plot(
                AR.time_base[spikeindex],
                AR.traces[trial][spikeindex] * 1e3,
                "ro",
                markersize=2.5,
            )

            sinds = np.array(SP.spikeIndices[trial]) * AR.sample_rate[trial]
            noutspikes += len(np.argwhere(sinds > deadtime))
            calx = 20.
            if protocol in ["AN", "runANSingles"]:
                if trial in list(d["Results"].keys()) and "inputSpikeTimes" in list(
                    d["Results"][trial].keys()
                ):
                    spkt = d["Results"][trial]["inputSpikeTimes"]
                elif "inputSpikeTimes" in list(d["Results"].keys()):
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
                    # ax.scatter(spkt[ian], vy, s=sz, marker="|", linewidths=0.35)
                    ninspikes += len(spkt[ian] > deadtime)

                ax.set_ylim(-80.0, 20.0)
                ax.set_xlim(50., np.max(AR.time_base))
                calx = 50.
            elif protocol in ["VC", "vc", "vclamp"]:
                pass  #
                # ax.set_ylim((-100.0, 100.0))
            else:
                ax.set_ylim(ymin, ymax)
                ax.set_xlim(0.080, np.max(AR.time_base))

        ftname = str(ivdatafile.name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = f"{ftname:s}"
        if protocol in ["IV"]:
            toptitle += (
                f"\nRin={RMA['Rin']:.1f} M$\Omega$  $\\tau_m$={RMA['taum']:.2f} ms"
            )

            if iax == 2:
                PH.calbar(
                    ax,
                    calbar=[120.0, -30.0, 10.0, 20.0],
                    scale=[1.0, 1.0],
                    axesoff=True,
                    orient="right",
                    unitNames={"x": "ms", "y": "mV"},
                    fontsize=11,
                    weight="normal",
                    color="k",
                    font="Arial",
                )
            else:
                PH.noaxes(ax)
            # insert IV curve
            secax = PLS.create_inset_axes([0.45, -0.05, 0.3, 0.3], ax, label=str(ax))
            secax.plot(
                RM.ivss_cmd_all * 1e12,
                RM.ivss_v_all * 1e3,
                "ks-",
                markersize=3,
                markerfacecolor="k",
                zorder=10,
                clip_on=False,
            )

            ltz = np.where(RM.ivss_cmd_all <= 0.0)[0]

            secax.plot(
                RM.ivpk_cmd_all[ltz] * 1e12,
                RM.ivpk_v_all[ltz] * 1e3,
                "ko-",
                markersize=3,
                markerfacecolor="w",
                zorder=10,
                clip_on=False,
            )
            secax.plot(
                RM.ivss_cmd[-1] * 1e12,
                RM.ivss_v[-1] * 1e3,
                "ro",
                markersize=3,
                markerfacecolor="r",
                zorder=100,
                clip_on=False,
            )
            PH.crossAxes(
                secax,
                xyzero=[0.0, -60.0],
                limits=[
                    np.min(RM.ivss_cmd_all) * 1e12,
                    -120,
                    np.max(RM.ivss_cmd_all) * 1e12,
                    -25.0,
                ],  #
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
            maxt = np.max(AR.time_base)
            tlen = 10.0  # ms
            PH.calbar(
                ax,
                calbar=[maxt - tlen, 2.0, tlen, 5],
                orient="right",
                unitNames={"x": "ms", "y": "nA"},
                fontsize=9,
            )
        else:
            if nax == 0:
                PH.calbar(
                    ax,
                    calbar=[calx, ymin, 10.0, 20.0],
                    unitNames={"x": "ms", "y": "mV"},
                    fontsize=9,
                )
            else:
                PH.noaxes(ax)
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
        if nax == 0:
            if figure is not None:
                figure.suptitle(toptitle, y=0.95, fontsize=9, verticalalignment="top")
            else:
                ax.set_title(toptitle, y=1.0, fontsize=8, verticalalignment="top")
        return (synno, noutspikes, ninspikes)

    def setup_VC_plots(self):
        sizer = OrderedDict(
            [
                ("A", {"pos": [0.125, 0.75, 0.6, 0.37]}),
                ("B", {"pos": [0.125, 0.75, 0.45, 0.12]}),
                ("C", {"pos": [0.125, 0.75, 0.06, 0.32]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.
        P = PH.arbitrary_grid(
            sizer, order="columnsfirst", label=True, figsize=(4.0, 6.0),
        )
        return P

    @winprint
    @time_func
    def analyzeVC(
        self, ax: object, fn: Union[Path, str], PD: dataclass, protocol: str,
    ) -> tuple:
        ax = ax[:, 0]
        changetimestamp = get_changetimestamp()
        x = self.get_data_file(fn, changetimestamp, PD)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if x is None:
            self.textappend("No simulation found that matched conditions", color="red")
            self.textappend(fn, color="red")
            return
        # unpack x
        inx = str(fn).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(fn[inx + 4 : inx + 7])
        if protocol in ["VC", "runVC"]:
            protocol = "VC"
        self.textappend(f"Protocol: {protocol:s}")
        par, stitle, ivdatafile, filemode, d = x
        AR, SP, RMA = self.analyze_data(ivdatafile, filemode, protocol)
        tss = [0, 0]
        sr = AR.sample_rate[0]
        tss[0] = int(
            AR.tstart / sr + (AR.tdur / 2.0) / sr
        )  # wrong duration stored in traces - need to fix.
        tss[1] = int(AR.tstart / sr + (AR.tdur / 1.0) / sr)
        I = AR.traces.asarray()
        V = AR.cmd_wave
        ntraces = np.shape(V)[0]

        # initialize all result arrays, lists and dicts
        vss = np.empty(ntraces)
        vmin = np.zeros(ntraces)
        vrmss = np.zeros(ntraces)
        vm = np.zeros(ntraces)
        i0 = np.zeros(ntraces)
        ic = np.zeros(ntraces)
        for j in range(0, ntraces):
            vss[j] = np.mean(V[j, tss[0] : tss[1]])  # steady-state voltage
            ic[j] = np.mean(I[j, tss[0] : tss[1]])  # corresponding currents
            vm[j] = np.mean(
                V[j, 0 : int((AR.tstart - 1.0) / AR.sample_rate[0])]
            )  # resting potential - for 1 msec prior to step
            i0[j] = np.mean(
                I[j, 0 : int((AR.tstart - 1.0) / AR.sample_rate[0])]
            )  # resting/holding current - for 1 msec prior to step
            ax[1].plot(AR.time, V[j, :], "k-")
            # ax[0].plot(AR.time[tss[0]:tss[1]], I[j,tss[0]:tss[1]], 'r--')

        # now fit traces to variation of g = gmax * I/(V-Vr)
        gmodel = Model(boltzG)
        Ek = -84.0 * 1e-3  # mV to V
        # transform to g to get conductance values and fit to Boltzmann
        gss = ic / (vss - Ek)
        gmodel.set_param_hint("gmax", value=20.0e-9, min=0.0, vary=True)
        gmodel.set_param_hint("vhalf", value=-38e-3, min=-90e-3, max=90e-3)
        gmodel.set_param_hint(
            "k", value=7, min=0.05 * 1e-3, max=200.0 * 1e-3, vary=True
        )
        gmodel.set_param_hint("E", value=Ek, vary=False)
        gparams = gmodel.make_params()
        weights = np.ones_like(vss)
        d_index = np.argwhere((vss >= -0.120) & (np.fabs(vss - Ek) > 0.008))
        weights_sel = weights[d_index]
        vss_sel = vss[d_index]
        gss_sel = gss[d_index]
        result = gmodel.fit(
            gss_sel, method="nedler", params=gparams, x=vss_sel, weights=weights_sel
        )

        # capacitance transient. Use the single trace nearest to -70 mV for the fit.
        mHypStep = np.argmin(np.fabs(vss - (-0.090)))
        sr = AR.sample_rate[mHypStep]
        t0 = int(AR.tstart / sr)

        pts = int(5.0 / sr)  # fit over first 5 msec
        tfit = AR.time_base[t0 + 1 : t0 + pts] - AR.time_base[t0]
        ifit = I[mHypStep, t0 + 1 : t0 + pts]
        expmodel = Model(exp2decay)
        # pars = expmodel.guess(ifit, x=tfit)
        expmodel.set_param_hint("tau0", value=0.1, min=0.005, max=10.0)
        expmodel.set_param_hint("tau1", value=0.3, min=0.005, max=10.0)
        expmodel.set_param_hint("a0", value=-1e-9, min=-50e-9, max=1e-12)
        expmodel.set_param_hint("a1", value=-1e-9, min=-50e-9, max=1e-12)
        expmodel.set_param_hint("offset", value=0.0, min=-1e-8, max=1e-8, vary=True)
        expparams = expmodel.make_params()
        exp_result = expmodel.fit(ifit, expparams, x=tfit)
        # print(exp_result.params)
        deltaV = vss[mHypStep] - vm[mHypStep]
        deltaI = I[mHypStep][t0] - i0[mHypStep]
        deltaI2 = exp_result.params["a0"].value + exp_result.params["a1"].value
        Q = np.trapz(I[mHypStep, t0 : t0 + 2 * pts], dx=sr * 1e-3)
        Cm = Q / deltaV
        self.textappend(
            f"Q: {Q*1e12:.3f} pC  deltaV: {deltaV*1e3:.1f} mV, Cm: {Cm*1e12:.1f} pF"
        )
        Rs_est = deltaV / deltaI2  # Estimate of Rs from peak current
        self.textappend(f"Estimated Rs: {Rs_est*1e-6:.1f} MOhm")
        # print(deltaV, deltaI, deltaI2)
        # other approach: use fastest tau in voltage clamp
        #
        tau0 = 1e-3 * exp_result.params["tau0"].value  # convert to seconds
        tau1 = 1e-3 * exp_result.params["tau1"].value
        a0 = exp_result.params["a0"].value
        a1 = exp_result.params["a1"].value
        #
        # Here we use the fastest time constant and the
        # associated current from the fit to estimate cm (perisomatic)
        # Note: do not use total input resistance!
        # See Golowasch et al., J. Neurophysiol. 2009 and references therein

        if tau0 < tau1:
            tau = tau0
            R0 = deltaV / a0
            cm = tau / R0
        else:
            tau = tau1
            R0 = deltaV / a1
            cm = tau / R0
        # for curiosity, weighted tau (as if compensation was done as
        # much as possible for the whole transient)
        #
        tauw = (a0 * tau0 + a1 * tau1) / (a0 + a1)
        R0w = deltaV / a0
        R1w = deltaV / a1
        cm1 = tau1 / R1w
        cmw = ((a0 * tau0 / R0w) + (a1 * tau1 / R1w)) / (a0 + a1)  # tauw/(R0w + R1w)
        self.textappend("By Coeffs: ")
        self.textappend(
            f"RCoeff0 = {R0w*1e-6:.2f} MOhm, tau0: {tau0*1e3:.3f} ms,  cm0: {cm*1e12:.1f} pF"
        )
        self.textappend(
            f"RCoeff1 = {R1w*1e-6:.2f} MOhm, tau1: {tau1*1e3:.3f} ms,  cm1: {cm1*1e12:.1f} pF"
        )
        self.textappend(
            f"Weighted: Rw={(R0w+R1w)*1e-6:.2f} tauw: {tauw*1e3:.3f} ms, Weighted cm: {cmw*1e12:.1f} pF"
        )
        tfit2 = AR.time_base[t0 : t0 + pts] - AR.time_base[t0]
        bfit = expmodel.eval(params=exp_result.params, x=tfit2)
        ax[0].plot(
            tfit2 + AR.time_base[t0], bfit * 1e9, "r-", dashes=[6, 2], linewidth=1
        )
        # ax[0].plot(tfit+AR.time_base[t0], exp_result.best_fit*1e9, 'r-')
        # ax[0].plot(tfit+AR.time_base[t0], ifit*1e9, 'g--')
        ax[2].plot(vss_sel, gss_sel * 1e9, "ko", markersize=2)
        ax[2].plot(vss_sel, result.best_fit * 1e9, "r-")
        # print('vhalf: ', result.params['vhalf'])
        textstr = r"g$_{max}$  = "
        textstr += f"{result.params['gmax'].value*1e9:.1f} nS\n"
        textstr += r"$V_{0.5}$  = "
        textstr += f"{result.params['vhalf'].value*1e3:.1f} mV\n"
        textstr += f"k  = {1e3*result.params['k'].value:.1f}\n"
        textstr += f"Cm = {cm*1e12:.1f} pF\n"
        textstr += r"$tau_{0}$ = "
        textstr += f"{tau*1e3:.3f} ms"
        props = dict(boxstyle="square", facecolor="None", alpha=0.5)
        ax[2].text(
            0.05,
            0.95,
            textstr,
            transform=ax[2].transAxes,
            fontsize=8,
            verticalalignment="top",
            bbox=props,
        )
        mpl.show()

    def trace_viewer(
        self,
        filename: Union[Path, str],
        PD: Union[object, None] = None,
        runProtocol: Union[str, None] = None,
    ):
        movie = False
        changetimestamp = get_changetimestamp()
        x = self.get_data_file(filename, changetimestamp, PD)
        mtime = Path(filename).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if x is None:
            print("No simulation found that matched conditions")
            print(filename)
            return
        # unpack x
        inx = str(filename).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(str(filename)[inx + 4 : inx + 7])
        if runProtocol in ["IV", "runIV"]:
            runProtocol = "IV"
        elif runProtocol in ["VC", "runVC"]:
            runProtocol = "VC"
        print("Protocol: ", runProtocol)
        par, stitle, ivdatafile, filemode, d = x
        AR, SP, RMA = self.analyze_data(ivdatafile, filemode, runProtocol)

        self.ntr = len(AR.traces)  # number of trials
        self.pwin = None
        self.pfull = None
        self.lines = [None]*self.parent.n_trace_sel
        self.punchtape = None
        self.pt_text = None
        self.inputtimes = None
        self.first = True
        self.parent.trace_selector.setValue(0)
        if self.allspikes is not None:
            self.nspikes = len(self.allspikes)
            for i, trn in enumerate(range(0, self.parent.n_trace_sel)):
                self.plot_spikes_withinputs(
                    int(i), n=trn, color=pg.intColor(i, hues=20), first=self.first
                )

            self.parent.trace_selector_plot.setXRange(0, self.nspikes)
            self.parent.trace_selector.setBounds((0, self.nspikes))
            self.parent.frameTicks.setXVals(
                range(0, self.nspikes, self.parent.n_trace_sel)
            )
            self.first = False
            QtGui.QApplication.processEvents()
        else:
            return

        self.parent.Dock_Traces.raiseDock()
        self.parent.trace_selector.sigPositionChanged.connect(self.timeLineChanged)

        nrep = 1
        self.parent.trace_plots.setXRange(-5.0, 2.5, padding=0.1)
        self.check_yscaling()
        self.parent.trace_plots.clear()
        print("movie state: ", self.parent.movie_state)

        if self.parent.movie_state:  # start with the movie button

            while self.parent.movie_state:
                for i in range(0, nrep * self.nspikes, self.parent.n_trace_sel):
                    # self.parent.trace_selector.setValue(int(np.mod(i, self.nspikes)))
                    self.parent.trace_plots.clear()

                    # QtGui.QApplication.processEvents()
                    for ix, trn in enumerate(range(i, i + self.parent.n_trace_sel)):
                        self.plot_spikes_withinputs(
                            ix,
                            trn,
                            color=pg.intColor(
                                np.mod(trn, self.parent.n_trace_sel),
                                hues=self.parent.n_trace_sel,
                            ),
                            first=self.first,
                        )

                    QtGui.QApplication.processEvents()
                    if not self.parent.movie_state:
                        break
                    time.sleep(self.parent.frame_interval)
                    QtGui.QApplication.processEvents()
                # print(self.parent.movie_state)
                QtGui.QApplication.processEvents()
            return
        else:
            i = 0
            # self.parent.trace_plots.clear()

            # QtGui.QApplication.processEvents()
            for ix, trn in enumerate(range(i, i + self.parent.n_trace_sel)):
                self.plot_spikes_withinputs(
                    ix,
                    trn,
                    color=pg.intColor(
                        np.mod(trn, self.parent.n_trace_sel),
                        hues=self.parent.n_trace_sel,
                    ),
                )
            QtGui.QApplication.processEvents()

    def check_yscaling(self):
        if self.parent.V_disp_sel == "dV/dt":
            self.parent.trace_plots.setYRange(-100.0, 200.0, padding=0.1)
        else:
            self.parent.trace_plots.setYRange(-75.0, 10.0, padding=0.1)
            pass

    def timeLineChanged(self):
        # cprint('r', 'Time line changed')
        if self.first:
            return
        t = int(self.parent.trace_selector.value())
        if t < 0:
            self.parent.trace_selector.setValue(0)
            t = 0
        if t >= self.nspikes - self.parent.n_trace_sel:
            self.parent.trace_selector.setValue(self.nspikes - self.parent.n_trace_sel)
            t = self.nspikes - self.parent.n_trace_sel
        if np.mod(t, self.parent.n_trace_sel) == 0:
            self.parent.trace_plots.clear()
            for ix, trn in enumerate(range(t, t + self.parent.n_trace_sel)):
                self.plot_spikes_withinputs(
                    ix,
                    trn,
                    color=pg.intColor(
                        np.mod(trn, self.parent.n_trace_sel),
                        hues=self.parent.n_trace_sel,
                    ),
                )

    def plot_spikes_withinputs(
        self, ix: int = 0, n: int = 0, color: object = None, first=False
    ):
        """
            plot a spike and indicate its inputs.
        ix : index into this run (counter): for plotting a block of spikes
        n : spike to plot within the block
        color: indexed color.
        """
        print(self.nspikes, color, first, ix, n)
        # self.parent.trace_plots.plot(np.arange(10), np.ones(10))
        if self.nspikes > n and n >= 0:
            self.check_yscaling()
            # print('n: ', n)
            spk = self.allspikes[n]
            if self.parent.V_disp_sel == "dV/dt":
                if first:
                    self.lines[ix] = self.parent.trace_plots.plot(
                        spk.dt * np.arange(0, len(spk.waveform))[:-1] + spk.start_win,
                        np.diff(spk.waveform / spk.dt),
                        pen=color,
                    )
                else:
                    # self.parent.trace_plots.plot(
                    #     spk.dt * np.arange(0, len(spk.waveform))[:-1] + spk.start_win,
                    #     np.diff(spk.waveform / spk.dt),
                    #     pen=color,
                    # )
                    self.lines[ix].setData(
                        spk.dt*np.arange(0, len(spk.waveform))[:-1]+spk.start_win, 
                        np.diff(spk.waveform/spk.dt), 
                        pen=color)
                papertape_yoff = 120.0
                spkin_off = -50.0
                dy = 5.0
            else:
                if first:
                    print('plotting first')
                    tx = spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win
                    # print(np.min(tx), np.max(tx), tx.shape)
 #                    print(np.min(spk.waveform), np.max(spk.waveform), spk.waveform.shape)
                    # print(dir(self.parent.trace_plots))
                    self.parent.trace_plots.setEnabled(True)
                    self.parent.trace_plots.plotItem.plot(np.arange(10), np.ones(10)*ix, pen="r")
                    self.lines[ix] = self.parent.trace_plots.plotItem.plot(
                        x=tx,
                        y=spk.waveform,
                        # pen=color,
                    )
                    self.lines[ix].curve.show()
                    #print(dir(self.lines[ix]))
                 
                else:
                    self.check_yscaling()
                    # self.parent.trace_plots.plot(
                    #     spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win,
                    #     spk.waveform,
                    #     pen=color,
                    # )
                    self.parent.trace_plots.setEnabled(True)
                    u = self.parent.trace_plots.plotItem.plot(np.arange(10), np.ones(10)*ix, pen='b')
                    tx = spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win
                    print(np.min(tx), np.max(tx), tx.shape)
                    print(np.min(spk.waveform), np.max(spk.waveform), spk.waveform.shape)
                    print(dir(self.lines[ix]))
                    self.lines[ix].curve.setData(
                        x=tx,
                        y=spk.waveform,
                        pen='w',# pen=color,
                        )
                    self.lines[ix].curve.show()
                papertape_yoff = 0.0
                spkin_off = -65.0
                dy = 2.0
            # QtGui.QApplication.processEvents()
            # print('first: ', first, 'lines: ', self.lines)
            # prespt = np.zeros(len(spk.prespikes))
            # prespv = np.nan * np.zeros(len(spk.prespikes))
            # prespk_time = np.zeros(len(spk.prespikes))
            # prespv2 = np.nan * np.zeros(len(spk.prespikes))
            # # print(spk.prespikes)
            # for i, prespk in enumerate(spk.prespikes):
            #     if prespk is not np.nan:
            #         prespt[i] = (
            #             0.0 - 2.5 + 0.07 * np.mod(n, self.parent.n_trace_sel)
            #         )  # prespk
            #         # prespk_time[i] = prespk
            #         prespv[i] = -i * dy + papertape_yoff
            #         prespv2[i] = spkin_off - i * dy
            # if ix == 0:
            #     for j in range(self.ninputs):
            #         text = pg.TextItem(f"{j:>2d}", anchor=(1.0, 0.5))
            #         if first:
            #             self.pttext = self.parent.trace_plots.addItem(text)
            #         text.setPos(0 - 2.7, -j * dy + papertape_yoff)
            # indices = np.logical_not(np.logical_or(np.isnan(prespt), np.isnan(prespv)))
            # indices = np.array(indices)
            # sysizes = np.zeros((len(prespt), len(prespv)))
            #
            # # prespt = prespt[indices]
            # # prespv = prespv[indices]
            # # print(prespt, prespv)
            # # print(prespt, prespv)
            # # if len(prespt) > 0:
            # if first:
            #     self.punchtape = self.parent.trace_plots.plot(
            #         prespt,
            #         prespv,
            #         pen=None,
            #         symbol="o",
            #         symbolBrush=color,
            #         symbolSize=8,
            #     )
            #     self.inputtimes = self.parent.trace_plots.plot(
            #         spk.prespikes,
            #         prespv2,
            #         pen=None,
            #         symbol="o",
            #         symbolBrush=color,
            #         symbolSize=8,
            #     )
            # else:
            #     # self.punchtape.setData(prespt, prespv, pen=None, symbol='o', symbolBrush=color, symbolSize=8)
            #     # self.inputtimes.setData(spk.prespikes, prespv2, pen=None, symbol='o', symbolBrush=color, symbolSize=8)
            #     self.parent.trace_plots.plot(
            #         prespt,
            #         prespv,
            #         pen=None,
            #         symbol="o",
            #         symbolBrush=color,
            #         symbolSize=8,
            #     )
            #     self.parent.trace_plots.plot(
            #         spk.prespikes,
            #         prespv2,
            #         pen=None,
            #         symbol="o",
            #         symbolBrush=color,
            #         symbolSize=8,
            #     )

    def plot_revcorr_map(
        self,
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

    def make_patch_spines_invisible(self, axn):
        axn.set_frame_on(True)
        axn.patch.set_visible(False)
        for sp in axn.spines.values():
            sp.set_visible(False)

    def plot_revcorr2(self, P: object, PD: dataclass, RCP: dataclass, RCD: dataclass):
        seaborn.set_style("ticks")
        # secax = twinax(P.figure_handle, ax, pos=maxwin)
        ax = P.axdict["B"]
        secax = P.axdict["A"]
        PH.noaxes(secax, "xy")

        secax.set_facecolor((1, 1, 1, 0))
        ax.spines["top"].set_visible(False)
        summarySiteTC = {}
        RCD.max_coin_rate = 0.0
        maxrevcorr = 0
        for isite in range(
            RCP.ninputs
        ):  # range(ninputs):  # for each ANF input (get from those on first trial)
            if RCD.sv_all.shape == ():
                continue
            stepsize = int(RCD.sv_all.shape[0] / 20)
            if stepsize > 0:
                sel = list(range(0, RCD.sv_all.shape[0], stepsize))
            else:
                sel = list(range(0, RCD.sv_all.shape[0], 1))
            sel = list(range(0, RCD.sv_all.shape[0], 1))
            refzero = int(RCP.minwin / RCP.binw)
            if RCD.C[isite] is not None:
                nc = int(len(RCD.C[isite]) / 2)
                RCD.TC = RCD.TC / len(RCD.st)
                summarySiteTC[isite] = RCD.TC
                color = mpl.cm.viridis(norm(RCD.sites, isite))
                maxrevcorr = np.max((maxrevcorr, np.max(RCD.CB[isite])))
                totalrevcorr = np.sum(RCD.CB[isite])

                if isite in range(RCP.ninputs):
                    if RCP.algorithm == "RevcorrSPKS":
                        ax.plot(
                            RCD.tx,
                            RCD.C[isite][:nc],
                            color=color,
                            label=(
                                "Input {0:2d} N={1:3d}".format(
                                    isite, int(RCD.sites[isite])
                                )
                            ),
                            linewidth=1.5,
                            zorder=5,
                        )
                        RCD.max_coin_rate = np.max(
                            (RCD.max_coin_rate, np.max(RCD.C[:nc]))
                        )

                    elif RCP.algorithm == "RevcorrSimple":
                        ax.plot(
                            RCD.CBT - RCP.maxwin,
                            RCD.CB[isite] / float(RCD.npost_spikes),
                            color=color,
                            label=(
                                "Input {0:2d} N={1:3d}".format(
                                    isite, int(RCD.sites[isite])
                                )
                            ),
                            linewidth=1,
                            zorder=5,
                            alpha=1.0,
                        )
                        RCD.max_coin_rate = np.max(
                            (
                                RCD.max_coin_rate,
                                np.max(RCD.CB[isite] / float(RCD.npost_spikes)),
                            )
                        )

                    elif RCP.algorithm == "RevcorrSTTC":  # use spike time tiling method
                        ax.plot(
                            RCD.tx,
                            RCD.STTC[isite][:nc],
                            color=color,
                            label=(
                                "Input {0:2d} N={1:3d}".format(
                                    isite, int(RCD.sites[isite])
                                )
                            ),
                            linewidth=1.5,
                            zorder=5,
                        )

            if (
                isite == 0
            ):  # only plot the first time through - the APs are the same no matter the input
                for t in sel:
                    secax.plot(
                        RCD.ti_avg,
                        RCD.sv_all[t],
                        color="#666666",
                        linewidth=0.2,
                        zorder=1,
                    )

                secax.plot(RCD.ti_avg, RCD.sv_avg, color="k", linewidth=0.75, zorder=2)
                secax.plot([0.0, 0.0], [-120.0, 10.0], "r", linewidth=0.5)
                PH.noaxes(secax)
                PH.calbar(
                    secax,
                    calbar=[1.35, -30, 1.0, 20.0],
                    axesoff=True,
                    orient="right",
                    unitNames={"x": "ms", "y": "mV"},
                )
                PH.referenceline(secax, -60.0)
                PH.noaxes(secax)
        print(f"Total spikes plotted: {RCD.nsp_avg:d}")

        seaborn.despine(ax=ax)
        ax.set_ylabel("Rate of coincidences/bin (Hz)", fontsize=10)
        ax.set_xlabel("T (ms)", fontsize=10)
        ax.set_xlim((RCP.minwin, RCP.maxwin))

        # print(RCD.max_coin_rate)
        if RCD.max_coin_rate > 0.0:
            ns = PH.NiceScale(0.0, RCD.max_coin_rate)
            ax.set_ylim(0, ns.niceMax)
        else:
            ax.set_ylim(0, 0.25)
        yls = ax.get_ylim()
        secax.set_ylim([-70.0, 10.0])
        secax.set_xlim((RCP.minwin, RCP.maxwin))
        secax.tick_params(direction="in", length=5.0, width=1.0, labelsize=9)
        ax.tick_params(direction="in", length=5.0, width=1.0, labelsize=9)
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            pointSize=10,
            tickPlacesAdd={"x": 1, "y": 1},
            floatAdd={"x": 1, "y": 1},
        )

        return summarySiteTC

    def get_synaptic_info(self, gbc: str) -> tuple:
        SC = cell_config.CellConfig()
        syninfo = SC.VCN_Inputs[gbc]
        return (SC, syninfo)

    def get_data(
        self, fn: Union[Path, str], PD: dataclass, changetimestamp, protocol
    ) -> Union[None, tuple]:

        X = self.get_data_file(fn, changetimestamp, PD)
        if X is None:
            print("No simulation found that matched conditions")
            print("Looking for file: ", fn)
            return None
        # unpack x
        par, stitle, ivdatafile, filemode, d = X
        if "time" in list(d["Results"].keys()):
            d["Results"] = self._data_flip(d["Results"])

        # 2. find spikes
        AR, SP, RMA = self.analyze_data(ivdatafile, filemode, protocol)
        # set up analysis parameters and result storage
        RCP = RevCorrPars()
        RCD = RevCorrData()

        RCD.npost = 0  # number of postsynaptic spikes
        RCD.npre = 0  # number of presynaptic spikes

        trials = range(len(d["Results"]))
        RCP.ntrials = len(trials)
        print("res keys: ", d["Results"].keys())
        for i, tr in enumerate(trials):
            trd = d["Results"][tr]  # trial data
            ti = trd["time"]
            for n in range(len(trd["inputSpikeTimes"])):  # for each sgc
                RCD.npre += len(trd["inputSpikeTimes"][n])
            RCD.npost += len(trd["spikeTimes"])
            RCD.st = SP.spikeIndices
            RCD.npost += len(RCD.st[tr])

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

    @time_func
    def compare_revcorrs(self):
        plabels = [f"VCN_c{int(self.parent.cellID):02d}"]
        pgbc = plabels[0]
        revcorrtype = "RevcorrSimple"
        PSum = PH.regular_grid(
            rows=1,
            cols=1,
            order="rowsfirst",
            figsize=(6, 6),
            showgrid=False,
            verticalspacing=0.1,
            horizontalspacing=0.1,
            margins={
                "bottommargin": 0.15,
                "leftmargin": 0.15,
                "rightmargin": 0.15,
                "topmargin": 0.15,
            },
            label=["A"],
            labelposition=(-0.05, 1.05),
        )

        for iax, index_row in enumerate(self.parent.selected_index_rows):
            selected = self.parent.table_manager.get_table_data(
                index_row
            )  # table_data[index_row]
            if selected is None:
                return
            sfi = Path(selected.simulation_path, selected.files[0])
            res = self.compute_revcorr(
                None, pgbc, sfi, PData(), selected.runProtocol, revcorrtype
            )
            # unpack
            ninputs, ynspike, sites, participation, nspikes = res
            PSum.axarr[0, 0].plot(
                np.arange(ninputs) + 1,
                ynspike,
                label=f"{selected.dendriteMode:s} {selected.synapseExperiment}",
            )
            PSum.axarr[0, 0].set_ylim(0, 1.0)

        PSum.axarr[0, 0].legend(fontsize=7)
        PSum.axarr[0, 0].set_xlabel("# of inputs active prior to spike", fontsize=10)
        PSum.axarr[0, 0].set_ylabel("Cumulative Fraction of bushy spikes", fontsize=10)

        PSum.figure_handle.show()

    def plot_revcorr_figure(self, selected, revcorrtype):
        PD = PData()

        plabels = [f"VCN_c{int(self.parent.cellID):02d}"]
        pgbc = plabels[0]

        sizer = {
            "A": {
                "pos": [0.05, 0.40, 0.52, 0.40],
                "labelpos": (0.02, 1.00),
                "noaxes": True,
            },
            "B": {"pos": [0.05, 0.40, 0.08, 0.35], "labelpos": (0.02, 1.00)},
            "C": {"pos": [0.52, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
            "D": {"pos": [0.52, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
            "E": {"pos": [0.78, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
            "F": {"pos": [0.78, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
        }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            figsize=(9, 5),
            label=True,
            # verticalspacing=0.12,
            # horizontalspacing=0.12,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.1,
                "rightmargin": 0.1,
                "topmargin": 0.1,
            },
            fontsize={"tick": 7, "label": 9, "panel": 12},
            fontweight={"tick": "normal", "label": "normal", "panel": "bold"},
        )

        dPD = PData()
        sfi = Path(selected.simulation_path, selected.files[0])
        res = self.compute_revcorr(P, pgbc, sfi, PD, selected.runProtocol, revcorrtype)

        P.figure_handle.show()

    # @time_func
    def revcorr(
        self,
        st1: Union[np.ndarray, List] = None,
        st2: Union[np.ndarray, List] = None,
        binwidth: float = 0.1,
        datawindow: Union[List, Tuple] = [
            0.0,
            100.0,
        ],  # time window for selecteing spikes
        corrwindow: Union[List, Tuple] = [
            -5.0,
            1.0,
        ],  # time window to examine correlation relative to st1
    ) -> np.ndarray:
        if st1 is None or st2 is None:
            raise ValueError(
                "coincident_spikes_correlation: reference and comparator must be defined"
            )
        refa = st1  # [a for a in st1 if (datawindow[0] <= a <= datawindow[1])]
        refb = st2  # [b for b in st2 if (datawindow[0] <= b <= datawindow[1])]
        # xds  = [None]*len(refa)
        xds = np.zeros(int((corrwindow[1] - corrwindow[0]) / binwidth))
        for i, sp in enumerate(refa):
            diff = refb - sp
            v = diff[np.where((corrwindow[0] <= diff) & (diff <= corrwindow[1]))]
            # print('v: ', v)
            if len(v) > 0:
                indxs = [int(vx / binwidth) for vx in v]
                xds[indxs] = xds[indxs] + 1

        return xds, len(st1)  # return the n postsynaptic spikes

    def _count_spikes_in_window(self, d, trial, site, s, RCP, pre_w):
        an_i = d["Results"][trial]["inputSpikeTimes"][
            site
        ]  # input spike times for one input
        an_i = an_i[
            (an_i > RCP.min_time) & (an_i < RCP.max_time)
        ]  # restrict to those only within the response window
        an_i = an_i - s  # get relative latency from spike to it's inputs
        # print('ani: ', ani)
        pre_indx = np.asarray((an_i >= pre_w[0]) & (an_i <= pre_w[1])).nonzero()[0]
        pre_times = [an_i[k] for k in pre_indx]
        if len(pre_times) > 0:
            # print("pre times, prewindow: ", pre_times, pre_w)

            npre_i = len(pre_times)
        else:
            npre_i = 0
            pre_times = []
        return npre_i, pre_times

    @winprint
    # @time_func
    def compute_revcorr(
        self,
        P: object,
        gbc: str,
        fn: Union[str, Path],
        PD: object,
        protocol: str,
        revcorrtype: str = "RevcorrSPKS",
        thr: float = -20.0,
        width: float = 4.0,
    ) -> Union[None, tuple]:

        changetimestamp = get_changetimestamp()
        self.allspikes = None

        #
        # 1. Gather data
        #
        print("Getting data")
        SC, syninfo = self.get_synaptic_info(gbc)

        res = self.get_data(fn, PD, changetimestamp, protocol)
        if res is None:
            return None
        (d, AR, SP, RMA, RCP, RCD) = res
        si = d["Params"]
        ri = d["runInfo"]
        RCP.algorithm = revcorrtype
        if RCP.algorithm == "RevcorrSTTC":
            sttccorr = STTC.STTC()

        print("Preparing for computation")
        RCP.ninputs = len(syninfo[1])
        self.ninputs = RCP.ninputs  # save for trace viewer.
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

        if isinstance(si, dict):
            min_time = (si["pip_start"] + 0.025) * 1000.0  # push onset out of the way
            max_time = (si["pip_start"] + si["pip_duration"]) * 1000.0
            F0 = si["F0"]
            dB = (si["dB"],)
        else:
            if isinstance(ri.pip_start, list):
                start = ri.pip_start[0]
            else:
                start = ri.pip_start
            min_time = (start + 0.025) * 1000.0
            max_time = (start + ri.pip_duration) * 1000.0
        expt = ri.Spirou
        print("expt: ", expt)

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
        nbins = int((max_time - min_time) / RCP.binw)

        RCD.CBT = np.arange(RCP.minwin, RCP.maxwin, RCP.binw)  # refcorr time base
        for isite in range(RCP.ninputs):
            RCD.CB[isite] = np.zeros_like(RCD.CBT)
            ncpts = len(np.arange(RCP.minwin, -RCP.minwin, RCP.binw))
            RCD.C[isite] = np.zeros(ncpts)
            RCD.STTC[isite] = np.zeros(ncpts)

        #
        # sum across trials, and sort by inputs
        #
        # print("starting loop")

        start_time = datetime.datetime.now()
        for trial in range(RCP.ntrials):  # sum across trials
            spikeindex = [int(t / si.dtIC) for t in d["Results"][trial]["spikeTimes"]]
            stx = AR.time_base[spikeindex]
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

                if RCD.nsp_avg == 0:  # init arrays
                    RCD.sv_avg = d["Results"][trial]["somaVoltage"][areltime]
                    RCD.nsp_avg = 1
                    RCD.ti_avg = RCD.ti[0 : len(areltime)] + RCP.minwin
                    RCD.sv_all = np.zeros(
                        (RCP.ntrials, RCD.ti_avg.shape[0])
                    )  # initialize the array
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
                RCD.sv_all[trial] = d["Results"][trial]["somaVoltage"][areltime]

                nspk_plot += RCD.nsp_avg
                RCD.sv_sites.append(RCD.sv_all)

            # Now get reverse  correlation for each input
            for isite in range(RCP.ninputs):  # for each ANF input
                anx = d["Results"][trial]["inputSpikeTimes"][
                    isite
                ]  # get input AN spikes and trim list to window
                anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
                if len(anx) == 0:
                    continue
                RCD.npre_spikes += len(anx)  # count up pre spikes.
                if revcorrtype == "RevcorrSPKS":
                    RCD.C[isite] += SPKS.correlogram(
                        stx, anx, width=-RCP.minwin, binwidth=RCP.binw, T=None
                    )

                elif revcorrtype == "RevcorrSimple":
                    refcorr, npost = self.revcorr(
                        stx,
                        anx,
                        binwidth=RCP.binw,
                        datawindow=[RCP.min_time, RCP.max_time],
                        corrwindow=[RCP.minwin, RCP.maxwin],
                    )
                    RCD.CB[isite] = RCD.CB[isite] + refcorr

                elif revcorrtype == "RevcorrSTTC":
                    sttccorr.set_spikes(RCP.binw, stx, anx, RCP.binw)
                    refcorr = sttccorr.calc_ccf_sttc(
                        corrwindow=[RCP.minwin, RCP.maxwin], binwidth=RCP.binw
                    )
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
        print("Time for calculation: ", elapsed_time)
        pre_w = [-2.7, -0.5]

        ################
        # now some pariwise and participation stats on input events prior to a spike
        ################
        spikedata = SpikeData()  # storage in a dataclass
        self.allspikes = []
        pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
        participation = np.zeros(RCP.ninputs)
        pre_spike_counts = np.zeros(
            RCP.ninputs + 1
        )  # there could be 0, or up to RCP.ninputs pre spikes
        nperspike = []
        nspikes = 0
        sellist = [True] * RCP.ninputs
        if ri.Spirou == "largestonly":
            for i in range(1, len(sellist)):
                sellist[i] = False
        elif ri.Spirou == "twolargest":
            for i in range(2, len(sellist)):
                sellist[i] = False
        elif ri.Spirou == "removelargest":
            sellist[0] = False

        for trial in range(RCP.ntrials):  # accumulate across all trials
            spikeindex = [int(t / si.dtIC) for t in d["Results"][trial]["spikeTimes"]]
            spks = AR.time_base[spikeindex]  # get postsynaptic spikes for the trial
            for n, s in enumerate(spks):  # for each postsynaptic spike
                if (
                    s < RCP.min_time or s > RCP.max_time
                ):  # restrict post spikes to those only in a response window
                    continue
                reltime = np.around(RCD.ti, 5) - np.around(s, 5)
                areltime = np.argwhere(
                    (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
                ).squeeze()
                spikedata = SpikeData()  # storage in a dataclass
                spikedata.trial = trial
                spikedata.waveform = d["Results"][trial]["somaVoltage"][areltime]

                spikedata.time_index = n
                spikedata.prespikes = [np.nan] * RCP.ninputs
                nspikes += 1  # number of post spikes evaluated
                npre_spikes = 0  # number of pre spikes associated with this post spike
                for isite in range(RCP.ninputs):  # examine each input
                    if not sellist[isite]:
                        continue
                    npre_i, pre_times = self._count_spikes_in_window(
                        d, trial, isite, s, RCP, pre_w
                    )
                    if npre_i > 0:
                        spikedata.prespikes[isite] = pre_times[
                            0
                        ]  # save all event times even if more thn=an one
                        # print('tiime of pre spikes: ', pre_index*AR.sample_rate[0]+pre_w[0])
                        participation[
                            isite
                        ] += 1  # any spikes in the window = participation (but only count as 1)
                        npre_spikes += 1
                        # print(' spk: ', s, 'isite: ', isite)
                        for jsite in range(
                            isite + 1, RCP.ninputs
                        ):  # now do for joint combinations with other remaining inputs
                            npre_j, pre_times = self._count_spikes_in_window(
                                d, trial, jsite, s, RCP, pre_w
                            )
                            # print(npre_j)
                            if npre_j > 0:  # accumulate if coincident for this pair
                                pairwise[isite, jsite] += 1
                pre_spike_counts[
                    npre_spikes
                ] += 1  # increment the number of times there were npre_spikes input to this post spike
                self.allspikes.append(spikedata)

        print(pairwise)
        # print('pre spike_count associated with a post spike: ', pre_spike_counts)
        # plot the position of the prespikes for every trial as determined by the
        # second trial loop above.
        # f, ax = mpl.subplots(1,1)
        # for i in range(len(self.allspikes)):
        #     y = i*np.ones(len(self.allspikes[i].prespikes))
        #     print(list(self.allspikes[i].prespikes))
        #     print(y)
        #     ax.plot(self.allspikes[i].prespikes, y)
        # mpl.show()
        # print(self.allspikes)
        npartipating = np.sum(participation)
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
        # nperspike = [n for n in nperspike if n != 0]
        # nperspike = scipy.stats.itemfreq(nperspike).T
        # print('nperspike counts: ', nperspike)
        # nperspike = np.array(np.unique(nperspike, return_counts=True))/nspikes
        # properly fill out output
        # xnspike = np.arange(RCP.ninputs)
        # ynspike = np.zeros(RCP.ninputs)
        # for j, i in enumerate(nperspike[0]):
        #     # print(i, j, nperspike[1,j])
        #     ynspike[i - 1] = nperspike[1, j]

        # ynspike = np.cumsum(ynspike / nspikes)
        ynspike = np.cumsum(pre_spike_counts[1:]) / np.sum(pre_spike_counts[1:])
        # print(RCD.sites)
        # print(pos)
        maxp = np.max(pairwise)

        if P is None:
            return (RCP.ninputs, ynspike, RCD.sites, participation, nspikes)

        ax = P.axdict["B"]
        summarySiteTC = self.plot_revcorr2(P, PD, RCP, RCD)

        # ax.set_title(
        #     f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}" +
        #     f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s}",
        #     fontsize=11,
        # )
        # return summarySiteTC, RCD.sites

        sax = P.axdict
        # f, sax = mpl.subplots(3,1)
        # f.set_size_inches( w=3.5, h=9)
        sax["C"].plot(np.arange(RCP.ninputs) + 1, RCD.sites, "bo")
        # print('pairwise: ', pairwise)
        colormap = "plasma"
        if s_pair > 0.0:
            pclip = np.clip(pairwise, np.min(np.where(pairwise > 0)), np.max(pairwise))
            pclip[np.where(pclip == 0)] = np.nan
            pclipped = pclip - np.nanmin(pclip)
            sax["D"].scatter(
                pos[:, :, 0],
                pos[:, :, 1],
                s=200 * pairwise / maxp,
                c=pclipped,
                cmap=colormap,
            )
            vmax = np.nanmax(pclip) * 100
            vmin = np.nanmin(pclip) * 100
            # print("vmin, vmax: ", vmin, vmax)
        else:
            vmin = 0
            vmax = 1
        # sax['B'].plot(np.arange(RCP.ninputs)+1, participation/nspikes, 'gx')
        sax["E"].plot(RCD.sites, participation / nspikes, "gx")

        axcbar = PLS.create_inset_axes(
            [0.8, 0.05, 0.05, 0.5], sax["D"], label=str(P.axdict["D"])
        )
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        ticks = np.linspace(vmin, vmax, num=4, endpoint=True)
        cm_sns = mpl.cm.get_cmap(colormap)
        c2 = matplotlib.colorbar.ColorbarBase(
            axcbar, cmap=cm_sns, ticks=ticks, norm=norm
        )

        # PH.nice_plot(sax['C'], position=-0.2)
        sax["F"].plot(np.arange(len(ynspike)) + 1, ynspike, "m^-")

        sax["C"].set_ylim(bottom=0)
        sax["E"].set_ylim((0, 1.0))
        sax["E"].set_xlim(left=0)
        sax["F"].set_ylim(0, 1.05)
        sax["C"].set_ylabel("# Release Sites")
        sax["C"].set_xlabel("Input #")
        sax["E"].set_xlabel("# Release Sites")
        sax["E"].set_ylabel("Participation")
        sax["D"].set_ylabel("Input #")
        sax["D"].set_xlabel("Input #")
        sax["F"].set_xlabel(
            f"# Inputs in [{pre_w[0]:.1f} to {pre_w[1]:.1f}] before spike"
        )

        PH.cleanAxes(P.axarr.ravel())
        # PH.talbotTicks(sax["C"])

        PH.talbotTicks(
            sax["B"],
            tickPlacesAdd={"x": 1, "y": 2},
            floatAdd={"x": 1, "y": 2},
            pointSize=7,
        )
        PH.talbotTicks(
            sax["C"], pointSize=7,
        )

        PH.talbotTicks(
            sax["E"],
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            pointSize=7,
        )
        PH.talbotTicks(
            sax["F"],
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            pointSize=7,
        )
        PH.talbotTicks(
            axcbar,
            tickPlacesAdd={"x": 0, "y": 2},
            floatAdd={"x": 0, "y": 2},
            pointSize=7,
        )
        P.figure_handle.suptitle(
            f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}"
            + f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s}",
            fontsize=11,
        )
        mpl.show()

        return (summarySiteTC, RCD.sites)

    @time_func
    def plot_tuning(self, args, filename=None, filenames=None):
        PD = PData()
        changetimestamp = get_changetimestamp()
        channel = "nacncoop"
        channel = "klt"
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

        for ic, chdist in enumerate(truefile.keys()):
            args.dendritemode = truefile[chdist]
            if filename is None:
                # fna = select_filenames(filenames, args)
                fna = filenames
                print("\n plot data from: ".join([str(f) for f in fna]))
            else:
                fna = filename  # just one file
            fna = [Path(fn) for fn in fna if str(fn).find("ASA=") < 0]

            for k, fn in enumerate(fna):
                if str(fn).find(args.dendritemode.lower()) > 0:
                    self.plot_traces(P.axarr[0, ic], fn, PD, args.protocol, nax=k)

                    P.axarr[0, ic].text(50, 1.0, chdist, horizontalalignment="center")
            basen = fn.parts[-5]

            pngfile = Path(
                PD.renderpath, f"{basen:s}_{channel:s}_{truefile[chdist]:s}.png"
            )
            # print(pngfile)
            imaged = mpl.imread(pngfile)
            P.axarr[1, ic].imshow(imaged)
        return P

    def setup_PSTH(self):
        sizer = OrderedDict(  # define figure layout
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

        P = PH.arbitrary_grid(
            sizer, order="columnsfirst", label=True, figsize=(8.0, 6.0),
        )

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

    def vector_strength(self, spikes, freq, dmod, dB, experiment):
        """
        Calculate vector strength and related parameters from a spike train, for the specified frequency
        :param spikes: Spike train, in sec.
        :param freq: Stimulus frequency in Hz
        :return: a dictionary containing:
    
            r: vector strength
            n: number of spikes
            R: Rayleigh coefficient
            p: p value (is distribution not flat?)
            ph: the circularized spike train over period of the stimulus freq, freq, in radians
            d: the "dispersion" computed according to Ashida et al., 2010, etc.
        """
    
        per = 1/freq # c
        ph = 2*np.pi*np.fmod(spikes, per)/(per) # convert to radians within a cycle
        sumcos = np.sum(np.cos(ph))
        sumsin = np.sum(np.sin(ph))
        mean_phase = np.arctan2(sumsin,sumcos)
        sumc2 = sumcos**2
        sums2 = sumsin**2
        n = len(spikes)

        SC, syninfo = self.get_synaptic_info(f"VCN_c{int(self.parent.cellID):02d}")
        ninputs = len(syninfo[1])
        sites = np.zeros(ninputs)
        amax = 0.
        for isite in range(ninputs):  # precompute areas
            area = syninfo[1][isite][0]
            if area > amax:
                amax = area
            sites[isite] = int(np.around(area * SC.synperum2))

        if n == 0:
            return{'r': np.nan, 'n': 0, 'R': np.nan, 'p': 1.0, 'ph': np.nan, 'd': np.nan, 'amax': amax, 'ninputs': ninputs}
        vs = (1./n)*np.sqrt(sumc2+sums2)  # standard vector strength computation
        R = 2*n*vs*vs  # Raleigh coefficient
        Rp = np.exp(-n*vs*vs)  # p value for n > 50 (see Ashida et al. 2010).
        dx = np.sqrt(2.*(1-vs))/(2*np.pi*freq)
        self.spikes = spikes
        index_row = self.parent.selected_index_rows[0]
        selected = self.parent.table_manager.get_table_data(
            index_row
        )


        d = {'r': vs, 'n': n, 'R': R, 'p': Rp, 'ph': ph, 'd': dx, 'amax': amax, 'ninputs': ninputs}
        return d
        
    @winprint_continuous
    def print_VS(self, d, freq, dmod, dB, experiment):
        colnames  = f"Cell,Configuration,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,maxArea,ninputs"
        if 'an_vs' in d.keys():
            anvs = d['an_vs']
        else:
            anvs = np.nan
        print('anvs: ', anvs)
        line = f"{int(self.parent.cellID):d},{experiment:s},"
        line += f"{freq:.1f},{dmod:.1f},{dB:.1f},"
        line += f"{d['r']:.4f},"
        line += f"{d['n']:d},"
        line += f"{scipy.stats.circmean(d['ph'])/(2*np.pi):.4f},"
        line += f"{scipy.stats.circstd(d['ph'])/(2.*np.pi):.4f},"
        line += f"{d['R']:.4f},"
        line += f"{d['p']:.4e},"
        line += f"{anvs:.4f},"
        line += f"{d['amax']:.4f},"
        line += f"{d['ninputs']:d},"
        if self.firstline:
            self.textappend(colnames)
        self.textappend(line)
        # put on clipboard
        pyperclip.copy(line)
        pyperclip.paste()
        
        return d
    
    def psth_vs(self):
        """
        Generate tables of vs measures for all cells
        across the frequencies listed
        """
        self.textclear() # just at start
        PD = PData()
        if self.parent.selected_index_rows is None:
            return
        # P = self.PLT.setup_PSTH()
        P = None
        PD = PData()
        selrows = self.parent.table.selectionModel().selectedRows()
        for i,  index_row in enumerate(selrows):
            selected = self.parent.table_manager.get_table_data(index_row)  # table_data[index_row]
            if selected is None:
                return
            sfi = Path(selected.simulation_path, selected.files[0])
            if i == 0:
                self.firstline = True
            else:
                self.firstline = False
            self.plot_AN_response(P, sfi, PD, selected.runProtocol)
        self.firstline = True

    @winprint_continuous
    def plot_AN_response(self, P:Union[object, None]=None, 
            fn:Union[str, Path, None]=None, 
            PD: object=None,
            protocol:str=''):
        if P is None:
            plotflag = False
        else:
            plotflag = True
        changetimestamp = get_changetimestamp()
        x = self.get_data_file(fn, changetimestamp, PD)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if x is None:
            print("No simulation found that matched conditions")
            print(fn)
            return
        # unpack x
        par, stitle, ivdatafile, filemode, d = x
        AR, SP, RMA = self.analyze_data(ivdatafile, filemode, protocol)
        ntr = len(AR.traces)  # number of trials
        v0 = -160.0
        trstep = 25.0 / ntr
        inpstep = 5.0 / ntr
        sz = 50.0 / ntr
        # print(dir(AR))
        si = d["Params"]
        ri = d["runInfo"]
        gbc = f"VCN_c{int(self.parent.cellID):02d}"
        if plotflag:
            P.figure_handle.suptitle(
                f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}"
                + f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s}",
                fontsize=11,
            )

        if isinstance(si, dict):
            totaldur = (
                si["pip_start"]
                + np.max(si["pip_start"])
                + si["pip_duration"]
                + si["pip_offduration"]
            )
            soundtype = si["soundtype"]
            pip_start = si["pip_start"]
            pip_duration = si["pip_duration"]
            F0 = si["F0"]
            dB = (si["dB"],)
            fmod = (si["fmod"],)
            dmod = (si["dmod"],)
        else:
            totaldur = (
                ri.pip_start
                + np.max(ri.pip_start)
                + ri.pip_duration
                + ri.pip_offduration
            )
            pip_start = ri.pip_start
            pip_duration = ri.pip_duration
            soundtype = ri.soundtype
            F0 = ri.F0
            dB = ri.dB
            fmod = ri.fmod
            dmod = ri.dmod

        timescale = 1000.0
        ntr = len(AR.traces)
        all_an_st = []
        all_bu_st = []
        for i in range(ntr):  # for all trials in the measure.
            vtrial = AR.traces[i] * 1e3
            trd = d["Results"][i]
            waveform = trd["stimWaveform"].tolist()
            stb = trd["stimTimebase"]
            all_bu_st.append(trd["spikeTimes"])
            if i < 1 and plotflag:
                P.axdict["A"].plot(AR.time_base / 1000.0, vtrial, "k-", linewidth=0.5)
                spikeindex = [int(t / si.dtIC) for t in trd["spikeTimes"]]
                P.axdict["A"].plot(
                    np.array(trd["spikeTimes"]) / timescale,
                    vtrial[spikeindex],
                    "ro",
                    markersize=1.5,
                )
            if i == 0 and waveform is not None and plotflag:
                P.axdict["C"].plot(
                    stb, waveform, "k-", linewidth=0.5
                )  # stimulus underneath
            if plotflag:
                P.axdict["D"].plot(
                    np.array(trd["spikeTimes"]) / timescale,
                    i * np.ones(len(trd["spikeTimes"])),
                    "|",
                    markersize=1.5,
                    color="b",
                )
            inputs = len(trd["inputSpikeTimes"])
            for k in range(inputs):
                tk = trd["inputSpikeTimes"][k]
                all_an_st.extend(tk)
                y = (i + 0.1 + k * 0.05) * np.ones(len(tk))
                if i % 10 == 0 and plotflag:
                    P.axdict["G"].plot(
                        tk / 1000.0, y, "|", markersize=2.5, color="k", linewidth=0.5
                    )
        # get first spike latency information
        fsl, ssl = SPKANA.CNfspike(all_bu_st, ri.pip_start*1e3, ntr)
        index_row = self.parent.selected_index_rows[0]
        selected = self.parent.table_manager.get_table_data(
            index_row
        )
        print(f"'{int(self.parent.cellID):d}','{selected.synapseExperiment:s}',", end="")
        print(f"{np.nanmean(fsl):.3f},{np.nanstd(fsl):.3f},", end="")
        print(f"{np.nanmean(ssl):.3f},{np.nanstd(ssl):.3f}")
        all_bu_st_flat = functools.reduce(operator.iconcat, all_bu_st, [])
        all_bu_st = np.array(all_bu_st_flat) / 1000.0
        all_an_st = np.sort(np.array(all_an_st)) / 1000.0

        if soundtype == "tonepip":  # use panel F for FSL/SSL distributions
            sl_hbins = np.arange(0.0, 25.0, 0.5)
            if (len(fsl)) > 0.0 and plotflag:
                P.axdict["F"].hist(
                    fsl, bins=sl_hbins, facecolor="b", edgecolor="b", alpha=0.6
                )
            if (len(ssl)) > 0.0 and plotflag:
                P.axdict["F"].hist(
                    ssl, bins=sl_hbins, facecolor="r", edgecolor="r", alpha=0.6
                )
            if plotflag:
                P.axdict["F"].set_xlim(0.0, 25.0)
                P.axdict["F"].text(
                    1.0,
                    1.0,
                    f"FSL: {np.nanmean(fsl):8.3f} (SD {np.nanstd(fsl):8.3f} N={np.count_nonzero(~np.isnan(fsl)):3d})",
                    fontsize=7,
                    color="b",
                    fontfamily="monospace",
                    transform=P.axdict["F"].transAxes,
                    horizontalalignment="right",
                    verticalalignment="top",
                )
                P.axdict["F"].text(
                    1.0,
                    0.9,
                    f"SSL: {np.nanmean(ssl):8.3f} (SD {np.nanstd(ssl):8.3f} N={np.count_nonzero(~np.isnan(ssl)):3d})",
                    fontsize=7,
                    color="r",
                    fontfamily="monospace",
                    transform=P.axdict["F"].transAxes,
                    horizontalalignment="right",
                    verticalalignment="top",
                )
            # the histogram of the data
            bu_hbins = np.arange(
                0.0, np.max(AR.time / timescale), 1.0 / timescale
            )  # 1 ms msec bins
            st_hbins = np.arange(0.0, np.max(AR.time / 1000.0), 1e-3)  # 1 ms msec bins

            if len(all_bu_st) > 0:
                P.axdict["B"].hist(
                    all_bu_st, bins=bu_hbins, facecolor="k", edgecolor="k", alpha=1
                )
            else:
                P.axdict["B"].text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=P.axdict["C"].transAxes,
                    horizontalalignment="center",
                )
            if len(all_an_st) > 0:
                P.axdict["H"].hist(
                    all_an_st, bins=st_hbins, facecolor="k", edgecolor="k", alpha=1
                )
            else:
                P.axdict["H"].text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=P.axdict["H"].transAxes,
                    horizontalalignment="center",
                )

            for a in ["A", "B", "C", "D", "G", "H"]:  # set some common layout scaling
                P.axdict[a].set_xlim((0.0, totaldur))
            if a in ["F"] and soundtype == "SAM":
                P.axdict[a].set_xlim((0.0, totaldur))

        if soundtype == "SAM":  # calculate vs and plot histogram
            if len(all_bu_st) == 0 and plotflag:
                P.axdict["H"].text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=P.axdict["H"].transAxes,
                    horizontalalignment="center",
                )
            else:

                phasewin = [
                    pip_start + 0.25 * pip_duration,
                    pip_start + pip_duration,
                ]

                spkin = all_bu_st[np.where(all_bu_st > phasewin[0])]  # bu spikes
                spikesinwin = spkin[np.where(spkin <= phasewin[1])]
                an_spkin = all_an_st[np.where(all_an_st > phasewin[0])]  # matching AN spikes
                an_spikesinwin = an_spkin[np.where(an_spkin <= phasewin[1])]
                # set freq for VS calculation
                if soundtype == "tone":
                    print("Tone: F0=%.3f at %3.1f dbSPL, cell CF=%.3f" % (F0, dB, F0))
                if soundtype == "SAM" and plotflag:
                    tstring = (
                        "SAM Tone: F0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f"
                        % (F0, dB, fmod, dmod, F0,)
                    )
                    # print(tstring)
                    P.figure_handle.suptitle(tstring, fontsize=10)
                # print('spikes: ', spikesinwin)
                # print('fmod: ', fmod, '# spikes: ', len(spikesinwin),' spikesinwin: ', spikesinwin)
                vs = self.vector_strength(
                    spikesinwin, fmod, dmod, dB, ri.Spirou, # time must be in msec
                )  # vs expects spikes in msec
                print(" Sound type: ", soundtype)
                print(
                    "CN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
                    % (F0, vs["r"], vs["d"] * 1e6, vs["R"], vs["p"], vs["n"])
                )
                vs_an = self.vector_strength(
                    an_spikesinwin, fmod, dmod, dB, ri.Spirou, # time must be in msec
                )  # vs expects spikes in msec
                print(" Sound type: ", soundtype)
                print(
                    "AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
                    % (F0, vs_an["r"], vs_an["d"] * 1e6, vs_an["R"], vs_an["p"], vs_an["n"])
                )
                vs['an_vs'] = vs_an["r"] # copy over to array
                self.print_VS(vs, fmod, dmod, dB, ri.Spirou)
                # print(vs['ph'])
                if plotflag:
                    P.axdict["E"].hist(
                        vs["ph"],
                        bins=2 * np.pi * np.arange(30) / 30.0,
                        facecolor="k",
                        edgecolor="k",
                    )
                    P.axdict["E"].set_xlim((0.0, 2 * np.pi))
                    P.axdict["E"].set_title(
                        "Phase (VS = {0:.3f})".format(vs["r"]),
                        fontsize=9,
                        horizontalalignment="center",
                    )

    def select_filenames(self, fng, args) -> Union[list, None]:
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
    PS = PlotSims(parent=None)
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

    if args.analysis == "tuning":

        fng = []
        for ig, gbc in enumerate(PD.gradeA):
            basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
            pgbc = f"VCN_c{gbc:02d}"
            allf = list(Path(basefn).glob("*"))
            allfiles = sorted([f for f in allf if f.is_dir()])
            # print(allfiles)
            fnlast3 = allfiles[-3:]
        for f in fnlast3:
            fn = list(f.glob("*"))
            fng.append(fn[0])
        P = PS.plot_tuning(args, filename=None, filenames=fng)

    else:
        for ig, gbc in enumerate(PD.gradeA):
            basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
            pgbc = f"VCN_c{gbc:02d}"
            if args.protocol == "IV":
                name_start = f"IV_Result_VCN_c{gbc:02d}_inp=self_{modelName:s}*.p"
                args.experiment = None
            elif args.protocol == "AN":
                name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"

            # print(f"Searching for:  {str(Path(basefn, name_start)):s}")
            fng = list(Path(basefn).glob(name_start))
             # print(f"Found: {len(fng):d} files.")
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
                for k, fn in enumerate(fng):
                    if args.analysis == "traces":
                        PS.plot_traces(ax, fn, PD, args.protocol, nax=k)
                    elif args.analysis == "revcorr":
                        res = PS.compute_revcorr(ax, pgbc, fn, PD, args.protocol)
                        if res is None:
                            return
            elif args.analysis == "singles":
                fna = PS.select_filenames(fng, args)
                print(fna)

            elif args.analysis == "tuning":
                P = PS.plot_tuning(args, filename=None, filenames=fng)
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
    # for key, value in vargs.items():
    #      if key in parnames:
    #          # print('key: ', key)
    #          # print(str(value))
    #          exec(f"PD.{key:s} = {value!r}")
    #      # elif key in runnames:
    #      #     exec(f"runinfo.{key:s} = {value!r}")
    cmdline_display(args, PD)
    return (args, PD)


def main():
    args, PD = getCommands()


if __name__ == "__main__":
    main()
