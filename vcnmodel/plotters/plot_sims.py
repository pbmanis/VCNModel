import argparse
import dataclasses
import datetime
import functools
import itertools
import operator
import pickle
import string
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from functools import wraps
from pathlib import Path
from typing import List, Tuple, Union


import lmfit
import matplotlib.colorbar  # type: ignore
import matplotlib.colors  # type: ignore
from matplotlib import ticker
import matplotlib.cm
import numpy as np  # type: ignore
import pandas as pd
import pyperclip
import pyqtgraph as pg  # type: ignore
import quantities as pq
import rich as RI
from rich.text import Text
from rich.console import Console
import scipy.stats  # type: ignore
import seaborn as sns

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
from vcnmodel import cell_config as cell_config
from vcnmodel.analyzers import analysis as SPKANA
from vcnmodel.analyzers import spikestatistics as SPKS
from vcnmodel.analyzers import sttc as STTC
from vcnmodel.analyzers import sac as SAC
import vcnmodel.analyzers.reverseCorrelation as REVCORR
import vcnmodel.util.fixpicklemodule as FPM
from vcnmodel.analyzers import vector_strength as VS
import toml

config = toml.load(open("wheres_my_data.toml", "r"))

cprint = CP.cprint
console = Console()
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


    Supported primarily by R01DC015901 (Spirou, Manis, Ellisman),
    Early development: R01 DC004551 (Manis, until 2019)
    Later development: R01 DC019053 (Manis, 2020-2025)

"""


# make a shortcut for each of the clases
AR = vcnmodel.util.readmodel.ReadModel()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()
rc("text", usetex=True)
rc("text.latex", preamble=r"\usepackage{xcolor}")
rc("mathtext", fontset="stixsans")

sns.set_style(rc={"pdf.fonttype": 42})
mpl.style.use("~/.matplotlib/figures.mplstyle")

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

    return [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]


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
    basepath: str = config["baseDataDirectory"]
    renderpath: str = str(Path(config["codeDirectory"], "Renderings"))
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


@dataclass
class SACPars:
    twin: float = 5.0
    binw: float = 0.05
    delay: float = 0.0
    dur: float = 1000.0
    displayDuration: float = 10.0
    ntestrep: int = 20
    baseper: float = 4 / 3.0
    stimdur: float = 1000


@dataclass
class Inflate:
    soma_inflation: float = 0.0
    dendrite_inflation: float = 0.0


def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


def def_empty_dict():
    return {}


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
    si: dict = field(default_factory=def_empty_dict)
    ri: dict = field(default_factory=def_empty_dict)


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
    mean_pre_intervals: np.array = field(default_factory=def_empty_np)
    mean_post_intervals: float = 0.0
    max_coin_rate: float = 0
    participation: np.array = field(default_factory=def_empty_np)
    s_pair: float = 0.0
    ynspike: np.array = field(default_factory=def_empty_np)
    pre_w: list = field(default_factory=def_empty_list)
    pre_st: list = field(default_factory=def_empty_list)


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


class TraceCalls(object):
    """ Use as a decorator on functions that should be traced. Several
        functions can be decorated - they will all be indented according
        to their call depth.
        from: https://eli.thegreenplace.net/2012/08/22/easy-tracing-of-nested-function-calls-in-python
    """

    def __init__(
        self, stream=sys.stdout, indent_step=2, show_args=False, show_ret=False,
                show=True):
        self.stream = stream
        self.indent_step = indent_step
        self.show_ret = show_ret
        self.show_args = show_args
        self.show = show

        # This is a class attribute since we want to share the indentation
        # level between different traced functions, in case they call
        # each other.
        TraceCalls.cur_indent = 0

    def __call__(self, fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            indent = " " * TraceCalls.cur_indent
            if self.show_args:
                argstr = ", ".join(
                    [repr(a) for a in args]
                    + ["%s=%s" % (a, repr(b)) for a, b in kwargs.items()]
                )
                if self.show:
                    self.stream.write("-->%s%s(%s)\n" % (indent, fn.__name__, argstr))
            else:
                # self.stream.write('-->%s%s\n' % (indent, fn.__name__))
                text = Text.assemble(
                    (f"* {indent:s}-->{fn.__name__:s}\n", "bold yellow")
                )
                if self.show:
                    console.print(text)

            TraceCalls.cur_indent += self.indent_step
            ret = fn(*args, **kwargs)
            TraceCalls.cur_indent -= self.indent_step

            if self.show_ret and self.show:
                self.stream.write("%s returns--> %s\n" % (indent, ret))
            return ret

        return wrapper


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
    DOES NOT clear the print area first,
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
        self.VS = VS.VectorStrength()

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

    def _convert_params(self, d, filemode):
        """
        Capture the parameters in different versions of the files,
        and return as a simple dictionary.
        
        Parameters
        ----------
        d : dict from the top level of the file
            Must have 'runInfo' and/or 'Params' as entries
        
        filemode : str (no default)
            A string corresponding to the file mode for the file
            Must be either "vcnmodel.v0" or "vcnmodel.v1"
        """
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
        else:
            raise ValueError("File mode must be either vcnmodel.v0 or vcnmodel.v1")

        return par

    def _get_scaling(self, fn, PD, par_in):
        """
        Get the dendrite/soma scaling information from this data file
        and return a string. Also prints out the scaling information as we go.
        
        Parameters
        ----------
        PD : Parameter Dataclass (no default)
        par : parameter list.
        Returns
        -------
        string with the type of scaling applied to the data.
        """
        stitle = "Bare Scaling"
        if isinstance(par_in, dict):
            par = Inflate(par_in["soma_inflation"], par_in["dendrite_inflation"],)
        if PD.soma_inflate and PD.dend_inflate:
            if par.soma_inflation > 0.0 and par.dendrite_inflation > 0.0:
                ivdatafile = Path(fn)
                stitle = "Soma and Dend scaled"
                print(stitle)
        elif PD.soma_inflate and not PD.dend_inflate:
            if par.soma_inflation > 0.0 and par.dendrite_inflation < 0.0:
                ivdatafile = Path(fn)
                stitle = "Soma only scaled"
                print(stitle)
        elif PD.dend_inflate and not PD.soma_inflate:
            if par.soma_inflation < 0.0 and par.dendrite_inflation > 0.0:
                ivdatafile = Path(fn)
                stitle = "Dend only scaled"
                print(stitle)
        elif not PD.soma_autoinflate and not PD.dendrite_autoinflate:
            print("\nConditions x: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)
            ivdatafile = Path(fn)
            stitle = "No scaling (S, D)"
            print(stitle)
        else:
            ivdatafile = Path(fn)
            self.textappend(
                f"Bare file: no identified soma/dendrite inflation conditions", "red"
            )
            stitle = "Bare Scaling"
        return "      " + stitle

    @TraceCalls(show=False)
    def get_data_file(
        self, fn: Union[str, Path], changetimestamp: object, PD: dataclass, verbose=False
    ) -> Union[None, tuple]:
        """
        Get a data file, and also parse information from the file
        for display,
    
        Parameters
        ----------
        fn : str or Path
            the file to read
        changetimestamp : object
            the timestamp associated with this file.
        PD : dataclass
            
        """
        fnp = Path(fn)
        fns = str(fn)
        ivdatafile = Path(fn)
        if self.firstline:
            if not fnp.is_file():
                cprint("r", f"   File: {str(fnp):s} NOT FOUND")
                self.textappend(f"   File: {str(fnp):s} NOT FOUND", color="red")
                return None
            else:
                if verbose:
                    self.textappend(f"   File: {str(fnp):s} OK")
        mtime = fnp.stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if verbose and self.firstline:
            self.textappend(
                f"pgbcivr2: Checking file: {fnp.name:s} [{timestamp_str:s}]"
            )
        if mtime > changetimestamp:
            filemode = "vcnmodel.v1"
        else:
            filemode = "vcnmodel.v0"
        if self.firstline and verbose:
            self.textappend(f"pgbcivr2: file mode: {filemode:s}")
        with (open(fnp, "rb")) as fh:
            d = FPM.pickle_load(fh)
        if verbose:
            self.textappend(f"   ...File read, mode={filemode:s}")
        par = self._convert_params(d, filemode)
        stitle = self._get_scaling(fn, PD, par)

        if verbose:
            print("Read file: ", ivdatafile)
        if ivdatafile is None or not ivdatafile.is_file():
            if self.firstline and verbose:
                self.textappend(f"no file matching conditions : {str(ivdatafile):s}")
            return None

        if self.firstline and verbose:
            self.textappend(f"\npgbcivr2: datafile to read: {str(ivdatafile):s}")
        if isinstance(d["Results"], dict):
            if "time" in list(d["Results"].keys()):
                d["Results"] = self._data_flip(d["Results"])
        return par, stitle, ivdatafile, filemode, d

    @TraceCalls(show=False)
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

        for i, tr in enumerate(list(d["Results"].keys())):
            trd = d["Results"][tr]  # trial data
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            if 'inputSpikeTimes' in list(trd.keys()):
                for n in range(len(trd["inputSpikeTimes"])):  # for each sgc
                    RCD.npre += len(trd["inputSpikeTimes"][n])
                RCD.npost += len(trd["spikeTimes"])
                RCD.npost += len(RCD.st[tr])
            else:
                RCD.npost += sum(SP.spikes[i])
            RCD.st = SP.spikeIndices

        RCD.ti = time_base
        # print(f"Detected {RCD.npost:d} Post spikes")
        # print(f"Detected {RCD.npre:d} Presynaptic spikes")
        # print("# trials: ", RCP.ntrials)

        # clip trace to avoid end effects
        RCP.max_time = (
            np.max(RCD.ti) - RCP.min_time
        )  # this window needs to be at least as long as maxwin
        RCD.tx = np.arange(RCP.minwin, 0, RCP.binw)
        return (d, AR, SP, RMA, RCP, RCD)

    def _data_flip(self, d):
        """
        Convert from old data file format to new
        
        Parameters
        ----------
        d : dict from the top of the data file,
            or from results (d['Results'])
        
        Returns
        -------
        d['Results'] with the data format adjusted.
        """
        # flip order to put trials first (this was an old format)
        if "runInfo" in list(d.keys()):  # marker for top level
            data_res = d["Results"]
        else:
            data_res = d
        # if d is None:
        #     trials = range(len(data_res["somaVoltage"]))
        # else:
        #     trials = range(d["runInfo"].nReps)

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
        # print("_data_flip: ", data_res.keys())
        return data_res

    @time_func
    @TraceCalls(show=False)
    def analyze_data(
        self,
        ivdatafile: Union[Path, str],
        filemode: str,
        protocol: str,
        spike_shape=False,
    ) -> tuple:
        """
        Provides basic spike detection, spike shape analysis, and
        IV analysis if appropriate.
        We use ephys.acq4read.read_pfile to read the pickled data
        file into acq4 format, then we can use the ephys
        analysis tools to analyze the data
        
        Parameters
        ----------
        ivdatafile : path or str path to file
        filemode : str file mode (.v0 or .v1)
        protocol : str 
            name of the protocol
        
        Returns
        -------
        AR : Acq4 Read object
        SP : Spike analysis object
        RMA : Analysis summary object
        """
        AR.read_pfile(ivdatafile, filemode=filemode)

        bridge_offset = 0.0
        threshold = -0.035  # V
        tgap = 0.0  # gap before fittoign taum
        RM.setup(AR.MC, SP, bridge_offset=bridge_offset)
        SP.setup(
            clamps=AR.MC,
            threshold=threshold,  # in units of V
            refractory=0.001,  # In units of sec
            peakwidth=0.001,  # in units of sec
            interpolate=True,
            verify=True,
            mode="peak",
            data_time_units="ms",  # pass units for the DATA
            data_volt_units="V",
        )

        SP.set_detector("Kalluri")  # spike detector: argrelmax, threshold, Kalluri
        SP.analyzeSpikes()
        if spike_shape:
            SP.analyzeSpikeShape()
        RMA = None
        print(SP.analysis_summary)
        SP.analysis_summary['pulseDuration'] = 0.1
        if protocol == "IV":
            SP.fitOne(function="fitOneOriginal")
            RM.analyze(
                rmpregion=[0.0, AR.MC.tstart - 0.001],
                tauregion=[
                    AR.MC.tstart,
                    AR.MC.tstart + (AR.MC.tend - AR.MC.tstart) / 5.0,
                ],
                to_peak=True,
                tgap=tgap,
            )

            RMA = RM.analysis_summary
        return AR, SP, RMA

    def simple_plot_traces(
        self, rows=1, cols=1, width=5.0, height=4.0, stack=True, ymin=-120.0, ymax=0.0
    ):
        self.P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            figsize=(width, height),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.01,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.03,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )

        PD = PData()
        for iax, index_row in enumerate(self.parent.selected_index_rows):
            selected = self.parent.table_manager.get_table_data(index_row)
            if selected is None:
                return
            sfi = Path(selected.simulation_path, selected.files[0])
            if stack:
                self.plot_traces(
                    self.P.axarr[iax, 0],
                    sfi,
                    PD,
                    protocol=selected.runProtocol,
                    ymin=ymin,
                    ymax=ymax,
                    iax=iax,
                    figure=self.P.figure_handle,
                )
            else:
                self.plot_traces(
                    self.P.axarr[0, iax],
                    sfi,
                    PD,
                    protocol=selected.runProtocol,
                    ymin=ymin,
                    ymax=ymax,
                    iax=iax,
                    figure=self.P.figure_handle,
                )
        self.P.figure_handle.show()

    @winprint
    def print_file_info(self, selected, mode="list"):
        if mode not in ["list", "dict"]:
            raise ValueError()
        self.textappend("For copy into figure_data.py: ")
        if mode == "dict":
            br = "{}"
            self.textappend(f"{int(self.parent.cellID):d}: {br[0]:s}")
        if mode == "list":
            br = "[]"
            self.textappend(f"    {int(self.parent.cellID):d}: {br[0]:s}")
        for sel in selected:
            data = self.parent.table_manager.get_table_data(sel)
            fn = Path(data.files[0])
            fnr = str(fn.parts[-2])
            fkey = data.dendriteExpt
            if mode == "dict":
                self.textappend(f'    "{fkey:s}": "{fnr:s}",')
            if mode == "list":
                self.textappend(f'        "{fnr:s}",')
        if mode == "dict":
            self.textappend(f"{br[1]:s},")
        if mode == "list":
            self.textappend(f"    {br[1]:s},")

    @winprint
    @TraceCalls()
    def plot_traces(
        self,
        ax: object,
        fn: Union[Path, str],
        PD: dataclass,
        protocol: str,
        ymin=-80.0,
        ymax=20.0,
        xmin=0.0,
        xmax=None,
        yoffset: float = 0.0,
        iax=None,
        nax=0,
        rep=None,  # which rep : none is for all.
        figure=None,
        longtitle=True,
        trace_color="k",
        ivaxis=None,
        ivcolor="k",
        iv_spike_color="r",
        spike_marker_size=2.5,
        spike_marker_color="r",
        calx=0.0,
        caly=0.0,
        calt=10.0,
        calv=20.0,
        calx2=20.0,
        caly2=20.0,
        calt2=10.0,
        calv2=10.0,
    ) -> tuple:
        changetimestamp = get_changetimestamp()
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        inx = str(fn).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(str(fn)[inx + 4 : inx + 7])
        if protocol in ["IV", "runIV"]:
            protocol = "IV"
        elif protocol in ["VC", "runVC"]:
            protocol = "VC"

        d, AR, SP, RMA, RCP, RCD = self.get_data(fn, PD, changetimestamp, protocol)
        si = d["Params"]
        ri = d["runInfo"]
        if figure is None:  # no figure... just analysis...
            return AR, SP, RMA

        cprint("c", "plot_traces: reparing for plot")
        ntr = len(AR.MC.traces)  # number of trials
        v0 = -160.0
        if isinstance(ri, dict):
            deadtime = ri["stimDelay"]
        else:
            deadtime = ri.stimDelay
        trstep = 25.0 / ntr
        inpstep = 2.0 / ntr
        sz = 50.0 / ntr
        noutspikes = 0
        ninspikes = 0
        ispikethr = None
        spike_rheobase = None
        if isinstance(ax, list):
            ax1 = ax[0]
            ax2 = ax[1]
        elif hasattr("ax", "len") and len(ax) == 2:
            ax1 = ax[0]
            ax2 = ax[1]
        elif hasattr("ax", "len") and len(ax) == 1:
            ax1 = ax
            ax2 = None
        elif not hasattr("ax", "len"):
            ax1 = ax
            ax2 = None
        for trial, icurr in enumerate(d["Results"]):
            if rep is not None and trial != rep:
                continue
            # cprint('c', f"trial: {trial:d}")
            AR.MC.traces[trial][0] = AR.MC.traces[trial][1]
            if protocol in ["VC", "vc", "vclamp"]:
                AR.MC.traces[trial] = AR.MC.traces[trial].asarray() * 1e9  # nA

                cmd = AR.MC.cmd_wave[trial] * 1e3  # from V to mV
            else:
                AR.MC.traces[trial] = AR.MC.traces[trial].asarray() * 1e3  # mV
                cmd = AR.MC.cmd_wave[trial] * 1e9  # from A to nA
            ax1.plot(
                AR.MC.time_base,
                AR.MC.traces[trial] + yoffset,
                linestyle="-",
                color=trace_color,
                linewidth=0.5,
            )
            if ax2 is not None:
                ax2.plot(AR.MC.time_base, cmd, linewidth=0.5)
            if "spikeTimes" in list(d["Results"][icurr].keys()):
                # cprint('r', 'spiketimes from results')
                # print(d["Results"][icurr]["spikeTimes"])
                #  print(si.dtIC)
                spikeindex = [
                    int(t * 1e3 / (si.dtIC)) for t in d["Results"][icurr]["spikeTimes"]
                ]
            else:
                # cprint('r', 'spikes from SP.spikeIndices')
                spikeindex = SP.spikeIndices[trial]
            # print(f"Trial: {trial:3d} Nspikes: {len(spikeindex):d}")
            ax1.plot(
                AR.MC.time_base[spikeindex],
                AR.MC.traces[trial][spikeindex] + yoffset,
                "o",
                color=spike_marker_color,
                markerfacecolor=spike_marker_color,
                markersize=spike_marker_size,
            )
            sinds = np.array(spikeindex) * AR.MC.sample_rate[trial]
            # print('sinds: ', sinds, deadtime, ri.stimDelay, ri.stimDur)
            nspk_in_trial = len(np.argwhere(sinds > deadtime))
            if (
                nspk_in_trial > 0
                and ispikethr is None
                and sinds[0] < (ri.stimDelay + ri.stimDur)
            ):
                # cprint('c', f"Found threshold spike:  {icurr:.2f}, {trial:d}")
                spike_rheobase = icurr
                ispikethr = trial
            noutspikes += nspk_in_trial
            if protocol in ["AN", "runANSingles"]:
                if trial in list(d["Results"].keys()) and "inputSpikeTimes" in list(
                    d["Results"][icurr].keys()
                ):
                    spkt = d["Results"][icurr]["inputSpikeTimes"]
                elif "inputSpikeTimes" in list(d["Results"].keys()):
                    spkt = d["Results"]["inputSpikeTimes"][trial]
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

                ax.set_ylim(ymin, ymax)
                if xmin is None:
                    xmin = 0.050
                if xmax is None:
                    xmax = np.max(AR.MC.time_base)
                ax.set_xlim(xmin, xmax)
            elif protocol in ["VC", "vc", "vclamp"]:
                pass  #
                # ax.set_ylim((-100.0, 100.0))
            else:
                ax.set_ylim(ymin, ymax)
                if xmin is None:
                    xmin = 0.050
                if xmax is None:
                    xmax = np.max(AR.MC.time_base)
                ax.set_xlim(xmin, xmax)
        ftname = str(Path(fn).name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = ""
        if longtitle:
            toptitle = f"{ftname:s}"
        else:
            toptitle = si.dendriteExpt
        if protocol in ["IV"]:
            cprint("r", f"RMA taum: {RMA['taum']:.2f}")
            toptitle += (
                f"\nRin={RMA['Rin']:.1f} M$\Omega$  $\\tau_m$={RMA['taum']:.2f} ms"
            )

            if iax == 2 and calx is not None:
                PH.calbar(
                    ax,
                    calbar=[calx, caly, calt, calv],
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
            if ivaxis is None:
                secax = PLS.create_inset_axes(
                    [0.45, -0.05, 0.3, 0.3], ax, label=str(ax)
                )
                color = "k"
                ticklabelsize = 6
            else:
                secax = ivaxis
                color = ivcolor
                ticklabelsize = 8
            secax.plot(
                RM.ivss_cmd_all * 1e9,
                RM.ivss_v_all * 1e3,
                f"{color:s}s-",
                markersize=3,
                markerfacecolor="k",
                zorder=10,
                clip_on=False,
            )

            ltz = np.where(RM.ivss_cmd_all <= 0.0)[0]

            secax.plot(
                RM.ivpk_cmd_all[ltz] * 1e9,
                RM.ivpk_v_all[ltz] * 1e3,
                f"{color:s}o-",
                markersize=3,
                markerfacecolor="w",
                zorder=10,
                clip_on=False,
            )
            if ispikethr is not None:  # decorate spike threshold point
                secax.plot(
                    RM.ivss_cmd_all[ispikethr] * 1e9,
                    RM.ivss_v_all[ispikethr] * 1e3,
                    "o",
                    markersize=3,
                    color=iv_spike_color,
                    markerfacecolor=iv_spike_color,
                    zorder=100,
                    clip_on=False,
                )
            PH.crossAxes(
                secax,
                xyzero=[0.0, -60.0],
                limits=[
                    np.min(RM.ivss_cmd_all) * 1e9,
                    -120,
                    np.max(RM.ivss_cmd_all) * 1e9,
                    -25.0,
                ],  #
            )
            PH.talbotTicks(
                secax,
                axes="xy",
                density=(1.0, 1.0),
                insideMargin=0.02,
                pointSize=ticklabelsize,
                tickPlacesAdd={"x": 1, "y": 0},
                floatAdd={"x": 1, "y": 0},
            )
        elif protocol in ["VC", "vc", "vclamp"]:
            maxt = np.max(AR.MC.time_base)
            if calx is None:
                tlen = 10  # ms
                calx = maxt - tlen

            PH.calbar(
                ax1,
                calbar=[calx, caly, calt, calv],
                orient="left",
                unitNames={"x": "ms", "y": "nA"},
                fontsize=9,
            )
            if ax2 is not None:
                PH.calbar(
                    ax2,
                    calbar=[calx2, caly2, calt2, calv2],
                    scale=[1.0, 1.0],
                    axesoff=True,
                    orient="left",
                    unitNames={"x": "ms", "y": "mV"},
                    fontsize=9,
                    weight="normal",
                    color="k",
                    font="Arial",
                )
        else:
            if nax == 0 or calx is not None:
                PH.calbar(
                    ax,
                    calbar=[calx, caly, calt, calv],
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
        cprint("y", toptitle)
        cprint("m", figure)
        cprint("m", ax)
        return (synno, noutspikes, ninspikes)

    @TraceCalls()
    def plot_VC(
        self,
        selected_index_rows=None,
        sfi=None,
        parent_figure=None,
        loc=None,
        show=True,
    ):
        if selected_index_rows is not None:
            n_columns = len(selected_index_rows)
            if n_columns == 0:
                return
            sfi = [None for i in range(n_columns)]
            protocol = [None for i in range(n_columns)]
            dendriteMode = [None for i in range(n_columns)]
            for i in range(n_columns):
                index_row = selected_index_rows[i]
                selected = self.parent.table_manager.get_table_data(index_row)
                sfi[i] = Path(selected.simulation_path, selected.files[0])
                protocol[i] = selected.runProtocol
                dendriteMode[i] = selected.dendriteMode

        else:
            n_columns = len(sfi)
            protocol = [None for i in range(n_columns)]
            dendriteMode = [None for i in range(n_columns)]
            for i in range(n_columns):
                with open(sfi[i], "rb") as fh:
                    sfidata = FPM.pickle_load(fh)  # , encoding='latin-1')
                    # print(sfidata['Params'])
                    dendriteMode[i] = sfidata["Params"].dendriteMode
                    protocol[i] = sfidata["runInfo"].runProtocol
        P = self.setup_VC_plots(n_columns, parent_figure=parent_figure, loc=loc)
        titlemap = {"normal": "Half-active", "passive": "Passive", "active": "Active"}

        for i in range(n_columns):
            PD = PData()
            trace_ax = P.axarr[i * 3, 0]
            cmd_ax = P.axarr[i * 3 + 1, 0]
            self.plot_traces(
                ax=[trace_ax, cmd_ax],
                fn=sfi[i],
                PD=PD,
                protocol=protocol[i],
                calx=120.0,
                caly=5.0,
                calt=10.0,
                calv=2.0,
                calx2=120.0,
                caly2=-40.0,
                calt2=10.0,
                calv2=20.0,
                figure=P.figure_handle,
            )
            trace_ax.set_ylim((-1, 15))
            self.analyzeVC(P.axarr[i * 3 + 2, 0], sfi[i], PD, protocol=protocol[i])
            trace_ax.text(
                0.5,
                1.1,
                titlemap[dendriteMode[i]],
                transform=trace_ax.transAxes,
                fontsize=9,
                verticalalignment="bottom",
                horizontalalignment="center",
            )
            vref = -80.0
            PH.referenceline(cmd_ax, vref)
            cmd_ax.text(
                -5.25,
                vref,
                f"{int(vref):d} mV",
                fontsize=9,
                color="k",
                # transform=secax.transAxes,
                horizontalalignment="right",
                verticalalignment="center",
            )
        if show:
            P.figure_handle.show()
        return P

    @TraceCalls()
    def setup_VC_plots(
        self, n_columns, parent_figure=None, loc: Union[None, tuple] = None
    ):
        sizer = OrderedDict()
        tmar = 0.0
        bmar = 0.0
        if loc is not None:
            tmar = loc[3]
            bmar = loc[2]
        left = 0.5
        right = 0.5
        hspc = 0.5
        if parent_figure is not None:
            figsize = parent_figure.figsize
        else:
            figsize = [8.0, 5.0]
        cwid = (figsize[0] - (left + right) - hspc * (n_columns - 1)) / n_columns
        for i, col in enumerate(range(n_columns)):
            panel_letter = string.ascii_uppercase[i]
            xl = left + i * (cwid + hspc)
            sizer[panel_letter + "1"] = {
                "pos": [xl, cwid, 2.65 + bmar, 1.8],
                "labelpos": (-0.15, 1.0),
            }
            sizer[panel_letter + "2"] = {
                "pos": [xl, cwid, 2.05 + bmar, 0.4],
                "labelpos": (-0.15, 1.0),
            }
            sizer[panel_letter + "3"] = {
                "pos": [xl, cwid, 0.4 + bmar, 1.35],
                "labelpos": (-0.15, 1.05),
            }

        P = PH.arbitrary_grid(
            sizer,
            units="in",
            order="columnsfirst",
            label=True,
            showgrid=False,
            figsize=figsize,
            parent_figure=parent_figure,
        )
        P.figure_handle.show()
        return P

    @winprint
    @time_func
    @TraceCalls()
    def analyzeVC(
        self, ax: object, fn: Union[Path, str], PD: dataclass, protocol: str,
    ) -> tuple:
        changetimestamp = get_changetimestamp()
        # x = self.get_data_file(fn, changetimestamp, PD)
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
        res = self.get_data(fn, PD, changetimestamp, protocol)
        if res is None:
            return None
        (d, AR, SP, RMA, RCP, RCD) = res
        # AR, SP, RMA = self.analyze_data(ivdatafile, filemode, protocol)
        tss = [0, 0]
        sr = AR.MC.sample_rate[0]
        tss[0] = int(
            AR.MC.tstart / sr + (AR.MC.tdur / 2.0) / sr
        )  # wrong duration stored in traces - need to fix.
        tss[1] = int(AR.MC.tstart / sr + (AR.MC.tdur / 1.0) / sr)
        I = AR.MC.traces.asarray()
        V = AR.MC.cmd_wave
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
                V[j, 0 : int((AR.MC.tstart - 1.0) / AR.MC.sample_rate[0])]
            )  # resting potential - for 1 msec prior to step
            i0[j] = np.mean(
                I[j, 0 : int((AR.MC.tstart - 1.0) / AR.MC.sample_rate[0])]
            )  # resting/holding current - for 1 msec prior to step
            # ax[1].plot(AR.MC.time, V[j, :], "k-")
            # ax[0].plot(AR.MC.time[tss[0]:tss[1]], I[j,tss[0]:tss[1]], 'r--')

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
        sr = AR.MC.sample_rate[mHypStep]
        t0 = int(AR.MC.tstart / sr)

        pts = int(5.0 / sr)  # fit over first 5 msec
        tfit = AR.MC.time_base[t0 + 1 : t0 + pts] - AR.MC.time_base[t0]
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
        tfit2 = AR.MC.time_base[t0 : t0 + pts] - AR.MC.time_base[t0]
        bfit = expmodel.eval(params=exp_result.params, x=tfit2)
        # ax[0].plot(
        #     tfit2 + AR.MC.time_base[t0], bfit * 1e9, "r-", dashes=[6, 2], linewidth=1
        # )
        # ax[0].plot(tfit+AR.MC.time_base[t0], exp_result.best_fit*1e9, 'r-')
        # ax[0].plot(tfit+AR.MC.time_base[t0], ifit*1e9, 'g--')
        ax.plot(vss_sel * 1e3, gss_sel * 1e9, "ko", markersize=2)
        ax.plot(vss_sel * 1e3, result.best_fit * 1e9, "r-")
        ax.set_ylim(0, 100.0)
        ax.set_xlim(-100.0, 40.0)
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            # pointSize=10,
            tickPlacesAdd={"x": 0, "y": 0},
            floatAdd={"x": 0, "y": 0},
        )
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
        cprint("c", textstr)
        ax.text(
            1.0,
            0.05,
            textstr,
            transform=ax.transAxes,
            size=8,
            verticalalignment="bottom",
            horizontalalignment="right",
            bbox=props,
        )
        ax.set_xlabel("V (mV)")
        ax.set_ylabel("I (nA)")
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

        self.ntr = len(AR.MC.traces)  # number of trials
        self.pwin = None
        self.pfull = None
        self.lines = [None] * self.parent.n_trace_sel
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
                        spk.dt * np.arange(0, len(spk.waveform))[:-1] + spk.start_win,
                        np.diff(spk.waveform / spk.dt),
                        pen=color,
                    )
                papertape_yoff = 120.0
                spkin_off = -50.0
                dy = 5.0
            else:
                if first:
                    print("plotting first")
                    tx = spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win
                    # print(np.min(tx), np.max(tx), tx.shape)
                    #                    print(np.min(spk.waveform), np.max(spk.waveform), spk.waveform.shape)
                    # print(dir(self.parent.trace_plots))
                    self.parent.trace_plots.setEnabled(True)
                    self.parent.trace_plots.plotItem.plot(
                        np.arange(10), np.ones(10) * ix, pen="r"
                    )
                    self.lines[ix] = self.parent.trace_plots.plotItem.plot(
                        x=tx,
                        y=spk.waveform,
                        # pen=color,
                    )
                    self.lines[ix].curve.show()
                    # print(dir(self.lines[ix]))

                else:
                    self.check_yscaling()
                    # self.parent.trace_plots.plot(
                    #     spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win,
                    #     spk.waveform,
                    #     pen=color,
                    # )
                    self.parent.trace_plots.setEnabled(True)
                    u = self.parent.trace_plots.plotItem.plot(
                        np.arange(10), np.ones(10) * ix, pen="b"
                    )
                    tx = spk.dt * np.arange(0, len(spk.waveform)) + spk.start_win
                    self.lines[ix].curve.setData(
                        x=tx, y=spk.waveform, pen="w",  # pen=color,
                    )
                    self.lines[ix].curve.show()
                papertape_yoff = 0.0
                spkin_off = -65.0
                dy = 2.0

    def make_patch_spines_invisible(self, axn):
        axn.set_frame_on(True)
        axn.patch.set_visible(False)
        for sp in axn.spines.values():
            sp.set_visible(False)

    def plot_revcorr2(
        self,
        P: object,
        PD: dataclass,
        RCP: dataclass,
        RCD: dataclass,
        axarray=None,
        calbar_show=True,
        calbar_fontsize=11,
        yaxis_label=True,
        start_letter="A",
        colormap="Set3",
        max_inputs=12,
        synlabel: Union[bool, None] = None,
        cbar_vmax: float = 300.0,
    ):
        sns.set_style("ticks")

        # secax = twinax(P.figure_handle, ax, pos=maxwin)
        str_a = string.ascii_uppercase
        p_labels = str_a[str_a.find(start_letter) : str_a.find(start_letter) + 4]
        # print('revcorr2: sax: ', sax)
        # print('plabels', p_labels)
        if axarray is None:
            sax = P.axdict
            sax0 = sax[p_labels[0]]
            sax1 = sax[p_labels[1]]
            sax2 = sax[p_labels[2]]
            sax3 = sax[p_labels[3]]
            ax = sax1
            secax = sax0
        else:
            ax = axarray[0]
            secax = axarray[1]

        secax.set_facecolor((1, 1, 1, 0))
        ax.spines["top"].set_visible(False)
        summarySiteTC = {}
        RCD.max_coin_rate = 0.0
        maxrevcorr = 0
        colors = [None] * RCP.ninputs
        linehandles = [None] * RCP.ninputs
        cmx = sns.color_palette(colormap, as_cmap=True)

        cell_n = int(RCP.si.cellID[-2:])
        SC, syninfo = self.get_synaptic_info(cell_n)
        syn_ASA = np.array([syninfo[1][isite][0] for isite in range(RCP.ninputs)])
        max_ASA = np.max(syn_ASA)

        for isite in reversed(range(RCP.ninputs)):
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
                # color = mpl.cm.viridis(norm(RCD.sites, isite))
                # cmx = mpl.cm.get_cmap(colormap)
                # color = mpl.cm.brg(norm(RCD.sites, isite))
                # print("index: ", int(cmx.N*syn_ASA[isite]/max_ASA))
                # print("cmx.N: ", cmx.N)
                color = cmx.colors[int(cmx.N * syn_ASA[isite] / cbar_vmax) - 1]
                colors[isite] = color
                maxrevcorr = np.max((maxrevcorr, np.max(RCD.CB[isite])))
                totalrevcorr = np.sum(RCD.CB[isite])
                # print("Using algorithm:", RCP.algorithm)
                label = "Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))
                label = None
                if isite in range(RCP.ninputs):
                    if RCP.algorithm == "RevcorrSPKS":
                        linehandles[isite] = ax.plot(
                            RCD.tx,
                            RCD.C[isite][:nc],
                            color=color,
                            label=label,
                            linewidth=0.9,
                            zorder=5,
                        )
                        RCD.max_coin_rate = np.max(
                            (RCD.max_coin_rate, np.max(RCD.C[:nc]))
                        )
                        # mpl.figure()
                        # mpl.plot(RCD.tx, RCD.C[isite][:nc])

                    elif RCP.algorithm == "RevcorrSimple":
                        linehandles[isite] = ax.plot(
                            RCD.CBT - RCP.maxwin,
                            RCD.CB[isite] / float(RCD.npost_spikes),
                            color=color,
                            label=label,
                            linewidth=0.9,
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
                        linehandles[isite] = ax.plot(
                            RCD.tx,
                            RCD.STTC[isite][:nc],
                            color=color,
                            label=label,
                            linewidth=0.9,
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
                if calbar_show:
                    PH.calbar(
                        secax,
                        calbar=[1.35, -30, 1.0, 20.0],
                        axesoff=True,
                        orient="left",
                        unitNames={"x": "ms", "y": "mV"},
                        fontsize=calbar_fontsize,
                    )
                    secax.text(
                        -5.25,
                        -60.0,
                        "-60 mV",
                        fontsize=7,
                        color="k",
                        # transform=secax.transAxes,
                        horizontalalignment="right",
                        verticalalignment="center",
                    )
                PH.referenceline(secax, -60.0)
                PH.noaxes(secax)
        print(f"Total spikes: {RCD.nsp_avg:d}")

        # finally, put up a legend that shows which inputs map to what colors
        # # the colors are in the legend
        PH.nice_plot(ax, position=-0.05)

        if synlabel is None or synlabel is True:
            self.axins = PH.make_colorbar(
                ax,
                bbox=[0.05, 0.8, 0.8, 0.2],
                vmin=0,
                vmax=300.0,
                nticks=4,
                palette=colormap,
            )
            self.axins.tick_params(axis="x", length=2, labelsize=6)

        if yaxis_label:
            ax.set_ylabel("Coinc. Rate (Hz)")
        ax.set_xlabel("T (ms)")
        if RCD.max_coin_rate > 0.0:
            ns = PH.NiceScale(0.0, RCD.max_coin_rate)
            ax.set_ylim((0, ns.niceMax))
        else:
            ax.set_ylim((0, 0.25))
        ax.set_xlim((RCP.minwin, RCP.maxwin))
        # PH.talbotTicks(
        #     ax,
        #     axes="xy",
        #     # density=(5.0, 2.0),
        #     insideMargin=0.05,
        #     # pointSize=10,
        #     tickPlacesAdd={"x": 1, "y": 1},
        #     floatAdd={"x": 1, "y": 1},
        # )
        secax.set_ylim([-70.0, 10.0])
        secax.set_xlim((RCP.minwin, RCP.maxwin))
        secax.tick_params(direction="in", length=3.0, width=1.0)
        return summarySiteTC

    def get_synaptic_info(self, gbc: str) -> tuple:
        SC = cell_config.CellConfig()
        # print(SC.VCN_Inputs)
        if isinstance(gbc, int):
            gbc_string = f"VCN_c{gbc:02d}"
        elif isinstance(gbc, str) and len(gbc) <= 2:
            gbc_string = f"VCN_c{int(gbc):02d}"
        elif isinstance(gbc, str) and gbc.startswith("VCN_"):
            gbc_string = gbc
        syninfo = SC.VCN_Inputs[gbc_string]
        return (SC, syninfo)

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
        self.parent.selected_index_rows = (
            self.parent.table.selectionModel().selectedRows()
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
                label=f"{selected.dendriteMode:s} {selected.dendriteExpt:s} {selected.synapseExperiment}",
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
        """
        Compute reverse correlation from data
        
        """
        changetimestamp = get_changetimestamp()
        self.allspikes = None

        #
        # 1. Gather data
        #
        print(f"Getting data for gbc: {gbc:s}")
        SC, syninfo = self.get_synaptic_info(gbc)
        res = self.get_data(fn, PD, changetimestamp, protocol)
        if res is None:
            return None
        (d, AR, SP, RMA, RCP, RCD) = res
        PD.thiscell = gbc
        si = d["Params"]
        ri = d["runInfo"]
        RCP.si = si
        RCP.ri = ri
        RCP.algorithm = revcorrtype
        if RCP.algorithm == "RevcorrSTTC":
            sttccorr = STTC.STTC()

        print("Preparing for computation", gbc)
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

        """
         2. set up parameters
        """

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
        RCD.mean_post_interval = 0.0
        RCD.mean_pre_intervals = np.zeros(RCP.ninputs)
        RCD.nsp_avg = 0
        nspk_plot = 0
        # spksplotted = False
        RCP.min_time = min_time  # driven window without onset
        RCP.max_time = max_time

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

        #
        # sum across trials, and sort by inputs
        #
        # print("starting loop")
        post_intervals = []
        pre_intervals = [[] for x in range(RCP.ninputs)]
        start_time = datetime.datetime.now()
        srate = (
            si.dtIC * 1e-3
        )  # this needs to be adjusted by the date of the run, somewhere...
        for trial in range(RCP.ntrials):  # sum across trials
            # spiketimes is in msec, so leave si.dtIC in msec
            spikeindex = [int(t / (srate)) for t in d["Results"][trial]["spikeTimes"]]
            stx = AR.MC.time_base[spikeindex]
            stx = stx[  # get postsynaptic spikes and trim to analysis window
                (stx > RCP.min_time) & (stx < RCP.max_time)
            ]
            if len(stx) == 0:
                continue
            RCD.npost_spikes += len(stx)
            post_intervals.extend(np.diff(stx))
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
                # print(len( d["Results"][trial]["inputSpikeTimes"]), isite)
                anx = d["Results"][trial]["inputSpikeTimes"][
                    isite
                ]  # get input AN spikes and trim list to window
                anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
                if len(anx) == 0:
                    continue
                RCD.npre_spikes += len(anx)  # count up pre spikes.
                pre_intervals[isite].extend(
                    np.diff(anx)
                )  # keep track of presynapit intervals by input
                if revcorrtype == "RevcorrSPKS":
                    RCD.C[isite] += SPKS.correlogram(
                        stx * pq.ms,
                        anx * pq.ms,
                        width=-RCP.minwin * pq.ms,
                        bin_width=RCP.binw * pq.ms,
                        T=None,
                    )

                elif revcorrtype == "RevcorrSimple":
                    refcorr, npost = REVCORR.reverse_correlation(
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
                    print(len(refcorr))
                    print(len(RCD.STTC[isite]))
                    print(RCD.STTC[isite])
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
        RCD.mean_pre_intervals = [0] * RCP.ninputs
        RCD.pre_st = [[]] * RCP.ninputs
        print(len(pre_intervals))
        # print(pre_intervals)
        for isite in range(RCP.ninputs):
            RCD.mean_pre_intervals[isite] = np.mean(pre_intervals[isite])
            RCD.pre_st[isite] = pre_intervals[isite]
        RCD.mean_post_intervals = np.mean(post_intervals)
        elapsed_time = datetime.datetime.now() - start_time
        print("Time for calculation: ", elapsed_time)
        RCD.pre_w = [-2.7, -0.5]
        # RCD.pre_w = [-3.75, -0.1]

        ################
        # now some pariwise and participation stats on input events prior to a spike
        ################
        spikedata = SpikeData()  # storage in a dataclass
        self.allspikes = []
        RCD.pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
        RCD.participation = np.zeros(RCP.ninputs)
        pre_spike_counts = np.zeros(
            RCP.ninputs + 1
        )  # there could be 0, or up to RCP.ninputs pre spikes
        pre_solo_spikes = np.zeros(RCP.ninputs + 1)

        nperspike = []
        nfilt_spikes = 0
        nfilt2_spikes = 0
        filttable = []
        filttable2 = []
        RCD.nspikes = 0
        sellist = [True] * RCP.ninputs
        # if ri.Spirou == "largestonly":
        #     for i in range(1, len(sellist)):
        #         sellist[i] = False
        # elif ri.Spirou == "twolargest":
        #     for i in range(2, len(sellist)):
        #         sellist[i] = False
        # elif ri.Spirou == "removelargest":
        #     sellist[0] = False
        # print('RCP.ninputs: ', RCP.ninputs, RCD.sites)
        for trial in range(RCP.ntrials):  # accumulate across all trials
            # spiketimes is in msec, so si.dtIC should be in msec
            spikeindex = [int(t / (srate)) for t in d["Results"][trial]["spikeTimes"]]
            spks = AR.MC.time_base[spikeindex]  # get postsynaptic spikes for the trial
            for n, s in enumerate(spks):  # for each postsynaptic spike
                if (
                    s < RCP.min_time or s > RCP.max_time
                ):  # restrict post spikes to those only in a response window
                    continue
                reltime = np.around(RCD.ti, 5) - np.around(s, 5)
                areltime = np.argwhere(
                    (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
                ).squeeze()
                spikedata = SpikeData()  # store spike waveform using a dataclass
                spikedata.trial = trial
                spikedata.waveform = d["Results"][trial]["somaVoltage"][areltime]

                spikedata.time_index = n
                spikedata.prespikes = [[np.nan] for x in range(RCP.ninputs)]
                RCD.nspikes += 1  # number of post spikes evaluated
                n_active_inputs = (
                    0  # number of active inputs associated with this post spike
                )
                spike_pattern = np.zeros(RCP.ninputs)
                solo = np.zeros(RCP.ninputs)
                lack_largest = np.zeros(RCP.ninputs)
                lack_two_largest = np.zeros(RCP.ninputs)

                for isite in range(RCP.ninputs):  # examine each input
                    if not sellist[isite]:
                        continue
                    npre_i, pre_times = self._count_spikes_in_window(
                        d, trial, isite, s, RCP, RCD.pre_w
                    )
                    n_active_inputs += min(1, npre_i)  # only count once
                    if npre_i > 0:

                        spikedata.prespikes[isite] = pre_times[
                            0
                        ]  # save all event times even if more than one
                        RCD.participation[
                            isite
                        ] += 1  # any spikes in the window = participation (but only count as 1)
                        spike_pattern[isite] += 1
                        # print(' spk: ', s, 'isite: ', isite)
                        for jsite in range(
                            isite + 1, RCP.ninputs
                        ):  # now do for joint combinations with other remaining inputs
                            npre_j, pre_times = self._count_spikes_in_window(
                                d, trial, jsite, s, RCP, RCD.pre_w
                            )
                            # print(npre_j)
                            if npre_j > 0:  # accumulate if coincident for this pair
                                RCD.pairwise[isite, jsite] += 1
                # increment the number of times there were npre_spikes input to this post spike
                pre_spike_counts[n_active_inputs] += 1
                if np.sum(spike_pattern) == 1:  # only one input was active
                    which_input = np.where(spike_pattern == 1)[0]
                    pre_solo_spikes[which_input] += 1
                if sum(spike_pattern[0:5]) == 0:
                    cprint(
                        "magenta",
                        f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                    )
                    nfilt_spikes += 1
                    filttable.append(spike_pattern)
                elif sum(spike_pattern[0:4]) == 0:
                    cprint(
                        "cyan",
                        f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                    )
                    nfilt_spikes += 1
                    filttable.append(spike_pattern)
                elif sum(spike_pattern[0:3]) == 0:
                    cprint(
                        "blue",
                        f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                    )
                    nfilt_spikes += 1
                    filttable.append(spike_pattern)
                elif sum(spike_pattern[0:2]) == 0:
                    cprint(
                        "green",
                        f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                    )
                    nfilt_spikes += 1
                    filttable.append(spike_pattern)
                    if spike_pattern[2] == 1:
                        filttable2.append(spike_pattern)
                        nfilt2_spikes += 1
                # elif spike_pattern[0]== 0:
                #     cprint('yellow', f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}")
                # nfilt_spikes += 1
                # filttable.append(spike_pattern)
                else:
                    pass
                    # cprint("w", f"{str(spike_pattern):s}, {which_input[0]:d}")

                self.allspikes.append(spikedata)

        print("\nPairwise matrix: \n", RCD.pairwise)
        print("nprespikes: \n", [int(p) for p in pre_spike_counts])
        print("\nTotal post spikes: ", RCD.nspikes)
        print("\nPre solo drive: ", pre_solo_spikes)
        print("\nFiltered Spikes: ", nfilt_spikes)

        # print('filttable: \n', filttable)
        filttable = np.array(filttable)
        filttable2 = np.array(filttable2)
        print("Counts: ", filttable.sum(axis=0))
        if filttable.shape[0] > 0:
            print("\nFilt Spike Proportions: ", filttable.sum(axis=0) / nfilt_spikes)
            fs1 = np.array(filttable)[:, 2:].sum(axis=0) / nfilt_spikes

        fsa = np.array(RCD.sites[2:]) / np.sum(RCD.sites[2:])
        print("Input Proportions: ", fsa)

        filttable2 = np.array(filttable2)
        # print('filttable 2 shape: ', filttable2.shape)
        if filttable2.shape[0] > 0:
            print(
                "\nFilt Spike Proportions on input #3: ",
                filttable2.sum(axis=0) / nfilt2_spikes,
            )
            fs2 = np.array(filttable2)[:, 2:].sum(axis=0) / nfilt2_spikes

        # if filttable.shape[0] > 0:
        #     f=mpl.figure()
        #     mpl.plot(fs1, fsa)
        #     mpl.show()
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
        npartipating = np.sum(RCD.participation)
        RCD.s_pair = np.sum(RCD.pairwise)
        if RCD.s_pair > 0.0:
            RCD.pairwise /= RCD.s_pair

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
        cprint("r", f"prespikecounts: {np.sum(pre_spike_counts[1:]):f}")
        if np.sum(pre_spike_counts[1:]) == 0:
            return (P, PD, RCP, RCD)

        RCD.ynspike = np.cumsum(pre_spike_counts[1:]) / np.sum(pre_spike_counts[1:])

        if P is None:
            return (P, PD, RCP, RCD)
        else:
            self.plot_revcorr_details(P, PD, RCP, RCD)

    def plot_revcorr_details(self, P, PD, RCP, RCD):
        ax = P.axdict["B"]
        summarySiteTC = self.plot_revcorr2(P, PD, RCP, RCD)

        # ax.set_title(
        #     f"Cell {gbc:s} {str(si.shortSimulationFilename):s}\n[{ri.runTime:s}] dB:{ri.dB:.1f} Prot: {ri.runProtocol:s}" +
        #     f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s}",
        #     fontsize=11,
        # )
        # return summarySiteTC, RCD.sites
        maxp = np.max(RCD.pairwise)
        psh = RCD.pairwise.shape
        pos = np.zeros((psh[0], psh[1], 2))
        for i in range(RCP.ninputs):
            for j in range(RCP.ninputs):
                # print(f"{pairwise[i,j]:.3f}", end='  ')
                pos[i, j, 0] = i + 1
                pos[i, j, 1] = j + 1

        sax = P.axdict
        # f, sax = mpl.subplots(3,1)
        # f.set_size_inches( w=3.5, h=9)
        sax["C"].plot(np.arange(RCP.ninputs) + 1, RCD.sites, "bo")
        # print('pairwise: ', pairwise)
        colormap = "plasma"
        if RCD.s_pair > 0.0:
            # pclip = np.clip(RCD.pairwise, np.min(np.where(RCD.pairwise > 0)), np.max(RCD.pairwise))
            pclip = np.clip(RCD.pairwise, 0, np.max(RCD.pairwise))
            pclip[np.where(pclip == 0)] = np.nan
            pclipped = pclip - np.nanmin(pclip)
            sax["D"].scatter(
                pos[:, :, 0],
                pos[:, :, 1],
                s=160 * RCD.pairwise / maxp,
                c=pclipped,
                cmap=colormap,
            )
            vmax = np.nanmax(pclip) * 100
            vmin = np.nanmin(pclip) * 100
            # print("vmin, vmax: ", vmin, vmax)
        else:
            vmin = 0
            vmax = 1
        vmin = 0
        # sax['B'].plot(np.arange(RCP.ninputs)+1, participation/nspikes, 'gx')
        sax["E"].plot(RCD.sites, RCD.participation / RCD.nspikes, "gx")

        axcbar = PLS.create_inset_axes(
            [0.8, 0.05, 0.05, 0.5], sax["D"], label=str(P.axdict["D"])
        )
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        # ticks = [int(p) for p in np.linspace(vmin, vmax, num=4, endpoint=True)]
        ticks = [int(p) for p in np.linspace(0, vmax, num=5, endpoint=True)]
        cm_sns = mpl.cm.get_cmap(colormap)
        colorbar = matplotlib.colorbar.ColorbarBase(
            axcbar, cmap=cm_sns, ticks=ticks, norm=norm
        )
        ticks = [int(p) for p in np.arange(1, RCP.ninputs + 0.5, 1)]

        # PH.nice_plot(sax['C'], position=-0.2)
        sax["F"].plot(np.arange(len(RCD.ynspike)) + 1, RCD.ynspike, "m^-")

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
        maxv = RCP.ninputs + 0.5
        sax["D"].set_xlim(0, maxv)
        sax["D"].set_ylim(0, maxv)
        sax["D"].xaxis.set_major_formatter(ticker.ScalarFormatter())
        sax["D"].yaxis.set_major_formatter(ticker.ScalarFormatter())
        sax["D"].set_xticks(ticks)
        sax["D"].set_yticks(ticks)

        sax["F"].set_xlabel(
            f"# Inputs in [{RCD.pre_w[0]:.1f} to {RCD.pre_w[1]:.1f}] before spike"
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
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 1},
            pointSize=7,
            colorbar=colorbar,  # use the colorbar tick settings instead
        )
        P.figure_handle.suptitle(
            f"Cell {PD.thiscell:s} {str(RCP.si.shortSimulationFilename):s}\n[{RCP.ri.runTime:s}] dB:{RCP.ri.dB:.1f} Prot: {RCP.ri.runProtocol:s}"
            + f"\nExpt: {RCP.ri.Spirou:s}  DendMode: {RCP.si.dendriteMode:s} DendExpt: {RCP.si.dendriteExpt:s}"
            + f"\nDepr: {RCP.si.ANSynapticDepression:d} ANInput: {RCP.si.SRType:s}",
            fontsize=11,
        )
        mpl.show()

        return (summarySiteTC, RCD.sites)

    def analyze_singles(self, index_row, selected):
        nfiles = len(selected.files)

        P = PH.regular_grid(
            nfiles,
            1,
            order="rowsfirst",
            figsize=(6.0, 8.0),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.01,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.15,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )

        PD = PData()
        sfi = sorted(selected.files)

        df = pd.DataFrame(
            index=np.arange(0, nfiles),
            columns=[
                "cell",
                "syn#",
                "nout",
                "nin",
                "efficacy",
                "ASA",
                "SDRatio",
                "nsites",
            ],
        )
        for i in range(nfiles):
            synno, nout, nin = self.plot_traces(
                P.axarr[i, 0],
                sfi[i],
                PD,
                selected.runProtocol,
                nax=i,
                figure=P.figure_handle,
            )
            SC, syninfo = self.get_synaptic_info(self.parent.cellID)
            r, ctype = SC.make_dict(f"VCN_c{int(self.parent.cellID):02d}")
            area = syninfo[1][synno][0]
            nsites = int(np.around(area * SC.synperum2))
            eff = f"{(float(nout) / nin):.5f}"
            area = f"{float(area):.2f}"
            cell_number = int(self.parent.cellID)
            soma_area = SC.SDSummary.loc[
                SC.SDSummary["Cell Number"] == cell_number
            ].iloc[0]["Somatic Surface Area"]
            dendrite_area = SC.SDSummary.loc[
                SC.SDSummary["Cell Number"] == cell_number
            ].iloc[0]["Dendritic Surface Area"]
            SDRatio = f"{dendrite_area/soma_area:.2f}"
            df.iloc[i] = [
                self.parent.cellID,
                synno,
                nout,
                nin,
                eff,
                area,
                SDRatio,
                nsites,
            ]
        u = df.head(n=nfiles)
        self.textappend(df.to_csv(sep="\t"))
        P.figure_handle.show()

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

    @TraceCalls()
    def setup_PSTH(self):
        sizer = OrderedDict(  # define figure layout
            [
                ("A", {"pos": [0.08, 0.4, 0.73, 0.16]}),
                ("B", {"pos": [0.08, 0.4, 0.50, 0.16]}),
                ("C", {"pos": [0.08, 0.4, 0.22, 0.24]}),
                ("D", {"pos": [0.08, 0.4, 0.06, 0.10]}),
                # rhs
                ("E", {"pos": [0.55, 0.18, 0.73, 0.16]}),  # vs and sac side by side
                ("F", {"pos": [0.78, 0.18, 0.73, 0.16]}),
                ("G", {"pos": [0.55, 0.4, 0.50, 0.16]}),
                ("H", {"pos": [0.55, 0.4, 0.22, 0.24]}),
                ("I", {"pos": [0.55, 0.4, 0.06, 0.10]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.

        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            label=True,
            figsize=(8.0, 6.0),
            # labelposition=(-0.05, 1.02),
        )

        P.axdict["A"].set_ylabel("mV", fontsize=8)

        P.axdict["B"].set_title("Bushy Spike Raster", fontsize=9)
        P.axdict["B"].set_ylabel("Trial")

        P.axdict["C"].set_title("PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["C"].set_ylabel("Spikes/second", fontsize=9)
        P.axdict["D"].set_title("Stimulus", fontsize=9)
        P.axdict["D"].set_ylabel("Amplitude (Pa)", fontsize=8)
        P.axdict["D"].set_xlabel("T (s)", fontsize=9)

        P.axdict["E"].set_title("Phase", fontsize=9)
        P.axdict["F"].set_title("SAC", fontsize=9)

        P.axdict["G"].set_title("ANF Spike Raster", fontsize=9)
        P.axdict["H"].set_title("ANF PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["I"].set_title("FSL/SSL", fontsize=9)
        for axl in [
            "B",
            "C",
            "D",
            "G",
            "H",
        ]:
            P.axdict[axl].sharex(P.axdict["A"])
        return P

    @winprint_continuous
    def print_VS(self, d, freq, dmod, dB, experiment):
        print("Getting data")
        SC, syninfo = self.get_synaptic_info(self.parent.cellID)
        amax = 0
        n_inputs = len(syninfo[1])
        sites = np.zeros(n_inputs)
        for isite in range(n_inputs):  # precompute areas
            area = syninfo[1][isite][0]
            if area > amax:
                amax = area
            sites[isite] = int(np.around(area * SC.synperum2))

        self.VS_colnames = f"Cell,Configuration,carrierfreq,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,AN_phase,AN_phasesd,SAC_AN,SAC_Bu,SAC_AN_HW,SAC_Bu_HW,maxArea,ninputs"
        line = f"{int(self.parent.cellID):d},{experiment:s},"
        line += f"{d.carrier_frequency:.1f},{freq:06.1f},{dmod:.1f},{dB:.1f},"
        line += f"{d.vs:.4f},"
        line += f"{d.n_spikes:d},"
        line += f"{d.circ_phaseMean:.4f},"
        line += f"{d.circ_phaseSD:.4f},"
        line += f"{d.Rayleigh:.4f},"
        line += f"{d.pRayleigh:.4e},"
        line += f"{d.an_vs:.4f},"
        line += f"{d.an_circ_phaseMean:.4f},"
        line += f"{d.an_circ_phaseSD:.4f},"
        line += f"{d.sac_an_CI:.4f},"
        line += f"{d.sac_bu_CI:.4f},"
        line += f"{d.sac_an_hw:.6f},"
        line += f"{d.sac_bu_hw:.6f},"
        line += f"{amax:.4f},"
        line += f"{d.n_inputs:d}"
        self.VS_line = line
        if self.firstline:
            self.textappend(self.VS_colnames)
        self.textappend(line)
        # put on clipboard
        pyperclip.copy(line)
        pyperclip.paste()

        return d

    @winprint_continuous
    @TraceCalls()
    def plot_SAC(self, selected=None):
        print("SAC analysis and plotting starting")

        if self.parent.selected_index_rows is None:
            return

        selrows = self.parent.table.selectionModel().selectedRows()

        P = PH.regular_grid(
            3,
            1,  # top: SAC, middle: raster, bottom: waveform
            order="rowsfirst",
            figsize=(6.0, 8.0),
            showgrid=False,
            verticalspacing=0.1,
            horizontalspacing=0.1,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.15,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )
        PD = PData()
        protocol = selected.runProtocol
        changetimestamp = get_changetimestamp()
        sfi = sorted(selected.files)
        plot_win = [0.1, 1.0]
        plot_dur = np.fabs(np.diff(plot_win))
        time_scale = 1.0
        for j, index_row in enumerate(selrows):
            selected = self.parent.table_manager.get_table_data(index_row)
            if selected is None:
                return
            nfiles = len(selected.files)
            fn = selected.files[0]
            gbc = f"VCN_c{int(self.parent.cellID):02d}"
    #         self.do_one_SAC(fn, gbc, PD, protocol)
    #
    # def do_one_SAC(self, fn, gbc, PD, protocol):
        
            print(f"Getting data for gbc: {gbc:s}")
            SC, syninfo = self.get_synaptic_info(gbc)
            mtime = Path(fn).stat().st_mtime
            timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
                "%Y-%m-%d-%H:%M"
            )

            changetimestamp = get_changetimestamp()
            
            res = self.get_data(fn, PD, changetimestamp, protocol)
            if res is None:
                return None
            (d, AR, SP, RMA, RCP, RCD) = res
            PD.thiscell = gbc
            si = d["Params"]
            ri = d["runInfo"]
            RCP.si = si
            RCP.ri = ri
            (
                totaldur,
                soundtype,
                pip_start,
                pip_duration,
                F0,
                dB,
                fmod,
                dmod,
            ) = self.get_stim_info(si, ri)
            sim_dir = Path(fn).parts[-2]
            tstring = f"{gbc:s}_{str(sim_dir):s}_{soundtype:s}"
            mpl.get_current_fig_manager().canvas.set_window_title(tstring)
            fstring = f"{gbc:s}_{soundtype:s}_{ri.dB} dBSPL"
            P.figure_handle.suptitle(fstring, fontsize=13, fontweight="bold")
            
            all_bu_st = []
            all_bu_st_trials = []
            ntr = len(AR.MC.traces)  # number of trials
            for i in range(ntr):  # for all trials in the measure.
                time_base = AR.MC.time_base / 1000.0  # convert to seconds
                trd = d["Results"][i]
                dt = si.dtIC / 1000.0  # convert from msec to seconds
                idx = (int(plot_win[0] / dt), int(plot_win[1] / dt))
                if i == 0:
                    waveform = trd["stimWaveform"]
                stb = trd["stimTimebase"]  # convert to seconds
                stim_dt = stb[1] - stb[0]
                if i == 0:
                    n_inputs = len(trd["inputSpikeTimes"])
                try:
                    if len(trd["spikeTimes"]) > 0 and np.nanmax(trd["spikeTimes"]) > 2.0:  # probably in msec
                        time_scale = 1e-3  # so scale to seconds
                except:
                    raise ValueError(trd["spikeTimes"])
                sptimes = np.array(trd["spikeTimes"]) * time_scale  # convert to seconds
                if not isinstance(trd["spikeTimes"], list) and not isinstance(
                    trd["spikeTimes"], np.ndarray
                ):
                    cprint("r", "spiketimes is not list")
                    cprint("r", f"    {type(trd['spikeTimes'])=}")
                    return
                all_bu_st.extend(sptimes)
                all_bu_st_trials.append(sptimes)
                ispt = [  # plot spike times in the SAC analysis window
                    i
                    for i in range(len(sptimes))
                    if sptimes[i] >= pip_start and sptimes[i] < pip_duration - pip_start
                ]
                P.axdict["B"].plot(
                    np.array(sptimes[ispt]),
                    i * np.ones(len(ispt)),
                    "|",
                    markersize=1.5,
                    color="b",
                )
                w_tb = np.linspace(0.0, stim_dt * len(waveform), num=len(waveform))

                i_wpt = np.where((w_tb > pip_start) & (w_tb <= pip_duration))[0]
                P.axdict["C"].plot(w_tb[i_wpt], waveform[i_wpt], linewidth=0.33)
            if ri.soundtype.endswith("Clicks"):
                print("Clickpars")
                pars = {
                    "twin": 0.002,
                    "binw": 3 * dt,
                    "delay": ri.clickStart + 0.2 * ri.clickTrainDuration,
                    "dur": 0.8 * ri.clickTrainDuration,
                    "displayDuration": 0.002,
                    "nrep": len(all_bu_st_trials),
                    "baseper": 1e-3 * 1.3333333333333333,
                    "stimdur": pip_duration * 0.8,
                    "maxn": 100000000,
                }

            else:
                print("Sam Pars")
                pars = {
                    "twin": 0.020,
                    "binw": 3 * dt,
                    "delay": pip_start + 0.2 * pip_duration,
                    "dur": 0.8 * pip_duration,
                    "displayDuration": 0.050,
                    "nrep": len(all_bu_st_trials),
                    "baseper": 1e-3 * 1.3333333333333333,
                    "stimdur": pip_duration * 0.8,
                    "maxn": 100000000,
                }
            print(ri.soundtype)
            sac = SAC.SAC()
            yh, bins = sac.SAC_with_histo(
                all_bu_st_trials, pars=pars, engine="cython", dither=dt / 2.0
            )
            sac_bu_CI, peaks, HW, FW = sac.SAC_measures(yh, bins)
            sac_bu_hw = HW[0][0]*pars["binw"]
            print("BU SAC Report: \n  ")
            print(f"    HW:    {sac_bu_hw:.6f}  CI: {sac_bu_CI:.2f}  Left IPS: {HW[2][0]:.2f}  Right IPS: {HW[3][0]:.2f}")
            # Fs: float = 100e3  # cochlea/zilany model rate
            # F0: float = 16000.0  # stimulus frequency
            # dB: float = 30.0  # in SPL
            # RF: float = 2.5e-3  # rise-fall time
            # fmod: float = 20  # hz, modulation if SAM
            # dmod: float = 0.0  # percent if SAM
            if ri.soundtype.endswith("Clicks"):
                sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL  HW={1e3*sac_bu_hw:.3f} ms  CI={sac_bu_CI:6.2f}"
            else:
                sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL Fmod={ri.fmod:5.1}fHz Dmod={ri.dmod:5.1f}\%"
            
            return 
            P.axdict["A"].plot(
                bins,
                yh,
                #   'k-',
                label=sac_label,
            )
            self.textappend(sac_label)
        
        P.axdict["C"].set_xlim(pip_start, pip_duration - pip_start)
        P.axdict["B"].set_xlim(pip_start, pip_duration - pip_start)
        P.axdict["A"].set_xlim(-pars['twin'], pars['twin'])
        P.axdict["A"].set_ylim(0, sac_bu_CI*1.05)
        P.axdict["A"].legend(prop={"family": "monospace"}, fontsize=7)

        P.axdict["B"].get_shared_x_axes().join(P.axdict["B"], P.axdict["C"])
        
        mpl.show()

    @TraceCalls()
    def psth_vs(self):
        """
        Generate tables of vs measures for all cells
        across the frequencies listed
        """
        self.textclear()  # just at start
        PD = PData()
        if self.parent.selected_index_rows is None:
            return
        # P = self.PLT.setup_PSTH()
        P = None
        PD = PData()
        selrows = self.parent.table.selectionModel().selectedRows()
        for i, index_row in enumerate(selrows):
            selected = self.parent.table_manager.get_table_data(
                index_row
            )  # table_data[index_row]
            if selected is None:
                return
            sfi = Path(selected.simulation_path, selected.files[0])
            if i == 0:
                self.firstline = True
            else:
                self.firstline = False
            self.plot_AN_response(P, sfi, PD, selected.runProtocol)
        self.firstline = True

    def _ntrials(self, spike_times: Union[list, np.ndarray]):
        if spike_times.ndim == 1:  # get trial count by looking at data
            num_trials = 1
        else:
            num_trials = spike_times.shape[0]
        return num_trials

    @TraceCalls()
    def plot_psth(
        self,
        spike_times: Union[list, np.ndarray],
        run_info: object,
        stim_onset: float = 0.0,
        zero_time: float = 0.0,
        max_time: float = 1.0,
        bin_width: float = 1e-3,
        ax: Union[object, None] = None,
        scale: float = 1.0,
    ):
        """
        Correctly plot PSTH with spike rate in spikes/second
        All values are in seconds (times, binwidths)
        """
        num_trials = run_info.nReps
        bins = np.arange(0.0, max_time - zero_time, bin_width)
        spf = []
        for x in spike_times:
            if isinstance(x, list):
                spf.extend(x)
            else:
                spf.extend([x])
        spike_times_flat = np.array(spf, dtype=object).ravel() - zero_time
        # h, b = np.histogram(spike_times_flat, bins=bins)
        bins = np.arange(0.0, max_time - zero_time, bin_width)
        if (not np.isnan(np.sum(spike_times_flat))) and (len(spike_times_flat) > 0) and (ax is not None):
            ax.hist(
                x=spike_times_flat,
                bins=bins,
                density=False,
                # weights=scale * np.ones_like(spike_times_flat),
                histtype="stepfilled",
                facecolor="k",
                edgecolor="k",
            )
        else:
            if ax is not None:
                ax.text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=ax.transAxes,
                    horizontalalignment="center",
                )
        return  # h, b

    @TraceCalls()
    def plot_fsl_ssl(
        self,
        spike_times: Union[list, np.ndarray],
        run_info: object,
        max_time: float,
        min_time: float = 0.0,
        fsl_win: Union[None, tuple] = None,
        bin_width: float = 1e-3,
        ax: Union[object, None] = None,
        zero_time: float = 0.0,
        offset: float = 0.0,
        cellID: Union[int, None, str] = None,
        show_values: bool = True,
        trange: Union[tuple, list] = [0, 25.0],
    ):
        # spike_times = np.array(spike_times, dtype=object)

        num_trials = run_info.nReps  # self._ntrials(spike_times)
        # for s in range(len(spike_times)):
        #     spike_times[s] = [t*1e-3 for t in spike_times[s]]
        # print('FSL SSL: ', zero_time, run_info.pip_start, fsl_win)
        fsl, ssl = SPKANA.CNfspike(
            spike_times,
            run_info.pip_start - zero_time,
            nReps=run_info.nReps,
            fsl_win=fsl_win,
        )

        fsl *= 1e3  # better to show FSL/SSL in milliseconds
        ssl *= 1e3
        if cellID is None:
            index_row = self.parent.selected_index_rows[0]
            selected = self.parent.table_manager.get_table_data(index_row)
            print(
                f"Cell: '{int(self.parent.cellID):d}','{selected.synapseExperiment:s}',",
                end="",
            )
        # else:
        #     print(f"Cell: '{int(cellID):d}','{selected.synapseExperiment:s}',", end="")

        if len(fsl) > 0:
            print(f"{np.nanmean(fsl):.3f},{np.nanstd(fsl):.3f},", end="")
        else:
            print("No FSL spike data. ", end="")
        if len(ssl) > 0:
            print(f"{np.nanmean(ssl):.3f},{np.nanstd(ssl):.3f}")
        else:
            print("No SSL spike data")
        if ax is not None:
            latency_hbins = np.arange(min_time, max_time, bin_width)  # use msec here
            if (len(fsl)) > 0.0:
                ax.hist(
                    fsl,
                    bins=latency_hbins,
                    facecolor="b",
                    edgecolor="b",
                    alpha=0.6,
                    bottom=offset,
                )
            if (len(ssl)) > 0.0:
                ax.hist(
                    ssl,
                    bins=latency_hbins,
                    facecolor="r",
                    edgecolor="r",
                    alpha=0.6,
                    bottom=offset,
                )
            ax.set_xlim(min_time, max_time)
            if show_values:
                if len(fsl) > 0:
                    fsl_text = f"FSL: {np.nanmean(fsl):.3f} (SD {np.nanstd(fsl):.3f})"
                else:
                    fsl_text = "No spikes for FSL measure. "
                if len(ssl) > 0:
                    fsl_text += (
                        f"\nSSL: {np.nanmean(ssl):.3f} (SD {np.nanstd(ssl):.3f})"
                    )
                else:
                    fsl_text = "No spikes for SSL measure. "
                ax.text(
                    0.45,
                    0.95,
                    fsl_text,
                    fontsize=7,
                    color="k",
                    # fontfamily="monospace",
                    transform=ax.transAxes,
                    horizontalalignment="left",
                    verticalalignment="top",
                )

        return (fsl, ssl)

    def get_stim_info(self, si, ri):
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
        return totaldur, soundtype, pip_start, pip_duration, F0, dB, fmod, dmod

    @winprint_continuous
    @TraceCalls()
    def plot_AN_response(
        self,
        P: Union[object, None] = None,
        fn: Union[str, Path, None] = None,
        PD: object = None,
        protocol: str = "",
        sac_flag: bool = True,  # set true to compute SAC as well as standard VS
        sac_engine: str = "cython",
    ):
        gbc = f"VCN_c{int(self.parent.cellID):02d}"
        if P is None:
            plotflag = False
        else:
            plotflag = True
        changetimestamp = get_changetimestamp()
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        d, AR, SP, RMA, RCP, RCD = self.get_data(fn, PD, changetimestamp, protocol)

        ntr = len(AR.MC.traces)  # number of trials
        v0 = -160.0
        trstep = 25.0 / ntr
        inpstep = 5.0 / ntr
        sz = 50.0 / ntr
        si = d["Params"]
        ri = d["runInfo"]
        (
            totaldur,
            soundtype,
            pip_start,
            pip_duration,
            F0,
            dB,
            fmod,
            dmod,
        ) = self.get_stim_info(si, ri)
        if plotflag:
            sim_dir = fn.parts[-2]
            tstring = f"{gbc:s}_{str(sim_dir):s}_{soundtype:s}.pdf"
            mpl.get_current_fig_manager().canvas.set_window_title(tstring)

        bu_voltage_panel = "A"
        bu_raster_panel = "B"
        bu_psth_panel = "C"
        stimulus_panel = "D"
        vs_panel = "E"
        sac_panel = "F"
        an_raster_panel = "G"
        an_psth_panel = "H"
        bu_fsl_panel = "I"

        # We must get the units correct here.
        # AR time_base is in milliseconds, so convert to seconds
        ntr = len(AR.MC.traces)
        all_an_st = []
        all_bu_st = []
        all_bu_st_trials = []
        waveform = []
        plot_win = (
            pip_start - 0.05,
            pip_start + pip_duration + 0.05,
        )  # manually change to shorten...
        plot_dur = np.fabs(np.diff(plot_win))

        for i in range(ntr):  # for all trials in the measure.
            vtrial = AR.MC.traces[i] * 1e3  # convert to mV
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            dt = si.dtIC / 1000.0  # convert from msec to seconds
            trd = d["Results"][i]
            idx = (int(plot_win[0] / dt), int(plot_win[1] / dt))
            waveform.append(trd["stimWaveform"].tolist())
            stb = trd["stimTimebase"]  # convert to seconds
            if i == 0:
                n_inputs = len(trd["inputSpikeTimes"])
            if len(trd["spikeTimes"]) == 0:
                trd["spikeTimes"] = np.array([np.nan])
                sf = 1.0
            elif np.nanmax(trd["spikeTimes"]) > 2.0:  # probably in miec...
                sf = 1e-3  # so scale to seconds
            else:
                sf = 1.0
            sptimes = np.array(trd["spikeTimes"]) * sf  # convert to seconds

            stim_dt = stb[1] - stb[0]
            idxs = (int(plot_win[0] / stim_dt), int(plot_win[1] / stim_dt))
            if not isinstance(trd["spikeTimes"], list) and not isinstance(
                trd["spikeTimes"], np.ndarray
            ):
                cprint("r", "spiketimes is not list")
                cprint("r", f"    {type(trd['spikeTimes'])=}")
                return
            all_bu_st.extend(sptimes)
            all_bu_st_trials.append(sptimes)

            if i == 0 and plotflag:
                P.axdict[bu_voltage_panel].plot(
                    time_base[idx[0] : idx[1]] - plot_win[0],
                    vtrial[idx[0] : idx[1]],
                    "k-",
                    linewidth=0.5,
                )
                spikeindex = [
                    int(t / dt) for t in sptimes if not np.isnan(t)
                ]  # dt is in sec, but spiketimes is in msec
                ispt = [
                    t
                    for t in range(len(sptimes))
                    if sptimes[t] >= plot_win[0] and sptimes[t] < plot_win[1]
                ]
                P.axdict[bu_voltage_panel].plot(
                    sptimes[ispt] - plot_win[0],
                    vtrial[spikeindex][ispt],
                    "ro",
                    markersize=1.5,
                )
                P.axdict[bu_voltage_panel].set_xlim(0.0, plot_dur)
            if waveform is not None and plotflag:  # could also do if i == 0
                dw = 0.5 * np.max(waveform[i])

                P.axdict[stimulus_panel].plot(
                    stb[idxs[0] : idxs[1]] - plot_win[0],
                    np.array(waveform[i][idxs[0] : idxs[1]]) + i * dw,
                    "-",
                    linewidth=0.5,
                )  # stimulus underneath

            if plotflag:  # plot the raster of spikes
                ispt = [
                    i
                    for i in range(len(sptimes))
                    if sptimes[i] >= plot_win[0] and sptimes[i] < plot_win[1]
                ]
                P.axdict[bu_raster_panel].plot(
                    np.array(sptimes[ispt]) - plot_win[0],
                    i * np.ones(len(ispt)),
                    "|",
                    markersize=1.5,
                    color="b",
                )
                P.axdict[bu_raster_panel].set_xlim(0, plot_dur)
                
                for k in range(n_inputs):  # raster of input spikes
                    if np.max(trd["inputSpikeTimes"][k]) > 2.0:  # probably in miec...
                        tk = np.array(trd["inputSpikeTimes"][k]) * 1e-3
                    else:
                        tk = np.array(trd["inputSpikeTimes"][k])

                    y = (i + 0.1 + k * 0.05) * np.ones(len(tk))
                    in_spt = [
                        i
                        for i in range(tk.shape[0])
                        if (tk[i] >= plot_win[0]) and (tk[i] < plot_win[1])
                    ]
                    y = (i + 0.1 + k * 0.05) * np.ones(len(in_spt))
                    if i % 10 == 0 and plotflag:
                        P.axdict[an_raster_panel].plot(
                            tk[in_spt] - plot_win[0],
                            y,
                            "|",
                            markersize=2.5,
                            color="k",
                            linewidth=0.5,
                        )
        an_st_by_input = [
            [] for x in range(n_inputs)
        ]  # accumlate by input across trials
        all_an_st = []  # collapsed completel = all spikes in one array
        an_st_grand = [
            [] for x in range(ntr)
        ]  # accumulate all input spike times for each trial

        for i in range(ntr):  # for all trials in the measure.
            trd = d["Results"][i]  # trial i data
            n_inputs = len(trd["inputSpikeTimes"])  # get AN spike time array (by input)
            if i == 0:
                an_st_by_input = [[] for x in range(n_inputs)]
            for k in range(n_inputs):  # raster of input spikes
                if np.max(trd["inputSpikeTimes"][k]) > 2.0:
                    tk = np.array(trd["inputSpikeTimes"][k]) * 1e-3
                else:
                    tk = np.array(trd["inputSpikeTimes"][k])
                all_an_st.extend(tk)  # strictly for the psth and VS calculation
                an_st_by_input[k].append(
                    tk
                )  # index by input, then trials for that input
                an_st_grand[i].extend(tk)

        cprint("r", f"{soundtype=}")

        if (
            soundtype
            in [
                "tonepip",
                "SAM",
                "noise",
                "stationaryNoise",
                "regularClicks",
                "poissonClicks",
            ]
            and plotflag
        ):  # use panel F for FSL/SSL distributions
            # the histograms of the data
            psth_binw = 0.2e-3
            self.plot_psth(
                all_bu_st,
                run_info=ri,
                zero_time=plot_win[0],
                max_time=plot_win[1],
                bin_width=psth_binw,
                ax=P.axdict[bu_psth_panel],
            )
            self.plot_psth(
                all_an_st,
                run_info=ri,
                zero_time=plot_win[0],
                max_time=plot_win[1],
                bin_width=psth_binw,
                ax=P.axdict[an_psth_panel],
            )
            self.plot_fsl_ssl(
                all_bu_st_trials,
                run_info=ri,
                fsl_win=(2.5e-3, 5e-3),
                max_time=25.0,
                bin_width=0.25,
                ax=P.axdict[bu_fsl_panel],
            )

        if soundtype in [
            "SAM",
            "stationaryNoise",
            "regularClicks",
            "poissonClicks",
        ]:  # also calculate vs and plot histogram
            if len(all_bu_st) == 0 and plotflag:
                P.axdict[vs_panel].text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=P.axdict[vs_panel].transAxes,
                    horizontalalignment="center",
                )
            else:
                # fold spike times on period, for the last part of the stimulus
                if soundtype.endswith("Clicks"):
                    phasewin = [
                        ri.clickStart + 0.25 * ri.clickTrainDuration,
                        ri.clickStart + ri.clickTrainDuration,
                    ]
                else:
                    phasewin = [
                        pip_start + 0.25 * pip_duration,
                        pip_start + pip_duration,
                    ]

                x = np.array(all_bu_st)
                v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
                bu_spikesinwin = x[v]
                x = np.array(all_an_st)
                v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
                an_spikesinwin = x[v]
                if soundtype == "tone":
                    print("Tone: F0=%.3f at %3.1f dbSPL, cell CF=%.3f" % (F0, dB, F0))
                if (
                    soundtype
                    in ["SAM", "stationaryNoise", "regularClicks", "poissonClicks"]
                    and plotflag
                ):
                    sim_dir = fn.parts[-2]
                    tstring = f"{gbc:s} {str(sim_dir):s}\nSAM Tone: F0={F0:.3f} at {dB:3.1f} dbSPL, fMod={fmod:3.1f}  dMod={dmod:5.2f}, cell CF={F0:.3f}"
                    # print('si: ', si)
                    tstring += f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s} DendExpt: {si.dendriteExpt:s}"
                    tstring += (
                        f" Depr: {si.ANSynapticDepression:d} ANInput: {si.SRType:s}"
                    )

                    P.figure_handle.suptitle(tstring, fontsize=8)
                vs = self.VS.vector_strength(
                    bu_spikesinwin, fmod
                )  # vs expects spikes in msec
                vs.carrier_frequency = F0
                vs.n_inputs = n_inputs
                if sac_flag:
                    S = SAC.SAC()
                    spars = SAC.SACPars()
                    spars.binw = 3 * si.dtIC * 1e-3
                    spars.maxn = 100000000
                    if soundtype.endswith("Clicks"):
                        spars.delay = ri.clickStart + 0.2 * ri.clickTrainDuration
                        spars.dur = 0.8 * ri.clickTrainDuration
                        spars.displayDuration = 2e-3  # 2.0/ri.clickRate
                        spars.twin = 10.0 / ri.clickRate
                    else:
                        spars.delay = ri.pip_start + 0.2 * ri.pip_duration
                        spars.dur = 0.8 * ri.pip_duration
                        spars.displayDuration = 2.0 / ri.fmod
                        spars.twin = 10.0 / ri.fmod
                    bu_sac, bu_sacbins = S.SAC_with_histo(
                        all_bu_st_trials,
                        spars,
                        engine=sac_engine,
                        dither=1e-3 * si.dtIC / 2.0,
                    )
                    vs.sac_bu_CI, peaks, HW, FW = S.SAC_measures(bu_sac, bu_sacbins)
                    vs.sac_bu_hw = HW[0][0]*spars.binw
                    print("BU SAC Report: \n  ")
                    print(f"    HW:    {vs.sac_bu_hw:.6f}  CI: {vs.sac_bu_CI:.2f}  Left IPS: {HW[2][0]:.2f}  Right IPS: {HW[3][0]:.2f}")
                else:
                    bu_sac = None
                    bu_sacbins = None
                    vs.sac_bu_CI = 0.0
                    vs.sac_bu_hw = 0.0
                print(
                    "CN Vector Strength at %.1f: %7.3f, dispersion=%.2f (us) Rayleigh: %7.3f  p = %.3e  n_spikes = %d"
                    % (
                        fmod,
                        vs.vs,
                        vs.circ_timeSD * 1e6,
                        vs.Rayleigh,
                        vs.pRayleigh,
                        vs.n_spikes,
                    )
                )
                vs_an = self.VS.vector_strength(an_spikesinwin, fmod,)
                vs_an.n_inputs = n_inputs
                if sac_flag:
                    SA = SAC.SAC()
                    spars = SAC.SACPars()
                    spars.binw = 3 * si.dtIC * 1e-3
                    if soundtype.endswith("Clicks"):
                        spars.delay = ri.clickStart + 0.2 * ri.clickTrainDuration
                        spars.dur = 0.8 * ri.clickTrainDuration
                        spars.displayDuration = 2e-3  # 2.0/ri.clickRate
                        spars.twin = 10.0 / ri.clickRate
                    else:
                        spars.delay = ri.pip_start + 0.2 * ri.pip_duration
                        spars.dur = 0.8 * ri.pip_duration
                    spars.maxn = 100000000
                    if soundtype in ["SAM"]:
                        spars.twin = 10.0 / ri.fmod
                        spars.displayDuration = 2.0 / ri.fmod

                    an_sac = [None] * n_inputs
                    an_hw = [None] * n_inputs
                    an_fw = [None] * n_inputs
                    for i_an in range(n_inputs):
                        an_sac[i_an], an_sacbins = S.SAC_with_histo(
                            an_st_by_input[i_an],
                            spars,
                            engine=sac_engine,
                            dither=1e-3 * si.dtIC / 2.0,
                        )
                        an_CI_windowed, peaks, an_hw[i_an], an_fw[i_an] = S.SAC_measures(an_sac[i_an], an_sacbins)
                    an_hw = np.array(an_hw).T[0]  # reorganize the data
                    vs.sac_an_hw = np.mean(an_hw[0,:])*spars.binw
                    CI = np.mean(an_hw[1,:])
                    vs.sac_an_CI = CI
                    lips = np.mean(an_hw[2,:])
                    rips = np.mean(an_hw[3,:])
                    print("AN SAC Report: ")
                    print(f"    HW:    {vs.sac_an_hw:.6f}  CI: {vs.sac_an_CI:.2f}  Left IPS: {lips:.2f}  Right IPS: {rips:.2f}")
                else:
                    an_sac = None
                    an_sacbins = None
                    vs.sac_an_CI = 0.0
                    vs.sac_an_hw = 0.0
                
                print(" Sound type: ", soundtype)
                print(
                    "AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
                    % (
                        fmod,
                        vs_an.vs,
                        vs_an.circ_timeSD * 1e6,
                        vs_an.Rayleigh,
                        vs_an.pRayleigh,
                        vs_an.n_spikes,
                    )
                )
                if sac_flag:
                    print(vs.sac_bu_CI)
                    print(np.min(bu_sac))
                    print(vs.sac_bu_hw)
                    print(
                        f"BU SAC: CI = {vs.sac_bu_CI:7.3f}, min = {np.min(bu_sac):7.3f}, hw = {1e3*vs.sac_bu_hw:7.3f} ms"
                    )
                    print(
                        f"AN SAC: CI = {vs.sac_an_CI:7.3f}, min = {np.min(an_sac):7.3f}, hw = {1e3*vs.sac_an_hw:7.3f} ms"
                    )
                vs.an_vs = vs_an.vs  # copy
                vs.an_circ_phaseMean = vs_an.circ_phaseMean
                vs.an_circ_phaseSD = vs_an.circ_phaseSD
                vs.an_circ_timeSD = vs_an.circ_timeSD

                self.print_VS(vs, fmod, dmod, dB, ri.Spirou)
                if plotflag:
                    # set bin width based on sampling rate.
                    per = 1.0 / fmod  # period, in seconds
                    est_binw1 = per / 30.0  # try this
                    dt = 1e-3 * si.dtIC
                    nints = np.floor(est_binw1 / dt)
                    if nints == 0:
                        nints = 1
                    if nints < est_binw1 / dt:
                        est_binw = nints * dt
                    else:
                        est_binw = est_binw1
                    # print('est_binw: ', est_binw, est_binw1, dt, nints, per)
                    self.plot_psth(
                        vs.circ_phase,
                        run_info=ri,
                        max_time=2 * np.pi,
                        bin_width=est_binw,
                        ax=P.axdict[vs_panel],
                    )

                    # P.axdict["E"].hist(
                    #     vs["ph"],
                    #     bins=2 * np.pi * np.arange(30) / 30.0,
                    #     facecolor="k",
                    #     edgecolor="k",
                    # )
                    P.axdict[vs_panel].set_xlim((0.0, 2 * np.pi))
                    P.axdict[vs_panel].set_title(
                        f"VS: AN = {vs.an_vs:.3f} BU = {vs.vs:.3f}",
                        fontsize=8,
                        horizontalalignment="center",
                    )
                    if sac_flag:
                        for i_an in range(n_inputs):
                            if i_an == 0:
                                P.axdict[sac_panel].plot(
                                an_sacbins, an_sac[i_an], "r-", label="AN"
                            )
                            else:
                                P.axdict[sac_panel].plot(
                                an_sacbins, an_sac[i_an], "r-",
                            )
                                
                        P.axdict[sac_panel].plot(
                            bu_sacbins, bu_sac, "b-", label="BU"
                        )
                        if soundtype in ["SAM"]:
                            P.axdict[sac_panel].set_xlim(
                                (-2.0 / ri.fmod, 2.0 / ri.fmod)
                            )
                        elif soundtype.endswith("Clicks"):
                            P.axdict[sac_panel].set_xlim(
                                (-spars.displayDuration, spars.displayDuration)
                            )
                        else:
                            print(soundtype)
                            P.axdict[sac_panel].set_xlim(
                                (-spars.displayDuration, spars.displayDuratoin)
                            )
                        P.axdict[sac_panel].set_ylim(
                            (0.0, np.max([vs.sac_an_CI, vs.sac_bu_CI]) * 1.05)
                        )
                        P.axdict[sac_panel].set_title(
                            f"SAC\nCI: AN={vs.sac_an_CI:.3f}  BU={vs.sac_bu_CI:.3f}\nHW: AN={1e3*vs.sac_an_hw:.3f}  BU={1e3*vs.sac_bu_hw:.3f} ms",
                            fontsize=8,
                            horizontalalignment="center",
                        )
                        P.axdict[sac_panel].legend(fontsize=7)

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

    parser.add_argument(
        "-v",
        "--vectorstrength",
        dest="vstest",
        action="store_true",
        help=("test the vs calculation"),
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
            basefn = f"{config['cellDataDirectory']:s}/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
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
            basefn = f"{config['cellDataDirectory']:s}/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
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
    if args.vstest:
        test_vs()
        exit()

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


def plot_ticks(ax, data, color):
    ax.plot(np.tile(data, (2, 1)), [np.zeros_like(data), np.ones_like(data)], color)


def compute_vs(freq=100.0, nsp=1000.0, sd=0.0, rg=None):
    if sd <= 1:
        sdn = (
            (2.0 / freq) + sd * rg.standard_normal(size=nsp) + np.pi
        )  # just express in temrs of time
    else:
        sdn = (1.0 / freq) * rg.uniform(size=nsp)  # uniform across the interval
    spikes = (
        np.cumsum(np.ones(nsp) * (1.0 / freq)) + sdn
    )  # locked spike train with jitter
    phsp = 2 * np.pi * freq * np.fmod(spikes, 1.0 / freq)
    return (spikes, phsp)


def test_vs():
    from numpy.random import default_rng

    rg = default_rng(12345)
    psi = PlotSims(parent=None)
    freq = 100.0
    nsp = 1000
    sd = 0.002  # in seconds, unles > 1 then is uniform
    x = np.array(
        [
            0.0,
            0.00005,
            0.0001,
            0.0002,
            0.0003,
            0.0005,
            0.00075,
            0.001,
            0.002,
            0.003,
            0.004,
            0.0050,
            0.0075,
            2,
        ]
    )
    y = np.zeros_like(x)
    ph = np.zeros_like(x)
    vs = np.zeros_like(x)
    fig, ax = mpl.subplots(3, 1)
    ax = ax.ravel()
    for i, sd in enumerate(x):
        spikes, phsp = compute_vs(freq, nsp, sd, rg)
        vsd = psi.VS.vector_strength(spikes, freq=freq)
        y[i] = vsd.circ_timeSD
        vs[i] = vsd.vs
        ph[i] = vsd.circ_phaseSD
    x[x == 2] = 0.020
    ax[0].plot(x, y, "o")
    ax[0].plot([0, np.max(x)], [0, np.max(x)], "--k", alpha=0.5)
    ax[1].plot(x, vs, "x-")
    ax[2].plot(x, ph, "s-")
    mpl.show()

    fig, ax = mpl.subplots(2, 1)
    ax = ax.ravel()
    plot_ticks(ax[0], phsp, "r-")
    plot_ticks(ax[1], spikes, "b-")
    mpl.show()
    print("mean isi: ", np.mean(np.diff(spikes)), "stim period: ", 1.0 / freq)
    vsd = psi.VS.vector_strength(spikes, freq=freq)
    print("VS: ", vsd.vs, "SD sec: ", vsd.circ_timeSD)
    print("circ_phase min max: ", np.min(vsd.circ_phase), np.max(vsd.circ_phase))
    bins = np.linspace(0, 2 * np.pi, 36)
    mpl.hist(vsd.circ_phase, bins)
    mpl.show()


if __name__ == "__main__":
    # main()
    test_vs()
