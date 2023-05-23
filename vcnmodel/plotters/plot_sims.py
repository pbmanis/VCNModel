"""
Functions to compute some results and plot the simulation results from
model_run2. Reads the new format filenames.

This is written as a class with a set of methods that can be called from
elsewhere. Typically, these will be called from within DataTablesVCN.py, which
provides graphical table-based access to the simulation results. Included at the
end is a parser so that the plots can be called from the command line as well,
but this is discouraged.

Some of the parameters must be instantiated by creating an instance of PData
that is passed into the routnines.

Wrapper for various analysis functions, handles multiple cells.

Command line usage (discouraged; use DataTables instead):

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


This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)


Copyright 2017-2023 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""

import re
import string
import time
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import matplotlib.colorbar  # type: ignore
import matplotlib.colors  # type: ignore
import mplcursors
import numpy as np  # type: ignore
import pandas as pd
import pint
import pyperclip
import pyqtgraph as pg  # type: ignore

# import quantities as pq
import seaborn as sns
import vcnmodel.cell_config as cell_config
import vcnmodel.group_defs as GRPDEF
import vcnmodel.util.readmodel as readmodel
from vcnmodel.util.get_data_paths import get_data_paths, update_disk
from vcnmodel.util import trace_calls
from vcnmodel.util import fixpicklemodule as FPM
from ephys.ephysanalysis import RmTauAnalysis, SpikeAnalysis
from lmfit import Model  # type: ignore
from matplotlib import pyplot as mpl  # type: ignore
from matplotlib import rc  # type: ignore
from matplotlib import ticker
from pylibrary.plotting import plothelpers as PH
from pylibrary.plotting import styler as PLS
from pylibrary.tools import cprint as CP
from . import plot_functions as PF
from pyqtgraph.Qt import QtGui
from rich.console import Console
from vcnmodel.analyzers import reverse_correlation as REVCORR
from vcnmodel.analyzers import sac as SAC
from vcnmodel.analyzers import vector_strength as VS
from vcnmodel.analyzers.analyzer_data_classes import PData, RevCorrData, RevCorrPars
from vcnmodel.plotters import (
    figure_data as FD,
)  # table of simulation runs used for plotting figures

TRC = trace_calls.TraceCalls
cprint = CP.cprint
UR = pint.UnitRegistry


@dataclass
class SACPars:
    """
    Shuffled autocorrelation parameter dataclass - copied from sac.py
    """

    twin: float = 5.0
    binw: float = 0.05
    delay: float = 0.0
    dur: float = 1000.0
    displayDuration: float = 10.0
    ntestrep: int = 20
    baseper: float = 4 / 3.0
    stimdur: float = 1000


def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


def def_empty_dict():
    return {}


def norm(p: Union[list, np.ndarray], n: int) -> np.ndarray:
    """
    Simple function to normalize the n'th point of p
    by the min and max
    """

    pmin = np.min(p)
    pmax = np.max(p)
    return (p[n] - pmin) / float(pmax - pmin)


def boltzI(x, gmax, vhalf, k, E):
    return gmax * (x - E) * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def boltzG(x, gmax, vhalf, k, E):
    return gmax * (1.0 / (1.0 + np.exp(-(x - vhalf) / k)))


def expdecay(x, decay, amplitude, offset):
    return offset + amplitude * np.exp(-x / decay)


def exp2decay(x, a0, a1, tau0, tau1, offset):
    return offset + a0 * np.exp(-x / tau0) + a1 * np.exp(-x / tau1)


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


class PlotSims:
    def __init__(self, parent):
        self.parent = (
            parent  # mostly to provide access to datatables elements (display)
        )
        self.datapaths = get_data_paths()
        self.ReadModel = readmodel.ReadModel()
        self.ReadModel.set_parent(
            parent=self.parent, my_parent=self
        )  # pass our parent and US to the reader
        self.firstline = True
        self.VS = VS.VectorStrength()
        self.axis_offset = -0.02
        self.config = get_data_paths()
        mpl.style.use(Path(Path.cwd(), "styles", "figures.mplstyle"))
        self.allspikes = None
        self.in_Parallel = (
            False  # flag to prevent graphics access during parallel processing.
        )
        self.RevCorrData = RevCorrData()
        self.RevCorrPars = RevCorrPars()

    def newPData(self):
        """
        Return Pdata with the paths set from self.config
        and the GradeA Cell list.

        Returns:
            PData: dataclass
        """

        return PData(
            gradeA=GRPDEF.gradeACells,
            basepath=self.config["baseDataDirectory"],
            renderpath=str(Path(self.config["codeDirectory"], "Renderings")),
            revcorrpath=self.config["RevCorrDataDirectory"],
        )

    def textclear(self):
        if self.parent is None or self.in_Parallel:
            return
        else:
            self.parent.textbox.clear()

    def textappend(self, text, color="white"):
        if self.parent is None or self.in_Parallel:
            cprint(color, text)  # just go straight to the terminal
        else:
            self.parent.textbox.setTextColor(self.parent.QColor(color))
            self.parent.textbox.append(text)
            self.parent.textbox.setTextColor(self.parent.QColor("white"))

    @trace_calls.time_func
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

        PD = self.newPData()
        for iax, index_row in enumerate(self.parent.selected_index_rows):
            selected = self.parent.table_manager.get_table_data(index_row)
            if selected is None:
                return
            sfi = Path(selected.simulation_path, selected.files[0])
            sfi = update_disk(sfi, self.datapaths)
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

    @trace_calls.winprint
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

    @trace_calls.winprint
    @TRC()
    def plot_traces(
        self,
        ax: object,
        fn: Union[Path, str],
        PD: dataclass,
        protocol: str,
        ymin: float = -80.0,
        ymax: float = 20.0,
        xmin: float = 0.0,
        xmax: Union[float, None] = None,
        yoffset: float = 0.0,
        iax: Union[int, None] = None,
        rep: Union[int, list, None] = None,  # which rep : none is for all.
        figure: object = None,
        show_title: bool = True,
        longtitle: bool = True,
        trace_color: str = "k",
        ivaxis: object = None,
        ivcolor: str = "k",
        iv_spike_color: str = "r",
        spike_marker_size: float = 2.5,
        spike_marker_color: str = "r",
        spike_marker_shape: str = "o",
        calx: Union[float, None] = 0.0,
        caly: Union[float, None] = 0.0,
        calt: Union[float, None] = 10.0,
        calv: Union[float, None] = 20.0,
        calx2: Union[float, None] = 20.0,
        caly2: Union[float, None] = 20.0,
        calt2: Union[float, None] = 10.0,
        calv2: Union[float, None] = 10.0,
        clipping: bool = False,
        axis_index: int = 0,  # index for axes, to prevent replotting text
    ) -> tuple:
        """Plot traces in a general way
        Yes, this should be broken up with fewer parameters,
        probably with dataclasses for the cals, etc... but plotting
        is just grunge, so here we go:

        Parameters
        ----------
        ax : object
            target matplotlib axis
        fn : Union[Path, str]
            filename of data to plot
        PD : dataclass
            PData object with info about the dataset
        protocol : str
            String name of the protocol
        ymin : float, optional
            Minimum for y axis, normally mV, by default -80.0
        ymax : float, optional
            Maximum for y axis, normally mV, by default 20.0
        xmin : float, optional
            Minimum for x axis, normally msec, but may be sec, by default 0.0
        xmax : Union[float, None], optional
            Maximum for x axis, normally msec, but may be sec, by default None
        yoffset : float, optional
            Value to offset Y by with an axis (for stacked traces), by default 0.0
        iax : Union[int, None], optional
            axis number, by default None
        rep : Union[int, list, None], optional
            repetition, by default None
        show_title : bool, optional
            display descriptive title, by default True
        longtitle : bool, optional
            display a long descriptive title, by default True
        trace_color : str, optional
            color of the traces, by default "k"
        ivaxis : object, optional
            optional axis for an IV plot, by default None
        ivcolor : str, optional
            color for an ivplot trace, by default "k"
        iv_spike_color : str, optional
            color to mark spikes in IV plot, by default "r"
        spike_marker_size : float, optional
            size of spike marker in traces, by default 2.5
        spike_marker_color : str, optional
            color of spike marker in traces, by default "r"
        spike_marker_shape : str, optional
            shape of spike marker in traces, by default "o"
        calx : Union[float, None], optional
            calibration bar x length, by default 0.0
        caly : Union[float, None], optional
           calibration bar y length, by default 0.0
        calt : Union[float, None], optional
            calibration bar position along x axis (time), by default 10.0
        calv : Union[float, None], optional
            calibration bar position along y axis (voltage), by default 20.0
        calx2 : Union[float, None], optional
            secondary cal bar, by default 20.0
        caly2 : Union[float, None], optional
            secondary cal bar, by default 20.0
        calt2 : Union[float, None], optional
            secondary cal bar, by default 10.0
        calv2 : Union[float, None], optional
            secondary cal bar, by default 10.0
        axis_index : int, optional
            index for axes, to prevent replotting text, by default 0

        Returns
        -------
        tuple :
            synno : number of synapses on this cell
            noutspikes : number of spikes from the cell
            ninspikes : number of input spikes to the cell

        """

        inx = str(fn).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(str(fn)[inx + 4 : inx + 7])
        if protocol in ["IV", "runIV"]:
            protocol = "IV"
        elif protocol in ["VC", "runVC"]:
            protocol = "VC"
        fn = update_disk(fn, self.datapaths)
        model_data = self.ReadModel.get_data(fn, PD=PD, protocol=protocol)
        data = model_data.data
        si = model_data.SI
        ri = model_data.RI
        if figure is None:  # no figure... just analysis...
            return model_data.AR, model_data.SP, model_data.RM
        AR = model_data.AR
        SP = model_data.SP
        RM = model_data.RM
        cprint("c", "plot_traces: preparing for plot")
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
        if xmax is None and protocol not in ["IV", "VC"]:
            xmax = 1e3 * (ri.pip_start + ri.pip_duration)
            xmin = 1e3 * ri.pip_start
        elif xmax is not None:
            pass
        elif xmax is None and protocol in ["IV"]:
            xmax = ri.stimDelay + ri.stimDur + ri.stimPost
        elif xmax is None and protocol in ["VC"]:
            xmax = ri.vstimDelay + ri.vstimDur + ri.vstimPost
        else:
            print("xmax: ", xmax)
            print("protocol: ", protocol)
            raise ValueError("Need to specificy xmax for plot")
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
        for trial, icurr in enumerate(data["Results"]):
            if rep is not None and trial != rep:
                continue
            AR.MC.traces[trial][0] = AR.MC.traces[trial][1]
            if protocol in ["VC", "vc", "vclamp"]:
                AR.MC.traces[trial] = AR.MC.traces[trial].asarray() * 1e9  # nA
                cmd = AR.MC.cmd_wave[trial] * 1e3  # from V to mV
            else:
                AR.MC.traces[trial] = AR.MC.traces[trial].asarray() * 1e3  # mV
                cmd = AR.MC.cmd_wave[trial] * 1e9  # from A to nA
            xclip = np.argwhere((AR.MC.time_base >= xmin) & (AR.MC.time_base <= xmax))
            # plot trace
            ax1.plot(
                AR.MC.time_base[xclip],
                AR.MC.traces[trial][xclip] + yoffset,
                linestyle="-",
                color=trace_color,
                linewidth=0.5,
                clip_on=clipping,
            )
            if ax2 is not None:
                ax2.plot(AR.MC.time_base[xclip], cmd[xclip], linewidth=0.5)
            if "spikeTimes" in list(data["Results"][icurr].keys()):
                # cprint('r', 'spiketimes from results')
                # print(data["Results"][icurr]["spikeTimes"])
                #  print(si.dtIC)
                spikeindex = [
                    int(t * 1e3 / (si.dtIC))
                    for t in data["Results"][icurr]["spikeTimes"]
                ]
            else:
                cprint("r", "spikes from SP.spikeIndices")
                spikeindex = SP.spikeIndices[trial]
            # print(f"Trial: {trial:3d} Nspikes: {len(spikeindex):d}")
            # plot spike peaks
            ax1.plot(
                AR.MC.time_base[spikeindex],
                AR.MC.traces[trial][spikeindex] + yoffset,
                marker=spike_marker_shape,  # "o",
                color=spike_marker_color,
                markerfacecolor=spike_marker_color,
                markersize=spike_marker_size,
                linestyle="none",
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
                if trial in list(data["Results"].keys()) and "inputSpikeTimes" in list(
                    data["Results"][icurr].keys()
                ):
                    spkt = data["Results"][icurr]["inputSpikeTimes"]
                elif "inputSpikeTimes" in list(data["Results"].keys()):
                    spkt = data["Results"]["inputSpikeTimes"][trial]
                tr_y = trial * (trstep + len(spkt) * inpstep)
                if synno is None:
                    for ian in range(len(spkt)):
                        vy = v0 + tr_y * np.ones(len(spkt[ian])) + inpstep * ian
                        ax1.scatter(spkt[ian], vy, s=sz, marker="|", linewidths=0.35)
                else:
                    ian = synno
                    vy = v0 + tr_y * np.ones(len(spkt[ian])) + inpstep * ian
                    # ax.scatter(spkt[ian], vy, s=sz, marker="|", linewidths=0.35)
                    ninspikes += len(spkt[ian] > deadtime)

                ax1.set_ylim(ymin, ymax)
                if xmin is None:
                    xmin = 0.050
                if xmax is None:
                    xmax = np.max(AR.MC.time_base)
                ax1.set_xlim(xmin, xmax)
            elif protocol in ["VC", "vc", "vclamp"]:
                pass  #
                # ax.set_ylim((-100.0, 100.0))
            else:
                ax1.set_ylim(ymin, ymax)
                if xmin is None:
                    cprint("r", "2. xmin is None")
                    xmin = 0.050
                if xmax is None:
                    xmax = np.max(AR.MC.time_base)
                ax1.set_xlim(xmin, xmax)
        ftname = str(Path(fn).name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = ""
        if longtitle and show_title:
            toptitle = f"{ftname:s}"
        else:
            if show_title:
                toptitle = si.dendriteExpt
        if protocol in ["IV"]:
            cprint("r", f"RM analysis taum: {RM.analysis_summary['taum']:.2f}")
            if show_title:
                omega = r"$\Omega$"
                tau_m = r"$\tau_m$"
                toptitle += f"\nRin={RM.analysis_summary['Rin']:.1f} M{omega:s}  {tau_m}{RM.analysis_summary['taum']:.2f} ms"

            if iax is not None and calx is not None:
                PH.calbar(
                    ax1,
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
                RM.ivss_v_all.asarray() * 1e3,
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
                    marker=spike_marker_shape,
                    markersize=spike_marker_size,
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
                    ymin,
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
            if axis_index == 0:
                secax.text(2.0, -70.0, "nA", ha="center", va="top", fontweight="normal")
                secax.text(
                    0.0, -40.0, "mV ", ha="right", va="center", fontweight="normal"
                )
            self.traces_ax = ax1
            self.crossed_iv_ax = secax

        elif protocol in ["VC", "vc", "vclamp"]:
            maxt = np.max(AR.MC.time_base)
            if calx is not None:
                PH.calbar(
                    ax1,
                    calbar=[calx, caly, calt, calv],
                    orient="left",
                    unitNames={"x": "ms", "y": "nA"},
                    fontsize=9,
                )
            else:
                PH.noaxes(ax1)
            if ax2 is not None and calx2 is not None:
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
                if ax2 is not None:
                    PH.noaxes(ax2)
        else:
            if calx is not None and iax is not None:
                cprint("r", "**** making cal bar")
                print("calx, y, calt, v: ", iax, calx, caly, calt, calv)
                PH.calbar(
                    ax,
                    calbar=[calx, caly, calt, calv],
                    unitNames={"x": "ms", "y": "mV"},
                    fontsize=8,
                )
            else:
                PH.noaxes(ax)
            if protocol in ["IV"]:  # only valid for an IV
                if RM.analysis_summary is not None:
                    PH.referenceline(ax, RM.analysis_summary["RMP"])
                ax.text(
                    -1.0,
                    RM.analysis_summary["RMP"],
                    f"{RM.analysis_summary['RMP']:.1f}",
                    verticalalignment="center",
                    horizontalalignment="right",
                    fontsize=9,
                )
        if show_title:
            cprint("r", "ShowTitle in plot_traces")
            toptitle += f"\n{model_data.timestamp:s}"
            figure.suptitle(toptitle, fontsize=9)

        return (synno, noutspikes, ninspikes)

    @TRC()
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
                sfi[i] = Path(str(sfi[i]).replace("/._", "/"))
                with open(sfi[i], "rb") as fh:
                    sfidata = FPM.pickle_load(fh)
                    # print(sfidata['Params'])
                    dendriteMode[i] = sfidata["Params"].dendriteMode
                    protocol[i] = sfidata["runInfo"].runProtocol
        P = self.setup_VC_plots(n_columns, parent_figure=parent_figure, loc=loc)
        titlemap = {"normal": "Half-active", "passive": "Passive", "active": "Active"}

        for i in range(n_columns):
            PD = self.newPData()
            trace_ax = P.axarr[i * 3, 0]
            cmd_ax = P.axarr[i * 3 + 1, 0]
            if i == 0:
                calx = 120.0
            else:
                calx = None
            self.plot_traces(
                ax=[trace_ax, cmd_ax],
                fn=sfi[i],
                PD=PD,
                protocol=protocol[i],
                calx=calx,
                caly=5.0,
                calt=10.0,
                calv=2.0,
                calx2=calx,
                caly2=-40.0,
                calt2=10.0,
                calv2=20.0,
                xmax=150.0,
                figure=P.figure_handle,
                clipping=True,
                show_title=False,
            )
            trace_ax.set_ylim((-1, 15))
            trace_ax.set_clip_on(True)
            self.analyzeVC(P.axarr[i * 3 + 2, 0], sfi[i], PD, protocol=protocol[i])
            PH.nice_plot(
                P.axarr[i * 3 + 2, 0], position=-0.03, direction="outward", ticklength=3
            )
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

    @TRC()
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

    @trace_calls.winprint
    @trace_calls.time_func
    @TRC()
    def analyzeVC(
        self,
        ax: object,
        fn: Union[Path, str],
        PD: dataclass,
        protocol: str,
    ) -> tuple:
        # if x is None:
        #     self.textappend("No simulation found that matched conditions", color="red")
        #     self.textappend(fn, color="red")
        #     return
        # # unpack x
        inx = str(fn).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(fn[inx + 4 : inx + 7])
        if protocol in ["VC", "runVC"]:
            protocol = "VC"
        self.textappend(f"Protocol: {protocol:s}")
        model_data = self.ReadModel.get_data(fn, PD=PD, protocol=protocol)
        if not model_data.success:
            return None
        data = model_data.data
        AR = model_data.AR

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
        boltz_result = gmodel.fit(
            gss_sel, method="nedler", params=gparams, x=vss_sel, weights=weights_sel
        )
        boltz_x = np.linspace(np.min(vss_sel), np.max(vss_sel), 100, endpoint=True)
        boltz_fit = gmodel.eval(params=boltz_result.params, x=boltz_x)

        # capacitance transient. Use the single trace nearest to -70 mV for the fit.
        mHypStep = np.argmin(np.fabs(vss - (-0.090)))
        sr = AR.MC.sample_rate[mHypStep]
        t0 = int(AR.MC.tstart / sr)
        print(AR.MC.tstart, AR.MC.tdur, AR.MC.tend)
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
        deltaV = (vss[mHypStep] - vm[mHypStep]) * UR.V
        deltaI = (I[mHypStep][t0] - i0[mHypStep]) * UR.A
        deltaI2 = (exp_result.params["a0"].value + exp_result.params["a1"].value) * UR.A
        qpts = int(20.0 / sr)  # 20 msec duration for integration
        qend = int(AR.MC.tend / sr)
        Q = (
            np.trapz(
                I[mHypStep, t0 : t0 + qpts]
                - np.mean(I[mHypStep, t0 + qpts : qend]),  # remove DC part
                dx=sr * 1e-3,
            )
            * UR.C
        )
        Cm_U = (Q / deltaV).to(UR.F)  # with units
        Aq = Cm_U / (0.9 * UR.uF / (UR.cm * UR.cm))
        Aq = Aq.to(UR.um * UR.um)
        self.textappend(
            f"Q: {Q:.3f#~P}  deltaV: {deltaV:.1f#~P}, Cm: {Cm_U:.1f#~P}  Area: {Aq:.1f#~P})"
        )
        Rs_est = (deltaV / deltaI2).to(UR.ohm)  # Estimate of Rs from peak current
        self.textappend(f"Estimated Rs: {Rs_est:.1f#~P}")
        # print(deltaV, deltaI, deltaI2)
        # other approach: use fastest tau in voltage clamp
        #
        tau0 = exp_result.params["tau0"].value * UR.ms
        tau1 = exp_result.params["tau1"].value * UR.ms
        a0 = exp_result.params["a0"].value * UR.A
        a1 = exp_result.params["a1"].value * UR.A
        #
        # Here we use the fastest time constant and the
        # associated current from the fit to estimate cm (perisomatic)
        # Note: do not use total input resistance!
        # See Golowasch et al., J. Neurophysiol. 2009 and references therein
        print(tau0, tau1)
        if tau0 < tau1:
            tau = tau0
            R0 = (deltaV / a0).to(UR.ohm)
            cm = (tau / R0).to(UR.F)
        else:
            tau = tau1
            R0 = (deltaV / a1).to(UR.ohm)
            cm = (tau / R0).to(UR.F)
        # for curiosity, weighted tau (as if compensation was done as
        # much as possible for the whole transient)
        #
        tauw = (a0 * tau0 + a1 * tau1) / (a0 + a1)
        R0w = (deltaV / a0).to(UR.ohm)
        R1w = (deltaV / a1).to(UR.ohm)
        cm1 = (tau1 / R1w).to(UR.F)
        cmw = ((a0 * tau0 / R0w) + (a1 * tau1 / R1w)) / (a0 + a1)  # tauw/(R0w + R1w)
        cmw = cmw.to(UR.F)
        Acm = cm / (0.9 * UR.uF / (UR.cm * UR.cm))
        Acmw = cmw / (0.9 * UR.uF / (UR.cm * UR.cm))
        self.textappend("By Coeffs: ")
        self.textappend(
            f"    RCoeff0 = {R0w:.2f#~P}, tau0: {tau0:6.1f#~P},  cm0: {cm:.1f#~P}, Area: {Acm:.1f#~P}"
        )
        self.textappend(
            f"    RCoeff1 = {R1w:.2f#~P}, tau1: {tau1:6.1f#~P},  cm1: {cm1:.1f#~P}"
        )
        self.textappend(
            f"Weighted: Rw={(R0w+R1w):.2f#~P} tauw: {tauw:6.1f#~P}, Weighted cm: {cmw:.1f#~P}, Area: {Acmw:.1f#~P}"
        )
        tfit2 = AR.MC.time_base[t0 : t0 + pts] - AR.MC.time_base[t0]
        expfit = expmodel.eval(params=exp_result.params, x=tfit2)
        # ax[0].plot(
        #     tfit2 + AR.MC.time_base[t0], bfit * 1e9, "r-", dashes=[6, 2], linewidth=1
        # )
        # ax[0].plot(tfit+AR.MC.time_base[t0], exp_result.best_fit*1e9, 'r-')
        # ax[0].plot(tfit+AR.MC.time_base[t0], ifit*1e9, 'g--')

        ax.plot(vss_sel * 1e3, gss_sel * 1e9, "ko", markersize=2, clip_on=False)
        ax.plot(boltz_x * 1e3, boltz_fit * 1e9, "r-", clip_on=False)
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
        # Align the parameters on the = sign by making 2 texts right and left justified
        textstr1 = (
            r"$\mathrm{g_{max}}$" + f" = {boltz_result.params['gmax'].value:.1f} nS"
        )
        textstr2 = (
            r"$\mathrm{V_{0.5}}$"
            + f" = {boltz_result.params['vhalf'].value*1e3:.1f} mV"
        )
        textstr3 = r"$\mathrm{k}$" + f" = {1e3*boltz_result.params['k'].value:.1f}"
        textstr4 = r"$\mathrm{C_m}$" + f" = {cm:.1f#~P}"
        textstr5 = r"$\mathrm{\tau_{0}}$" + f" = {tau:6.1f#~P}"
        texts = [textstr1, textstr2, textstr3, textstr4, textstr5]
        props = dict(boxstyle="square", facecolor="None", alpha=0.5)

        y = 0.40
        x = 0.60
        for i, txt in enumerate(texts):
            tx1, tx2 = txt.split("=")
            ax.text(
                x,
                y - i * 0.1,
                tx1,
                transform=ax.transAxes,
                #  fontfamily="monospace",
                size=8,
                verticalalignment="bottom",
                horizontalalignment="right",
                # bbox=props,
            )
            ax.text(
                x,
                y - i * 0.1,
                " = " + tx2,
                transform=ax.transAxes,
                # fontfamily="monospace",
                size=8,
                verticalalignment="bottom",
                horizontalalignment="left",
                # bbox=props,
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
        model_data = self.ReadModel.get_data(filename, PD=PD)
        if model_data is None or not model_data.success:
            print("No simulation found that matched conditions")
            print(filename)
            return

        inx = str(filename).find("_Syn")
        synno = None
        if inx > 0:
            synno = int(str(filename)[inx + 4 : inx + 7])
        if runProtocol in ["IV", "runIV"]:
            runProtocol = "IV"
        elif runProtocol in ["VC", "runVC"]:
            runProtocol = "VC"
        # par, stitle, ivdatafile, filemode, d = x
        AR = model_data.AR
        RM = model_data.RM
        SP = model_data.SP

        self.ntr = len(AR.MC.traces)  # number of trials
        self.pwin = None
        self.pfull = None
        self.lines = [None] * self.parent.n_trace_sel
        self.punchtape = None
        self.pt_text = None
        self.inputtimes = None
        self.first = True
        self.parent.trace_selector.setValue(0)
        model_data.SP.analyzeSpikes()
        self.sample_interval = model_data.SP.Clamps.sample_interval
        print("dt: ", self.sample_interval)
        if self.allspikes is not None:
            self.nspikes = len(self.allspikes)
            for i, trn in enumerate(range(0, self.parent.n_trace_sel)):
                self.plot_spikes_withinputs(
                    int(i),
                    n=trn,
                    color=pg.intColor(i, hues=20),
                    first=self.first,
                    dt=self.sample_interval,
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
                            dt=self.sample_interval,
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
                    dt=self.sample_interval,
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
                    dt=self.sample_interval,
                )

    def plot_spikes_withinputs(
        self,
        ix: int = 0,
        n: int = 0,
        color: object = None,
        first=False,
        dt: float = 1.0,
    ):
        """
        plot a spike and indicate its inputs.
        
        Params
        ------
            ix : index into this run (counter): for plotting a block of spikes
            n : spike to plot within the block
            color: indexed color.
            first: bool: if first time plotting
            dt : float - not used
        """

        # self.parent.trace_plots.plot(np.arange(10), np.ones(10))
        if self.nspikes > n and n >= 0:
            self.check_yscaling()
            # print('n: ', n)
            spk = self.allspikes[n]
            if self.parent.V_disp_sel == "dV/dt":
                if first:
                    self.lines[ix] = self.parent.trace_plots.plot(
                        self.sample_interval * np.arange(0, len(spk.waveform))[:-1]
                        + spk.start_win,
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
                        self.sample_interval * np.arange(0, len(spk.waveform))[:-1]
                        + spk.start_win,
                        np.diff(spk.waveform / spk.dt),
                        pen=color,
                    )
                papertape_yoff = 120.0
                spkin_off = -50.0
                dy = 5.0
            else:
                if first:
                    print("plotting first")
                    tx = (
                        self.sample_interval * np.arange(0, len(spk.waveform))
                        + spk.start_win
                    )
                    # print(np.min(tx), np.max(tx), tx.shape)
                    #                    print(np.min(spk.waveform), np.max(spk.waveform), spk.waveform.shape)
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
                    tx = (
                        self.sample_interval * np.arange(0, len(spk.waveform))
                        + spk.start_win
                    )
                    self.lines[ix].curve.setData(
                        x=tx,
                        y=spk.waveform,
                        pen="w",  # pen=color,
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
        cell_number: Union[int, None] = None,
        axarray=None,
        calbar_show=True,
        calbar_fontsize=11,
        yaxis_label=True,
        reflevel_show=False,
        start_letter="A",
        colormap="Set3",
        show_average: bool = False,
        synlabel: bool = True,
        cbar_vmax: float = 300.0,
    ):
        sns.set_style("ticks")
        max_example_spikes = 10
        # secax = PF.twinax(P.figure_handle, ax, pos=maxwin)
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

        if cell_number is not None:
            cell_n = cell_number
        elif isinstance(self.parent.cellID, int):
            cell_n = self.parent.cellID
        else:
            cell_n = int(self.parent.cellID[-2:])

        SC, syninfo = self.get_synaptic_info(cell_n)
        syn_ASA = np.array([syninfo[1][isite][0] for isite in range(RCP.ninputs)])
        max_ASA = np.max(syn_ASA)
        for isite in reversed(range(RCP.ninputs)):
            # print(dataclasses.fields(RCD))
            cprint("y", f"len rcd sv trials: {len(RCD.sv_trials):d}")
            if len(RCD.sv_trials) == 0:
                continue
            # if not 'sv_trials' in dataclasses.fields(RCD) or len(RCD.sv_trials) == 0:
            #     cprint("r", "SV_TRIALS not found in the fields for RCD")
            #     continue
            # stepsize = int(RCD.sv_all.shape[0] / 20)
            if isite == 0:  # all sites have the same voltage, but here
                # we select some spikes waveforms
                prewin = np.where(RCD.ti_avg < -0.5)
                nprewin = len(prewin)
                ex_spikes = []
                for trial in range(len(RCD.sv_trials)):
                    # if isite == 0:
                    #     print(np.max(RCD.sv_all[trial][prewin]))
                    for spk in RCD.sv_trials[trial]:
                        if np.max(spk[prewin]) < -40.0:
                            if len(ex_spikes) < max_example_spikes:
                                ex_spikes.append(spk)
            # if stepsize > 0:
            #     sel = list(range(0, RCD.sv_all.shape[0], stepsize))
            # else:
            #     sel = list(range(0, RCD.sv_all.shape[0], 1))
            # sel = list(range(0, RCD.sv_all.shape[0], 1))
            refzero = int(RCP.minwin / RCP.binw)

            if RCD.C[isite] is not None:
                nc = int(len(RCD.C[isite]) / 2)
                RCD.TC = RCD.TC / len(RCD.st)
                summarySiteTC[isite] = RCD.TC
                color = cmx.colors[int(cmx.N * syn_ASA[isite] / cbar_vmax) - 1]
                colors[isite] = color
                maxrevcorr = np.max((maxrevcorr, np.max(RCD.CB[isite])))
                totalrevcorr = np.sum(RCD.CB[isite])
                label = "Input {0:2d} N={1:3d}".format(isite, int(RCD.sites[isite]))
                label = None
                print("************************* RCP ALGORITHM::: ", RCP.algorithm)
                if isite in range(RCP.ninputs):
                    if RCP.algorithm == "RevcorrEleph":
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
                            RCD.CBT,  # - RCP.maxwin,
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
            ax.set_ylim(0, 0.8)
            ax.set_xlim(-5, 2.5)
            PH.talbotTicks(  # set ticks for revcorrs
                ax,
                tickPlacesAdd={"x": 1, "y": 1},
                floatAdd={"x": 1, "y": 1},
                # pointSize=7,
            )
            if (
                isite == 0
            ):  # only plot the first time through - the APs are the same no matter the input
                for t in range(len(ex_spikes)):  # this plots the individual traces
                    secax.plot(
                        RCD.ti_avg,
                        ex_spikes[t],
                        color="#666666",
                        linewidth=0.2,
                        zorder=1,
                    )
                # on top of that plot the average trace
                if show_average:
                    secax.plot(
                        RCD.ti_avg, RCD.sv_avg, color="k", linewidth=0.75, zorder=2
                    )
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
                if reflevel_show:
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
                PH.nice_plot(secax, position=self.axis_offset, direction="outward")
                PH.referenceline(secax, -60.0)
                PH.talbotTicks(
                    secax,
                    density=(1.0, 0.5),
                    insideMargin=0,
                    tickPlacesAdd={"x": 1, "y": 0},
                    floatAdd={"x": 1, "y": 0},
                    axrange={"x": (-5.0, 2.5), "y": (-60, 0)},
                    pointSize=None,
                )
                PH.showaxes(secax)
        if yaxis_label:
            secax.set_ylabel("Postsynaptic Voltage", fontsize=9, x=0.5)
        print(f"Total spikes in voltage trace: {RCD.nsp_avg:d}")

        # finally, put up a legend that shows which inputs map to what colors
        # # the colors are in the legend
        PH.nice_plot(ax, position=self.axis_offset, direction="outward")

        if yaxis_label:
            ax.set_ylabel("Pre-post\nCoinc. Rate (Hz)")
        # ax.set_xlabel("T (ms)")
        secax.set_xlabel("T (ms)")
        if RCD.max_coin_rate > 0.0:
            ns = PH.NiceScale(0.0, RCD.max_coin_rate)
            ax.set_ylim((0, ns.niceMax))
        else:
            ax.set_ylim((0, 0.25))
        ax.set_xlim((RCP.minwin, RCP.maxwin))
        secax.set_ylim([-70.0, 10.0])
        secax.set_xlim((RCP.minwin, RCP.maxwin))
        secax.tick_params(direction="in", length=3.0, width=1.0)
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            pointSize=None,
            tickPlacesAdd={"x": 1, "y": 1},
            floatAdd={"x": 1, "y": 1},
            axrange={"x": (-5.0, 2.5), "y": (0, 0.8)},
        )
        if synlabel:
            self.axins = PH.make_colorbar(
                ax,
                bbox=[0.125, 0.7, 0.75, 0.2],
                vmin=0,
                vmax=300.0,
                nticks=4,
                palette=colormap,
            )
            self.axins.tick_params(axis="x", length=2, labelsize=6)
            self.axins.text(
                x=0.5,
                y=1.25,
                s=r"Input ASA (${\mu m^2}$)",
                horizontalalignment="center",
                fontsize=8,
            )
            # transform=self.axins.transAxes)

        return summarySiteTC

    def get_synaptic_info(self, gbc: str, add_inputs="none") -> tuple:
        SC = cell_config.CellConfig(add_inputs=add_inputs)
        # print(SC.VCN_Inputs)
        if isinstance(gbc, int):
            gbc_string = f"VCN_c{gbc:02d}"
        elif isinstance(gbc, str) and len(gbc) <= 2:
            if gbc.isnumeric():
                gbc_string = f"VCN_c{int(gbc):02d}"
            else:
                clean_gbc = re.sub(r"\D", "", gbc)  # remove non-numeric parts
                gbc_string = f"VCN_c{int(clean_gbc):02d}"
        elif isinstance(gbc, str) and gbc.startswith("VCN_"):
            gbc_string = gbc
        elif isinstance(gbc, str) and gbc.startswith("BC"):
            gbc_string = gbc.replace("BC", "VCN_c")
        syninfo = SC.VCN_Inputs[gbc_string]
        return (SC, syninfo)

    @trace_calls.time_func
    def compare_revcorrs(self):
        plabels = [f"{FD.BC_name:s}{int(self.parent.cellID):02d}"]
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
                None, pgbc, sfi, self.newPData(), selected.runProtocol, revcorrtype
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
        PD = self.newPData()

        plabels = [f"{FD.BC_name:s}{int(self.parent.cellID):02d}"]
        pgbc = plabels[0]
        sizer = {
            "A": {
                "pos": [0.05, 0.35, 0.52, 0.40],
                "labelpos": (0.02, 1.00),
                "noaxes": True,
            },
            "B": {"pos": [0.05, 0.35, 0.08, 0.35], "labelpos": (0.02, 1.00)},
            "C": {"pos": [0.45, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
            "D": {"pos": [0.45, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
            "E": {"pos": [0.70, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
            "F": {"pos": [0.70, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
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

        dPD = self.newPData()
        sfi = Path(selected.simulation_path, selected.files[0])
        res = self.compute_revcorr(P, pgbc, sfi, PD, selected.runProtocol, revcorrtype)

        P.figure_handle.show()

    @trace_calls.winprint
    # @trace_calls.time_func
    def compute_revcorr(
        self,
        P: Union[object, None] = None,
        gbc: str = "",
        fn: Union[str, Path] = "",
        PD: object = None,
        protocol: str = "",
        revcorrtype: str = "RevcorrSPKS",
        width: float = 4.0,
    ) -> Union[None, tuple]:
        """

        Compute reverse correlation from data
        This is a wrapper that sets up the call to the routine in
        reverse_correlation.py

        Parameters
        ----------
        P : object
            The figure object from plothelpers
        gbc : str
            The name of the gbc ('VCN_cnn')
        fn : Path or str
            name of the source data file
        PD : dataclass
            The dataclass of the Data.
        protocol : str
            name of the acq4 protocol
        revcorrtype: str
            revcorr algorithm to use. The only valid (tested)
            protocol is RevcorrSPKS
        width : width of the revcorr calculation, in msec.

        Returns
        -------
        RCP: updated parameters with results.
        RCD: updated data structure with results.
        """

        self.allspikes = None
        #
        # 1. Gather data
        #
        print(f"    compute_revcorr  Getting data for gbc: {gbc:s}")
        fn = update_disk(fn, self.datapaths)
        MD = self.ReadModel.get_data(fn, PD=PD, protocol=protocol)
        if not MD.success:
            return None
        # (d, AR, SP, RM, RCP, RCD) = res
        PD.thiscell = gbc
        RCP = MD.RCP
        RCD = MD.RCD
        RCP.algorithm = revcorrtype
        SC, syninfo = self.get_synaptic_info(gbc)

        print(f"    compute_revcorr: Preparing for computation for: {str(gbc):s}")
        RCP.ninputs = len(syninfo[1])
        print("ninputs from syninfo: ", RCP.ninputs)
        self.ninputs = RCP.ninputs  # save for trace viewer.
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
            min_time = (
                MD.SI["pip_start"] + 0.025
            ) * 1000.0  # push onset out of the way
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
        print("    compute_revcorr experiment type: ", expt)

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
        if self.parent is not None:
            revcorr_params = {
                "deselect": self.parent.deselect_flag,
                "threshold": self.parent.deselect_threshold,
                "window": self.parent.revcorr_window,
                "ASAs": syn_ASA,
            }
        else:
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
            MD, printflag=True
        )  # perfrom an analysis of input spike patters
        self.allspikes = PAT.all_spikes
        if P is not None:
            self.plot_revcorr_details(P, PD, MD2, RCP, RCD)
        return P, PD, RCP, RCD

    def plot_revcorr_details(self, P, PD, MD, RCP, RCD):
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
        E_plot = sax["E"].plot(RCD.sites, RCD.participation, "gx")

        axcbar = PLS.create_inset_axes(
            [0.8, 0.05, 0.05, 0.5], sax["D"], label=str(P.axdict["D"])
        )
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        # ticks = [int(p) for p in np.linspace(vmin, vmax, num=4, endpoint=True)]
        ticks = [int(p) for p in np.linspace(0, vmax, num=5, endpoint=True)]
        cm_sns = mpl.colormaps[colormap]
        colorbar = matplotlib.colorbar.ColorbarBase(
            axcbar, cmap=cm_sns, ticks=ticks, norm=norm
        )
        ticks = [int(p) for p in np.arange(1, RCP.ninputs + 0.5, 1)]

        # PH.nice_plot(sax['C'], position=-0.2)
        F_plot = sax["F"].plot(
            np.arange(len(RCD.ynspike)) + 1, np.cumsum(RCD.ynspike), "m^-"
        )

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
            sax["C"],
            pointSize=7,
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
            f"Cell {PD.thiscell:s} {str(MD.SI.shortSimulationFilename):s}\n[{MD.RI.runTime:s}] dB:{MD.RI.dB:.1f} Prot: {MD.RI.runProtocol:s}"
            + f"\nExpt: {MD.RI.Spirou:s}  DendMode: {MD.SI.dendriteMode:s} DendExpt: {MD.SI.dendriteExpt:s}"
            + f"\nDepr: {MD.SI.ANSynapticDepression:d} ANInput: {MD.SI.SRType:s}",
            fontsize=11,
        )
        self.F_cursor = mplcursors.cursor(F_plot, hover=2)
        self.E_cursor = mplcursors.cursor(E_plot, hover=2)

        @self.F_cursor.connect("add")
        def on_add_F(select):
            select.annotation.set(text=f"Y={RCD.ynspike[int(select.index)]:6.3f}")

        @self.E_cursor.connect("add")
        def on_add(select):
            select.annotation.set(text=f"Y={RCD.participation[int(select.index)]:6.3f}")

        mpl.show()

        return summarySiteTC, RCD.sites, sax

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
                "leftmargin": 0.125,
                "rightmargin": 0.05,
                "topmargin": 0.15,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )

        PD = self.newPData()
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
                iax=i,
                figure=P.figure_handle,
            )
            sfi[i] = update_disk(sfi[i], self.datapaths)
            model_data = self.ReadModel.get_data(
                sfi[i],
                PD=PD,
                protocol=selected.runProtocol,
            )
            if not model_data.success:
                return None

            PD.thiscell = self.parent.cellID
            si = model_data.SI
            ri = model_data.RI
            try:
                print("Ri: ", ri.SpirouSubs)
            except:
                ri.SpirouSubs = "none"
            SC, syninfo = self.get_synaptic_info(
                self.parent.cellID, add_inputs=ri.SpirouSubs
            )
            r, ctype = SC.make_dict(f"VCN_c{int(self.parent.cellID):02d}")
            if synno is None:
                raise NotImplementedError(
                    f"Plotting singles is not permitted for protocols with synaptic input (detected protocol: {selected.runProtocol:s})"
                )
            isite = i
            area = syninfo[1][isite][0]
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
        self.textappend(df.to_csv(sep="\t"))
        P.figure_handle.show()

    @trace_calls.time_func
    def plot_tuning(self, args, filename=None, filenames=None):
        PD = self.newPData()
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

    @TRC()
    def setup_PSTH(self):
        sizer = OrderedDict(  # define figure layout
            [
                ("A", {"pos": [0.08, 0.4, 0.78, 0.16]}),
                ("B", {"pos": [0.08, 0.4, 0.58, 0.16]}),
                ("C", {"pos": [0.08, 0.4, 0.34, 0.20]}),
                ("D", {"pos": [0.08, 0.4, 0.18, 0.12]}),
                ("E", {"pos": [0.08, 0.4, 0.05, 0.09]}),
                # rhs
                ("F", {"pos": [0.55, 0.18, 0.78, 0.16]}),  # vs and sac side by side
                ("G", {"pos": [0.78, 0.18, 0.78, 0.16]}),
                ("H", {"pos": [0.55, 0.4, 0.58, 0.16]}),
                ("I", {"pos": [0.55, 0.4, 0.34, 0.20]}),
                ("J", {"pos": [0.55, 0.4, 0.18, 0.12]}),
                ("K", {"pos": [0.55, 0.4, 0.05, 0.09]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.

        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            label=True,
            figsize=(8.0, 8.0),
            # labelposition=(-0.05, 1.02),
        )
        # PH.show_figure_grid(P)  # use to configure the plot.
        # mpl.show()
        # return
        P.axdict["A"].set_ylabel("mV", fontsize=8)

        P.axdict["B"].set_title("Bushy Spike Raster", fontsize=9)
        P.axdict["B"].set_ylabel("Trial")

        P.axdict["C"].set_title("PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["C"].set_ylabel("Spikes/second", fontsize=9)
        P.axdict["D"].set_title("ISI", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["D"].set_ylabel("# spikes")

        P.axdict["E"].set_title("Stimulus", fontsize=9)
        P.axdict["E"].set_ylabel("Amplitude (Pa)", fontsize=8)
        P.axdict["E"].set_xlabel("T (s)", fontsize=9)

        P.axdict["F"].set_title("Phase", fontsize=9)
        P.axdict["G"].set_title("SAC", fontsize=9)

        P.axdict["H"].set_title("ANF Spike Raster", fontsize=9)
        P.axdict["I"].set_title("ANF PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["J"].set_title("ANF ISI", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["K"].set_title("BU FSL/SSL", fontsize=9)
        for axl in [
            "B",
            "C",
            "E",
            "H",
            "I",
        ]:
            P.axdict[axl].sharex(P.axdict["A"])
        for axl in ["J"]:
            P.axdict[axl].sharex(P.axdict["D"])
        return P

    @trace_calls.winprint_continuous
    def print_VS(
        self,
        cell_number,
        d,
        freq,
        dmod,
        dB,
        experiment,
        filename: str = "",
        save_to_clipboard: bool = False,
    ):
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

        VS_colnames = f"Cell,Filename,Configuration,carrierfreq,frequency,dmod,dB,"
        line = f"{str(cell_number):s},{filename:s},{experiment:s},"
        line += f"{d.carrier_frequency:.1f},{freq:06.1f},{dmod:.1f},{dB:.1f},"

        VS_colnames += f"VectorStrength,SpikeCount,rMTF,entrainment,phase,phasesd,Rayleigh,RayleighP,VS_mean,VS_SD,VS_Ns,VS_groups,"
        line += f"{d.vs:.4f},"
        line += f"{d.n_spikes:d},"
        line += f"{d.rMTF:.2f},"
        line += f"{d.entrainment:.4f},"
        line += f"{d.circ_phaseMean:.4f},"
        line += f"{d.circ_phaseSD:.4f},"
        line += f"{d.Rayleigh:.4f},"
        line += f"{d.pRayleigh:.4e},"
        line += f"{d.vs_mean:.4f},"
        line += f"{d.vs_sd:.4f},"
        line += f"{int(d.vs_Ns):d},"
        line += f"{int(d.vs_groups):d},"

        VS_colnames += f"AN_VS,AN_rMTF,AN_entrainment,AN_phase,AN_phasesd,SAC_AN,SAC_Bu,SAC_AN_HW,SAC_Bu_HW,maxArea,ninputs"
        line += f"{d.an_vs:.4f},"
        line += f"{d.an_rMTF/d.n_inputs:.2f},"  # normalize across all inputs
        line += f"{d.an_entrainment:.4f},"
        line += f"{d.an_circ_phaseMean:.4f},"
        line += f"{d.an_circ_phaseSD:.4f},"
        line += f"{d.sac_an_CI:.4f},"
        line += f"{d.sac_bu_CI:.4f},"
        line += f"{d.sac_an_hw:.6f},"
        line += f"{d.sac_bu_hw:.6f},"
        line += f"{amax:.4f},"
        line += f"{d.n_inputs:d}"
        if self.firstline:
            self.textappend(VS_colnames)
        self.textappend(line)
        # put on clipboard
        if save_to_clipboard:
            pyperclip.copy(line)
            pyperclip.paste()

        return VS_colnames, line

    @trace_calls.winprint_continuous
    @TRC()
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
        PD = self.newPData()
        protocol = selected.runProtocol
        sfi = sorted(selected.files)
        data_window = [0.1, 1.0]
        plot_dur = np.fabs(np.diff(data_window))
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
            fn = update_disk(fn, self.datapaths)
            model_data = self.ReadModel.get_data(fn, PD=PD, protocol=protocol)
            if not model_data.success:
                return None
            data = model_data.data
            AR = model_data.AR
            SP = model_data.SP
            RM = model_data.RM
            RCP = model_data.RCP
            RCD = model_data.RCD
            PD.thiscell = gbc
            si = model_data.SI
            ri = model_data.RI
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
            mpl.get_current_fig_manager().canvas.setWindowTitle(tstring)
            fstring = f"{gbc:s}_{soundtype:s}_{ri.dB} dBSPL"
            P.figure_handle.suptitle(fstring, fontsize=13, fontweight="bold")

            all_bu_st = []
            all_bu_st_trials = []
            ntr = len(AR.MC.traces)  # number of trials
            for i in range(ntr):  # for all trials in the measure.
                time_base = AR.MC.time_base / 1000.0  # convert to seconds
                trd = data["Results"][i]
                dt = si.dtIC / 1000.0  # convert from msec to seconds
                idx = (int(data_window[0] / dt), int(data_window[1] / dt))
                if i == 0:
                    waveform = trd["stimWaveform"]
                stb = trd["stimTimebase"]  # convert to seconds
                stim_dt = stb[1] - stb[0]
                if i == 0:
                    n_inputs = len(trd["inputSpikeTimes"])
                try:
                    if (
                        len(trd["spikeTimes"]) > 0
                        and np.nanmax(trd["spikeTimes"]) > 2.0
                    ):  # probably in msec
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
            sac = SAC.SAC()
            yh, bins = sac.SAC_with_histo(
                all_bu_st_trials, pars=pars, engine="cython", dither=dt / 2.0
            )
            sac_bu_CI, peaks, HW, FW = sac.SAC_measures(yh, bins)
            sac_bu_hw = HW[0][0] * pars["binw"]
            print("BU SAC Report: \n  ")
            print(
                f"    HW:    {sac_bu_hw:.6f}  CI: {sac_bu_CI:.2f}  Left IPS: {HW[2][0]:.2f}  Right IPS: {HW[3][0]:.2f}"
            )
            # Fs: float = 100e3  # cochlea/zilany model rate
            # F0: float = 16000.0  # stimulus frequency
            # dB: float = 30.0  # in SPL
            # RF: float = 2.5e-3  # rise-fall time
            # fmod: float = 20  # hz, modulation if SAM
            # dmod: float = 0.0  # percent if SAM
            if ri.soundtype.endswith("Clicks"):
                sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL  HW={1e3*sac_bu_hw:.3f} ms  CI={sac_bu_CI:6.2f}"
            else:
                sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL Fmod={ri.fmod:5.1}fHz Dmod={ri.dmod:5.1f}%"

            # return
            P.axdict["A"].plot(
                bins[:-1],
                yh,
                #   'k-',
                label=sac_label,
            )
            self.textappend(sac_label)

        P.axdict["C"].set_xlim(pip_start, pip_duration - pip_start)
        P.axdict["B"].set_xlim(pip_start, pip_duration - pip_start)
        P.axdict["A"].set_xlim(-pars["twin"], pars["twin"])
        P.axdict["A"].set_ylim(0, sac_bu_CI * 1.05)
        P.axdict["A"].legend(prop={"family": "monospace"}, fontsize=7)

        P.axdict["B"].get_shared_x_axes().join(P.axdict["B"], P.axdict["C"])

        mpl.show()

    @TRC()
    def psth_vs(self):
        """
        Generate tables of vs measures for all cells
        across the frequencies listed
        """

        self.textclear()  # just at start
        PD = self.newPData()
        if self.parent.selected_index_rows is None:
            return
        # P = self.PLT.setup_PSTH()
        P = None
        PD = self.newPData()
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
            self.plot_AN_response(P, sfi, PD, fn=selected.runProtocol)
        self.firstline = True

    def _ntrials(self, spike_times: Union[list, np.ndarray]):
        if spike_times.ndim == 1:  # get trial count by looking at data
            num_trials = 1
        else:
            num_trials = spike_times.shape[0]
        return num_trials

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

    @trace_calls.winprint_continuous
    @TRC()
    def plot_AN_response(
        self,
        P: Union[object, None] = None,
        fn: Union[str, Path, None] = None,
        PD: object = None,
        protocol: str = "",
        sac_flag: bool = True,  # set true to compute SAC as well as standard VS
        sac_engine: str = "cython",
        filename: str = "",
        make_VS_raw: bool = False,
        cell_ID: Union[str, int, None] = None,
    ):
        if cell_ID is None:
            cell_ID = self.parent.cellID
        if isinstance(cell_ID, str):
            if ("U" in cell_ID) or ("I" in cell_ID):
                gbc = f"VCN_c{int(cell_ID[0:1]):02d}"
                cellN = int(cell_ID[0])
            else:
                gbc = f"VCN_c{int(cell_ID[0:2]):02d}"
                cellN = int(cell_ID[0:2])
        else:
            gbc = f"VCN_c{int(cell_ID):02d}"
            cellN = cell_ID

        plotflag = False
        SC, syninfo = self.get_synaptic_info(cellN)
        fn = update_disk(fn, self.datapaths)
        model_data = self.ReadModel.get_data(fn, PD, protocol)
        AR = model_data.AR
        data = model_data.data
        ntr = len(AR.MC.traces)  # number of trials
        nreps = ntr
        n_inputs = len(data["Results"][0]["inputSpikeTimes"])
        show_vtrial = 0  # which trial to show for votage
        v0 = -160.0
        trstep = 25.0 / ntr
        inpstep = 5.0 / ntr
        sz = 50.0 / ntr
        si = model_data.SI
        ri = model_data.RI
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

        if P is not None:
            plotflag = True
            sim_dir = fn.parts[-2]
            tstring = f"{gbc:s}_{str(sim_dir):s}_{ri.soundtype:s}.pdf"
            mpl.get_current_fig_manager().set_window_title(tstring)
            bu_voltage_panel = "A"
            bu_raster_panel = "B"
            bu_psth_panel = "C"
            bu_isi_panel = "D"
            stimulus_panel = "E"
            vs_panel = "F"
            sac_panel = "G"
            an_raster_panel = "H"
            an_psth_panel = "I"
            an_isi_panel = "J"
            bu_fsl_panel = "K"
        # We must get the units correct here.
        # AR time_base is in milliseconds, so convert to seconds
        all_an_st = []
        all_bu_st = []
        all_bu_st_trials = []  # bushy spike times, but per trial

        an_st_by_input = [
            [] for x in range(n_inputs)
        ]  # accumlate by input across trials
        an_st_grand = [
            [] for x in range(ntr)
        ]  # accumulate all input spike times for each trial

        waveform = []
        data_window = (
            pip_start - 0.05,
            pip_start + pip_duration + 0.05,
        )  # manually change to shorten...
        plot_dur = np.fabs(np.diff(data_window))

        for i_trial in range(ntr):  # for all trials in the measure.
            if i_trial == show_vtrial:
                v_show_trial = (
                    AR.MC.traces[i_trial] * 1e3
                )  # save for display, and convert to mV
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            dt = si.dtIC / 1000.0  # convert from msec to seconds
            trd = data["Results"][i_trial]
            idx = (int(data_window[0] / dt), int(data_window[1] / dt))
            waveform.append(trd["stimWaveform"].tolist())
            stb = trd["stimTimebase"]  # convert to seconds

            if len(trd["spikeTimes"]) == 0:
                trd["spikeTimes"] = np.array([np.nan])
                sf = 1.0
            elif np.nanmax(trd["spikeTimes"]) > 2.0:  # probably in miec...
                sf = 1e-3  # so scale to seconds
            else:
                sf = 1.0
            sptimes = np.array(trd["spikeTimes"]) * sf  # convert to seconds

            stim_dt = stb[1] - stb[0]
            idxs = (int(data_window[0] / stim_dt), int(data_window[1] / stim_dt))
            if not isinstance(trd["spikeTimes"], list) and not isinstance(
                trd["spikeTimes"], np.ndarray
            ):
                cprint("r", "spiketimes is not list")
                cprint("r", f"    {type(trd['spikeTimes'])=}")
                return
            all_bu_st.extend(sptimes)
            all_bu_st_trials.append(sptimes)

            n_inputs = len(trd["inputSpikeTimes"])
            if i_trial == 0:
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
                an_st_grand[i_trial].extend(tk)

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

        # prepare spike arrays
        x = np.array(all_bu_st)
        v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
        # window_duration = phasewin[1]-phasewin[0]
        bu_spikesinwin = x[v]

        x = np.array(all_an_st)
        v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
        an_spikesinwin = x[v]
        cprint("r", len(all_bu_st_trials))

        if soundtype in [
            "SAM",
            "stationaryNoise",
            "regularClicks",
            "poissonClicks",
        ]:  # also calculate vs and plot histogram
            print(soundtype, fmod)

            # break into 10 groups of 10 trials (if 100 total) and
            # compute vs for each group, then get mean and sd to report
            n_vs_groups = int(len(all_bu_st_trials) / 10)
            cprint("c", f"n vs groups:  {n_vs_groups:d}")
            cprint("c", f"len bust trials: {len(all_bu_st_trials):d}")
            vs_n = np.zeros(n_vs_groups)
            vs_nspikes = np.zeros(n_vs_groups)
            bu_vs_data = []
            for i in range(n_vs_groups):  # divide trials into trial groups (of 10)
                bu_vs_data = []
                for j in range(i * n_vs_groups, (i + 1) * n_vs_groups):  # within each
                    v = np.where(
                        (phasewin[0] <= all_bu_st_trials[j])
                        & (all_bu_st_trials[j] < phasewin[1])
                    )[0]
                    bu_vs_data.extend(all_bu_st_trials[j][v])
                vs_calc = self.VS.vector_strength(
                    bu_vs_data, fmod, time_window=phasewin, nreps=nreps
                )  # vs expects spikes in msec
                vs_n[i] = vs_calc.vs
                vs_nspikes[i] = int(vs_calc.n_spikes)

                if make_VS_raw:
                    VS_file_raw = Path(f"VS_raw_SAM_{cellN:02d}_{int(dB):02d}")
                    if (
                        not VS_file_raw.is_file()
                    ):  # if file does not exist, write header
                        with open(VS_file_raw, "w") as fh:
                            fh.write(f"celln,experiment,trial,fmod,VS,nspikes\n")
                    with open(VS_file_raw, "a") as fh:  # add the data to this file
                        fh.write(
                            f"{cellN:02d},{ri.Spirou:s},{i:d},{int(fmod):d},{vs_calc.vs:.5f},{vs_calc.n_spikes:d}\n"
                        )

            print("vs_nspikes: ", vs_nspikes)
            print("vs by group: ", vs_n)
            cprint(
                "c",
                f"mean VS: {np.mean(vs_n):8.4f}  SD={np.std(vs_n):8.4f}, total spikes: {int(np.sum(vs_nspikes)):d}",
            )

            vs = self.VS.vector_strength(
                bu_spikesinwin,
                fmod,
                time_window=phasewin,
                nreps=nreps,
            )  # vs expects spikes in msec
            cprint("c", f"Grand mean VS: {vs.vs:8.4f}  N={vs.n_spikes:6d}")
            vs.vs_mean = np.mean(vs_n)
            vs.vs_sd = np.std(vs_n)
            vs.vs_Ns = np.sum(vs_nspikes)
            vs.vs_groups = n_vs_groups

            vs.carrier_frequency = F0
            vs.n_inputs = n_inputs
            if sac_flag and len(all_bu_st) > 50:
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
                print("fmod, 10/ri.fmod: ", ri.fmod, 10.0 / ri.fmod)
                vs.sac_bu_CI, peaks, HW, FW = S.SAC_measures(
                    bu_sac, bu_sacbins, twin=spars.twin
                )
                if HW is not None:
                    vs.sac_bu_hw = HW[0][0] * spars.binw
                    left_ips = HW[2][0]
                    right_ips = HW[3][0]
                else:
                    vs.sac_bu_hw = np.nan
                    left_ips = np.nan
                    right_ips = np.nan
                print("BU SAC Report: \n  ")
                print(
                    f"    HW:    {vs.sac_bu_hw:.6f}  CI: {vs.sac_bu_CI:.2f}  Left IPS: {left_ips:.2f}  Right IPS: {right_ips:.2f}, Binw: {spars.binw:.6f}"
                )
            else:
                bu_sac = None
                bu_sacbins = None
                vs.sac_bu_CI = 0.0
                vs.sac_bu_hw = np.nan
            print(
                "CN Vector Strength at %.1f: %7.3f, rMTF=%.2f, sp/s entrainment:%.4f,  dispersion=%.2f (us) Rayleigh: %7.3f  p = %.3e  n_spikes = %d"
                % (
                    fmod,
                    vs.vs,
                    vs.rMTF,
                    vs.entrainment,
                    vs.circ_timeSD * 1e6,
                    vs.Rayleigh,
                    vs.pRayleigh,
                    vs.n_spikes,
                )
            )
            vs_an = self.VS.vector_strength(
                an_spikesinwin,
                fmod,
                time_window=phasewin,
                nreps=nreps,
            )
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
                    (
                        an_CI_windowed,
                        peaks,
                        an_hw[i_an],
                        an_fw[i_an],
                    ) = S.SAC_measures(an_sac[i_an], an_sacbins, twin=spars.twin)
                an_hw = np.array(an_hw).T[0]  # reorganize the data
                vs.sac_an_hw = np.mean(an_hw[0, :]) * spars.binw
                CI = np.mean(an_hw[1, :])
                vs.sac_an_CI = CI
                lips = np.mean(an_hw[2, :])
                rips = np.mean(an_hw[3, :])
                print("AN SAC Report: ")
                print(
                    f"    HW:    {vs.sac_an_hw:.6f}  CI: {vs.sac_an_CI:.2f}  Left IPS: {lips:.2f}  Right IPS: {rips:.2f}"
                )
            else:
                an_sac = None
                an_sacbins = None
                vs.sac_an_CI = 0.0
                vs.sac_an_hw = 0.0

            print(" Sound type: ", soundtype)
            print(
                "AN Vector Strength at %.1f: %7.3f, rMTF=%.2f sp/s entrainment:%.4f  dispersion=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
                % (
                    fmod,
                    vs_an.vs,
                    vs_an.rMTF,
                    vs_an.entrainment,
                    vs_an.circ_timeSD * 1e6,
                    vs_an.Rayleigh,
                    vs_an.pRayleigh,
                    vs_an.n_spikes,
                )
            )
            if sac_flag:
                print(
                    f"BU SAC: CI = {vs.sac_bu_CI:7.3f}, min = {np.min(bu_sac):7.3f}, hw = {1e3*vs.sac_bu_hw:7.3f} ms"
                )
                print(
                    f"AN SAC: CI = {vs.sac_an_CI:7.3f}, min = {np.min(an_sac):7.3f}, hw = {1e3*vs.sac_an_hw:7.3f} ms"
                )
            vs.an_vs = vs_an.vs  # copy
            vs.an_rMTF = vs_an.rMTF
            vs.an_entrainment = vs_an.entrainment
            vs.an_circ_phaseMean = vs_an.circ_phaseMean
            vs.an_circ_phaseSD = vs_an.circ_phaseSD
            vs.an_circ_timeSD = vs_an.circ_timeSD

            colnames, line = self.print_VS(
                cell_ID, vs, fmod, dmod, dB, ri.Spirou, filename=str(filename)
            )

        #########################

        if not plotflag or P is None:
            return colnames, line

        # Otherwise, lets do some plotting:
        # plot the raster of spikes
        PF.plot_spike_raster(
            P,
            mode="postsynaptic",
            spike_times=all_bu_st_trials,
            data_window=data_window,
            panel=bu_raster_panel,
        )
        PF.plot_stacked_spiketrain_rasters(
            spike_times_by_input=an_st_by_input,
            ax=P.axdict[an_raster_panel],
            si=si,
            syninfo=syninfo,
            plot_win=data_window,
            max_trials=5,
        )
        # PF.plot_spike_raster(P, mode='presynaptic', n_inputs=1, spike_times=trd["inputSpikeTimes"],
        #     data_window=data_window, panel=an_raster_panel)

        # be sure we are using the right data, not overwriting
        time_base = AR.MC.time_base / 1000.0  # convert to seconds
        dt = si.dtIC / 1000.0  # convert from msec to seconds
        i_trial = 0
        P.axdict[bu_voltage_panel].plot(
            np.array(time_base[idx[0] : idx[1]] - data_window[0]),
            np.array(v_show_trial[idx[0] : idx[1]]),
            "k-",
            linewidth=0.5,
        )
        # get the spike indices for the trial that will be shown at the top
        spikeindex = [
            int(t / dt)
            for t in all_bu_st_trials[show_vtrial]
            if not np.isnan(t) and t >= data_window[0] and t < data_window[1]
        ]
        P.axdict[bu_voltage_panel].plot(
            np.array(time_base[spikeindex] - data_window[0]),
            v_show_trial[spikeindex],
            "ro",
            markersize=1.5,
        )
        P.axdict[bu_voltage_panel].set_xlim(0.0, plot_dur)

        dw = 0.5 * np.max(waveform[i_trial])
        P.axdict[stimulus_panel].plot(
            stb[idxs[0] : idxs[1]] - data_window[0],
            np.array(waveform[i_trial][idxs[0] : idxs[1]]) + i_trial * dw,
            "-",
            linewidth=0.5,
        )

        if soundtype in [
            "tonepip",
            "SAM",
            "noise",
            "stationaryNoise",
            "regularClicks",
            "poissonClicks",
        ]:  # use panel F for FSL/SSL distributions
            # the histograms of the data
            import vcnmodel.util.make_sound_waveforms as msw

            MSW = msw.Make_Sound_Waveform()

            psth_binw = 0.2e-3
            isi_binw = 0.1e-3
            PF.plot_psth(
                all_bu_st_trials,
                run_info=ri,
                zero_time=data_window[0],
                max_time=data_window[1],
                bin_width=psth_binw,
                ax=P.axdict[bu_psth_panel],
                stimbar={"sound_params": None, "waveform": [stb, waveform[0]]},
            )
            PF.plot_isi(
                all_bu_st_trials,
                run_info=ri,
                zero_time=data_window[0],
                max_time=data_window[1],
                bin_width=isi_binw,
                ax=P.axdict[bu_isi_panel],
            )
            P.axdict[bu_isi_panel].set_xlim(0, 5.0 / ri.fmod)
            PF.plot_psth(
                # all_an_st_grand,
                an_st_by_input[0],
                run_info=ri,
                zero_time=data_window[0],
                max_time=data_window[1],
                bin_width=psth_binw,
                ax=P.axdict[an_psth_panel],
            )
            PF.plot_isi(
                # all_an_st_grand,
                an_st_by_input[0],
                run_info=ri,
                zero_time=data_window[0],
                max_time=data_window[1],
                bin_width=psth_binw,
                ax=P.axdict[an_isi_panel],
            )
            P.axdict[an_isi_panel].set_xlim(0, 5.0 / ri.fmod)
            # just print the AN rates
            print("an st by input: ", len(an_st_by_input))
            PF.print_AN_rates(
                an_st_by_input,
                run_info=ri,
            )

            PF.plot_fsl_ssl(
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
        ]:  # also calculate vs and plot phase locking, CI, etc. histogram
            if len(all_bu_st) == 0:
                P.axdict[vs_panel].text(
                    0.5,
                    0.5,
                    "No Spikes",
                    fontsize=14,
                    color="r",
                    transform=P.axdict[vs_panel].transAxes,
                    horizontalalignment="center",
                )
            sim_dir = fn.parts[-2]
            tstring = f"{gbc:s} {str(sim_dir):s}\nSAM Tone: F0={F0:.3f} at {dB:3.1f} dbSPL, fMod={fmod:3.1f}  dMod={dmod:5.2f}, cell CF={F0:.3f}"
            # print('si: ', si)
            tstring += f"\nExpt: {ri.Spirou:s}  DendMode: {si.dendriteMode:s} DendExpt: {si.dendriteExpt:s}"
            tstring += f" Depr: {si.ANSynapticDepression:d} ANInput: {si.SRType:s}"
            P.figure_handle.suptitle(tstring, fontsize=8)

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
            PF.plot_psth(
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
                            an_sacbins[:-1], an_sac[i_an], "r-", label="AN"
                        )
                    else:
                        P.axdict[sac_panel].plot(
                            an_sacbins[:-1],
                            an_sac[i_an],
                            "r-",
                        )

                P.axdict[sac_panel].plot(bu_sacbins[:-1], bu_sac, "b-", label="BU")
                if soundtype in ["SAM"]:
                    P.axdict[sac_panel].set_xlim((-2.0 / ri.fmod, 2.0 / ri.fmod))
                elif soundtype.endswith("Clicks"):
                    P.axdict[sac_panel].set_xlim(
                        (-spars.displayDuration, spars.displayDuration)
                    )
                else:
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
