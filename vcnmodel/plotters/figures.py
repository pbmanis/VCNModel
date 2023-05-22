"""figures.py - generate figures for the paper
This class generates final figures for the SBEM manuscript by reaching back to
the original simulation data, including, in some cases, refitting. Both primary
"examplar" figures, and supplemental figures, are generated. The figures are
made consistent by using both sns.set_style and mpl.style for figures.mplstyle,
which overrides some defaults in mpl. The resulting figures are editable in
Illustrator without any missing fonts.

This is just "plotting code". No apologies for the length of the file.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import datetime
import multiprocessing as MPROC
import pickle
import string
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from ssl import SSL_ERROR_EOF
from typing import List, Union

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import image as mpimg
from matplotlib.lines import Line2D
from pylibrary.plotting import plothelpers as PH
from pylibrary.tools import cprint as CP
from pyqtgraph import multiprocess as MP

import vcnmodel.util.fixpicklemodule as FPM
import vcnmodel.util.readmodel as readmodel
from vcnmodel import group_defs as GRPDEF
from vcnmodel.analyzers import analyze_data
from vcnmodel.analyzers import isi_cv as ISI
from vcnmodel.analyzers import pattern_summary as PATSUM
from vcnmodel.analyzers import sac as SAC
from vcnmodel.plotters import SAC_plots as SACP
from vcnmodel.plotters import SAM_VS_vplots
from vcnmodel.plotters import efficacy_plot as EF
from vcnmodel.plotters import \
    figure_data as FD  # table of simulation runs used for plotting figures
from vcnmodel.plotters import morphology_thr_correlations
from vcnmodel.plotters import plot_functions as PF
from vcnmodel.plotters import plot_z as PZ
from vcnmodel.util.get_data_paths import get_data_paths
from vcnmodel.util.set_figure_path import set_figure_path

cprint = CP.cprint


syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x"]
colors = ["c", "k", "m", "r"]
# seaborn default palette, first 10 colors
sns_colors = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
]
title_text = {
    "passive": "Passive",
    "normal": "Half-active",
    "active": "Active",
    "pasdend": "Passive",
    "actdend": "Active",
}
font_colors = {"passive": "c", "normal": "k", "active": "m"}
spike_colors = {"passive": "r", "normal": "r", "active": "r"}


def get_changetimestamp():
    # trip filemode based on date of simulation
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate, "%Y-%m-%d-%H:%M")
    changetimestamp = datetime.datetime.timestamp(dts)
    return changetimestamp


@dataclass
class PData:
    """
    data class for some parameters that control what we read
    """

    gradeA: list = field(default_factory=GRPDEF.gradeACells)
    default_modelName: str = "XM13_nacncoop"
    soma_inflate: bool = True
    dend_inflate: bool = True
    basepath: str = ""  # config["baseDataDirectory"]
    renderpath: str = ""  # str(Path(config["codeDirectory"], "Renderings"))
    revcorrpath: str = ""  # config["revcorrDataDirectory"]
    thiscell: str = ""


def title_data():
    return {"title": "", "x": None, "y": None}


@dataclass
class FigInfo:
    """
    data class to hold figure information for savefig
    """

    P: object = None
    show_figure_name: bool = False  # if True, figure filename appears in top right corner
    filename: Union[str, Path] = ""
    title: dict = field(default_factory=title_data)
    title2: dict = field(default_factory=title_data)


class Figures(object):
    """
    This class generates final figures for the SBEM manuscript by reaching back
    to the original simulation data, including, in some cases, refitting. Both
    primary "exemplar" figures, and supplemental figures, are generated. The
    figures are made consistent by using both sns.set_style and mpl.style for
    figures.mplstyle, which overrides some defaults in mpl. The resulting
    figures are editable in Illustrator without any missing fonts.

    This is overall ugly code, but it gets the job done. Note that some plots
    are modifications of what is present in plot_sims, and also that plot_sims
    is used here, accessed through self.parent.PLT. Some plotting routines have
    been moved into plot_functions.py.

    This is part of the DataTablesVCN interactive analysis tool for the
    simulations in the SBEM project.

    """

    def __init__(self, parent):
        self.parent = parent  # point back to caller's space
        self.config = get_data_paths()
        self.axis_offset = -0.02
        self.ReadModel = readmodel.ReadModel()
        self.ReadModel.set_parent(
            parent=self.parent, my_parent=self
        )  # pass our parent and US to the reader

    def newPData(self):
        """
        Return Pdata with the paths set from self.config
        and the GradeA Cell list.

        Returns:
            PData: dataclass
        """
        return PData(
            gradeA=GRPDEF.gradeACells,
            basepath=str(Path(self.config["disk"], self.config["baseDataDirectory"])),
            renderpath=str(Path(self.config["codeDirectory"], "Renderings")),
            revcorrpath=self.config["revcorrDataDirectory"],
        )

    def reset_style(self):
        sns.set_style(rc={"pdf.fonttype": 42})
        mpl.style.use("~/.matplotlib/figures.mplstyle")

    def make_figure(self, figure_name: Union[str, None] = None):
        self.reset_style()
        # dispatch
        dispatch_table = {
            "Figure4-Ephys_1_Main": self.Figure4_Main,
            "Figure4-Supplemental2_VC": self.Figure4_Supplemental2_VC,
            "Figure4-Supplemental3_CC": self.Figure4_Supplemental3_CC,
            "Figure4-Supplemental4_Zin_removed": self.Figure4_Supplemental4_Zin_removed,
            "Figure4-Supplemental4_PSTH": self.Figure4_Supplemental4_PSTH,
            "Figure5-Ephys_2_Main": self.Figure5_Main,
            
            "Figure5-Ephys_2_Supplemental1": self.Figure5_Supplemental1,
            "Figure5-Ephys_2_Supplemental2": self.Figure5_Supplemental2_removed,
            "Figure5-Ephys_2_Supplemental2": self.Figure5_Supplemental2,
            
            "Figure6-Ephys_3_Main": self.Figure6_Main,
            "Figure6-Ephys_3_Supplemental2 (VS, rMTF)": self.Figure6_Supplemental2,
            "Figure6-Ephys_3_Supplemental3 (Entrainment)": self.Figure6_Supplemental3,
            "Figure6-Ephys_3_Supplemental4 (SAC)": self.Figure6_Supplemental4,
            
            "Figure8-Ephys_4": self.Figure8_Panels_IJKLM,
            
            # Misc figures follow
            "Figure: IV Figure": self.plotIV,
            "Figure: All_IVs": self.allIVs,
            "Figure: CombinedEffRevCorr": self.plot_combined_Eff_Rev,
            "Figure: Efficacy": self.plot_efficacy,
            "Figure: Efficacy Supplement": self.plot_efficacy_supplement,
            "Figure: Revcorr Example": self.plot_revcorr,
            "Figure: All Revcors": self.plot_all_revcorr,
            "Figure: Revcorr at Spont": self.plot_revcorr_supplement_spont,
            "Figure: Revcorr at 30dB": self.plot_revcorr_supplement_30dB,
            "Figure: Revcorr at 40dB": self.plot_revcorr_supplement_40dB,
            "Figure: Compare Revcorrs": self.plot_revcorr_compare,
            "Figure: PSTHs": self.plot_PSTH,
            "Analyze VS-SAM table @ 15dBSPL": self.plot_VS_SAM_15,
            "Analyze VS-SAM table @ 30dBSPL": self.plot_VS_SAM_30,
            "Analyze VS-SAM BC09 table @ 15dBSPL": self.plot_VS_SAM_15_BC09,
        }
        # print(figure_name)
        # print(dispatch_table.keys())
        if figure_name in list(dispatch_table.keys()):
            fig = dispatch_table[figure_name]()
            if fig is not None:
                mpl.show()
                self.save_figure(fig)
        else:
            cprint("r", f"Figure name '{figure_name:s}' was not in dispatch table.")

    def save_figure(self, fig, show_figure_name:bool=False):
        """
        Save a figure to a disk file.
        This routine adds metadata to the figure, along with
        a title if specified.

        Parameters
        ----------
        fig : figinfo dataclass object. This should be created and
            populated by the calling routine.

        Returns
        -------
        Nothing
        """

        if show_figure_name:
            fig.P.figure_handle.text(
                0.98,
                1.0,
                str(fig.filename.name),  # .replace('_', '\_'),
                transform=fig.P.figure_handle.transFigure,
                fontdict = {"fontsize": 7, "fontweight": "normal",
                    "horizontalalignment": "right",
                    "verticalalignment": "top",
                }
            )

            if hasattr(
                fig, "title2"
            ):  # ggplot figures do not have a title or title2 attribute
                if fig.title2["title"] is not None or len(fig.title2["title"]) > 0:
                    fig.P.figure_handle.text(
                        fig.title2["x"],
                        fig.title2["y"],
                        fig.title2["title"],  # .replace('_', '\_'),
                        transform=fig.P.figure_handle.transFigure,
                        horizontalalignment="right",
                        verticalalignment="top",
                )
        if isinstance(fig.filename, str):
            fig.filename = Path(fig.filename)
        if not fig.filename.is_absolute():  # we were given a filename without the full path
            out_file = Path(self.config["baseDataDirectory"], # make a full path
                            self.config["figureDirectory"],
                            fig.filename)
        else:
            out_file = fig.filename
        Path(out_file).parent.mkdir(parents=True, exist_ok=True)
        cprint(
            "g",
            f"Saving figure to: {str(out_file):s}",
        )
        mpl.savefig(
            out_file,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": fig.title["title"],
            },
        )
        # make png file as well
        fn2 = Path(out_file).with_suffix(".png")
        mpl.savefig(
            Path(self.config["baseDataDirectory"], "Figures", fn2),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": fig.title["title"],
            },
        )
        fig.P.figure_handle.show()

    def force_log_ticks(self, ax):
        """
        Set up log spacing ticks for the specified axis
        """

        locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0,), numticks=100)
        ax.xaxis.set_major_locator(locmaj)

        locmin = matplotlib.ticker.LogLocator(
            base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100
        )
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    def allIVs(self):
        """
        Plot an IV for every cell in the grade A cell list
        """
        for cell in GRPDEF.grAList():
            self.plotIV(cell)

    def get_dendmode(self, dendmode: str = ""):
        """
        Provide a remapping for the dendrite mode names
        """
        if dendmode == "passive":
            dendm = "pasdend"
        elif dendmode == "active":
            dendm = "actdend"
        elif dendmode == "normal":
            dendm = "normal"
        else:
            raise ValueError(f"figures.py: dendmode not recognized: {dendmode:s}")
        return dendm

    def plotIV(
        self,
        cell=None,
        parent_figure=None,
        loc: Union[None, tuple] = None,
        toponly: bool = False,
        show_title: bool = False,
    ):
        """
        Plot the traces and current-voltage relationship for the specified cell.

        Parameters
        ----------
        cell : int (default None)
            Specify the cell number to plot.
            If no cell number is specified, this will plot the "default" cell specified
            in figure_data.py, figure_IV dictionary.
        parent_figure : figure object.
            Sets up the arrangement/order of panel labels based on whether this is being plotted
            as a "sub" figure of a parent figure or not.
        loc: location of the "sub" figure relative to the parent figure. ranges 0-1 in each axis

        toponly: just for replotting the raw traces.
        parts: str
            either "all" (respecting top only) or "DE"

        """
        if cell is None:
            cellN = list(FD.figure_IV.keys())[0]  # there will be only one
            d1 = FD.figure_IV[cellN]["normal"]
        else:
            cellN = cell
            d1 = FD.figure_AllIVs[cellN]["normal"]
        cellpath = Path(
            self.config["cellDataDirectory"], f"VCN_c{cellN:02d}", "Simulations", "IV"
        )
        rows = 1  # basic layout for figure, but applies only to the top
        cols = 3
        height = 4.0
        width = 8.0
        if parent_figure is not None:
            fsize = parent_figure.figsize
            p0panels = ["D", "E", "F"]
            p2panels = ["G", "H", "I", "J"]
            p3panels = ["H1", "H2"]
            phaseplot1 = "H1"
            phaseplot2 = "H2"
            rinplot = "I"
            tauplot = "J"
        else:
            fsize = (width, height)
            p0panels = ["A", "B", "C"]
            p2panels = ["D", "E", "F", "G"]
            p3panels = ["E1", "E2"]
            phaseplot1 = "E1"
            phaseplot2 = "E2"
            rinplot = "F"
            tauplot = "G"
        if loc is not None:
            bmar = loc[2]
            tmar = fsize[1] - loc[3]
        else:
            bmar = 0.0
            tmar = 0.0

        PD = self.newPData()
        ymin = -125.0
        ymax = 20.0
        # We generate subplot setups here by defining
        # multiple grids in the figure.
        self.P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            units="in",
            figsize=fsize,
            showgrid=False,
            verticalspacing=0.25,
            horizontalspacing=0.35,
            margins={
                "bottommargin": 2.00 + bmar,
                "leftmargin": 0.5,
                "rightmargin": 0.5,
                "topmargin": 0.25 + tmar,
            },
            labelposition=(-0.05, 1.05),
            parent_figure=parent_figure,
            panel_labels=p0panels,
        )

        self.P2 = PH.regular_grid(  # lower row of 4
            1,
            4,
            order="rowsfirst",
            units="in",
            showgrid=False,
            figsize=fsize,
            verticalspacing=0.25,
            horizontalspacing=0.15,
            margins={
                "bottommargin": 0.5 + bmar,
                "leftmargin": 0.2,
                "rightmargin": 0.2,
                "topmargin": 1.8 + tmar,
            },
            labelposition=(-0.05, 1.05),
            parent_figure=self.P,
            panel_labels=p2panels,
        )

        self.P3 = PH.regular_grid(  # replace with stacked plots in E
            2,
            1,
            order="rowsfirst",
            units="in",
            showgrid=False,
            figsize=fsize,
            verticalspacing=0.1,
            horizontalspacing=0.35,
            margins={
                "bottommargin": 0.5 + bmar,
                "leftmargin": 2.35,
                "rightmargin": 4.15,
                "topmargin": 1.8 + tmar,
            },
            labelposition=(-0.25, 1.05),
            parent_figure=self.P,
            panel_labels=p3panels,
        )
        ivaxis = self.P2.axarr[0, 0]
        for iax, iv in enumerate(["pasdend", "normal", "actdend"]):
            dfile = FD.figure_AllIVs[cellN][iv]
            sfi = Path(cellpath, Path(dfile).name)
            if not sfi.is_dir():
                print("Missing file: sfi 1: ", sfi)
                return None
            fn = list(sfi.glob("*"))
            if len(fn) == 0:
                print("sfi: ", sfi)
                raise ValueError("no files found")
            sfi = Path(sfi, fn[0])
            self.parent.PLT.plot_traces(
                self.P.axarr[0, iax],
                sfi,
                PD,
                protocol="IV",
                ymin=ymin,
                ymax=ymax,
                xmax=150.0,
                iax=iax,
                figure=self.P.figure_handle,
                ivaxis=ivaxis,  # accumulate IV's in bottom left plot
                ivcolor=colors[iax],
                calx=120.0,
                caly=-10.0,
                show_title=False,
            )
            # self.P.axarr[0, iax].set_title(
        #      title_text[iv], color="k", fontweight="normal"
        #  )

        ls = ["-", "-", "-", ""]
        mc = ["w", "w", "w", "r"]
        lines = [
            Line2D(
                [0],
                [0],
                color=c,
                linewidth=0.75,
                linestyle=ls[ic],
                marker="o",
                mfc=mc[ic],
                mew=0.5,
            )
            for ic, c in enumerate(colors)
        ]
        # self.P2.axarr[0, 0].line.set_label(iv)
        # self.P2.axarr[0,0].legend(["passive", "normal", "active", "spike thresholds"])

        labels = ["Passive", "Half-active", "Active", "Spk Thresh"]
        self.P2.axarr[0, 0].legend(
            lines, labels, bbox_to_anchor=(0.32, 0.35), fontsize=7
        )
        res_label = r"$\mMTCegular{R_{in} (M\Omega)}$"
        phase_label = r"$\mMTCegular{\phi (radians)}$"

        # plot overlays of all cell z/phase
        for iax, mode in enumerate(["Z_passive", "Z_normal", "Z_active"]):
            if mode not in FD.figure_IV.keys():
                continue
            sfi = Path(
                self.config["cellDataDirectory"],
                self.config["impedanceDirectory"],
                FD.figure_IV[mode],
            )
            if not sfi.is_file():
                cprint("r", f"File not found!!!!!!\n->>> {str(sfi):s}")
                return None
            ax = self.P3.axdict[phaseplot1]
            label = sfi.name  # .replace("_", "\_")
            with open(sfi, "rb") as fh:
                d = FPM.pickle_load(fh)
            pz = ax.plot(
                d["f"],
                d["zin"],
                color=colors[iax],
                marker=syms[iax],
                markersize=3,
                label=label[:8],
            )
            ax.set_ylim(0, 60.0)
            ax.set_xlim(1.0, 10000.0)
            ax.set_ylabel(res_label)
            PH.noaxeslabels(ax, whichaxes="x")
            # ax.set_xlabel("Frequency (Hz)", fontsize=7)
            ax.set_xscale("log")
            self.force_log_ticks(ax)

            # secax = PLS.create_inset_axes([0.0, 0.0, 1, 0.5], ax, label=str(ax))
            secax = self.P3.axdict[phaseplot2]
            secax.set_xscale("log")
            secax.set_xlim(1.0, 10000.0)
            self.force_log_ticks(secax)
            pp = secax.plot(
                d["f"],
                d["phase"],
                color=colors[iax],
                marker=syms[iax],
                markersize=3,
                # label=filename
            )
            secax.set_ylim(-0.6, 0.2)
            secax.set_ylabel(phase_label)
            secax.set_xlabel("Frequency (Hz)")
            PH.nice_plot(secax, position=self.axis_offset, direction="outward", ticklength=3)

        # PH.show_figure_grid(self.P.figure_handle)

        # get Rin and RM from all the examples and make a summary distribution plot
        rins, taus = self.get_Rin_Tau()

        self.plot_Rin_Tau(rins=rins, taus=taus, 
            ax_rin=self.P2.axdict[rinplot], ax_tau=self.P2.axdict[tauplot])
        
        fig = FigInfo()
        fig.P = self.P
        if cell is None:
            fig.filename = Path(
                self.config["figureIntermediateDirectory"], "Fig_M1.pdf"
            )
        else:
            fig.filename = Path(
                self.config["figureIntermediateDirectory"],
                f"Fig_IV/IV_cell_VCN_c{cellN:02d}.pdf",
            )
        fig.title["title"] = "SBEM Project Figure 1 Modeling (Main)"
        # self.save_figure(self.P, save_file, title)
        return fig
    
    def get_Rin_Tau(self):
                # get Rin and RM from all the examples and make a summary distribution plot
        rins = {}
        taus = {}
        PD = self.newPData()
        k = 0
        for rax, iv in enumerate(FD.figure_AllIVs.keys()):
            # cellpath = Path(self.config['cellDataDirectory'], f"VCN_c{iv:02d}", 'Simulations', 'IV')
            for iax, dendmode in enumerate(["passive", "normal", "active"]):
                dendm = self.get_dendmode(dendmode)

                sfi = Path(
                    self.config["cellDataDirectory"],
                    f"VCN_c{iv:02d}",
                    "Simulations",
                    "IV",
                    FD.figure_AllIVs[iv][dendm],
                )
                if not sfi.is_dir():
                    print(f"File '{str(sfi):s}' is not a directory")
                    continue
                fn = list(sfi.glob("*"))
                sfi = Path(sfi, fn[0])
                AR, SP, RM = self.parent.PLT.plot_traces(
                    None,  # self.P.axarr[rax, iax+1],
                    sfi,
                    PD,
                    protocol="IV",
                    ymin=-100,
                    ymax=20,
                    iax=iax,
                    figure=None,  # just returns values when figure is None.
                )
                rins[k] = {
                    "Cell": iv,
                    "Rin": RM.analysis_summary["Rin"],
                    "dendrites": dendmode,
                }
                taus[k] = {
                    "Cell": iv,
                    "taum": RM.analysis_summary["taum"],
                    "dendrites": dendmode,
                }
                k += 1
        return rins, taus

    def plot_Rin_Tau(self, rins:dict, taus:dict, ax_tau=None, ax_rin=None):

        df_rin = pd.DataFrame.from_dict(rins, orient="index")  # , orient='index')
        df_tau = pd.DataFrame.from_dict(taus, orient="index")
        res_label = r"$\mathrm{R_{in} (M\Omega)}$"
        tau_label = r"$\mathrm{\tau_{m} (ms)}$"
        sns.boxplot(
            data=df_rin,
            x="dendrites",
            y="Rin",
            ax=ax_rin,
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_rin,
            x="dendrites",
            y="Rin",
            ax=ax_rin,
            color="0.6",
            size=4.0,
            edgecolor="k",
        )
        ax_rin.set_ylim(0, 30.0)
        ax_rin.set_ylabel(res_label)
        ax_rin.set_xlabel("Dendrite Decoration")
        ax_rin.set_xticklabels(["Passive", "Half-active", "Active"])
        PH.nice_plot(ax_rin, position=-0.03, direction="outward", ticklength=3)
        sns.boxplot(
            data=df_tau,
            x="dendrites",
            y="taum",
            ax=ax_tau,
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_tau,
            x="dendrites",
            y="taum",
            ax=ax_tau,
            color="0.6",
            size=4.0,
            edgecolor="k",
        )
        ax_tau.set_ylim(0, 2.5)
        ax_tau.set_ylabel(tau_label)
        ax_tau.set_xticklabels(["Passive", "Half-active", "Active"])
        PH.nice_plot(ax_tau, position=-0.03, direction="outward", ticklength=3)

        ax_tau.set_xlabel("Dendrite Decoration")


    def Figure4_Main(self):
        message = """
        This Figure was made by using Illustrator to combine parts of other figures/files as follows:
        Panel A came from a syglass rendition of the obj (mesh).
        Panel B (and the inset in C) came from a vispy rendition (neuronvis) of the HOC file.
        Panel C is a combination of Illustrator text, and cnmodel output (see the notebook folder,
        nb/Figure3_main_PanelC.ipynb, for the generation code; the matplotlib windows were then pulled
        into Illustrator to make the figure)
        
        Panels D and E are taken from Figure4_Supplemental2_CC.pdf, for BC17.
        Panels F, G and H are from Figure4_Supplemental5_PSTH.pdf for BC17, with the stimuli underneath

        """
        cprint("y", message)

    def Figure4_Supplemental2_VC(self):
        """
        Figure 3, Supplemental Figure 2
        Combined voltage clamp traces,  Rin, taum plots
        """
        P0 = PH.regular_grid(  # dummy figure space
            1,
            1,
            figsize=(8.0, 7.0),
            order="rowsfirst",
            units="in",
            showgrid=True,
            parent_figure=None,
        )

        figp1 = self.plot_VC_gKLT(parent_figure=P0, loc=(0, 8, 2., 5.))
        P2 = PH.regular_grid(  # lower row of 2
            1,
            2,
            order="rowsfirst",
            units="in",
            showgrid=False,
            figsize=(8.0, 6.0),
            verticalspacing=0.25,
            horizontalspacing=0.5,
            margins={
                "bottommargin": 0.4,
                "leftmargin": 2.0,
                "rightmargin": 2.0,
                "topmargin": 4.6,
            },
            labelposition=(-0.05, 1.05),
            parent_figure=P0,
            panel_labels=["D", "E"],
        )
        rins, taus = self.get_Rin_Tau()
        self.plot_Rin_Tau(rins, taus, P2.axdict["D"], P2.axdict["E"])
        # P2 = self.plotIV(parent_figure=figp1.P, loc=(0, 8, 0.0, 4.0))
        fig = FigInfo()
        fig.P = P0
        fig.filename = set_figure_path(fignum=4, filedescriptor="VC_Rin_Taum_V2",
            suppnum=2)
        fig.title["title"] = "SBEM Project Figure 4 Supplemental Figure 2 VC_Rin_Taum"
        return fig

    def Figure4_Supplemental3_CC(self, parent_figure=None, show_pngs=False):
        """
        Plot all of the IVS, for a supplemental figure
        Passive, normal, active, plus the crossed IV
        Also put the PNG for the cell on the left.
        """
        nivs = len(FD.figure_AllIVs)
        cprint("c", "Plotting Figure4_Supplemental3_CC")
        rows = nivs
        if show_pngs:
            cols = 5
        else:
            cols = 4
        height = 1.5 * nivs
        width = 8.5
        PD = self.newPData()
        ymin = -125.0
        ymax = 40.0
        calx = 120.0
        caly = -110.0
        self.P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            figsize=(width, height),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.05,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.08,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=parent_figure,
            # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        cellpath = self.config["cellDataDirectory"]
        png_path = Path(self.config["baseDataDirectory"], self.config["pngDirectory"])
        column_titles = {'Passive': 'Passive', 'Normal': "Half-active", 'Active': "Active"}

        for rax, iv in enumerate(FD.figure_AllIVs.keys()):
            # if iv not in [9, 10]:
            #      continue
            cprint(
                "c", f"    Doing Cell {FD.BC_name:s} {iv:02d} -----------------------------------"
            )
            celln = Path(png_path, f"VCN_c{iv:02d}.png")
            if celln.is_file() and show_pngs:  # add images from png files
                img = mpimg.imread(str(celln))
                self.P.axarr[rax, 0].imshow(img, aspect="equal")
                ylim = self.P.axarr[rax, 0].get_ylim()
                if iv != 10:
                    self.P.axarr[rax, 0].set_xlim(900, 1500)
                else:
                    self.P.axarr[rax, 0].set_xlim(-1200, 1500)
                PH.noaxes(self.P.axarr[rax, 0])
            # plot 3 dendrite decorations

            for iax, dendmode in enumerate(["passive", "normal", "active"]):
                dendm = self.get_dendmode(dendmode)
                sfi = Path(
                    cellpath,
                    f"VCN_c{iv:02d}",
                    "Simulations",
                    "IV",
                    FD.figure_AllIVs[iv][dendm],
                )
                if not sfi.is_dir():
                    cprint("r", f"Unable to find dir: {str(sfi):s}")
                    continue
                fn = list(sfi.glob("*"))
                sfi = Path(sfi, fn[0])
                if rax > 0 or iax > 0:
                    calx = None  # only one cal bar on this plot, top row.

                if show_pngs:
                    jax = iax + 1
                    jrax = 4
                else:
                    jax = iax
                    jrax = 3
                self.parent.PLT.plot_traces(
                    self.P.axarr[rax, jax],
                    sfi,
                    PD,
                    protocol="IV",
                    ymin=ymin,
                    ymax=ymax,
                    xmax=150.0,
                    iax=iax,
                    figure=self.P.figure_handle,
                    show_title=False,
                    ivaxis=self.P.axarr[rax, jrax],  # accumulate IV's in right side
                    ivcolor=colors[iax],
                    iv_spike_color=spike_colors[dendmode],
                    spike_marker_size=1.5,
                    spike_marker_color=spike_colors[dendmode],
                    calx=calx,
                    caly=caly,
                    axis_index=iax,
                )
                if rax == 0:
                    self.P.axarr[rax, jax].set_title(column_titles[dendmode.title()], 
                        x=55.0, y=20., transform=self.P.axarr[rax, jax].transData,
                        ha="center")
                if iax == 0 and not show_pngs:
                    self.P.axarr[rax, 0].text(
                        -0.02, 0.5, f"{FD.BC_name:s}{iv:02d}", 
                        fontdict={"fontweight": "normal", "fontsize":11,
                        "ha": "right", "va": "center"},
                        transform = self.P.axarr[rax, 0].transAxes,
                    )
        if parent_figure is None:
            fig = FigInfo()
            fig.P = self.P
            fig.filename = set_figure_path(fignum=4, filedescriptor="CC", suppnum=3)
            timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
            fig.title[
                "title"
            ] = f"SBEM Project Figure 4 Modeling (Supplemental 3) ({timestamp_str:s})"
            return fig
        else:
            return self.P

    def Figure4_Supplemental4_Zin_removed(self):
        """
        All of the Zins for a supplemental figure
        """
        PZ.PlotZ()  # in plot_z.py

    def make_eff_fig(self):
        """
        Set up plotter grid for efficacy figure.
        """
        row1y = 5.5
        row1h = 3.0
        row2y = 3.0
        row2h = 2.0
        row3h = 2.0
        row3y = 0.5
        col1 = 1.0
        col2 = 4.5
        sizer = {
            "A": {
                "pos": [col1, 2.75, row1y, row1h],
                "labelpos": (-0.15, 1.02),
                "noaxes": True,
            },
            "B": {"pos": [col2, 2.5, row1y, row1h], "labelpos": (-0.15, 1.02)},
        }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) f
        sizer["C"] = {
            "pos": [col1, 2.5, row2y, row2h],
            "labelpos": (-0.15, 1.02),
        }
        sizer["D"] = {
            "pos": [col1, 2.5, row3y, row2h],
            "labelpos": (-0.15, 1.02),
        }
        sizer["E"] = {
            "pos": [col2, 2.5, row2y, row2h],
            "labelpos": (-0.15, 1.02),
        }
        sizer["F"] = {
            "pos": [col2, 2.5, row3y, row3h],
            "labelpos": (-0.15, 1.02),
        }
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            units="in",
            figsize=(8.0, 9.0),
            label=True,
            # showgrid=True,
        )
        return P

    def plot_combined_Eff_Rev(self):
        P = self.make_eff_fig()
        cell_number = 17
        cmap = "Dark2"
        figP2 = self.plot_revcorr(
            parent_figure=P, loc=(0, 8, 0.0, 5.5), start_letter="C", colormap=cmap
        )
        figP1 = self.plot_efficacy(
            parent_figure=P, loc=(0, 8, 6.0, 10.0), colormap=cmap
        )
        save_file = f"Fig_M2.pdf"
        fig = FigInfo()
        fig.P = P
        fig.filename = Path(
                self.config["figureIntermediateDirectory"], save_file)
        fig.title = {
            "title": "SBEM Project Figure Modeling Singles and Revcorr",
            "x": 0.99,
            "y": 0.95,
        }
        fig.title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.05}
        return fig

    def plot_efficacy(
        self,
        dendrite_mode="Full",
        parent_figure: Union[None, object] = None,
        loc: tuple = (0, 0, 0, 0),
        colormap="tab10",
    ):
        cell_number = 17
        example = FD.figure_efficacy_supplement_30dB[cell_number]

        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )

        sfi = Path(cellpath, Path(example[dendrite_mode]).name)
        if not sfi.is_dir():
            return
        fn = sorted(list(sfi.glob("*")))
        PD = self.newPData()
        calx = 800.0
        if parent_figure is None:
            parent_figure = self.make_eff_fig()
        EFP = EF.EfficacyPlots(parent_figure=parent_figure)
        EFP.plot_efficacy(
            datasetname="Full", ax=EFP.parent_figure.axdict["B"], loc=loc, clean=True
        )

        # stacked IV in first column:
        self.P_Eff_SingleInputs = PH.regular_grid(
            len(fn),
            1,
            order="rowsfirst",
            figsize=(8, 10),
            units="in",
            showgrid=False,
            verticalspacing=0.0,
            horizontalspacing=0.5,
            margins={
                "bottommargin": 6.0,
                "leftmargin": 1,
                "rightmargin": 4.5,
                "topmargin": 0.5,
                # "width": 0.17,
                # "height": 0.8,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=EFP.P,
        )
        nfiles = len(fn)
        cmx = mpl.colormaps[colormap]
        colors = [cmx(float(i) / nfiles) for i in range(nfiles)]

        calx = 800.0
        for n in range(nfiles):
            sfile = Path(sfi, fn[n])
            if n == nfiles - 1:
                calxv = calx
            else:
                calxv = None

            self.parent.PLT.plot_traces(
                ax=self.P_Eff_SingleInputs.axarr[n, 0],
                fn=sfile,
                PD=PD,
                protocol="runANSingles",
                ymin=-70.0,
                ymax=0.0,
                xmin=400.0,
                xmax=900.0,
                iax=n,
                rep=0,
                figure=self.P_Eff_SingleInputs.figure_handle,
                longtitle=True,
                trace_color=colors[n],
                ivaxis=None,
                ivcolor="k",
                iv_spike_color="k",
                spike_marker_size=1.5,
                spike_marker_color="r",
                calx=calxv,
                caly=-20.0,
                calt=50.0,
                calv=20.0,
            )

            self.P_Eff_SingleInputs.axarr[n, 0].text(
                -0.02,
                0.5,
                f"Syn {n+1:d}",
                fontsize=8,
                horizontalalignment="right",
                verticalalignment="center",
                transform=self.P_Eff_SingleInputs.axarr[n, 0].transAxes,
            )

        EFP.P.figure_handle.text(
            0.99,
            0.01,
            f"Cell {cell_number:d}",
            fontsize=10,
            horizontalalignment="right",
            verticalalignment="bottom",
        )
        save_file = Path(
            self.config["figureIntermediateDirectory"], "Fig_M2_Efficacy_Revcorr.pdf"
        )
        fig = FigInfo()
        fig.P = EFP.P
        fig.filename = save_file
        fig.title["title"] = "SBEM Project Figure 5 Modeling: Efficacy and Revcorr"
        return fig

    def plot_stacked_traces(
        self,
        cells=None,
        figure=None,
        dBSPL:str="30dB",
        simulation_experiment:str = "Full",
        show_title=True,
        axes: Union[list, None] = None,
        calxp: float = 800.0,
        calv: float = 20.0,
        maxstack: int = 9,
        maxtraces: int = 20, # up to 20 in a stack; fewer if desired
        cal_pos: int = 0,
    ):
        """
        What it says: plot traces in a stacked format

        Parameters
        ----------
        cells : a list of ints (if none, default is to use all of the cells in the
            GRPDEF.gradeAList)
        figure: Figure handle for the plot passed to plot_traces in plot_sims.py)
        axes : a list
             of axes passed to plot_traces in plot_sims.py
        calxp : float
            the position for the calibration bar in x (time), in msec
        calv : float
            the size of the calibration bar in voltage (in mV)
        maxstack: int (default 9)
            determines the height of the stack,  elative to the offset (see code)
        maxtraces: int (default 20)
            How many of the incoming traces to plot in a stack (first maxtraces are taken)

        Returns
        -------
        Nothing
        """
        assert simulation_experiment in ["Full", "NoUninnervated2"]
        if cells is None:
            cells = GRPDEF.grAList()
        trace_ht = 80  # mV
        if maxtraces >= maxstack:
            maxtraces = maxstack
        for ic, cellN in enumerate(cells):
            # if ic > 0:  # for quick tests
            #               continue
            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cellN:02d}",
                "Simulations",
                "AN",
            )
            sfiles = None
            if simulation_experiment == "NoUninnervated2":
                dBSPL = "30dB"
                sfiles = Path(
                    cellpath,
                    Path(
                        FD.figure_cell9_nouninnervated2[cellN][simulation_experiment]
                    ).name,
                )                
            elif dBSPL == "30dB":
                sfiles = Path(
                    cellpath,
                    Path(
                        FD.figure_efficacy_supplement_30dB[cellN][simulation_experiment]
                    ).name,
                )
            elif dBSPL == "Spont":
                sfiles = Path(
                    cellpath,
                    Path(
                        FD.figure_efficacy_supplement_Spont[cellN][
                            simulation_experiment
                        ]
                    ).name,
                )

            if not sfiles.is_dir():
                return
            fn = sorted(list(sfiles.glob("*")))
            PD = self.newPData()

            nfn = len(fn)
            ymax = 20.0
            if simulation_experiment != "Full":
                ymax = 40.0
            yoffset = -90.0
            ymin = maxtraces * yoffset
            ymax = 20.0
            if simulation_experiment != "Full":
                ymax = 40.0
            xmin = 400.0
            xmax = 900.0
            nmax = min(maxtraces, len(fn))
            for n, filename in enumerate(fn):
                if n >= maxtraces:
                    break
                if (n == (nmax - 1)) and (
                    ic == cal_pos
                ):  # (len(cells)-1): # cal bar on first axis
                    calxv = calxp
                    calyv = -120.0 + n * yoffset
                    iax = n
                else:
                    iax = None
                    calxv = None
                    calyv = -50.0
                y0 = n * yoffset
                cprint("c", f"iax: {str(iax):s}, calxv = {str(calxv):s}  calyv = {str(calyv):s}")
                cprint("c", f"xmin: {xmin:.1f}, xmax: {xmax:.1f}  ymin: {ymin:.1f}  ymax: {ymax:.1f}")
                cprint("c", f"y0: {y0:.1f}")
                self.parent.PLT.plot_traces(
                    ax=axes[ic],
                    fn=filename,
                    PD=PD,
                    protocol="runANSingles",
                    ymin=ymin,
                    ymax=ymax,
                    xmin=xmin,  # in msec
                    xmax=xmax,
                    yoffset=n * yoffset,
                    iax=n,
                    # nax=len(fn),
                    rep=0,
                    figure=figure,
                    longtitle=True,
                    show_title=show_title,
                    ivaxis=None,
                    ivcolor="k",
                    iv_spike_color="r",
                    spike_marker_size=1.5,
                    spike_marker_color="r",
                    calx=calxv,
                    caly=calyv,
                    calt=50.0,
                    calv=calv,
                )
                axes[ic].annotate(
                    text=f"{n+1:d} ",
                    xy=(xmin, y0 - 60.0),  # ms
                    fontsize=8,
                    horizontalalignment="right",
                    verticalalignment="center",
                    # transform=axes[ic].transAxes,
                )

                if n == 0:
                    axes[ic].set_title(
                        f"{FD.BC_name:s}{cellN:02d}",
                        loc="center",
                        fontdict={
                            "verticalalignment": "baseline",
                            "fontsize": 9,
                            "fontweight": "normal",
                        },
                    )
        return

    def plot_efficacy_supplement(self, cells=None, parent_figure=None, traces=True):
        if cells is None:
            cells = GRPDEF.grAList()
        print("Efficacy Supplement Plot")
        effp = []
        simulation_experiment = "Full"
        hspc = 0.04
        vspc = 0.05
        rowspc = 0.1
        ncol = 5
        nrow = int(len(cells) / 5)
        lmar = 0.05
        rmar = 0.05
        tmar = 0.05
        bmar = 0.05
        width = 0.14
        height = 0.38
        xp = lmar + np.arange(ncol) * (width + hspc)
        yp = bmar + np.arange(nrow) * (height + rowspc)
        if parent_figure is None:
            self.P = PH.regular_grid(  # dummy figure space
                1,
                1,
                figsize=(8.0, 8.0),
                order="rowsfirst",
                units="in",
                showgrid=False,
                parent_figure=None,
            )
            EFP = EF.EfficacyPlots(parent_figure=self.P)  # , cols=1, rows=1)
            EFP.parent_figure.figure_handle.set_size_inches((8.0, 8.0))
        else:
            EFP = EF.EfficacyPlots(parent_figure=parent_figure)

        for ic, cellN in enumerate(cells):
            # if ic > 0:  # for quick tests
            #               continue
            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cellN:02d}",
                "Simulations",
                "AN",
            )
            sfiles = Path(
                cellpath,
                Path(FD.figure_efficacy_supplement_30dB[cellN][simulation_experiment]).name,
            )
            if not sfiles.is_dir():
                return
            fn = sorted(list(sfiles.glob("*")))
            PD = self.newPData()

            row_n = 1 - int(ic / 5)
            col_n = ic % 5
            calxp = 800.0
            if traces:
                P_Eff = PH.regular_grid(
                    len(fn),
                    1,
                    order="rowsfirst",
                    # figsize=(8., 4.5),
                    showgrid=False,
                    verticalspacing=0.01,
                    # horizontalspacing=0.1,
                    margins={
                        "bottommargin": yp[row_n],
                        "leftmargin": xp[col_n],
                        "width": width,
                        "height": height,
                    },
                    labelposition=(-0.05, 1.06),
                    parent_figure=EFP.parent_figure,
                    # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
                )
                effp.append(P_Eff)
                nfn = len(fn)
                ymax = 20.0
                if simulation_experiment != "Full":
                    ymax = 40.0
                yoffset = -60.0
                ymin = nfn * -60.0
                for n in range(len(fn)):
                    if n == nfn - 1:
                        calxv = calxp
                    else:
                        calxv = None
                    self.parent.PLT.plot_traces(
                        ax=P_Eff.axarr[n, 0],
                        fn=fn[n],
                        PD=PD,
                        protocol="runANSingles",
                        ymin=ymin,
                        ymax=ymax,
                        xmin=400.0,
                        xmax=900.0,
                        yoffset=n * yoffset,
                        iax=n,
                        rep=0,
                        figure=P_Eff.figure_handle,
                        longtitle=True,
                        ivaxis=None,
                        ivcolor="k",
                        iv_spike_color="r",
                        spike_marker_size=1.5,
                        spike_marker_color="r",
                        calx=calxv,
                        caly=-110.0,
                        calt=50.0,
                        calv=20.0,
                    )
                    P_Eff.axarr[n, 0].text(
                        -0.02,
                        0.5,
                        f"{n+1:d}",
                        fontsize=8,
                        horizontalalignment="right",
                        verticalalignment="center",
                        transform=P_Eff.axarr[n, 0].transAxes,
                    )

                    if n == 0:
                        P_Eff.axarr[n, 0].set_title(
                            f"VCN_c{cellN:02d}",
                            loc="center",
                            fontdict={
                                "verticalalignment": "baseline",
                                "fontsize": 9,
                                "fontweight": "normal",
                            },
                        )

        title = "SBEM Project Figure 5 Modeling : Efficacy, Supplemental"
        save_file = Path(
            self.config["figureIntermediateDirectory"],
            f"Fig_M2_Supplemental_{simulation_experiment:s}.pdf",
        )
        fig = FigInfo()
        if traces:
            fig.P = self.P
        else:
            fig.P = None
        fig.filename = save_file
        fig.title["title"] = title
        return fig

    def _load_rcdata(self, dBSPL):
        rc_datafile = Path(
            self.config["baseDataDirectory"],
            self.config["revcorrDataDirectory"],
            f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl",
        )
        with open(rc_datafile, "rb") as fh:
            d = FPM.pickle_load(fh)
        return d
    



    def Figure5_Main(self, supplemental1=False, final_plot=True):
        """
        Generate Figure 5 for the paper. Combined bits from various other plots
        and revcorrs/singles
        
        Parameters
        ----------
        supplemental1 : set True to plot the associated supplement

        final_plot: set True to remove the "dbSPL" from the upper left corner of the plot
        """

        print(f"Figure 5 main: supplemental1={str(supplemental1):s}")
        if supplemental1:
            example_cells = [
                10,
                6,
                2,
                13,
                18,
                11,
            ]  # in order of spikes for largest input
        else:
            example_cells = [5, 30, 9, 17]  # in order of spikes for largest input

        Figure5_stim_level = "30dB"  # "30dB" or "Spont" are valid settings.
        participation_dB = 30  # the stimulus level at which the participation is
        # compared versus participation during spont

        start_letter = "A"
        parent_figure = None
        yrow0 = 0.75
        yrow1 = 2.50+0.75
        yh = 1.75
        #xw = 0.85 * yh
        xw = 1.15 * yh
        xl = 1.05

        xlp = [xl + (xw + 0.72) * i for i in range(4)]
        lpos = (-0.12, 1.06)
        if not supplemental1:
            sizer = {
                # "B": {"pos": [6.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                # "C": {"pos": [9.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                # "F": {"pos": [6.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
                # "G": {"pos": [9.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
                "D": {
                    "pos": [xlp[0], xw, yrow1, yh],
                    "labelpos": lpos},
                "E": {
                    "pos": [xlp[1], xw, yrow1, yh],
                    "labelpos": lpos,
                },
                "F": {
                    "pos": [xlp[2], xw, yrow1, yh],
                    "labelpos": lpos,
                },

                "G": {
                    "pos": [xlp[3], xw, yrow1, yh],
                    "labelpos": lpos,
                },
                "H": {
                    "pos": [xlp[0], xw, yrow0, yh],
                    "labelpos": lpos,
                },
                "I": {
                    "pos": [xlp[1], xw, yrow0, yh],
                    "labelpos": lpos,

                },
                "J": {
                    "pos": [xlp[2], xw, yrow0, yh],
                    "labelpos": lpos,

                },
                "K": {
                    "pos": [xlp[3], xw, yrow0, yh],
                    "labelpos": lpos,

                },
            }
            figsize = (12, 8+2.5)
            cal_pos = 0
        else:
            sizer = {}
            figsize = (9, 8)
            cal_pos = 1

        xw = 1.1
        xw2 = 1.0
        trace_axes = []
        if supplemental1:
            yA1 = 4.125
            yh1 = 3.5
            yB2 = 2.5
            yC2 = 0.5
            yh2 = 1.2
            grid = False
        else: # main figure
            yh2 = 1.2
            yB2 = 3.5 + 2.7 - 0.5 + 2.5
            yC2 = 3.5 + 0.6 - 0.5 + 2.5
            yA1 = 3.25 - 0.5 +2.5 - 0.105  # to align with B
            yh1 = 4.25
            grid = False
        for j in range(len(example_cells)):
            i = j + 1
            pan_rev, pan_vm = self.Figure5_assign_panel(supplemental1, i)
            if not supplemental1:
                xl1 = j * 1.25 + 0.75
                xl2 = j * 1.25 + 6.5  # set panels to the right
            else:
                xl1 = j * 1.25 + 0.75
                xl2 = j * 1.25 + 0.75  # set panels on the bottom matching columns
            axn = f"A{i:d}"
            trace_axes.append(axn)
            sizer[axn] = {
                "pos": [xl1, xw, yA1, yh1],
                "labelpos": (-0.15, 1.03),
                "noaxes": True,
            }
            sizer[pan_rev] = {  # reverse correlatoin
                "pos": [xl2, xw2, yB2, yh2],
                "labelpos": (-0.25, 1.1),
                # "noaxes": True,
            }
            sizer[pan_vm] = {
                "pos": [xl2, xw2, yC2, yh2],
                "labelpos": (-0.25, 1.1),
                "noaxes": True,
            }
        # dict pos elements are [left, width, bottom, height] for the axes in the plot.
        # gr = [(a, a+1, 0, 1) for a in range(0, 8)]
        # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            units="in",
            figsize=figsize,
            label=True,
            fontsize={"tick": 8, "label": 10, "panel": 13},
            showgrid=grid,
            parent_figure=parent_figure,
        )
        # PH.show_figure_grid(P)
        if not final_plot:
            P.figure_handle.text(x=0.01, y=figsize[1]-0.1, s=f"{Figure5_stim_level:s}",
            fontdict={"fontsize": 11, "fontweight": "bold", "verticalalignment": "top"},
            transform=P.figure_handle.dpi_scale_trans)

        synperum2 = 0.7686  # taken from cell_config.py, line 127 (11/15/2021)

        def plot_participation(ax, n, a, b, dB=0, color=None):
            ap = a[n][0].participation / a[n][0].npost_spikes
            bp = b[n][0].participation / b[n][0].npost_spikes
            ax.plot(
                [a[n][0].sites / synperum2, a[n][0].sites / synperum2],
                [ap, bp],
                "-",
                color=color,

            )
            ax.scatter(a[n][0].sites / synperum2, ap, marker="o", color=color, alpha=0.6)
            ax.scatter(a[n][0].sites / synperum2, bp, marker="x", color=color, alpha=0.6)
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 300)
            ax.set_ylim(0, 1.0)
            ax.set_ylabel(f"Participation at 0 and {dB:2d} dBSPL")
            PH.talbotTicks(ax, floatAdd={"x": 0, "y": 2})

        def plot_diff_participation(ax, n, a, b, dB=0, color=None, legend=True):
            ap = a[n][0].participation / a[n][0].npost_spikes
            bp = b[n][0].participation / b[n][0].npost_spikes
            mark = GRPDEF.get_group_symbol(n)
            ax.scatter(
                a[n][0].sites / synperum2,
                bp / ap,
                marker=mark,
                color=color,
                edgecolors='w',
                linewidths=0.5,
                label=f"{FD.BC_name:s}{n:02d}",
                clip_on=False,
                s=36,
                alpha=0.6
            )
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 350)
            ax.set_ylim(0, 3)
            ax.set_ylabel(f"Participation ratio {dB:2d}/{0:2d} dBSPL")
            PH.nice_plot(ax, position=self.axis_offset, direction="outward")
            PH.set_axes_ticks(ax=ax,
                xticks = [0, 100, 200, 300],
                xticks_str = ['0', '100', '200', '300'],
                # xticks_pad:Union[List, None]=None,
                x_minor = [50, 150, 250, 350],
                major_length = 3.0,
                minor_length =1.5,
                yticks = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
                yticks_str=["0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"], 
                yticks_pad=[1]*7,
                y_minor=None,
                fontsize=8,
            )
            PH.referenceline(ax, 1.0)
            if legend:
                ax.legend(fontsize=8, loc="upper right", ncol=2)

        def plot_clustering(ax):
            EF.EffClusters([ax], clip_on=False)
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 350)
            ax.set_ylim(0, 1.0)
            ax.set_ylabel(f"Efficacy (Bushy spikes/input spikes)")
            PH.nice_plot(ax, position=self.axis_offset, direction="outward")

            PH.talbotTicks(
                ax,
                density=(1.0, 1.0),
                insideMargin=0,
                tickPlacesAdd={"x": 0, "y": 2},
                floatAdd={"x": 0, "y": 2},
                axrange={"x": (0, 300.0), "y": (0, 1)},
                pointSize=None,
            )

        def plot_single_input(ax, legend: bool = True):
            EF.eff_one_input(ax=ax, legend=legend)

        # Plot the raw traces
        axl = [P.axdict[axi] for axi in trace_axes]
        self.plot_stacked_traces(
            cells=example_cells,
            dBSPL=Figure5_stim_level,
            figure=P.figure_handle,
            axes=axl,
            maxstack=10,
            show_title=False,
            calxp=600.0,
            cal_pos=cal_pos,
        )
        for ax in axl:
            ax.set_zorder(0)


        if not supplemental1:
            # plot all the different analyses.

            # Efficacy plot vs. input size (all)
            for s in ["D", "E", "F", "G", "H", "I", "J", "K"]:
                P.axdict[s].set_zorder(100)
            EFP = EF.EfficacyPlots(parent_figure=P)
            EFP.plot_efficacy(
                "Full", datasetname_added="Added", ax=P.axdict["D"], clean=True
            )

            # efficacy for a single sized input vs. dendritic area
            plot_single_input(P.axdict["H"], legend=False)

            # input pattern plot
            # PATSUM.Figure4F_pattern_plot(axin=P.axdict["E"], dataset="Spont", mode="multi")
            PATSUM.Figure5E_pattern_plot(axin=P.axdict["E"], dataset="Spont", mode='mmcd', cell_legend = True)  
            # participation
            ds = self._load_rcdata("Spont")
            drc = self._load_rcdata(f"{participation_dB:2d}dB")
            sns.set_palette(palette="tab10", n_colors=10)
            palette = sns.color_palette(palette="tab10", n_colors=len(ds.keys()))
            for i, c in enumerate(ds.keys()):
                # plot_participation(P.axdictax[0], c, ds, drc, dB=dB, color=palette[i])
                plot_diff_participation(
                    P.axdict["F"],
                    c,
                    ds,
                    drc,
                    dB=participation_dB,
                    color=palette[i],
                    legend=False,
                )


            # Cumulative plots
            self.plot_revcorr_compare(
                parent_figure=P,
                axlist=[P.axdict["G"]],
                dBSPLs=["Spont", "30dB"],
                legend=False,
            )
            # the following is set in plot_revcorr compare
            PH.set_axes_ticks(ax=P.axdict["G"],
                xticks =     [0,    4,   8,   12],
                xticks_str = ["0", "4", "8", "12"],
                xticks_pad=None,
                x_minor = [2, 6, 10],
                major_length = 3.0,
                minor_length =1.5,
                yticks =   [0,    0.2,   0.4,   0.6,   0.8,   1.0],
                yticks_str=["0", "0.2", "0.4", "0.6", "0.8", "1.0"], 
                yticks_pad=[1, 1, 1, 1, 1, 1,],
                y_minor=None,
                fontsize=8,
            )
            PH.nice_plot(P.axdict["G"], position=self.axis_offset, direction="outward", ticklength=3)
            MTC = morphology_thr_correlations
            MTC.AIS().MorphvsThr(ax=P.axdict["I"], compartment="dendrite")
            MTC.AIS().MorphvsThr(ax=P.axdict["J"], compartment="soma")
            MTC.AIS().AISLengthThr(ax=P.axdict["K"])

            synlabel_num = 5  # this sets which cell the scale bar will be plotted with
        else:
            synlabel_num = 10

        # revcorrs and traces
        self.plot_revcorr_panels(
            cells=example_cells,
            parent_figure=P,
            supplemental1=supplemental1,
            # dBSPL="30dB",
            dBSPL=Figure5_stim_level,
            synlabel_num=synlabel_num,
            show_title=False,
        )
        # self.plot_efficacy_supplement(cells=example_cells, parent_figure=P, traces=False)
        arrow_xy = {
            5: [-0.95, -45.0],
            30: [-0.98, -45.0],
            9: [-0.8, -46.0],
            17: [-0.660, -45.0],
            10: [-0.95, -45.0],
            6: [-0.95, -45.0],
            2: [-0.95, -45.0],
            13: [-0.75, -45.0],
            18: [-0.6, -45.0],
            11: [-0.5, -40.0],
        }
        # Revcorr axes cleanup
        for j, celln in enumerate(example_cells):
            pan_rev, pan_vm = self.Figure5_assign_panel(supplemental1, j + 1)
            ax = P.axdict[pan_rev]
            ax.set_ylim(0, 0.8)
            ax.set_xlim(-5.0, 2.5)
            ax2 = P.axdict[pan_vm]
            ax2.set_xlim(-5.0, 2.5)
            ax2.set_ylim(-70, 0)
            if celln in list(arrow_xy.keys()):
                ax2.annotate(
                    "",
                    xy=arrow_xy[celln],
                    xycoords="data",
                    xytext=(arrow_xy[celln][0] - 0.75, arrow_xy[celln][1] + 4),
                    arrowprops=dict(
                        facecolor="black",
                        edgecolor="black",
                        linewidth=0.5,
                        width=5,
                        headwidth=5.0,
                        headlength=8.0,
                        shrink=0.05,
                    ),
                )

            if j == 0:
                ax.set_ylabel("Pre-Post\nCoinc. Rate (Hz)", ha="center", fontsize=10)
                ax2.set_ylabel("Vm (mV)", ha="center", fontsize=10)
            else:
                PH.noaxes(ax, whichaxes="y")
                PH.noaxes(ax2, whichaxes="y")
            # ax.xaxis.set_minor_locator(MultipleLocator(2))

            ax.tick_params(which="major", length=4, direction="out")
            ax.tick_params(which="minor", length=2, direction="out")

            # ax2.xaxis.set_minor_locator(MultipleLocator(2))
            ax2.tick_params(which="major", length=4, direction="out")
            ax2.tick_params(which="minor", length=2, direction="out")

        fig = FigInfo()
        if parent_figure is not None:
            fig.P = parent_figure
        else:
            fig.P = P
        if not supplemental1:
            fig.filename = set_figure_path(fignum=5, filedescriptor="Ephys_2_main_v19")
            fig.title[
                "title"
            ] = "SBEM Project Figure 5 (main) Modeling: singles inputs, efficacy and revcorr, revised version 9"
        else:
            fig.filename = set_figure_path(fignum=5, filedescriptor="Revcorr_V5", suppnum=1)
            #"Figure5/Figure5_supp/Figure5_Supplemental1_Revcorr_V4.pdf"
            fig.title[
                "title"
            ] = "SBEM Project Figure 5 Modeling: Supplemental 1: other cells single inputs and revcorr"

        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        # fig.title2 = title2
        return fig

    def Figure5_Supplemental1(self):
        fig = self.Figure5_Main(supplemental1=True)
        return fig

    def Figure5_Supplemental2_removed(self):
        """Plot AIS and efficacy

        Returns:
            _type_: _description_

            NOTE: The rest of this figure is made in PRISM..... 

        """
        fig = EF.eff_ais(EF.data_Full, save_fig=True, figinfo=FigInfo(show_figure_name=False))
        return fig

    def Figure5_assign_panel(self, supplemental1: bool = False, index: int = 0):
        if not supplemental1:
            revcorr_panel = f"B{index:d}"
            vm_panel = f"C{index:d}"
        else:
            revcorr_panel = f"B{index:d}"
            vm_panel = f"C{index:d}"
        return revcorr_panel, vm_panel 

    def Figure5_Supplemental2(self):
        PATSUM.Figure5_Supplemental2_Patterns()  # writes its own figure to the directory

    def plot_all_revcorr(self):
        for cell in GRPDEF.grAList():
            fig = self.plot_revcorr(cell)
            if fig is not None:
                self.save_figure(fig)
        return None

    def _get_revcorr(
        self,
        cell_number: int,
        dBSPL="Spont",
        parent_figure: Union[object, None] = None,
        recompute=True,  # if you need to recompute the revcorrs for all the grade A cells, just set this True
    ) -> tuple:
        """
        Get the revcorr data associated with the cell number
        and the stimulus level (dbSPL, as a string)
        """
        cell_revcorr = FD.figure_revcorr[cell_number]
        run_calcs = False
        rc_datafile = Path(
            self.newPData().revcorrpath, f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl"
        )
        cprint("c", f"RevCorr datafile: {str(rc_datafile):s}")
        if rc_datafile.is_file() and not recompute:
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = FPM.pickle_load(fh)

            RCD = all_RCD_RCP[cell_number][0]
            RCP = all_RCD_RCP[cell_number][1]
            PD = self.newPData()
            P = None
        else:
            cprint(
                "r",
                "You must run revcorr_supplement plot first, then the data file we need will be present",
            )
            run_calcs = True
            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, Path(cell_revcorr[dBSPL]).name)
            if not sfi.is_dir():
                return None, None, None, None
            fn = sorted(list(sfi.glob("*")))
            P, PD, RCP, RCD = self.parent.PLT.compute_revcorr(
                P=None,  # no plotting (and return P is None)
                gbc=str(cell_number),
                fn=fn[0],
                PD=self.newPData(),
                protocol="runANPSTH",
                revcorrtype="RevcorrSPKS",
                # thr=-20.0,
                width=4.0,
            )
        return P, PD, RCP, RCD

    def plot_revcorr(
        self,
        cellN=None,
        parent_figure: Union[object, None] = None,
        loc: tuple = (0.0, 0.0, 0.0, 0.0),
        start_letter="B",  # default if no parent figure
        colormap="tab10",
    ):
        if cellN == None:
            cell_number = 17
            example = FD.figure_revcorr_example[cell_number]
        else:
            cell_number = cellN
            example = FD.figure_revcorr[cell_number]

        P, PD, RCP, RCD = self._get_revcorr(cell_number=cell_number, dBSPL="Spont")
        dBSPL = RCP.ri.dB
        if PD is None:
            return  # unable to get the revcorr
        str_a = string.ascii_uppercase
        p_labels = str_a[
            str_a.find(start_letter) : str_a.find(start_letter) + 4
        ]  # will fail if you have > 26
        if parent_figure is None:
            box_sizex = 2.25  # inches
            box_sizey = 2.2  # inches
            row3y = 0.5
            row2y = 3.0

            sizer = {
                p_labels[0]: {
                    "pos": [1, box_sizex, row2y, box_sizey],
                    "labelpos": (-0.15, 1.02),
                    "noaxes": True,
                },
                p_labels[1]: {
                    "pos": [1, box_sizex, row3y, box_sizey],
                    "labelpos": (-0.15, 1.02),
                },
                # "C": {"pos": [0.52, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
                # "D": {"pos": [0.52, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
                p_labels[2]: {
                    "pos": [4.25, box_sizex, row2y, box_sizey],
                    "labelpos": (-0.15, 1.02),
                },
                p_labels[3]: {
                    "pos": [4.25, box_sizex, row3y, box_sizey],
                    "labelpos": (-0.15, 1.02),
                },
            }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
            P = PH.arbitrary_grid(
                sizer,
                order="columnsfirst",
                units="in",
                figsize=(8, 10),
                label=True,
                fontsize={"tick": 9, "label": 10, "panel": 14},
                # margins={
                #     "bottommargin": 0.1,
                #     "leftmargin": 0.1,
                #     "rightmargin": 0.1,
                #     "topmargin": 0.1+loc[2],
                # },
                parent_figure=parent_figure,
            )
        else:
            P = parent_figure  # PH.show_figure_grid(P, 6, 7)
        # mpl.show()
        # return
        sax = P.axdict
        sax0 = sax[p_labels[0]]
        sax1 = sax[p_labels[1]]
        sax2 = sax[p_labels[2]]
        sax3 = sax[p_labels[3]]

        summarySiteTC = self.parent.PLT.plot_revcorr2(
            P,
            PD,
            RCP,
            RCD,
            cell_number=cell_number,
            start_letter=p_labels[0],
            colormap=colormap,
        )
        # sax1.set_xlim(-5, 0)

        maxp = np.max(RCD.pairwise)
        psh = RCD.pairwise.shape
        pos = np.zeros((psh[0], psh[1], 2))
        for i in range(RCP.ninputs):
            for j in range(RCP.ninputs):
                # print(f"{pairwise[i,j]:.3f}", end='  ')
                pos[i, j, 0] = i + 1
                pos[i, j, 1] = j + 1

        cmx = mpl.colormaps[colormap]
        clist = [cmx(float(isite) / RCP.ninputs) for isite in range(RCP.ninputs)]
        sax2.scatter(
            RCD.sites,
            RCD.participation / RCD.nspikes,
            marker="o",
            color=clist,
            sizes=[42],
        )
        sax2.set_ylim((0, 1.0))
        sax2.set_xlim(left=0)
        sax2.set_xlabel("# Release Sites")
        sax2.set_ylabel("Participation")

        sax3.plot(
            np.arange(len(RCD.ynspike)) + 1,
            RCD.ynspike,
            "k^-",
            markersize=5,
            clip_on=False,
        )
        sax3.set_ylim(0, 1.05)
        sax3.set_xlabel(
            f"# Inputs in [{RCD.pre_w[0]:.1f} to {RCD.pre_w[1]:.1f}] before spike"
        )
        sax3.set_ylabel(f"Cumulative Bushy Spikes\nwith N AN inputs")
        sax3.set_xlim(1, RCP.ninputs)
        PH.nice_plot(sax3, position=self.axis_offset, direction="outward")
        # PH.cleanAxes(P.axarr.ravel())
        # PH.noaxes(sax0)

        # PH.talbotTicks(
        #     sax1,
        #     tickPlacesAdd={"x": 1, "y": 2},
        #     floatAdd={"x": 1, "y": 2},
        #     # pointSize=7,
        # )
        PH.talbotTicks(
            sax2,
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            # pointSize=7,
        )
        PH.talbotTicks(
            sax3,
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            # pointSize=7,
        )

        if dBSPL in ["Spont", 0.0]:
            if cellN is None:
                save_file = Path(
                    self.config["figureIntermediateDirectory"], "Fig_M3.pdf"
                )
            else:
                save_file = Path(
                    self.config["figureIntermediateDirectory"],
                    f"Fig_Revcorr/Revcorr_VCN_c{cell_number:02d}.pdf",
                )
        else:
            if isinstance(dBSPL, float):
                db = f"{int(dBSPL):d}"
            else:
                db = dBSPL
            save_file = Path(
                self.config["figureIntermediateDirectory"], f"Fig_M3_{db:s}.pdf"
            )
        title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.01}
        # save_file = "Fig_M2_Efficacy_Revcorr.pdf"
        fig = FigInfo()
        if parent_figure is not None:
            fig.P = parent_figure
        else:
            fig.P = P
        fig.filename = save_file
        fig.title["title"] = "SBEM Project Figure 2 Modeling: Efficacy and Revcorr"
        fig.title2 = title2
        return fig

    def plot_revcorr_supplement_spont(self):
        fig = self.plot_revcorr_supplement("Spont", supplemental1=True)
        return fig

    def plot_revcorr_supplement_30dB(self):
        fig = self.plot_revcorr_supplement("30dB", supplemental1=True)
        return fig

    def plot_revcorr_supplement_40dB(self):
        fig = self.plot_revcorr_supplement("40dB", supplemental1=True)
        return fig

    def plot_revcorr_panels(
        self,
        dBSPL: str,
        cells=None,
        parent_figure=None,
        supplemental1=False,
        show_title=True,
        synlabel_num: int = 0,
        colormap: str = "magma",  # default color map
        save_calcs: bool = False,  # set to True if need to update.
    ):
        if cells is None:
            cells = GRPDEF.grAList()
        ncells = len(cells)

        if parent_figure is None:
            panel_labels = []
            chars = string.ascii_uppercase
            for cn in range(4):
                for r in range(len(cells)):
                    panel_labels.append(f"{chars[cn]:s}{r+1:d}")

            P = PH.regular_grid(
                rows=4,
                cols=len(cells),
                order="rowsfirst",
                figsize=(14, 11),
                # showgrid=True,
                margins={
                    "bottommargin": 0.1,
                    "leftmargin": 0.08,
                    "rightmargin": 0.05,
                    "topmargin": 0.1,
                },
                verticalspacing=0.03,
                horizontalspacing=0.03,
                panel_labels=panel_labels,
            )
            P.figure_handle.text(
                0.175,
                0.93,
                f"Coinc. Rate (Hz)",
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="bottom",
            )
            P.figure_handle.text(
                0.39,
                0.93,
                f"Membrane Potential",
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="bottom",
            )
            P.figure_handle.text(
                0.61,
                0.93,
                f"Participation",
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="bottom",
            )
            P.figure_handle.text(
                0.83,
                0.93,
                f"Cumulative Inputs to Spike",
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="bottom",
            )
            PH.cleanAxes(P.axarr.ravel())

        else:
            P = parent_figure

        run_calcs = False
        # P, PD, RCP, RCD = self._get_revcorr(cell_number=cellN, dBSPL = "Spont")
        rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
        if rc_datafile.is_file() and not run_calcs:
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = FPM.pickle_load(fh)
        # else:
        #     print(
        #         "Must run revcorr_supplement plot first, then the file we need will be present"
        #     )
        #     run_calcs = False
        #     all_RCD_RCP = {}  # all revcorr data

        i_plot = 0
        if cells is None:
            cells = GRPDEF.grAList()

        for i, cell_number in enumerate(cells):
            cprint("m", f"CELL: {cell_number:d}")
            PR, PD, RCP, RCD = self._get_revcorr(cell_number=cell_number, dBSPL=dBSPL, recompute=False)
            if PD is None:
                cprint("r", "PD is none in plot_revcorr_supplement")
                continue

            all_RCD_RCP[cell_number] = [RCD, RCP]
            # rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
            # if rc_datafile.is_file() and not run_calcs:
            #     with open(rc_datafile, "rb") as fh:
            #         all_RCD_RCP = FPM.pickle_load(fh)
            #
            dataset = FD.figure_revcorr[cell_number]

            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, Path(dataset[dBSPL]).name)
            if not sfi.is_dir():
                return
            fn = sorted(list(sfi.glob("*")))
            print("Revcorr supplement file: ", fn)
            print("Revcorr dbspl: ", dBSPL)
            if i == 0:
                tcal = False  # point font for cal bar
            else:
                tcal = False  # no cal bar
            if synlabel_num == 0:  # labels all of them
                synlabel = True
            else:
                if cell_number == synlabel_num:
                    synlabel = True
                else:
                    synlabel = False

            ax_top_row_name, ax_bot_row_name = self.Figure5_assign_panel(
                supplemental1, index=i_plot + 1
            )
            ax_top_row = P.axdict[ax_top_row_name]
            ax_bot_row = P.axdict[ax_bot_row_name]

            if parent_figure is None:
                axlist = P.axarr[0:2, i_plot]
            else:
                axlist = [
                    ax_top_row,
                    ax_bot_row,
                    None,
                ]
            summarySiteTC = self.parent.PLT.plot_revcorr2(
                P,
                PD,
                RCP,
                RCD,
                cell_number=cell_number,
                axarray=axlist,
                calbar_show=tcal,
                calbar_fontsize=7,
                yaxis_label=False,
                synlabel=synlabel,
                colormap=colormap,
            )

            cprint(
                "g",
                f"  Mean Pre:  {np.mean(RCD.mean_pre_intervals):7.3f} ms  ({1e3/np.mean(RCD.mean_pre_intervals):7.1f} Hz)",
            )
            cprint(
                "g",
                f"  Mean Post: {RCD.mean_post_intervals:7.3f} ms  ({1e3/RCD.mean_post_intervals:7.1f} Hz)",
            )
            # if parent_figure is None:
            ax_top_row.text(
                0.5,
                1.0,
                f"{FD.BC_name:s}{cell_number:02d}",
                fontsize=9,
                color="k",
                transform=ax_top_row.transAxes,
                horizontalalignment="center",
            )
            i_plot += 1

            maxp = np.max(RCD.pairwise)
            psh = RCD.pairwise.shape
            pos = np.zeros((psh[0], psh[1], 2))
            for j in range(RCP.ninputs):
                for k in range(RCP.ninputs):
                    # print(f"{pairwise[i,j]:.3f}", end='  ')
                    pos[j, k, 0] = j + 1
                    pos[j, k, 1] = k + 1

            if parent_figure is None:
                sax = P.axarr[:, i]
                if dBSPL == "Spont":
                    sax[0].set_ylim(0, 0.5)
                else:
                    sax[0].set_ylim(0, 0.8)
                PH.talbotTicks(  # revcorr
                    sax[0],
                    tickPlacesAdd={"x": 0, "y": 0},
                    floatAdd={"x": 1, "y": 1},
                    # pointSize=7,
                )

            if parent_figure is None:
                PH.noaxes(P.axarr[1, i])
                s2 = sax[2]

                s2.plot(
                    RCD.sites,
                    RCD.participation / RCD.nspikes,
                    "kx",
                    markersize=5,
                    clip_on=False,
                )
                s2.set_ylim((0, 1.0))
                s2.set_xlim(0, 225)

                if i == ncells - 1:
                    s2.set_xlabel("# Release Sites")
                # sax[2].set_ylabel("Participation")
                PH.talbotTicks(  # Participation
                    s2,
                    tickPlacesAdd={"x": 0, "y": 2},
                    floatAdd={"x": 0, "y": 2},
                    # pointSize=7,
                )
                s2.tick_params(direction="in", length=3.0, width=1.0)

            if parent_figure is None:
                s3 = sax[3]
                s3.plot(
                    np.arange(len(RCD.ynspike)) + 1,
                    RCD.ynspike,
                    "k^-",
                    markersize=5,
                    clip_on=False,
                )
                s3.set_ylim(0, 1.05)
                if i == ncells - 1:
                    s3.set_xlabel(
                        f"# Inputs in [{RCD.pre_w[0]:.1f} to {RCD.pre_w[1]:.1f}] before spike"
                    )
                PH.talbotTicks(  # cumulative inputs before spike
                    s3,
                    tickPlacesAdd={"x": 0, "y": 1},
                    floatAdd={"x": 0, "y": 2},
                    # pointSize=7,
                )
                s3.tick_params(direction="in", length=3.0, width=1.0)

        # save the accumulated RCD data
        if save_calcs:
            with open(rc_datafile, "wb") as fh:
                pickle.dump(all_RCD_RCP, fh)
        P.figure_handle.suptitle("")  # remove title.
        title = (
            "SBEM Project Supplemental Figure 4 Modeling : Reverse Correlation Summary",
        )
        save_file = set_figure_path(fignum=4, filedescriptor=f"Revcorr_Summary_Full_{dBSPL:s}", suppnum=99)
#       save_file = Path("Figure4_supp", f"Fig_M4_supplemental_Full_{dBSPL:s}.pdf")
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        # fig.title["title"] = title
        return fig

    def plot_revcorr_compare(  # x)
        self,
        parent_figure=None,
        axlist: Union[list, None] = None,
        dBSPLs: list = ["Spont", "30dB"],
        legend=True,
    ):
        if parent_figure is None:
            P = PH.regular_grid(
                1,
                2,
                figsize=(7, 4),
                margins={
                    "bottommargin": 0.15,
                    "leftmargin": 0.15,
                    "rightmargin": 0.15,
                    "topmargin": 0.15,
                },
                verticalspacing=0.03,
                horizontalspacing=0.1,
            )
            axes = P.axarr[0, :]
        else:
            if axlist is None:
                raise ValueError(
                    "plot_revcorr_compare: If plotting to a parent figure, axlist must be specified"
                )
            else:
                axes = axlist

        syms = ["^-", "o-"]
        for i, dBSPL in enumerate(dBSPLs):
            if len(axes) > 1:
                ax = axes[i]
            else:
                ax = axes[0]
            with open(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl", "rb") as fh:
                try:
                    R = FPM.pickle_load(fh)
                except:
                    print(
                        "You must run revcorr_supplement plot first, then the file we need will be present"
                    )
                    return
            hpoints = []
            if i == 0:
                fillstyle = "none"
            else:
                fillstyle = "full"
            r_keys = list(R.keys())
            sns.set_palette("tab10", len(r_keys))
            palette = sns.color_palette(palette="tab10", n_colors=len(r_keys))
            for c, cell_number in enumerate(r_keys):
                RCD = R[cell_number][0]
                RCP = R[cell_number][1]
                xnspike = np.arange(len(RCD.ynspike)) + 1
                ax.plot(
                    xnspike,
                    RCD.ynspike,
                    syms[i],
                    color=palette[c],  # force palette colors
                    markersize=3.5,
                    fillstyle=fillstyle,
                    clip_on=False,
                    label=f"{FD.BC_name:s}{cell_number:02d}",
                )
                if i == 0:
                    hpoints.extend(
                        [[xnspike[k], RCD.ynspike[k]] for k in range(len(xnspike))]
                    )
            # ax.set_xlim(0, 12)
            # ax.set_ylim(0, 1.0)
            ax.set_xlabel("Number of inputs prior to spike")
            ax.set_ylabel("Cumulative Fraction of Spikes")

            if i == 0 and legend:
                ax.legend()
            if len(axes) == 1:  # superimposed, so just show the symbols
                import alphashape
                import matplotlib.patches
                from descartes import PolygonPatch
                from matplotlib.lines import Line2D

                hidden_lines = [
                    Line2D(
                        [0],
                        [0],
                        color="gray",
                        marker="^",
                        markerfacecolor="none",
                        markersize=4,
                        lw=1,
                    ),
                    Line2D([0], [0], color="gray", marker="o", markersize=4, lw=1),
                ]
                ax.legend(hidden_lines, dBSPLs)
                hpoints = np.reshape(hpoints, (len(hpoints), 2))
                # ax.plot(hpoints[:,0], hpoints[:,1], 'o')
                # mpl.show()

                if i == 0:  # draw a alphashape hull around the data points
                    alphafactor = 0.4
                    poly = alphashape.alphashape(hpoints, alphafactor)
                    ax.add_patch(
                        PolygonPatch(poly, fc="grey", ec="darkgrey", alpha=0.3)
                    )

        if parent_figure is None:
            save_file = Path(
                self.config["figureIntermediateDirectory"],
                f"Fig_M4_Revcorr_Compare.pdf",
            )
            fig = FigInfo()
            fig.P = P
            fig.filename = save_file
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 Modeling : Reverse Correlation Comparison"
            return fig
        else:
            return None

    def Figure6_Main(self, parent_figure=None):
        example_cell_number = 17
        lh = 0.5
        hsp = 0.7
        xw = 2.25
        rh = lh + xw + hsp
        col3 = rh + xw + hsp
        col4 = col3 + xw + hsp
        yht = 0.8
        ypos = [0.5+0.75+0.15 + i*(yht+0.4) for i in range(6)]
        yhtr = 2.0
        yposr = [0.5 + i*(yhtr+0.82) for i in range(3)]
        sizer = OrderedDict(  # define figure layout
            [
                ("A", {"pos": [lh, xw, ypos[5], yht], "labelpos": [-0.1, 1.05]}),
                ("B", {"pos": [lh, xw, ypos[4], yht], "labelpos": [-0.1, 1.05]}),
                ("C", {"pos": [lh, xw, ypos[3], yht], "labelpos": [-0.1, 1.05]}),
                ("D", {"pos": [lh, xw, ypos[2], yht], "labelpos": [-0.1, 1.05]}),
                ("E", {"pos": [lh, xw, ypos[1], yht], "labelpos": [-0.1, 1.05]}),
                ("F", {"pos": [lh, xw, ypos[0], yht], "labelpos": [-0.1, 1.05]}),
                ("G", {"pos": [lh, xw, 0.5, 0.5], "labelpos": [-0.1, 1.05]}),
                ("H", {"pos": [rh, xw, ypos[5], yht], "labelpos": [-0.1, 1.05]}),
                ("I", {"pos": [rh, xw, ypos[4], yht], "labelpos": [-0.1, 1.05]}),
                ("J", {"pos": [rh, xw, ypos[3], yht], "labelpos": [-0.1, 1.05]}),
                ("K", {"pos": [rh, xw, ypos[2], yht], "labelpos": [-0.1, 1.05]}),
                ("L", {"pos": [rh, xw, ypos[1], yht], "labelpos": [-0.1, 1.05]}),
                ("M", {"pos": [rh, xw, ypos[0], yht], "labelpos": [-0.1, 1.05]}),
                ("N", {"pos": [rh, xw, 0.5, 0.5], "labelpos": [-0.1, 1.05]}),
                # right side, summary plots
                ("O1", {"pos": [col3, xw, yposr[2], yhtr], "labelpos": [-0.1, 1.05]}),
                ("O2", {"pos": [col4, xw, yposr[2], yhtr], "labelpos": [-0.1, 1.05]}),
                ("O3", {"pos": [col3, xw, yposr[1], yhtr], "labelpos": [-0.1, 1.05]}),
                ("O4", {"pos": [col4, xw, yposr[1], yhtr], "labelpos": [-0.1, 1.05]}),
                ("P1", {"pos": [col3, xw, yposr[0], yhtr], "labelpos": [-0.1, 1.05]}),
                ("P2", {"pos": [col4, xw, yposr[0], yhtr], "labelpos": [-0.1, 1.05]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.

        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            label=True,
            figsize=(12.0, 9.0),
            units="inches",
            showgrid=False, 
            # labelposition=(-0.05, 1.02),
        )
        
        # title the 2 left columns
        mpl.text(x=1.5, y=8.3, s=f"{FD.BC_name:s}{example_cell_number:02d}  200Hz 100% SAM", fontdict={
                "fontsize": 10, "fontweight": "bold", "ha": "center"},
                transform=P.figure_handle.dpi_scale_trans)
        mpl.text(x=4.5, y=8.3, s=f"{FD.BC_name:s}{example_cell_number:02d}  60 Hz Click Train", fontdict={
                "fontsize": 10, "fontweight": "bold", "ha": "center"},
                transform=P.figure_handle.dpi_scale_trans)

        # Column 1: SAM

        label_font = {"fontsize": 8, "fontweight": "normal"}
        title_font = {"fontsize": 9, "fontweight": "normal"}
        P.axdict["A"].set_ylabel("mV", fontdict=label_font)

        P.axdict["B"].set_title(f"{FD.BC_name:s} Spike Raster", fontdict=title_font)
        P.axdict["B"].set_ylabel("Trial", fontdict=label_font)

        P.axdict["C"].set_title(
            "Bushy PSTH", fontdict=title_font, verticalalignment="bottom", y=0.95
        )
        P.axdict["C"].set_ylabel("Spikes/second", fontdict=label_font)

        P.axdict["D"].set_title("ANF Spike Raster", fontdict=title_font)
        P.axdict["D"].set_ylabel("Trial", fontdict=label_font)

        P.axdict["E"].set_title(
            "ANF PSTH", fontdict=title_font, verticalalignment="top", y=0.95
        )
        P.axdict["E"].set_ylabel("Spikes/second", fontdict=label_font)

        P.axdict["F"].set_xlabel("Angle (radians)", fontdict=label_font)
        P.axdict["F"].set_ylabel("Spike Count", fontdict=label_font)
    
        # P.axdict["G"].set_title("Stimulus", fontdict=title_font)
        P.axdict["G"].set_ylabel("Stimulus", fontdict=label_font)
        P.axdict["G"].set_xlabel("Time (s)", fontdict=label_font)

        # column 2 Click train
        P.axdict["H"].set_ylabel("mV", fontdict=label_font)

        P.axdict["I"].set_title(f"{FD.BC_name:s} Spike Raster", fontdict=title_font)
        P.axdict["I"].set_ylabel("Trial", fontdict=label_font)

        P.axdict["J"].set_title(
            f"{FD.BC_name:s} PSTH", fontdict=title_font, verticalalignment="bottom", y=0.95
        )
        P.axdict["J"].set_ylabel("Spikes/second", fontdict=label_font)

        P.axdict["K"].set_title("ANF Spike Raster", fontdict=title_font)
        P.axdict["K"].set_ylabel("Trial", fontdict=label_font)

        P.axdict["L"].set_title(
            "ANF PSTH", fontdict=title_font, verticalalignment="top", y=0.95
        )
        P.axdict["L"].set_ylabel("Spikes/second", fontdict=label_font)

        P.axdict["M"].set_ylabel("CI", fontdict=label_font)
        P.axdict["M"].set_title("SAC", fontdict=title_font, verticalalignment="top", y=0.95)
        P.axdict["M"].set_ylim(0, 25)
        PH.talbotTicks(P.axdict["M"], axes='y', density=(1.0, 1.0), 
            tickPlacesAdd= {'x':2, "y": 0})
        # P.axdict["N"].set_title("Stimulus", fontdict=title_font)
        P.axdict["N"].set_ylabel("Stimulus", fontdict=label_font)
        P.axdict["N"].set_xlabel("Time (s)", fontdict=label_font)

        for axl in [
            "B",
            "C",
            "D",
            "E",
            "G",
        ]:
            P.axdict[axl].sharex(P.axdict["A"])

        for axl in [
            "I",
            "J",
            "K",
            "L",
            "N",
        ]:
            P.axdict[axl].sharex(P.axdict["H"])


        self.Figure6_one_column(
            mode = "SAM",
            cell_number = example_cell_number,
            dataset = FD.figure_SAM_SAC,
            P = P,
            pan = ["A", "B", "C", "D", "E", "F", "G"],
        )
        self.Figure6_one_column(
            mode = "SAC",
            cell_number = example_cell_number,
            dataset = FD.figure_SAM_SAC,
            P = P,
            pan = ["H", "I", "J", "K", "L", "M", "N"],
        )
        
        """ Now do the right side with the VS plots and "V" or line plots

        """
        inset_type = "rMTF"
        VSP15 = SAM_VS_vplots.VS_Plots(dBSPL=15)
        VSP15.plot_VS_Data(axin=P.axdict["P1"], legendflag=True)
        VSP30 = SAM_VS_vplots.VS_Plots(dBSPL=30)
        VSP30.plot_VS_Data(axin=P.axdict["P2"], legendflag=False)

        VSP = SAM_VS_vplots.VS_Plots(dBSPL=15)
        VSP.plot_VS_summary(2, axin=P.axdict["O1"], legendflag=False, inset_type=inset_type, keep_inset=False)
        VSP.plot_VS_summary(30, axin=P.axdict["O2"], legendflag=False, inset_type=inset_type, keep_inset=False)
        VSP.plot_VS_summary(9, axin=P.axdict["O3"], legendflag=False, inset_type=inset_type, keep_inset=False)
        VSP.plot_VS_summary(17, axin=P.axdict["O4"], legendflag=False, inset_type=inset_type, keep_inset=False)


        fig = FigInfo()
        if parent_figure is not None:
            fig.P = parent_figure
        else:
            fig.P = P
        filedescriptor = "Ephys_4_main_v5"

        filedescriptor += f"_{inset_type:s}"
        fig.filename = set_figure_path(fignum=6, filedescriptor=filedescriptor)
        fig.title["title"] = "SBEM Project Figure 6 Modeling: SAM, SAC, rMTF, Entrainment"
        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        fig.title2 = title2
        return fig

    def Figure6_one_column(
        self, mode: str, cell_number: int, dataset: dict, P: object, pan: object
    ):

        PD = self.newPData()
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )
        # print('dataset: ', dataset)
        # print('keys: ', dataset.keys())
        # print('mode: ', mode)
        # print('cellpath: ', cellpath)
        # print("cell_number: ", cell_number)
        # print(dataset.keys())
        # print('dataset[mode]: ', dataset[cell_number][mode])

        sfi = Path(cellpath, Path(dataset[cell_number][mode]).name)
        if not sfi.is_dir():
            print(f"File is not a directory: {str(sfi):s}")
            return

        fn = sorted(list(sfi.glob("*")))[0]
        MD = self.ReadModel.get_data(fn, PD)
        # mtime = Path(fn).stat().st_mtime
        # timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
        #     "%Y-%m-%d-%H:%M"
        # )
        if not MD.success:
            print("No simulation found that matched conditions")
            print("File was: ", fn)
            return
        # unpack x
        # par, stitle, ivdatafile, filemode, d = X
        axon_name = MD.SI.axonExpt
        AR = MD.AR
        SP = MD.SP
        RM = MD.RM
        protocol = "runANPSTH"
        # AR, SP, RM = analyze_data.analyze_data(ivdatafile, filemode, protocol)
        i = 0
        plot_win = (0.4, 0.450)
        psth_binw = 0.0001
        i = 0  # Voltage for first trial
        self.plot_voltage(
            ax=P.axdict[pan[0]], ntrace=i, d=MD.data, AR=AR, time_win=plot_win
        )
        self.plot_stim_waveform(
            ax=P.axdict[pan[6]], ntrace=i, d=MD.data, AR=AR, stim_win=plot_win
        )
        all_bu_st = self.get_bu_spikearray(AR, MD.data)
        PF.plot_spiketrain_raster(
            all_bu_st, ax=P.axdict[pan[1]], max_raster=10, plot_win=plot_win
        )
        self.plot_psth_psth(
            ax=P.axdict[pan[2]],
            data=all_bu_st,
            ri=MD.RI,
            psth_binw=psth_binw,
            ntr=len(all_bu_st),
            psth_win=plot_win,
        )

        an_st_by_input, all_an_st, an_st_grand = self.get_an_spikearray(AR, MD.data)
        ninputs = len(an_st_by_input)
        SC, syninfo = self.parent.PLT.get_synaptic_info(cell_number)
        PF.plot_stacked_spiketrain_rasters(
            an_st_by_input, ax=P.axdict[pan[3]], si=MD.SI, syninfo=syninfo, plot_win=plot_win,
            linewidth=2.0,
        )
        # P.axdict[pan[3]].set_xlabel("Time (s)")

        self.plot_psth_psth(
            ax=P.axdict[pan[4]],
            data=all_an_st,
            ri=MD.RI,
            psth_binw=psth_binw,
            psth_win=plot_win,
            ntr=len(all_an_st),
            ninputs=ninputs,
        )
        P.axdict[pan[4]].set_xlabel("Time (s)")
        P.axdict[pan[4]].set_title(
            "AN PSTH",
            y=1.0,
            fontdict={
                "fontsize": 9,
                "fontweight": "normal",
                "verticalalignment": "top",
            },
        )
        ri = MD.RI
        si = MD.SI
        (
            totaldur,
            soundtype,
            pip_start,
            pip_duration,
            F0,
            dB,
            fmod,
            dmod,
        ) = self.parent.PLT.get_stim_info(si, ri)
        pip_duration = ri.pip_duration
        print("Soundtype: ", soundtype)
        if soundtype.endswith("Clicks"):
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
            S = SAC.SAC()
            sac_engine = "cython"
            bu_sac, bu_sacbins = S.SAC_with_histo(
                all_bu_st,
                pars=spars,
                engine=sac_engine,
                dither=1e-3 * si.dtIC / 2.0,
            )
            # plot bu as red
            P.axdict[pan[5]].plot(
                bu_sacbins[:-1] * 1e3,
                bu_sac,
                "r-",
                # label=sac_label,
            )
            an_sac, an_sacbins = S.SAC_with_histo(
                an_st_grand,
                pars=spars,
                engine=sac_engine,
                dither=1e-3 * si.dtIC / 2.0,
            )
            P.axdict[pan[5]].plot(
                an_sacbins[:-1] * 1e3,
                an_sac,
                "k-",
                # label=sac_label,
            )
            custom_legend = [Line2D([0], [0], marker=None, color="k", lw=1, label='AN'),
                             Line2D([0], [0], marker=None, color="r", lw=1, label=f"{FD.BC_name:s}"),
                    ]
            P.axdict[pan[5]].legend(handles=custom_legend, handlelength=1, loc="upper right", fontsize=7, labelspacing=0.33, markerscale=0.5)
            P.axdict[pan[5]].set_xlabel("Time (ms)")
        else:
            phasewin = [
                pip_start + 0.25 * pip_duration,
                pip_start + pip_duration,
            ]
            y = []
            for z in all_bu_st:
                y.extend(z)
            window_duration = phasewin[1] - phasewin[0]
            x = np.array(y)
            v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
            nreps = len(AR.MC.traces)
            bu_spikesinwin = x[v]
            vs_bu = self.parent.PLT.VS.vector_strength(bu_spikesinwin, fmod, time_window=phasewin, nreps=nreps)
            vs_an = self.parent.PLT.VS.vector_strength(all_an_st, fmod, time_window=phasewin, nreps=nreps)
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
            # plot AN in red
            print("sync AN")
            PF.plot_psth(
                vs_an.circ_phase,
                run_info=ri,
                zero_time=0.0,
                max_time=2 * np.pi,
                bin_width=est_binw,
                ax=P.axdict[pan[5]],
                bin_fill=True,
                xunits="radians",
                edge_color="k",
                alpha=0.6,
            )
            # plot BU black
            print("sync BU")
            PF.plot_psth(
                vs_bu.circ_phase,
                run_info=ri,
                zero_time=0.0,
                max_time=2 * np.pi,
                bin_width=est_binw,
                bin_fill=True,
                ax=P.axdict[pan[5]],
                xunits="radians",
                edge_color='r',
                alpha=0.6,
            )
            # P.axdict["E"].hist(
            #     vs["ph"],
            #     bins=2 * np.pi * np.arange(30) / 30.0,
            #     facecolor="k",
            #     edgecolor="k",
            # )
            P.axdict[pan[5]].set_xlim((0.0, 2 * np.pi))
            custom_legend = [Line2D([0], [0], marker=None, color="k", lw=3, alpha=0.5, label=f"AN VS = {vs_an.vs:5.3f}"),
                             Line2D([0], [0], marker=None, color="r", lw=3, alpha=0.5, label=f"{FD.BC_name:s} VS = {vs_bu.vs:5.3f}"),
                    ]
            P.axdict[pan[5]].legend(handles=custom_legend, handlelength=1, loc="upper left", fontsize=7, labelspacing=0.33, markerscale=0.5)
            # P.axdict[pan[5]].text(x=0.05, y=1.0,
            #     s=f"VS: AN = {vs_an.vs:5.3f}\n    BU = {vs_bu.vs:5.3f}",
            #     fontdict={
            #         "fontsize": 8,
            #         "fontweight": "normal",
            #         "verticalalignment": "top",
            #     },
            #     transform=P.axdict[pan[5]].transAxes,
            # )

    def Figure6_Supplemental2(self):
        V = SAM_VS_vplots.VS_Plots()
        #fig, P = V.make_figure()
        fig, P = V.Figure6_Supplemental2()
        return fig
    
    def Figure6_Supplemental3(self):
        V = SAM_VS_vplots.VS_Plots()
        #fig, P = V.make_figure()
        fig, P = V.Figure6_Supplemental3()
        return fig
    
    def Figure6_Supplemental4(self):
        fig = SACP.plot_sacs(figinfo=FigInfo(show_figure_name=False))
        return fig

    def plot_psth_psth(
        self,
        ax: object,
        data: object,
        ri: dict,
        psth_binw: float = 0.001,
        psth_win: Union[list, np.array] = [0.0, 1.0],
        ntr: int = 1,
        ninputs: int = 1,
    ):
        PF.plot_psth(
            data,
            run_info=ri,
            zero_time=psth_win[0],
            max_time=psth_win[1],
            bin_width=psth_binw,
            edge_color=None,
            face_color="k",
            ax=ax,
            scale=1.0 / ninputs,
        )
        ax.set_xlim(0, np.fabs(np.diff(psth_win)))
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(2, 1.0),
            insideMargin=0.02,
            # pointSize=ticklabelsize,
            tickPlacesAdd={"x": 2, "y": 0},
            floatAdd={"x": 2, "y": 0},
        )
        ax.set_ylabel("Sp/s")

    def plot_psth_ANFSL(
        self,
        ax: object,
        wstart: float,
        cell_number: int,
        ri: dict,
        an_st_grand: object,
        fsl_win: Union[None, tuple] = None,
        label_x_axis=True,
    ):

        grand_fsl, grand_ssl = PF.plot_fsl_ssl(
            an_st_grand,
            run_info=ri,
            max_time=25.0,
            fsl_win=fsl_win,
            bin_width=0.25,
            ax=ax,
            zero_time=wstart * 1e-3,
            offset=0,  # (ninputs)*dy,
            cellID=cell_number,
            show_values=False,
        )
        ax.set_ylabel("Spikes")
        if label_x_axis:
            ax.set_xlabel("Latency (ms)")
        fsl_text = f"FSL: {np.nanmean(grand_fsl):.3f} (SD {np.nanstd(grand_fsl):.3f})"
        fsl_text += (
            f"\nSSL: {np.nanmean(grand_ssl):.3f} (SD {np.nanstd(grand_ssl):.3f})"
        )
        ax.text(
            0.18,
            0.95,
            fsl_text,
            # N={np.count_nonzero(~np.isnan(fsl)):3d})",
            fontsize=5,
            color="b",
            # fontfamily="monospace",
            transform=ax.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
        )

        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.02,
            # pointSize=ticklabelsize,
            tickPlacesAdd={"x": 0, "y": 0},
            floatAdd={"x": 0, "y": 0},
        )
        ax.set_title("AN")

    def plot_All_PSTH(self):
        for cell in GRPDEF.grAList():
            fig = self.plot_PSTH(cellN=cell)
        return fig

    def plot_PSTH(self, cellN=None):
        dBSPL = "30dB"
        if cellN is None:
            cell_number = 17
        else:
            cell_number = cellN
        print(f"Plotting PSTH for {FD.BC_name:s}{str(cell_number):s}")
        box_size = 0.32
        box_size_x = 0.45
        sizer = {
            "A": {  # trace
                "pos": [0.1, box_size_x, 0.7, box_size],
                "labelpos": (-0.15, 0.7),
                "noaxes": True,
            },
            "B": {
                "pos": [0.1, box_size_x, 0.65, 0.05],
                "labelpos": (-0.15, 1.02),
            },  # stimulus
            "C": {"pos": [0.1, box_size_x, 0.36, 0.24], "labelpos": (-0.15, 1.02)},
            "D": {"pos": [0.1, box_size_x, 0.1, 0.18], "labelpos": (-0.15, 1.02)},
            "E": {"pos": [0.65, 0.25, 0.55, box_size], "labelpos": (-0.15, 1.02)},
            "F": {"pos": [0.65, 0.25, 0.10, box_size], "labelpos": (-0.15, 1.02)},
        }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))

        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            figsize=(6, 6),
            label=True,
            # showgrid=True,
            # verticalspacing=0.12,
            # horizontalspacing=0.12,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.1,
                "rightmargin": 0.1,
                "topmargin": 0.1,
            },
            # fontsize={"tick": 7, "label": 9, "panel": 12},
            # fontweight={"tick": "normal", "label": "normal", "panel": "bold"},
        )
        tr_ax = P.axdict["A"]
        st_ax = P.axdict["B"]
        bupsth_ax = P.axdict["C"]
        anpsth_ax = P.axdict["D"]
        bufsl_ax = P.axdict["E"]
        anfsl_ax = P.axdict["F"]

        PH.cleanAxes(P.axarr.ravel())

        bu_fsl_win = (2.7e-3, 4.5e-3)
        an_fsl_win = (1.5e-3, 3.5e-3)

        self.plot_one_PSTH(
            cell_number=cell_number,
            dBSPL=dBSPL,
            tr_ax=tr_ax,
            st_ax=st_ax,
            bupsth_ax=bupsth_ax,
            anpsth_ax=anpsth_ax,
            psth_win=(0.15, 0.4),  # default is 0.15, 0.4 (start, end)
            bufsl_ax=bufsl_ax,
            anfsl_ax=anfsl_ax,
            bu_fsl_win=bu_fsl_win,
            an_fsl_win=an_fsl_win,
        )

        if cellN is None:
            save_file = Path(self.config["figureIntermediateDirectory"], f"Fig_M5.pdf")
        else:
            save_file = Path(
                self.config["figureIntermediateDirectory"],
                f"All_PSTH/PSTH_VCN_c{cell_number:02d}.png",
            )
        title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.01}
        title = "SBEM Project Figure 3 Modeling : PSTH Summary"
        # save_file = f"Fig_M5_supplemental_Full_{dBSPL:s}.pdf"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        fig.title2 = title2
        return fig

    def plot_voltage(
        self,
        ax: object,
        d: object = None,
        ntrace: int = 0,
        AR: object = None,
        time_win: tuple = (0, 0.25),
        cal_x_axis: bool = False,
    ):
        """
        Plot the voltage of a trace into an axis, with the spikes marked
        """
        vtrial = AR.MC.traces[ntrace] * 1e3
        time_base = AR.MC.time_base / 1000.0  # convert to seconds
        dt = d["Params"].dtIC / 1000.0  # convert from msec to seconds
        trd = d["Results"][ntrace]
        tb_beg = int(time_win[0] / dt)
        tb_end = int(time_win[1] / dt)

        # cprint('r', f"i is 0, psth: {psth_win=}")
        ax.plot((time_base - time_win[0]), np.array(vtrial), "k-", linewidth=0.5)
        spiketimes = np.array(trd["spikeTimes"])
        # trim spike mark array so we don't get spots at the edge of the plot
        spikeindex = [
            int(t / dt) for t in spiketimes if (t >= time_win[0] and t < (time_win[1]))
        ]
        ax.plot(
            (time_base[spikeindex] - time_win[0]),
            vtrial[spikeindex],
            "ro",
            markersize=1.5,
        )
        PH.referenceline(ax, -60.0)
        # if cal_x_axis:
        #     ax.text(
        #         0.0,
        #         -60.0,
        #         "-60 mV",
        #         fontsize=7,
        #         horizontalalignment="right",
        #         verticalalignment="center",
        #     )
        ax.set_xlim(0.0, time_win[1] - time_win[0])
        if cal_x_axis:
            PH.noaxes(ax)
            PH.calbar(
                ax,
                calbar=[0.180, -20, 0.020, 10.0],
                scale=[1, 1.0],
                axesoff=True,
                orient="right",
                unitNames={"x": "ms", "y": "mV"},
                fontsize=8,
                weight="normal",
                color="k",
                font="Arial",
            )

    def plot_stim_waveform(
        self,
        ax,
        ntrace: int = 0,
        d: object = None,
        AR: object = None,
        si: object = None,
        stim_win: tuple = (0, 0.25),
        color='b',
        scale:str='sec', # or "ms"
        y_label:bool=False,  # turn off mPa
    ):
        # stimulus waveform
        trd = d["Results"][ntrace]
        waveform = np.array(trd["stimWaveform"].tolist())
        stb = trd["stimTimebase"]
        timescale = 1.0
        if scale in ['ms', 'msec']:
            timescale = 1e3  # rescale so we are in milliseconds

        stimdt = np.mean(np.diff(stb))
        sttb_beg = int(stim_win[0] / stimdt)
        sttb_end = int(stim_win[1] / stimdt)
        # print('STIMWIN: ', stim_win)
        # print(stb[sttb_beg], stb[sttb_end])
        ax.plot(
            np.array(stb[sttb_beg:sttb_end] - stim_win[0])*timescale,
            waveform[sttb_beg:sttb_end] * 1e3,
            linestyle='-',
            color=color,
            linewidth=0.5,
        )  # stimulus underneath
        ax.set_xlim(0, timescale*(stim_win[1] - stim_win[0]))
        PH.nice_plot(ax, position={'bottom': -0.03}, direction='outward', ticklength=3)

        if y_label:
            ax.set_ylabel("mPa")
        else:
            tl = ax.get_yticklabels()
            ax.set_yticklabels([""]*len(tl))

    
    def get_bu_spikearray(self, AR, d):
        """
        Get the arrays of auditory nerve inputs, sorted by input, all, and grand...

        Parameters
        ----------
        AR : object
         acq4 read object
        d : object
         data object from runmodel
        """
        all_bu_st = []
        si = d["Params"]
        ntr = len(AR.MC.traces)
        for i in range(ntr):  # for all trials in the measure.
            trd = d["Results"][i]
            if not isinstance(trd["spikeTimes"], list):
                cprint("r", "spiketimes is not a list")
                return
            all_bu_st.append(trd["spikeTimes"])

        return all_bu_st

    def get_an_spikearray(self, AR, d):
        """
        Get the arrays of auditory nerve inputs, sorted by input, all, and grand...

        Parameters
        ----------
        AR : object
            acq4 read object
        d : object
            data object from runmodel
        """
        # Collapse the input data
        si = d["Params"]
        ntr = len(AR.MC.traces)
        all_an_st = []  # collapsed in trial order only
        an_st_grand = [[] for x in range(ntr)]  # accumulate across trials

        for i in range(ntr):  # for all trials in the measure.
            trd = d["Results"][i]  # trial i data
            ninputs = len(trd["inputSpikeTimes"])  # get AN spike time array (by input)
            if i == 0:
                an_st_by_input = [[] for x in range(ninputs)]
            for k in range(ninputs):  # raster of input spikes
                tk = trd["inputSpikeTimes"][k]
                all_an_st.extend(np.array(tk) * 1e-3)  # strictly for the psth
                an_st_by_input[k].append(
                    tk * 1e-3
                )  # index by input, then trials for that input
                an_st_grand[i].extend(tk * 1e-3)
        return an_st_by_input, all_an_st, an_st_grand


    def get_cell_PSTH_data(self, dataset, cell_number, dBSPL:float):

        PD = self.newPData()
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )
        sfi = Path(cellpath, Path(dataset[dBSPL]).name)
        if not sfi.is_dir():
            print("file not found: ", str(sfi))
            return None, None, None

        fn = sorted(list(sfi.glob("*")))[0]
        X = self.ReadModel.get_data_file(fn, PD)
        # mtime = Path(fn).stat().st_mtime
        # timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
        #     "%Y-%m-%d-%H:%M"
        # )
        if X is None:
            print("No simulation found that matched conditions")
            print(fn)
            return None, None, None
        # unpack x
        par, stitle, ivdatafile, filemode, d = X
        axon_name = par["axonExpt"]
        protocol = "runANPSTH"
        AR, SP, RM = analyze_data.analyze_data(
            ivdatafile=ivdatafile, filemode=filemode, protocol=protocol
        )
        self.all_bu_st = self.get_bu_spikearray(AR, d)
        return d, AR, axon_name

    def plot_one_PSTH(
        self,
        cell_number: int,
        dBSPL: float,
        tr_ax: object,
        st_ax: object,
        bupsth_ax: object,
        anpsth_ax: object,
        bufsl_ax: object,
        anfsl_ax: object,
        psth_win: tuple = (0, 0.25),  # (0.15, 0.3), # [start, end]
        plot_win: tuple = (0.4, 0.550),
        bu_fsl_win: Union[None, tuple] = None,
        an_fsl_win: Union[None, tuple] = None,
        stim_tr_ax:Union[None, object] = None,
        stim_psth_ax:Union[None, object] = None,
        stim_fsl_ax:Union[None, object] = None,
        label_x_axis=True,
    ):
        dataset = FD.figure_psth[cell_number]
        d, AR, axon_name = self.get_cell_PSTH_data(dataset, cell_number, dBSPL = dBSPL)
        if d is None:
            return None
        ntr = len(AR.MC.traces)  # number of trials
        waveform = None
        v0 = -160.0
        trstep = 25.0 / ntr
        inpstep = 5.0 / ntr
        sz = 50.0 / ntr
        psth_dur = psth_win[1] - psth_win[0]
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
        ) = self.parent.PLT.get_stim_info(si, ri)

        ntr = len(AR.MC.traces)
        #
        # Traces
 
        for i in range(ntr):  # for all trials in the measure.
            vtrial = AR.MC.traces[i] * 1e3
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            dt = si.dtIC / 1000.0  # convert from msec to seconds
            trd = d["Results"][i]
            tb_beg = int(psth_win[0] / dt)
            tb_end = int(psth_win[1] / dt)
            if i == 0:  # Voltage for first trial
                self.plot_voltage(ax=tr_ax, ntrace=i, d=d, AR=AR, time_win=psth_win)
            if i == 0 and waveform is not None and st_ax is not None:
                self.plot_stim_waveform(
                    ax=st_ax, ntrace=i, d=d, AR=AR, stim_win=plot_win
                )
        PH.nice_plot(tr_ax, direction="outward", ticklength=3.0)
        if label_x_axis:
            xticks_str = ["0.0", "0.1", "0.2"]
        else:
            xticks_str = []
        PH.set_axes_ticks(tr_ax,
                yticks = [-60, -30, 0],
                yticks_str = ["-60", "-30", "0"],
                y_minor=[-70, -50, -40, -20, -10],
                xticks = [0, 0.1, 0.2],
                xticks_str = xticks_str,
                x_minor = [0.05, 0.15, 0.25],
            )   
        tr_ax.set_ylim(-75, 5)
        self.stim_data = {'ntrace': i, 'd': d, 'AR': AR, 'stim_win': psth_win}
        if stim_tr_ax is not None:
            self.plot_stim_waveform(
                    ax=stim_tr_ax, ntrace=i, d=d, AR=AR, stim_win=psth_win
                )
            PH.set_axes_ticks(stim_tr_ax,
                yticks = [-1, -0, 1],
                yticks_str = [], # ["-1", "", "1"],
                xticks = [0, 0.1, 0.2],
                xticks_str = ["0.0", "0.1", "0.2"],
                x_minor = [0.05, 0.15, 0.25],
            )
            stim_tr_ax.set_xlabel("Time (s)")
        if label_x_axis:
            tr_ax.set_xlabel("Time (s)")    

        #
        # PSTH


        psth_binw = 0.5e-3
        ninputs = 1

        if bupsth_ax is not None:
            self.plot_psth_psth(
                ax=bupsth_ax,
                data=self.all_bu_st,
                ri=ri,
                psth_binw=psth_binw,
                psth_win=psth_win,
                ntr=ntr,
                ninputs=1,
            )
            PH.nice_plot(bupsth_ax, direction="outward", ticklength=3.0)
            if label_x_axis:
                xticks_str = ["0", "0.1", "0.2"]
            else:
                xticks_str = []
            PH.set_axes_ticks(bupsth_ax,
                yticks = [0,  500,  1000, 1500],
                yticks_str = ["0", "500",  "1000", "1500"],
                xticks = [0, 0.1, 0.2],
                xticks_str = xticks_str,
                x_minor = [0.05, 0.15, 0.25],
            )            
            if stim_psth_ax is not None:
                self.plot_stim_waveform(
                        ax=stim_psth_ax, ntrace=i, d=d, AR=AR, stim_win=psth_win
                    )
                PH.set_axes_ticks(stim_psth_ax,
                    yticks = [-1, -0, 1],
                    yticks_str = [], # ["-1", "", "1"],
                    xticks = [0, 0.1, 0.2],
                    xticks_str = ["0", "0.1", "0.2"],
                    x_minor = [0.05, 0.15, 0.25],
                )
                stim_psth_ax.set_xlabel("Time (s)")
            if label_x_axis:
                bupsth_ax.set_xlabel("Time (s)")

        an_st_by_input, all_an_st, an_st_grand = self.get_an_spikearray(AR, d)

        #
        # ANPSTH
        if anpsth_ax is not None:
            print("anpsth to plotpsth")
            self.plot_psth_psth(
                ax=anpsth_ax,
                data=all_an_st,
                ri=ri,
                psth_binw=psth_binw,
                psth_win=psth_win,
                ntr=ntr,
                ninputs=ninputs,
            )
            if label_x_axis:
                anpsth_ax.set_xlabel("Time (sec)")
            anpsth_ax.set_title("AN")
    #
    # FSL
    #
        if bufsl_ax is not None:
            PF.plot_fsl_ssl(
                self.all_bu_st,
                run_info=ri,
                max_time=25.0,
                bin_width=0.25,
                fsl_win=bu_fsl_win,
                ax=bufsl_ax,
                zero_time=psth_win[0] * 1e-3,
                cellID=cell_number,
            )
            PH.nice_plot(bufsl_ax, direction="outward", ticklength=3.0)
            if label_x_axis:
                xticks_str = ["0", "5", "10", "15", "20", "25"]
            else:
                xticks_str = []
            PH.set_axes_ticks(bufsl_ax,
                yticks = [0, 20, 40, 60],
                yticks_str = ["0", "20", "40", "60"],
                y_minor = [10, 30, 50],
                xticks = [0, 5, 10, 15, 20, 25],
                xticks_str = xticks_str,
            )
            bufsl_ax.set_ylim(0, 60)
            bufsl_ax.set_xlim(0, 25)
            bufsl_ax.set_ylabel("Spikes")

            if stim_fsl_ax is not None:
                self.plot_stim_waveform(
                        ax=stim_fsl_ax, ntrace=i, d=d, AR=AR, stim_win=[0.2, 0.225], scale='ms', color='b',
                    )
                stim_fsl_ax.set_xlabel("Time (ms)")
                PH.set_axes_ticks(stim_fsl_ax,
                    yticks = [-1, -0, 1],
                    yticks_str = [], # ["-1", "", "1"],
                    xticks = [0, 5, 10, 15, 20, 25],
                    xticks_str = ["0", "5", "10", "15", "20", "25"],
                )
                stim_fsl_ax.set_xlim(0, 25)


            if label_x_axis:
                bufsl_ax.set_xlabel("Time (ms)")

        # a raster plot
        # y = 0.
        # cyc = cycler(color='bgrcmyk')
        # anfsl_ax.set_prop_cycle(cyc)
        # for i in range(ntr):
        #     trd = d["Results"][i]
        #     ninputs = len(trd["inputSpikeTimes"])
        #     for k in range(ninputs):  # raster of input spikes
        #         tk = trd["inputSpikeTimes"][k]
        #         y += 1.0
        #         anfsl_ax.plot(
        #                 tk*1e-3, y+np.zeros_like(tk), "|", markersize=2., linewidth=0.5
        #             )
        #     y += 10
        # anfsl_ax.plot([0.200,0.200], [0, y], 'r-')
        # anfsl_ax.set_xlim(0, 1)

        dy = 50
        # for k in range(ninputs):
        #
        #     fsl, ssl = PF.plot_fsl_ssl(an_st_by_input[k],
        #          run_info=ri,
        #          max_time = 25.0,
        #          min_fsl = 1.0e-3,
        #          bin_width= 0.25,
        #          ax = ax,
        #          zero_time =  wstart*1e-3,
        #          offset = k*dy,
        #          cellID = cell_number,
        #          show_values = False,
        #          )

        if anfsl_ax is not None:
            print("figures:anfslwin: ", an_fsl_win)
            self.plot_psth_ANFSL(
                ax=anfsl_ax,
                wstart=psth_win[0],
                cell_number=cell_number,
                ri=ri,
                an_st_grand=an_st_grand,
                label_x_axis=label_x_axis,
                fsl_win=an_fsl_win,
            )

        return axon_name

    def plot_ISI(
            self,
            cell_number: int,
            isi_win: tuple = (0.05, 0.2),
            binwidth: float=1e-4,
            isi_ax: object=None,
        ):
        if self.all_bu_st is None or isi_ax is None:
            return
        cvisit, cvisi, cvt, cvm, cvs = ISI.isi_cv(
            self.all_bu_st,
            binwidth=binwidth,
            reftime=0.0,
            t0=isi_win[0],
            t1=isi_win[1],
        )
        

    def plot_one_CV(
        self,
        cell_number: int,
        cv_win: Union[list, tuple] = (0.220, 0.3),
        reftime: float = 0.0,
        t_grace: float = 0.0,
        cv_binw: float = 0.001,
        cv_ax: object = None,
        label_x_axis: bool = False,
        stim_ax: Union[object, None] = None,
    ):
        print("CVwin, ref_time: ", cv_win, reftime)# use the data from PSTH. If not present, skip the CV plot
        if self.all_bu_st is None or cv_ax is None:
            return
        cvisit, cvisi, cvt, cvm, cvs = ISI.isi_cv(
            self.all_bu_st,
            binwidth=cv_binw,
            reftime=0.0,
            t0=cv_win[0],
            t1=cv_win[1],
            tgrace=t_grace,
        )
        print('cvt min/max: ', np.min(cvt), np.max(cvt))
        cv_ax.plot((cvt - reftime) * 1e3, cvs / cvm, "k-")

        PH.nice_plot(cv_ax, direction="outward", ticklength=3.0)
        cv_ax.set_ylim(0, 1.2)
        cv_ax.set_xlim(0, 1e3 * (cv_win[1] - cv_win[0]) - 20.0)
        if stim_ax is not None:
            # print("Plotting CV waveform")
            self.plot_stim_waveform(
                ax=stim_ax, ntrace=self.stim_data["ntrace"], 
                d=self.stim_data["d"], 
                AR=self.stim_data["AR"],
                stim_win=(0.200, 0.3),
                scale='ms',
            )
            stim_ax.set_ylim(-1, 1)
            stim_ax.set_xlim(0, 1e3 * (cv_win[1] - cv_win[0]) - 20.0)
            stim_ax.set_xlabel("Time (ms)")
            # PH.nice_plot(stim_ax, direction="outward", ticklength=3)
            PH.set_axes_ticks(stim_ax,
                yticks = [-1, -0, 1],
                yticks_str = [], # ["-1", "", "1"],
                xticks = [0, 20, 40, 60, 80],
                xticks_str = ["0", "20", "40", "60", "80"],
                x_minor=[10, 30, 50, 70],
            )
            stim_ax.set_xlabel("Time (s)")

        # compute the mean CV from 10-60 msec for display on the plot
        itwin = np.where(((cvt - reftime) > 0.01) & ((cvt - reftime) < 0.06))[0]
        CVp = np.nanmean(cvs[itwin] / cvm[itwin])
        if label_x_axis:
            xticks_str = ["0", "20", "40", "60", "80"]
        else:
            xticks_str = []       
        PH.set_axes_ticks(cv_ax,
                yticks = [0, 0.4, 0.8, 1.2],
                yticks_str = ["0.0", "0.4", "0.8","1.2"],
                y_minor=[0.2, 0.6, 1.0],
                xticks = [0, 20, 40, 60, 80],
                xticks_str = xticks_str,
                x_minor=[10, 30, 50, 70],
            ) 
        cv_ax.set_ylim(0, 1.2)

        cvlabel = "" # r"$CV\prime$: "
        cv_ax.text(
            1.02,
            CVp,
            f"{cvlabel:s}{CVp:4.2f}",
            transform=cv_ax.transAxes,
            fontsize=7,
            horizontalalignment="left",
            verticalalignment="center",
        )
        cv_ax.set_ylabel("CV")
        if label_x_axis:
            cv_ax.set_xlabel("Time (ms)")

    def Figure4_Supplemental4_PSTH(self):
        print("Plotting Figure 4 Supplement 4 PSTH")
        dBSPL = "30dB"
        lmar = 0.125
        rmar = 0.1
        hspc = 0.08
        pght = 10

        P = PH.regular_grid(
            rows=len(GRPDEF.grAList()),
            cols=4,
            order="rowsfirst",
            figsize=(8, pght),
            # showgrid=True,
            margins={
                "bottommargin": 0.08,
                "leftmargin": lmar,
                "rightmargin": rmar,
                "topmargin": 0.1,
            },
            verticalspacing=0.03,
            horizontalspacing=hspc,
        )
        P2 = PH.regular_grid(
            rows=1,
            cols=4,
            order="rowsfirst",
            figsize=(8, pght),
            margins = {
                "bottommargin": 0.06,
                "leftmargin": lmar,
                "rightmargin": rmar,
                "topmargin": 0.93,
            },
            horizontalspacing=hspc,
            verticalspacing=0.0,
            parent_figure=P,
        )
        for c in range(4):
            P.axarr[-1, c].xaxis.set_ticklabels([])
            P2.axarr[0, c].yaxis.set_ticklabels([])

        PH.cleanAxes(P.axarr.ravel())
        ncells = len(GRPDEF.grAList())
        # center labels over columns
        cdata = {
            0: f"Soma Voltage",
            1: f"PSTH",
            2: f"FSL/SSL",
            3: f"CV",
        }
        for i in list(cdata.keys()):
            P.axarr[0, i].set_title(
                cdata[i],
                y=1.1,  # raise up
                transform=P.axarr[0, i].transAxes,
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="bottom",
            )

        tr_ax = P.axdict["A"]
        st_ax = P.axdict["B"]
        bupsth_ax = P.axdict["C"]
        anpsth_ax = P.axdict["D"]
        bufsl_ax = P.axdict["E"]
        anfsl_ax = P.axdict["F"]

        PH.cleanAxes(P.axarr.ravel())

        for i, cell_number in enumerate(GRPDEF.grAList()):
            print("Plotting psth for cell: ", cell_number)
            show_label = False
            if i == ncells - 1:
                # show_label = True  # put label into lower axis
                stim_ax = P2.axarr[0,:]
                stim_ax_cv = P2.axarr[0,3]
            else:
                stim_ax = [None, None, None]
                stim_ax_cv = None
            axon_name = self.plot_one_PSTH(
                cell_number=cell_number,
                dBSPL=dBSPL,
                psth_win=(0.15, 0.4),
                tr_ax=P.axarr[i, 0],
                st_ax=None,
                bupsth_ax=P.axarr[i, 1],
                anpsth_ax=None,
                bufsl_ax=P.axarr[i, 2],
                anfsl_ax=None,
                label_x_axis=show_label,
                bu_fsl_win=(2.2e-3, 4e-3),  # FSL window after onset, in msec
                stim_tr_ax=stim_ax[0],
                stim_psth_ax=stim_ax[1],
                stim_fsl_ax=stim_ax[2]

            )
            self.plot_one_CV(
                cell_number=cell_number,
                cv_win=(0.2, 0.300),
                reftime=0.2,
                cv_ax=P.axarr[i, 3],
                label_x_axis=show_label,
                t_grace=0.0,
                stim_ax=stim_ax_cv,
            )
           # Shim CV slightly to the right so that label does not overlap text in FSL?SSL_ERROR_EOF
            pos = P.axarr[i, 3].get_position()
            pos.x0 = pos.x0 + 0.02
            P.axarr[i, 3].set_position(pos)
            # also shi the stim plot
            if stim_ax_cv is not None:
                pos = stim_ax_cv.get_position()
                pos.x0 = pos.x0 + 0.02
                stim_ax_cv.set_position(pos)

            P.axarr[i, 0].text(
                -0.40,
                0.5,
                f"{FD.BC_name:s}{cell_number:02d}",
                fontsize=9,
                color="k",
                transform=P.axarr[i, 0].transAxes,
                horizontalalignment="right",
            )
            if axon_name == "standardized":
                P.axarr[i, 0].text(
                    -0.40,
                    0.3,
                    "Sub. axon",
                    fontsize=7,
                    color="k",
                    transform=P.axarr[i, 0].transAxes,
                    horizontalalignment="right",
                )

        save_file = set_figure_path(fignum=4, filedescriptor="PSTH_V2", suppnum=4)
        title = "SBEM Project Figure 4 Modeling Supplemental4 : PSTH and FSL, all grade A cells"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        return fig

    def plot_VS_SAM_15_BC09(self):
        self.generate_VS_data_file(dB=15, bc09=True)
    
    def plot_VS_SAM_15(self):
        self.generate_VS_data_file(dB=15)
    
    def plot_VS_SAM_30(self):
        self.generate_VS_data_file(dB=30)
    
    def analyze_VS_data1(
        self, VS_data, cell_number, fout, firstline=False, sac_flag=False, testmode=False,
        dBSPL:int=0, make_VS_raw:bool=True,
    ):
        """
        Generate tables of Vector Strength measures for all cells
        across the frequencies listed
        Include the rate MTF as well

        """
        
        # self.parent.PLT.textclear()  # just at start
        PD = self.newPData()
        P = None
        linesout = ""
        print("Cell number: ", cell_number)
        if isinstance(cell_number, str):
            cell_n = int(cell_number[0])
        else:
            cell_n = cell_number
        # self.parent.cellID = cell_number
    
        for i, filename in enumerate(VS_data.samdata[cell_number]):
            newline = self.analyze_VS_data_single(filename, cell_number, fout, 
                                        PD=PD, linesout=linesout, firstline=firstline,
                                        sac_flag=sac_flag, testmode=testmode,
                                        dBSPL = dBSPL, make_VS_raw=make_VS_raw)
            linesout += newline
        return linesout
    
    def analyze_VS_data_single(self, filename, cell_number, fout, PD, linesout:str, firstline=False,
             sac_flag=False, testmode=False,
                dBSPL:int=0, make_VS_raw:bool=True,
             ):
        print(f"Cell: {str(cell_number):s}  Filename: {filename:s}")
        i_cell_number = int(cell_number[0])
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{i_cell_number:02d}",
            "Simulations",
            "AN",
        )
        sfi = Path(cellpath, filename + ".pkl")
        print(f"Opening pkl file: {str(sfi):s}")
        if not testmode:
            with open(sfi, "rb") as fh:
                d = FPM.pickle_load(fh)
            VS_colnames, VS_line = self.parent.PLT.plot_AN_response(
                None, d.files[0], PD, "runANPSTH", sac_flag=sac_flag,
                filename=filename, make_VS_raw=make_VS_raw,
                cell_ID=cell_number,
            )
            # note that VSline has results of VS computation to put in the table
            if firstline:
                linesout = VS_colnames # self.parent.PLT.VS_colnames
                linesout += "\n"
                firstline = False
            linesout += VS_line # self.parent.PLT.VS_line
            # linesout += "\n"
            # with open(fout, "a") as fth:
            #     if firstline:
            #         fth.write(VS_colnames)
            #         fth.write("\n")
            #         firstline = False
            #     fth.write(VS_line)
            #     fth.write("\n")
        else:
            print(f"Testing only... , file found: {str(sfi.is_file()):s}")
            if firstline:
                linesout += f"Cell: {str(cell_number):s}  Filename: {filename:s}\n"
                linesout += "Data would be here\n"
                firstline = False
            else:
                linesout +=  f"Cell: {str(cell_number):s}  Filename: {filename:s}\n"
        
        print(f"Analyze VS data completed for cell: {str(cell_number):s}  file: {filename:s}")
        return linesout
    
    def _write_VS_Header(self, fout:Path, timestamp_str:str):
        with open(fout, "w") as fh:
            fh.write(f'"""\n')
            fh.write(
                "    Vector strength for models with SAM tones, different input configurations.\n"
            )
            fh.write("    17 Aug 2021 version (added rMTF data 4 April 2023).\n")
            fh.write(
                "    Results are a printout from DataTablesVCN after selecting the data runs.\n"
            )
            fh.write(f"Run started at: {timestamp_str:s}\n")
            fh.write(
                "WARNING: This table is automatically written by figures.py generate_VS_data_file\n"
            )
            fh.write("       and should not be directly edited.\n")
            fh.write(f"To Regenerate:\n   After running the simulations, enter the filenames into VS_datasets_xxdB.py\n"
            )
            fh.write(f"   You can select the datasets, and then click the 'Print File Info' button for each cell.\n")
            fh.write(
                f"    Copy the text in the 'Reports' dock in DataTablesVCN\n"
            )
            fh.write(
                f"  into a 'VS_datasets_xxdB.py' file, where xx is the sound pressure level.\n"
            )
            fh.write(
                f"Then select 'Analyze VS-SAM Table @ xxdB DPL' in DataTables Figures dropdown, and click 'Create Figure.\n"
            )
            fh.write(
                f"  No figure will be generated, but vector strength will be calculated\n"
                
            )
            fh.write(f"   and the VS_data_xxdB (timestanp).py file will be created.\n")
            fh.write(
                f"The VS_data_xxdB.py file holds all of the vector-strength information, in a text format,\n "
            )
            fh.write(f"  and is read by the plotting programs.\n")
            fh.write(f'--pbm 2014-2023\n"""\n')  # end of the comment
            fh.write('\ndata = """')  # start of the data

    def _write_VS_tail(self, fout:Path, end_timestamp:str, run_time_seconds: float):
        with open(fout, "a") as fh:
            fh.write('\n"""\n')
            fh.write(f"Run ended at: {end_timestamp:s}\n")
            fh.write(f"Run duration (seconds): {run_time_seconds:9.1f}\n")
            fh.write('\n"""')

    def generate_VS_data_file(self, testmode=False, dB=15, bc09=False, append=False):
        """
        Write the VS_data file from the selected datasets.
        This routine reanalyzes all of the data in the table
        to get the vector strength information.

        The datasets that to be analyzed are listed in VS_datasets_nndB.py 
        or VS_datasets_15dB_BC09_NoUninnervated.py (if bc09 is True) 
        """
        parallel = True
        testmode = False

        start_timestamp = datetime.datetime.now()
        timestamp_str = start_timestamp.strftime("%Y-%m-%d-%H.%M.%S")

        if dB == 30:
            if f"VS_datasets_{dB:d}dB" not in list(dir()):
                import VS_datasets_30dB as VS_datasets
                outfile = f"VS_data_{dB:d}dB_{timestamp_str:s}.py"
                TASKS = [(c, d) for c in GRPDEF.grAList() for d in VS_datasets.samdata[c]]

        if dB == 15 and not bc09:
            if f"VS_datasets_{dB:d}dB" not in list(dir()):
                import VS_datasets_15dB as VS_datasets
                outfile = f"VS_data_{dB:d}dB_{timestamp_str:s}.py"
                TASKS = [(c, d) for c in GRPDEF.grAList() for d in VS_datasets.samdata[c]]
                # TASKS = [s for s in GRPDEF.grAList()]

        elif dB == 15 and bc09:
            if f"VS_datasets_{dB:d}dB_BC09_NoUninnervated" not in list(dir()):
                import VS_datasets_15dB_BC09_NoUninnervated as VS_datasets
                outfile = f"VS_data_{dB:d}dB_BC09_{timestamp_str:s}.py"
                TASKS = [(c, d) for c in ["9I", "9U"] for d in VS_datasets.samdata[c]]
        # importlib.reload(VS_datasets)  # make sure have the current one
        # print("TASKS: ", TASKS)
        # print(f"Data set keys found: {str(list(VS_datasets.samdata.keys())):s}")
        
        """
        Generate the table in VS_data.py by analyzing the data from 
        VS_datasets.py
        """
        cprint("g", f"Generate VS Data for {dB:d} dB")

        fout = Path(outfile)  # we will generate this automatically, but time stamp it
        write_flag = True
        if not append and write_flag:
            self._write_VS_Header(fout, timestamp_str)

        fl = True
        linesout = ""
        tresults = [None] * len(TASKS)
        results = {}
        # run using pyqtgraph's parallel support
        nWorkers = MPROC.cpu_count()-2
        PD = self.newPData()
        if parallel:
            print("Tasks parallel: \n", TASKS, "-"*80, "\n")
            cprint("m", f"VS_DataAnalysis : Parallel with {nWorkers:d} processes")
            self.parent.PLT.in_Parallel = True  # notify caller
            with MP.Parallelize(
                    enumerate(TASKS), results=tresults, workers=nWorkers
                ) as tasker:
                    for j, t in tasker:
                        if j > 0:
                            fl = False
                        celln = t[0]
                        filename = t[1]
                        print("celln: ", celln)
                        cprint("m", f"Cell: {str(celln):s}  j={j:d}")
                    #    VS_file_raw = f"VS_raw_SAM_{dB:02d}_{celln:02d}.txt"
                        tresults[j] = self.analyze_VS_data_single(filename, celln, fout, PD=PD,
                                                               linesout=linesout,
                            firstline=fl, sac_flag=True, testmode=testmode, dBSPL=dB, make_VS_raw=True,
                            )
                        tasker.results[j] = tresults[j]
            self.parent.PLT.in_Parallel = False

            if write_flag:
                for j in range(len(TASKS)):
                    with open(fout, "a") as fh:
                        if tresults[j] is None:
                            continue
                        lines = tresults[j].split('\n')
                        for l in lines:
                            l = l.replace("\n", "")
                            fh.write(l)
                            fh.write("\n")
        else:
            print("Tasks noparallel: \n", TASKS, "-"*80, "\n")
            for j, t in enumerate(TASKS):
                celln = t[0]
                filename = t[1]
                tresults = self.analyze_VS_data_single(filename, cell_number = celln, fout=fout, PD=PD,
                                                       linesout=linesout,
                    firstline=fl, sac_flag=True, testmode=testmode, dBSPL=dB, make_VS_raw=True,
                    )
                if j > 0:
                    fl = False
                if write_flag:
                    with open(fout, "a") as fh:
                        fh.write(tresults)
                    # fh.write("\n")
        if write_flag:
            with open(fout, "a") as fh:
                fh.write(f'"""\n')  # close the data text.

        # save running time for calculations
        end_timestamp_ts = datetime.datetime.now()
        end_timestamp_str = end_timestamp_ts.strftime("%Y-%m-%d-%H.%M.%S")
        run_time = end_timestamp_ts - start_timestamp
        if write_flag:
            self._write_VS_tail(fout, end_timestamp=end_timestamp_str,
                run_time_seconds = run_time.total_seconds())
            cprint("g", f"The VS_data file {str(fout):s} has been generated.")
        else:
            cprint("r", "NO VS Data file was generated (writeflag is false in Figures.py in generate_VS_data_file )")

    def Figure8_Panels_IJKLM(self):
        """Make the lower panels for Figure 8

        Returns:
            FigData: The figure object
        """        
        rows = 1
        cols = 5
        bmar = 0.5
        tmar = 0.5
        fsize = (8, 2.5)
        panels = ["I", "J", "K", "L", "M"]
        self.P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            units="in",
            figsize=fsize,
            showgrid=False,
            verticalspacing=0.25,
            horizontalspacing=0.35,
            margins={
                "bottommargin": bmar,
                "leftmargin": 0.25,
                "rightmargin": 0.25,
                "topmargin": tmar,
            },

            labelposition=(-0.05, 1.05),
            parent_figure=None,
            panel_labels=panels,
        )
        PD = self.newPData()
        # plot the CCIV response in the first panel
        cellN = 9
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cellN:02d}",
            "Simulations",
            "IV",
        )
        # get the filenames that will be overplotted
        fns = []
        for fd in FD.figure_No_Dend[cellN].keys():
            dfile = FD.figure_No_Dend[cellN][fd]
            sfi = Path(cellpath, Path(dfile).name)
            if not sfi.is_dir():
                print("Missing file: sfi 1: ", sfi)
                return None
            fng = list(sfi.glob("*"))
            if len(fng) == 0:
                raise ValueError("no files found")
            fns.append(fng[0])
        ymin = -140.0
        ymax = 20.0
        iax = 0
        spkmarkcolors = ['c', 'r']
        spkmarksizes = [2.5, 3.5]
        spkmarkshape = ["^", "o"]
        ivaxis = self.P.axdict["J"]
        for i, fn in enumerate(fns):
            sfi = Path(sfi, fn)

            self.parent.PLT.plot_traces(
                self.P.axdict["I"],
                sfi,
                PD,
                protocol="IV",
                ymin=ymin,
                ymax=ymax,
                xmax=150.0,
                iax=i,
                figure=self.P.figure_handle,
                ivaxis=ivaxis,  # accumulate IV's
                ivcolor=colors[i],
                trace_color=colors[i],
                iv_spike_color = spkmarkcolors[i],
                spike_marker_size = spkmarksizes[i],
                spike_marker_color= spkmarkcolors[i],
                spike_marker_shape = spkmarkshape[i],
                calx=100.0,
                caly=-110.0,
                show_title=False,
                axis_index = 1,
            )
        pos = self.P.axdict["I"].get_position()
        pos.y0 = pos.y0 - 0.15
        self.P.axdict["I"].set_position(pos)

        PH.set_axes_ticks(self.parent.PLT.crossed_iv_ax,
                yticks =     [-140,    -120, -100, -80, -60, -40,  -20],
                yticks_str = ["-140", "-120", "-100",   "-80",  "", "-40", "-20 mV"],
                yticks_pad = [-20,      -20,      -20.,     -16.,   4.,   -16.,   -30.],
                xticks = [-1, 0, 1, 2],
                xticks_str = ["-1", "", "1", "2 nA"],
                x_minor = [-0.5, 0.5, 1.5],
            )
        self.parent.PLT.crossed_iv_ax.tick_params(axis='both', which="both", direction="inout")
        pos = ivaxis.get_position()
        pos = pos.translated(-0.03, 0)
        pos.y0 = pos.y0-0.1
        ivaxis.set_position(pos)
        # plot the singles responses in the third panel

        cal_pos = 0
        axl = [self.P.axdict["K"]]
        self.plot_stacked_traces(
            cells=[9],
            dBSPL=30,
            simulation_experiment="NoUninnervated2",
            figure=self.P.figure_handle,
            axes=axl,
            maxstack=4,
            maxtraces=4,
            show_title=False,
            calxp=600.0,
            cal_pos=cal_pos,
        )
        for ax in axl:
            ax.set_zorder(0)
            pos = ax.get_position()
            pos = pos.translated(-0.03, 0)
            pos.y0 = pos.y0 - 0.08
            ax.set_position(pos)

        # plot the efficacy curves and points in the 4th panel
        intact_sym = 'o'
        intact_color = sns_colors[0]
        pruned_sym = '^'
        pruned_color = sns_colors[1]
        EFP = EF.EfficacyPlots(parent_figure=self.P.figure_handle)
        EFP.plot_efficacy(datasetname="NoUninnervated2", ax=self.P.axdict["L"], clean=True, legend=False, marker=pruned_sym, symbol_color=pruned_color)
        EFP.plot_efficacy(datasetname="NoUninnervated2_ctl", ax=self.P.axdict["L"], clean=True, legend=False, marker=intact_sym, symbol_color=intact_color)
        EFP.plot_fits("Full", ax=self.P.axdict["L"])
        PH.set_axes_ticks(ax=self.P.axdict["L"],
            xticks = [0, 100, 200, 300],
            xticks_str = ['0', '100', '200', '300'],
            # xticks_pad:Union[List, None]=None,
            x_minor = [50, 150, 250, 350],
            major_length = 3.0,
            minor_length =1.5,
            yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            yticks_str=["0", "0.2", "0.4", "0.6", "0.8", "1.0"], 
            yticks_pad=[1]*6,
            y_minor=None,
            fontsize=8,
        )
        custom_legend = [
                Line2D([0], [0], marker=intact_sym, color='w', markerfacecolor=intact_color, markersize=5, label=f"{FD.BC_name:s}09 Intact"),
                Line2D([0], [0], marker=pruned_sym, color="w", markerfacecolor=pruned_color, markersize=6, label=f"{FD.BC_name:s}09 Pruned"),
                Line2D([0], [0], color="#FF0000", lw=1, label='Group 1'),
                Line2D([0], [0], color="#94c8ff", lw=1, label='Group 2'),
                ]
        self.P.axdict["L"].legend(handles=custom_legend, handlelength=1, 
            loc="upper left", bbox_to_anchor=(-0.07, 1.0),
            fontsize=7, labelspacing=0.33)
        
        # panel M : compare 9I(intact) and 9U (NoUninnervated) VS across frequencies
        # Plot the Vector Strength for SAM at different frequencies in the 5th panel
        cells = ["9I", "9U"]
        VSP = SAM_VS_vplots.VS_Plots(sels= cells, dBSPL=15, dends="9I9U")
        for i, cell in enumerate(cells):
            legend = False
            if i == 0:
                legend = True
            print("Cell: ", cell)
            VSP.plot_VS_summary(cell, axin=self.P.axdict["M"], legendflag=legend, legend_type="Dendrites", show_cell=True,
                barwidth=180, figure8_xscale=True, inset_type="rMTF")

#         VSP.plot_VS_summary(axin=self.P.axdict["M"], cell=9, 
#                             legendflag=True, legend_type="Dendrites", show_cell=True,
#                 barwidth=75, figure8_xscale=True, inset_type="rMTF", keep_inset=True)
# #                barwidth=75, legendflag=True, xscale="log", figure8_xscale=True)
        self.P.adjust_panel_labels(fontsize=28, fontweight="light", fontname="myriad")
        custom_legend = [Line2D([0], [0], marker="_", markersize=2, color="firebrick", lw=2, label='AN'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor=sns_colors[0], markersize=5, label="Intact"),
                Line2D([0], [0], marker='o', color='w', markerfacecolor=sns_colors[1], markersize=5, label="Pruned"),
                ]
        self.P.axdict["M"].legend(handles=custom_legend, handlelength=1, loc="lower left", fontsize=7, labelspacing=0.33)
        # PH.set_axes_ticks(ax=self.P.axdict["M"],
        #     xticks = [0, 1, 2, 3, 4, 5, 6, 7],
        #     xticks_str = ['50', '100', '200', '300', '400', '500', '750', '1000'],
        #     # xticks_pad:Union[List, None]=None,
        #     x_minor = None,
        #     x_rotation=60.,
        #     major_length = 3.0,
        #     minor_length =1.5,
        #     fontsize=8,
        # )        

        fig = FigInfo()
        fig.P = self.P
        fig.filename = set_figure_path(fignum=8, filedescriptor="Ephys_Pruning_IJKLM_V3")
        fig.title[
            "title"
        ] = "SBEM Project Figure 8 (main) Modeling: Dendrite pruning effects"

        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        # fig.title2 = title2

        return fig


    def plot_VC_gKLT(
        self, parent_figure=None, loc: Union[None, tuple] = None
    ) -> object:
        cell_number = 17
        dataset = FD.figure_VClamp[cell_number]

        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "VC",
        )
        sfi = []
        for i, ds in enumerate(dataset):
            sfd = Path(cellpath, Path(dataset[ds]).name)
            if not sfd.is_dir():
                print("Directory not found: ", str(sfd))
                return
            fn = sorted(list(sfd.glob("*")))[0]
            sfi.append(fn)
        P = self.parent.PLT.plot_VC(
            sfi=sfi, show=False, parent_figure=parent_figure, loc=loc
        )
        save_file = Path(
            self.config["figureIntermediateDirectory"], f"Fig_M0_VC_Adjustment.pdf"
        )
        fig = FigInfo()
        fig.P = P

        if cell_number == None:
            fig.filename = Path(self.config["figureIntermediateDirectory"], "Fig_M0.pdf")
        else:
            fig.filename = save_file
        fig.title["title"] = "SBEM Project Figure Modeling Supplemental : VC"
        fig.title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.05}
        return fig
