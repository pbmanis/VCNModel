import datetime
import importlib
import pickle
from dataclasses import dataclass, field
from pathlib import Path
import string
from typing import Union

import matplotlib
import numpy as np
import pandas as pd
import rich as RI
import seaborn as sns
from cycler import cycler
from matplotlib import image as mpimg
from matplotlib import pyplot as mpl
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from pylibrary.plotting import plothelpers as PH
from pylibrary.tools import cprint as CP
from rich.console import Console
from rich.text import Text

import toml
from vcnmodel.plotters import efficacy_plot as EF
from vcnmodel.plotters import plot_z as PZ
import vcnmodel.util.fixpicklemodule as FPM
from vcnmodel.plotters import (
    figure_data as FD,
)  # table of simulation runs for plotting

config = toml.load(open("wheres_my_data.toml", "r"))
cprint = CP.cprint


def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """
    return [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]


syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x"]
colors = ["c", "k", "m", "r"]

title_text = {
    "passive": "Passive",
    "normal": "Half-active",
    "active": "Active",
    "pasdend": "Passive",
    "actdend": "Active",
}
font_colors = {"passive": "c", "normal": "k", "active": "m"}
spike_colors = {"passive": "c", "normal": "k", "active": "m"}

print(FD)


def get_changetimestamp():
    # trip filemode based on date of simulatoin
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate, "%Y-%m-%d-%H:%M")
    changetimestamp = datetime.datetime.timestamp(dts)
    return changetimestamp


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
    revcorrpath: str = config["revcorrDataDirectory"]
    thiscell: str = ""


def title_data():
    return {"title": "", "x": None, "y": None}


@dataclass
class FigInfo:
    """
    data class to hold figure information for savefig
    """

    P: object = None
    filename: Union[str, Path] = ""
    title: dict = field(default_factory=title_data)
    title2: dict = field(default_factory=title_data)


class Figures(object):
    """
    Entry point.
    This generates final figures for the SBEM manuscript by reaching back to the 
    original simulation data, including, in some cases, refitting.
    Both primary "examplar" figures, and supplemental figures, are generated.
    The figures are made consistent by using both sns.set_style and 
    mpl.style for figures.mplstyle, which overrides some defaults in mpl.
    The resulting figures are editable in Illustrator without any missing fonts.
    
    This is ugly code, but it gets the job done.
    Note that some plots are modifications of what is present in plot_sims,
    and also that plot_sims is used here. 
    This is part of the DataTables interactive analysis tool for the simulations
    in the SBEM project.
    Supported primarily by R01DC015901 (Spirou, Manis, Ellisman),
    Early development: R01 DC004551 (Manis, until 2019)
    Later development: R01 DC019053 (Manis, 2020-2025)
    
    """

    def __init__(self, parent):
        self.parent = parent  # point back to caller's space
        self.config = toml.load(
            open("wheres_my_data.toml", "r")
        )  # sorry, have to reload it here.

    def reset_style(self):
        sns.set_style(rc={"pdf.fonttype": 42})
        mpl.style.use("~/.matplotlib/figures.mplstyle")

    def make_figure(self, figure_name: Union[str, None] = None):
        self.reset_style()
        print("make_figure:", figure_name)
        # dispatch
        dispatch_table = {
            # "Fig 3 Ephys-1 main": self.Figure3Main,
            "Fig 3 Supplemental1ABC_VC-KLTCalibration": self.plot_VC_gKLT,
            "Fig 3 Supplemental1DE_RinTaum": self.plot_combined_VC_IV,
            "Fig 3 Supplemental2_CC": self.plot_IVS,
            "Fig 3 Supplemental3_Zin": self.plot_ZinS,
            "Fig 3 Suppplemental4_PSTH": self.plot_All_PSTH,
            "Fig 4 Ephys-2 main": self.Figure4Main,
            "Fig 4 Ephys-2 Supplement1": self.Figure4_Supplmental1,
            "Fig M1: IV Figure (Fig_M1)": self.plotIV,
            "Fig_IV/ All IVs (Fig_IV/IV_cell_VCN_c##)": self.allIVs,
            "Fig M2: CombinedEffRevCorr (Fig_M2)": self.plot_combined_Eff_Rev,
            "Fig M2B: Efficacy (Fig_M2_Efficacy_Revcorr)": self.plot_efficacy,
            "Fig M2 Supp: Efficacy Supplement (Fig_M2_Supplemental_[experiment])": self.plot_efficacy_supplement,
            "Fig M3: Revcorr Ex (Fig_M3_spont|dBSPL)": self.plot_revcorr,
            "Fig_Revcorr/ All Revcorrs (Revcorr_VCN_c##)": self.plot_all_revcorr,
            "Fig M3 Supp1: Revcorr Supplement Spont (Fig_M3_supplemental_Full_Spont)": self.plot_revcorr_supplement_spont,
            "Fig M3 Supp3: Revcorr Supplement (Fig_M3_supplemental_Full_30dB)": self.plot_revcorr_supplement_30dB,
            "Fig M3 Supp2: Revcorr Supplement (Fig_M3_supplemental_Full_40dB)": self.plot_revcorr_supplement_40dB,
            "Fig M4: Revcorr Compare (Fig_M4_Revcorr_Compare)": self.plot_revcorr_compare,
            "Fig M5: PSTH-FSL (Fig_M5)": self.plot_PSTH,
            "Fig M5 Supp: PSTH-FSL Supplement (Fig_M5_supplemental_Full_dBSPL)": self.plot_PSTH_supplement,
            "VS-SAM Tone (no figure file - analysis only)": self.plot_VS_SAM,
        }
        if figure_name in list(dispatch_table.keys()):
            fig = dispatch_table[figure_name]()
            if fig is not None:
                self.save_figure(fig)

    def save_figure(self, fig):
        # P, save_file, title, title2={"title": None, 'x': 0.0, 'y': 0.0}):

        fig.P.figure_handle.text(
            0.98,
            0.98,
            fig.filename,  # .replace('_', '\_'),
            transform=fig.P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )

        if fig.title2["title"] is not None or len(fig.title2["title"]) > 0:
            fig.P.figure_handle.text(
                fig.title2["x"],
                fig.title2["y"],
                fig.title2["title"],  # .replace('_', '\_'),
                transform=fig.P.figure_handle.transFigure,
                horizontalalignment="right",
                verticalalignment="top",
            )
        ofile = Path(config["baseDataDirectory"], "Figures", fig.filename)
        ofile.parent.mkdir(exist_ok=True)
        print(
            "Saving to: ",
            str(Path(config["baseDataDirectory"], "Figures", fig.filename)),
        )
        print("fig title: ", fig.title["title"])
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", fig.filename),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": fig.title["title"],
            },
        )
        fig.P.figure_handle.show()

    def force_log_ticks(self, ax):
        locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0,), numticks=100)
        ax.xaxis.set_major_locator(locmaj)

        locmin = matplotlib.ticker.LogLocator(
            base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100
        )
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    def allIVs(self):
        for cell in grAList():
            self.plotIV(cell)

    def get_dendmode(self, dendmode: str = ""):
        if dendmode == "passive":
            dendm = "pasdend"
        elif dendmode == "active":
            dendm = "actdend"
        elif dendmode == "normal":
            dendm = "normal"
        else:
            raise ValueError(f"figures.py: dendmode not recognized: {dendmode:s}")
        return dendm

    def plotIV(self, cell=None, parent_figure=None, loc: Union[None, tuple] = None):
        if cell is None:
            cellN = FD.figure_IV["Cell"]
            d1 = FD.figure_IV["normal"]
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

        PD = PData()
        ymin = -125.0
        ymax = 20.0
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
                "bottommargin": 2.25 + bmar,
                "leftmargin": 0.5,
                "rightmargin": 0.5,
                "topmargin": 0.5 + tmar,
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
            horizontalspacing=0.65,
            margins={
                "bottommargin": 0.5 + bmar,
                "leftmargin": 0.5,
                "rightmargin": 0.5,
                "topmargin": 2.0 + tmar,
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
                "topmargin": 2.0 + tmar,
            },
            labelposition=(-0.25, 1.05),
            parent_figure=self.P,
            panel_labels=p3panels,
        )

        for iax, iv in enumerate(["pasdend", "normal", "actdend"]):
            if cell is None:
                dfile = FD.figure_IV[iv]
            else:
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
                iax=iax,
                figure=self.P.figure_handle,
                ivaxis=self.P2.axarr[0, 0],  # accumulate IV's in bottom left plot
                ivcolor=colors[iax],
                calx=120.0,
                caly=-10.0,
            )
            self.P.axarr[0, iax].set_title(
                title_text[iv], color="k", fontweight="normal"
            )

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
        labels = ["Passive", "Half-active", "Active", "Spk Thresh"]
        self.P2.axarr[0, 0].legend(
            lines, labels, bbox_to_anchor=(0.32, 0.35), fontsize=7
        )
        # self.P2.axarr[0, 0].line.set_label(iv)
        # self.P2.axarr[0,0].legend(["passive", "normal", "active", "spike thresholds"])

        res_label = r"$\mathregular{R_{in} (M\Omega)}$"
        tau_label = r"$\mathregular{\tau_{m} (ms)}$"
        phase_label = r"$\mathregular{\phi (radians)}$"

        # plot overlays of all cell z/phase
        for iax, mode in enumerate(["Z_passive", "Z_normal", "Z_active"]):
            if mode not in FD.figure_IV.keys():
                continue
            # cprint('r', f"doing iv: {str(mode):s}")
            sfi = Path(
                config["cellDataDirectory"],
                config["impedanceDirectory"],
                FD.figure_IV[mode],
            )
            if not sfi.is_file():
                cprint("r", f"File not found!!!!!!\n->>> {str(sfi):s}")
                return None
                continue
            # cprint('c', f"sfi Z: {str(sfi):s}")
            # ax = self.P2.axarr[0, 1]
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
            PH.nice_plot(secax)

        # PH.show_figure_grid(self.P.figure_handle)

        # get Rin and RMA from all the examples and make a summary distribution plot
        rins = {}
        taus = {}

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
                    print("sfi is not a dir")
                    continue
                fn = list(sfi.glob("*"))
                sfi = Path(sfi, fn[0])
                AR, SP, RMA = self.parent.PLT.plot_traces(
                    None,  # self.P.axarr[rax, iax+1],
                    sfi,
                    PD,
                    protocol="IV",
                    ymin=ymin,
                    ymax=ymax,
                    iax=iax,
                    figure=None,
                )
                rins[k] = {"Cell": iv, "Rin": RMA["Rin"], "dendrites": dendmode}
                taus[k] = {"Cell": iv, "taum": RMA["taum"], "dendrites": dendmode}
                k += 1

        df_rin = pd.DataFrame.from_dict(rins, orient="index")  # , orient='index')
        df_tau = pd.DataFrame.from_dict(taus, orient="index")

        sns.boxplot(
            data=df_rin,
            x="dendrites",
            y="Rin",
            ax=self.P2.axdict[rinplot],
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_rin,
            x="dendrites",
            y="Rin",
            ax=self.P2.axdict[rinplot],
            color="0.6",
            size=4.0,
            edgecolor="k",
        )
        self.P2.axdict[rinplot].set_ylim(0, 30.0)
        self.P2.axdict[rinplot].set_ylabel(res_label)
        self.P2.axdict[rinplot].set_xlabel("Dendrite Decoration")
        sns.boxplot(
            data=df_tau,
            x="dendrites",
            y="taum",
            ax=self.P2.axdict[tauplot],
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_tau,
            x="dendrites",
            y="taum",
            ax=self.P2.axdict[tauplot],
            color="0.6",
            size=4.0,
            edgecolor="k",
        )
        self.P2.axdict[tauplot].set_ylim(0, 2.5)
        self.P2.axdict[tauplot].set_ylabel(tau_label)
        self.P2.axdict[tauplot].set_xlabel("Dendrite Decoration")

        fig = FigInfo()
        fig.P = self.P
        if cell is None:
            fig.filename = "Fig_M1.pdf"
        else:
            fig.filename = f"Fig_IV/IV_cell_VCN_c{cellN:02d}.pdf"
        fig.title["title"] = "SBEM Project Figure 1 Modeling (Main)"
        # self.save_figure(self.P, save_file, title)
        return fig

    def plot_combined_VC_IV(self):
        P0 = PH.regular_grid(  # dummy figure space
            1,
            1,
            figsize=(8.0, 9.0),
            order="rowsfirst",
            units="in",
            showgrid=False,
            parent_figure=None,
        )

        figp1 = self.plot_VC_gKLT(parent_figure=P0, loc=(0, 8, 4.0, 9.0))
        P2 = self.plotIV(parent_figure=figp1.P, loc=(0, 8, 0.0, 4.0))
        fig = FigInfo()
        fig.P = self.P2
        fig.filename = f"Figure_M0-Combined_Supplemental.pdf"
        fig.title[
            "title"
        ] = "SBEM Project Figure Suupplemental Figure 1 Modeling (Main)"
        return fig

    def plot_IVS(self, parent_figure=None):
        """
        All of the IVS, for a supplemental figure
        Passive, normal, active, plus the crossed IV
        Also put the PNG for the cell on the left.
        """
        nivs = len(FD.figure_AllIVs)
        cprint("c", "plot_IVS.")
        rows = nivs
        cols = 5
        height = 1.5 * nivs
        width = 8.5
        PD = PData()
        ymin = -125.0
        ymax = 40.0
        calx = 120.0

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
        cellpath = config["cellDataDirectory"]
        png_path = Path(config["baseDataDirectory"], config["pngDirectory"])
        cprint("c", "prepping fo run")

        for rax, iv in enumerate(FD.figure_AllIVs.keys()):
            cprint("r", f"Doing Cell VCN_c{iv:02d} -----------------------------------")
            celln = Path(png_path, f"VCN_c{iv:02d}.png")
            if celln.is_file():  # add images from png files
                img = mpimg.imread(str(celln))
                self.P.axarr[rax, 0].imshow(img, aspect="equal")
                ylim = self.P.axarr[rax, 0].get_ylim()
                self.P.axarr[rax, 0].set_xlim(900, 1500)
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
                if rax > 0:
                    calx = None  # only one cal bar on this plot, top row.
                self.parent.PLT.plot_traces(
                    self.P.axarr[rax, iax + 1],
                    sfi,
                    PD,
                    protocol="IV",
                    ymin=ymin,
                    ymax=ymax,
                    iax=iax,
                    figure=self.P.figure_handle,
                    ivaxis=self.P.axarr[rax, 4],  # accumulate IV's in right side
                    ivcolor=colors[iax],
                    iv_spike_color=spike_colors[dendmode],
                    spike_marker_size=1.5,
                    spike_marker_color=spike_colors[dendmode],
                    calx=calx,
                    caly=-10.0,
                )
                if rax == 0:
                    self.P.axarr[rax, iax + 1].set_title(dendmode)
                if iax == 0:
                    self.P.axarr[rax, 0].text(-0.1, 0.5, str(iv))
        if parent_figure is None:
            fig = FigInfo()
            fig.P = self.P
            fig.filename = f"Fig_M1A_Supplemental.pdf"
            timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
            fig.title[
                "title"
            ] = f"SBEM Project Figure 1 Modeling (Supplemental A) ({timestamp_str:s})"
            return fig
        else:
            return self.P

    def plot_ZinS(self):
        """
        All of the Zins for a supplemental figure
        """
        PZ.PlotZ()

    def make_eff_fig(self):
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
        fig.filename = save_file
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
        example = FD.figure_efficacy_supplement[cell_number]

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
        print("plot_efficacy fn: ", fn)
        PD = PData()
        calx = 800.0
        if parent_figure is None:
            parent_figure = self.make_eff_fig()
        EFP = EF.EfficacyPlots(parent_figure=parent_figure)
        EFP.plot_efficacy(datasetname="Full", ax=EFP.parent_figure.axdict["B"], loc=loc)
        # EFP.P.figure_handle.set_size_inches(figsize[0], figsize[1])
        #        return
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
        cmx = mpl.cm.get_cmap(colormap)
        colors = [cmx(float(i) / nfiles) for i in range(nfiles)]

        calx = 800.0
        for n in range(nfiles):
            sfile = Path(sfi, fn[n])
            print(n, sfile)
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
                nax=len(fn),
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
        save_file = "Fig_M2_Efficacy_Revcorr.pdf"
        fig = FigInfo()
        fig.P = EFP.P
        fig.filename = save_file
        fig.title["title"] = "SBEM Project Figure 2 Modeling: Efficacy and Revcorr"
        return fig

    def plot_stacked_traces(
        self,
        cells=None,
        figure=None,
        axes: Union[list, None] = None,
        calxp: float = 800.0,
        calv: float = 20.0,
        maxstack: int = 9,
    ):
        if cells is None:
            cells = grAList()
        trace_ht = 80  # mV

        simulation_experiment = "Full"
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
                Path(FD.figure_efficacy_supplement[cellN][simulation_experiment]).name,
            )
            if not sfiles.is_dir():
                return
            fn = sorted(list(sfiles.glob("*")))
            # print("files to plot: ", fn)
            PD = PData()

            nfn = len(fn)
            ymax = 20.0
            if simulation_experiment != "Full":
                ymax = 40.0
            yoffset = -90.0
            ymin = maxstack * yoffset
            ymax = 20.0
            if simulation_experiment != "Full":
                ymax = 40.0
            for n in range(len(fn)):
                if n == 0 and ic == 0:
                    calxv = calxp
                else:
                    calxv = None
                y0 = n * yoffset
                self.parent.PLT.plot_traces(
                    ax=axes[ic],
                    fn=fn[n],
                    PD=PD,
                    protocol="runANSingles",
                    ymin=ymin,
                    ymax=ymax,
                    xmin=400.0,
                    xmax=900.0,
                    yoffset=n * yoffset,
                    iax=n,
                    nax=len(fn),
                    rep=0,
                    figure=figure,
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
                axes[ic].annotate(
                    text=f"{n+1:d} ",
                    xy=(400.0, y0-60.),  # ms
                    fontsize=8,
                    horizontalalignment="right",
                    verticalalignment="center",
                    # transform=axes[ic].transAxes,
                )

                if n == 0:
                    axes[ic].set_title(
                        f"VCN_c{cellN:02d}",
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
            cells = grAList()
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
            EFP.parent_figure.figure_handle.set_size_inches((8.0, 8.0))
            EFP = EF.EfficacyPlots(parent_figure=self.P)  # , cols=1, rows=1)
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
                Path(FD.figure_efficacy_supplement[cellN][simulation_experiment]).name,
            )
            if not sfiles.is_dir():
                return
            fn = sorted(list(sfiles.glob("*")))
            PD = PData()

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
                print("ymin, ymax: ", ymin, ymax)
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
                        nax=len(fn),
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

        save_file = f"Fig_M2_Supplemental_{simulation_experiment:s}.pdf"
        title = "SBEM Project Figure 2 Modeling : Efficacy, Supplemental"
        save_file = "Fig_M2_Efficacy_Supplement.pdf"
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
            config["baseDataDirectory"],
            config["revcorrDataDirectory"],
            f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl",
        )
        with open(rc_datafile, "rb") as fh:
            d = pickle.load(fh)
        # print(d.keys())
        return d

    def Figure4_Supplmental1(self):
        fig = self.Figure4Main(supplemental1=True)
        return fig

    def Figure4Main(self, supplemental1=False):
        """
        Generate Figure 4 for the paper. Combined bits from various other plots
        and revcorrs/singles
        """
        if not supplemental1:
            example_cells = [5, 9, 17, 30]
        else:
            example_cells = [2, 6, 10, 11, 13, 18]

        start_letter = "A"
        parent_figure = None

        if not supplemental1:
            sizer = {
                "D": {"pos": [6.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                "E": {"pos": [9.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                "F": {"pos": [6.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
                "G": {"pos": [9.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
            }
            figsize = (12,8)
        else:
            sizer = {}
            figsize = (9, 8)
        xw = 1.1
        trace_axes = []
        for j in range(len(example_cells)):
            i = j + 1
            xl = j * 1.25 + 0.75
            axn = f"A{i:d}"
            trace_axes.append(axn)
            sizer[axn] = {
                "pos": [xl, xw, 3.25, 4.25],
                "labelpos": (-0.15, 1.02),
                "noaxes": True,
            }
            sizer[f"B{i:d}"] = {
                "pos": [xl, xw, 2.0, 1.0],
                "labelpos": (-0.15, 1.02),
                # "noaxes": True,
            }
            sizer[f"C{i:d}"] = {
                "pos": [xl, xw, 0.5, 1.0],
                "labelpos": (-0.15, 0.9),
                "noaxes": True,
            }
        # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            units="in",
            figsize=figsize,
            label=True,
            showgrid=False,
            parent_figure=parent_figure,
        )
        # Efficacy plot
        if not supplemental1:
            EFP = EF.EfficacyPlots(parent_figure=P)
            EFP.plot_efficacy("Full", ax=P.axdict["D"], figuremode="clean")
        # participation plots
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
            ax.scatter(a[n][0].sites / synperum2, ap, marker="o", color=color)
            ax.scatter(a[n][0].sites / synperum2, bp, marker="x", color=color)
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 300)
            ax.set_ylim(0, 1.0)
            ax.set_ylabel(f"Participation at 0 and {dB:2d} dBSPL")
            PH.talbotTicks(ax, floatAdd={"x": 0, "y": 2})

        def plot_diff_participation(ax, n, a, b, dB=0, color=None, legend=True):
            ap = a[n][0].participation / a[n][0].npost_spikes
            bp = b[n][0].participation / b[n][0].npost_spikes
            ax.scatter(
                a[n][0].sites / synperum2,
                bp / ap,
                marker="o",
                color=color,
                label=f"VCN_c{n:02d}",
            )
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 300)
            ax.set_ylim(0, 3)
            ax.set_ylabel(f"Participation ratio {dB:2d}/{0:2d} dBSPL")
            PH.talbotTicks(ax, floatAdd={"x": 0, "y": 2})
            if legend:
                ax.legend(fontsize=8, loc="upper right", ncol=2)

        dB = 30
        if not supplemental1:
            ds = self._load_rcdata("Spont")
            drc = self._load_rcdata(f"{dB:2d}dB")
            palette = sns.color_palette(None, len(ds.keys()))
            for i, c in enumerate(ds.keys()):
                # plot_participation(P.axdictax[0], c, ds, drc, dB=dB, color=palette[i])
                plot_diff_participation(
                    P.axdict["E"], c, ds, drc, dB=dB, color=palette[i], legend=False
                )

        axl = [P.axdict[axi] for axi in trace_axes]
        self.plot_stacked_traces(cells=example_cells, figure=P, axes=axl, maxstack=10)
        if not supplemental1:
            self.plot_revcorr_compare(
                parent_figure=P,
                axlist=[P.axdict["F"], P.axdict["G"]],
                dBSPLs=["Spont", "30dB"],
                legend=False,
            )
            synlabel_num = 5
        else:
            synlabel_num = 2
        self.plot_revcorr_supplement(cells=example_cells, parent_figure=P, dBSPL="30dB", synlabel_num=synlabel_num)
        # self.plot_efficacy_supplement(cells=example_cells, parent_figure=P, traces=False)

        for j in range(len(example_cells)):
            ax = P.axdict[f"B{j+1:d}"]
            ax.set_ylim(0, 0.8)
            ax.set_xlim(-5.0, 0.0)

            if j > 0:
                PH.noaxes(ax, whichaxes="y")
            else:
                ax.set_ylabel("Coinc. Rate (Hz)")
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.tick_params(which="major", length=4, direction="in")
            ax.tick_params(which="minor", length=2, direction="in")
        fig = FigInfo()
        if parent_figure is not None:
            fig.P = parent_figure
        else:
            fig.P = P
        if not supplemental1:
            fig.filename = "Figure4_Ephys2_main_v4.pdf"
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 Modeling: Singles, Efficacy and Revcorr"
        else:
            fig.filename = "Figure4-Supplemental1_Revcorr.pdf"
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 Modeling: other cells Singles and Revcorr"

        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        fig.title2 = title2
        print("returnin fig: ", fig)
        return fig

    def plot_all_revcorr(self):
        for cell in grAList():
            fig = self.plot_revcorr(cell)
            if fig is not None:
                self.save_figure(fig)
        return None

    def _get_revcorr(
        self,
        cell_number: int,
        dBSPL="Spont",
        parent_figure: Union[object, None] = None,
        recompute=False,
    ) -> tuple:
        """
        Get the revcorr data associated with the cell number
        and the stimulus level (dbSPL, as a string)
        """
        cprint("c", "Calling _get_revcorr")
        cell_revcorr = FD.figure_revcorr[cell_number]
        run_calcs = False
        rc_datafile = Path(
            PData().revcorrpath, f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl"
        )
        cprint("c", f"Rcdatafile: {str(rc_datafile):s}")
        if rc_datafile.is_file() and not recompute:
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = FPM.pickle_load(fh)
            RCD = all_RCD_RCP[cell_number][0]
            RCP = all_RCD_RCP[cell_number][1]
            PD = PData()
            P = None
        else:
            print(
                "You must run revcorr_supplement plot first, then the data file we need will be present"
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
            (P, PD, RCP, RCD) = self.parent.PLT.compute_revcorr(
                P=None,  # no plotting (and return P is None)
                gbc=str(cell_number),
                fn=fn[0],
                PD=PData(),
                protocol="runANPSTH",
                revcorrtype="RevcorrSPKS",
                thr=-20.0,
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

        P, PD, RCP, RCD = self._get_revcorr(cell_number=cellN, dBSPL="Spont")
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
            P, PD, RCP, RCD, start_letter=p_labels[0], colormap=colormap
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

        cmx = mpl.cm.get_cmap(colormap)
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
        sax3.set_ylabel(f"Cumulative Bushy Spikes\nwitn N AN inputs")
        sax3.set_xlim(1, RCP.ninputs)
        PH.nice_plot(sax3, position=-0.05)
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

        if dBSPL == "Spont":
            if cellN is None:
                save_file = "Fig_M3.pdf"
            else:
                save_file = f"Fig_Revcorr/Revcorr_VCN_c{cell_number:02d}.pdf"
        else:
            save_file = f"Fig_M3_{dBSPL:s}.pdf"
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
        fig = self.plot_revcorr_supplement("Spont")
        return fig

    def plot_revcorr_supplement_30dB(self):
        fig = self.plot_revcorr_supplement("30dB")
        return fig

    def plot_revcorr_supplement_40dB(self):
        fig = self.plot_revcorr_supplement("40dB")
        return fig

    def plot_revcorr_supplement(
        self,
        dBSPL: str,
        cells=None,
        parent_figure=None,
        rate_ax: Union[list, None] = None,
        vm_ax: Union[list, None] = None,
        cumulative: Union[object, None] = None,
        synlabel_num: int=0,
        colormap: str="magma",  # default color map
    ):
        if cells is None:
            cells = grAList()
        ncells = len(cells)

        if parent_figure is None:
            P = PH.regular_grid(
                rows=4,
                cols=len(cells),
                order="columnsfirst",
                figsize=(8, 10),
                # showgrid=True,
                margins={
                    "bottommargin": 0.1,
                    "leftmargin": 0.125,
                    "rightmargin": 0.1,
                    "topmargin": 0.1,
                },
                verticalspacing=0.03,
                horizontalspacing=0.08,
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

        # run_calcs = True
        # P, PD, RCP, RCD = self._get_revcorr(cell_number=cellN, dBSPL = "Spont")
        # rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
        # if rc_datafile.is_file() and not run_calcs:
        #     with open(rc_datafile, "rb") as fh:
        #         all_RCD_RCP = FPM.pickle_load(fh)
        # else:
        #     print(
        #         "Must run revcorr_supplement plot first, then the file we need will be present"
        #     )
        #     run_calcs = False
        all_RCD_RCP = {}  # all revcorr data
        save_calcs = False
        i_plot = 0
        for i, cell_number in enumerate(grAList()):

            PR, PD, RCP, RCD = self._get_revcorr(cell_number=cell_number, dBSPL=dBSPL)
            if PD is None:
                cprint("r", "PD is none in plot_revcorr_supplement")
                continue
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
            print(sfi)
            fn = sorted(list(sfi.glob("*")))
            print("revcorr supplement fn: ", fn)

            # PD = PData()
            # RCD = all_RCD_RCP[cell_number][0]
            # RCP = all_RCD_RCP[cell_number][1]

            if i == 0:
                tcal = True  # point font for cal bar
            else:
                tcal = False  # no cal bar
            if synlabel_num == 0:
                synlabel= True
            else:
                if cell_number == synlabel_num:
                    synlabel = True
                else:
                    synlabel = False

            all_RCD_RCP[cell_number] = [RCD, RCP]
            if cell_number in cells:
                if parent_figure is None:
                    axlist = P.axarr[i_plot, 0:2]
                else:
                    axlist = [
                        P.axdict[f"B{i_plot+1:d}"],
                        P.axdict[f"C{i_plot+1:d}"],
                        None,
                    ]
                summarySiteTC = self.parent.PLT.plot_revcorr2(
                    P,
                    PD,
                    RCP,
                    RCD,
                    axarray=axlist,
                    calbar_show=tcal,
                    calbar_fontsize=7,
                    yaxis_label=False,
                    synlabel=synlabel,
                    colormap=colormap,
                )

                print(f"  Mean pre: {RCD.mean_pre_intervals=}")
                print(f"  Mean Post: {RCD.mean_post_intervals:.3f}")
                if parent_figure is None:
                    P.axarr[i, 0].text(
                        -0.25,
                        0.5,
                        f"VCN_c{cell_number:02d}",
                        fontsize=9,
                        color="k",
                        transform=axlist.transAxes,
                        horizontalalignment="right",
                    )
                else:
                    pass
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
                sax = P.axarr[i, :]
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
                PH.noaxes(P.axarr[i, 1])
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

        title = (
            "SBEM Project Supplemental Figure 3 Modeling : Reverse Correlation Summary",
        )
        save_file = f"Fig_M3_supplemental_Full_{dBSPL:s}.pdf"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        return fig

    def plot_revcorr_compare(
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
            axes = P.axarray[0, :]
        else:
            if axlist is None:
                raise ValueError(
                    "plot_revcorr_compare: If plotting to a parent figure, axlist must be specified"
                )
            else:
                axes = axlist

        for i, dBSPL in enumerate(dBSPLs):
            with open(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl", "rb") as fh:
                try:
                    R = FPM.pickle_load(fh)
                except:
                    print(
                        "You must run revcorr_supplement plot first, then the file we need will be present"
                    )
                    return

            for cell_number in R.keys():
                RCD = R[cell_number][0]
                RCP = R[cell_number][1]
                axes[i].plot(
                    np.arange(len(RCD.ynspike)) + 1,
                    RCD.ynspike,
                    "^-",
                    markersize=4,
                    clip_on=False,
                    label=f"VCN_c{cell_number:02d}",
                )
            axes[i].set_xlim(0, 12)
            axes[i].set_ylim(0, 1.0)
            axes[i].set_xlabel("Number of inputs prior to spike")
            axes[i].set_ylabel("Cumulative Fraction of Spikes")
            if i == 0 and legend:
                axes[i].legend()
            axes[i].set_title(dBSPL)

        if parent_figure is None:
            save_file = f"Fig_M4_Revcorr_Compare.pdf"
            fig = FigInfo()
            fig.P = P
            fig.filename = save_file
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 Modeling : Reverse Correlation Comparison"
            return fig
        else:
            return None

    def plot_psth_psth(
        self,
        ax: object,
        data: object,
        ri: dict,
        time_base: np.ndarray,
        psth_binw: float = 0.5,
        psth_win: Union[list, np.array] = [0.0, 1.0],
        ntr: int = 1,
        ninputs: int = 1,
    ):
        self.parent.PLT.plot_psth(
            data,
            run_info=ri,
            zero_time=psth_win[0],
            max_time=psth_win[1],
            bin_width=psth_binw,
            ax=ax,
            scale=1.0 / ntr / psth_binw / ninputs,
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

        grand_fsl, grand_ssl = self.parent.PLT.plot_fsl_ssl(
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
            0.45,
            0.95,
            fsl_text,
            # N={np.count_nonzero(~np.isnan(fsl)):3d})",
            fontsize=7,
            color="k",
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
        for cell in grAList():
            self.plot_PSTH(cellN=cell)

    def plot_PSTH(self, cellN=None):
        print("PSTH")
        dBSPL = "30dB"
        if cellN is None:
            cell_number = 17
        else:
            cell_number = cellN
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
            save_file = f"Fig_M5.pdf"
        else:
            save_file = f"All_PSTH/PSTH_VCN_c{cell_number:02d}.png"
        title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.01}
        title = ("SBEM Project Figure 5 Modeling : PSTH Summary",)
        # save_file = f"Fig_M5_supplemental_Full_{dBSPL:s}.pdf"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        fig.title2 = title2
        return fig

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
        bu_fsl_win: Union[None, tuple] = None,
        an_fsl_win: Union[None, tuple] = None,
        label_x_axis=True,
    ):
        dataset = FD.figure_psth[cell_number]
        PD = PData()
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )
        sfi = Path(cellpath, Path(dataset[dBSPL]).name)
        if not sfi.is_dir():
            print("file not found: ", str(sfi))
            return

        fn = sorted(list(sfi.glob("*")))[0]
        changetimestamp = get_changetimestamp()
        x = self.parent.PLT.get_data_file(fn, changetimestamp, PD)
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
        axon_name = par["axonExpt"]
        protocol = "runANPSTH"
        AR, SP, RMA = self.parent.PLT.analyze_data(ivdatafile, filemode, protocol)
        ntr = len(AR.MC.traces)  # number of trials
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
        all_bu_st = []

        for i in range(ntr):  # for all trials in the measure.
            vtrial = AR.MC.traces[i] * 1e3
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            dt = si.dtIC / 1000.0  # convert from msec to seconds
            trd = d["Results"][i]
            tb_beg = int(psth_win[0] / dt)
            tb_end = int(psth_win[1] / dt)
            # cprint('r', f"{psth_win=}")
            # cprint('r', f"{tb_beg=}, {tb_end=}")
            ninputs = len(trd["inputSpikeTimes"])
            if i == 0:
                all_an_st = [[] for x in range(ninputs)]  # by input, then trial
                an_st = [[] for x in range(ntr)]  # by trial, collapsed across inputs

            waveform = trd["stimWaveform"].tolist()
            stb = trd["stimTimebase"]
            stimdt = np.mean(np.diff(stb))
            sttb_beg = int(psth_win[0] / stimdt)
            sttb_end = int(psth_win[1] / stimdt)
            if not isinstance(trd["spikeTimes"], list):
                cprint("r", "spiketimes is not a list")
                return
            all_bu_st.append(trd["spikeTimes"])

            if i == 0:  # Voltage for first trial
                # cprint('r', f"i is 0, psth: {psth_win=}")
                tr_ax.plot((time_base - psth_win[0]) * 1e3, vtrial, "k-", linewidth=0.5)
                spiketimes = np.array(trd["spikeTimes"])
                # trim spike mark array so we don't get spots at the edge of the plot
                spikeindex = [
                    int(t / dt)
                    for t in spiketimes
                    if (t >= psth_win[0] and t < (psth_win[1]))
                ]
                tr_ax.plot(
                    (time_base[spikeindex] - psth_win[0]) * 1e3,
                    vtrial[spikeindex],
                    "ro",
                    markersize=1.5,
                )
                PH.referenceline(tr_ax, -60.0)
                if label_x_axis:
                    tr_ax.text(
                        0.0,
                        -60.0,
                        "-60 mV",
                        fontsize=7,
                        horizontalalignment="right",
                        verticalalignment="center",
                    )
                tr_ax.set_xlim(0.0, psth_dur * 1e3)
                PH.noaxes(tr_ax)
                if label_x_axis:
                    PH.calbar(
                        tr_ax,
                        calbar=[180.0, -20, 20.0, 10.0],
                        scale=[1, 1.0],
                        axesoff=True,
                        orient="right",
                        unitNames={"x": "ms", "y": "mV"},
                        fontsize=8,
                        weight="normal",
                        color="k",
                        font="Arial",
                    )

            if i == 0 and waveform is not None and st_ax is not None:
                # stimulus waveform
                st_ax.plot(
                    (stb[sttb_beg:sttb_end] - psth_win[0]),
                    np.array(waveform)[sttb_beg:sttb_end] * 1e3,
                    "k-",
                    linewidth=0.5,
                )  # stimulus underneath
                st_ax.set_xlim(0, psth_dur)
                PH.talbotTicks(
                    st_ax,
                    axes="xy",
                    density=(1.0, 0.5),
                    insideMargin=0.02,
                    # pointSize=ticklabelsize,
                    tickPlacesAdd={"x": 2, "y": 1},
                    floatAdd={"x": 2, "y": 1},
                )
                st_ax.set_ylabel("mPa")

        psth_binw = 0.5e-3
        ninputs = 1

        if bupsth_ax is not None:
            self.plot_psth_psth(
                ax=bupsth_ax,
                data=all_bu_st,
                ri=ri,
                time_base=time_base,
                psth_binw=psth_binw,
                psth_win=psth_win,
                ntr=ntr,
                ninputs=1,
            )
            if label_x_axis:
                bupsth_ax.set_xlabel("Time (s)")
        # Collapse the input data
        an_st_by_input = [[] for x in range(ninputs)]
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

        if anpsth_ax is not None:
            print("anpsth to plotpsth")
            self.plot_psth_psth(
                ax=anpsth_ax,
                data=all_an_st,
                ri=ri,
                time_base=time_base,
                psth_binw=psth_binw,
                psth_win=psth_win,
                ntr=ntr,
                ninputs=ninputs,
            )
            if label_x_axis:
                anpsth_ax.set_xlabel("Time (sec)")
            anpsth_ax.set_title("AN")

        if bufsl_ax is not None:
            self.parent.PLT.plot_fsl_ssl(
                all_bu_st,
                run_info=ri,
                max_time=25.0,
                bin_width=0.25,
                fsl_win=bu_fsl_win,
                ax=bufsl_ax,
                zero_time=psth_win[0] * 1e-3,
                cellID=cell_number,
            )
            PH.talbotTicks(
                bufsl_ax,
                axes="xy",
                density=(1.0, 1.0),
                insideMargin=0.02,
                # pointSize=ticklabelsize,
                tickPlacesAdd={"x": 0, "y": 0},
                floatAdd={"x": 0, "y": 0},
            )
            bufsl_ax.set_ylabel("Spikes")
            if label_x_axis:
                bufsl_ax.set_xlabel("Latency (ms)")

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
        #     fsl, ssl = self.parent.PLT.plot_fsl_ssl(an_st_by_input[k],
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

    def plot_PSTH_supplement(self):
        print("PSTH supplement")
        dBSPL = "30dB"
        lmar = 0.125
        rmar = 0.1
        hspc = 0.08
        P = PH.regular_grid(
            rows=len(grAList()),
            cols=3,
            order="rowsfirst",
            figsize=(8, 10),
            # showgrid=True,
            margins={
                "bottommargin": 0.1,
                "leftmargin": lmar,
                "rightmargin": rmar,
                "topmargin": 0.1,
            },
            verticalspacing=0.03,
            horizontalspacing=hspc,
        )
        PH.cleanAxes(P.axarr.ravel())
        ncells = len(grAList())
        # center labels over columns
        cdata = {
            0: f"Soma Voltage",
            1: f"PSTH",
            2: f"FSL/SSL",
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

        for i, cell_number in enumerate(grAList()):
            print("doing cell: ", cell_number)
            show_label = False
            if i == ncells - 1:
                show_label = True
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
            )

            P.axarr[i, 0].text(
                -0.60,
                0.5,
                f"VCN_c{cell_number:02d}",
                fontsize=9,
                color="k",
                transform=P.axarr[i, 0].transAxes,
                horizontalalignment="left",
            )
            P.axarr[i, 0].text(
                -0.60,
                0.3,
                "\nAxon: " + axon_name,
                fontsize=7,
                color="k",
                transform=P.axarr[i, 0].transAxes,
                horizontalalignment="left",
            )

        save_file = f"Fig_M5_Supplmental_PSTHs.pdf"
        title = "SBEM Project Figure 5 Modeling Supplemental : PSTH and FSL, All cells"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        return fig

    def plot_VS_SAM(self):
        self.generate_VS_data_file()
        pass

    def analyze_VS_data(self, VS_data, cell_number, fout, firstline=False):
        """
        Generate tables of vs measures for all cells
        across the frequencies listed
        """
        self.parent.PLT.textclear()  # just at start
        PD = PData()
        P = None
        self.parent.cellID = cell_number
        for i, filename in enumerate(VS_data.samdata[cell_number]):
            print(f"Cell: {cell_number:d}  Filename: {filename:s}")
            cellpath = Path(
                config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, filename + ".pkl")
            print(f"Opening pkl file: {str(sfi):s}")
            with open(sfi, "rb") as fh:
                d = FPM.pickle_load(fh)

            self.parent.PLT.plot_AN_response(P, d.files[0], PD, "runANPSTH")
            with open(fout, "a") as fth:
                if firstline:
                    fth.write(self.parent.PLT.VS_colnames)
                    fth.write("\n")
                    firstline = False
                fth.write(self.parent.PLT.VS_line)
                fth.write("\n")

        print("Analyze VS data completed for cell: ", cell_number)

    #
    def generate_VS_data_file(self, testmode=False, dB=15, append=True):
        """
        Write the VS_data file from the selected datasets.
        The datasets are in VS_datasets_nndB.py
        """

        if dB == 30:
            if f"vcnmodel.VS_datasets_{db:d}dB" not in list(dir()):
                from vcnmodel import VS_datasets_30dB as VS_datasets
            # print(dir())
            # importlib.reload(VS_datasets)
        if dB == 15:
            if f"vcnmodel.VS_datasets_{dB:d}dB" not in list(dir()):
                from vcnmodel import VS_datasets_15dB as VS_datasets
        importlib.reload(VS_datasets)
        print("VS_datasets: ", VS_datasets)
        print(f"Data set keys found: {str(list(VS_datasets.samdata.keys())):s}")

        config = toml.load(open("wheres_my_data.toml", "r"))

        """
        Generate the table in VS_data.py by analyzing the data from 
        VS_datasets.py
        """
        cprint("r", f"Generate VS Data for {dB:d} dB")

        fout = f"VS_data_{dB:d}dB.py"  # we will generate this automatically
        if not append:
            with open(fout, "w") as fh:
                fh.write(f'"""\n')
                fh.write(
                    "    Vector strength for models with SAM tones, different input configurations.\n"
                )
                fh.write("    17 Aug 2021 version.\n")
                fh.write(
                    "    Results are printout from DataTablesVCN after selecting the data runs.\n"
                )
                fh.write(
                    "NOTE: This table is automatically written by figures.py and should not be\n"
                )
                fh.write("      directly edited.")
                fh.write(f"To Regenerate:\n  After running the simulations, use 'Print File Info' for each cell, and copy the text into a 'VS_datasets_xxdB.py' file")
                fh.write(f"  Then select 'VS-SAMTone-no figure' in DataTables, and 'Create Figure.")
                fh.write(f"  No figure will be generated, but the VS_data_xxdB.py file will be created.")
                fh.write(f"  The VS_data file has all of the vector-strength information, and is read by the plotting programs.")
                fh.write(f'    pbm\n"""\n')
                fh.write('\ndata = """')

        fl = True
        for i, celln in enumerate(grAList()):
            if testmode and i != 1:
                continue
            self.analyze_VS_data(VS_datasets, celln, fout, firstline=fl)
            fl = False
        with open(fout, "a") as fh:
            fh.write(f'"""\n')
        print("VS_Finis")

    def plot_VC_gKLT(self, parent_figure=None, loc: Union[None, tuple] = None):
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
            # print("figures:klt:i, ds: ", i, ds)
            sfd = Path(cellpath, Path(ds).name)
            # print(sfd)
            if not sfd.is_dir():
                print("Directory not found: ", str(sfd))
                return
            fn = sorted(list(sfd.glob("*")))[0]
            sfi.append(fn)
        P = self.parent.PLT.plot_VC(
            sfi=sfi, show=False, parent_figure=parent_figure, loc=loc
        )
        save_file = f"Fig_M0_VC_Adjustment.pdf"
        fig = FigInfo()
        fig.P = P

        if cell_number == None:
            fig.filename = "Fig_M0.pdf"
        else:
            fig.filename = save_file
        fig.title["title"] = "SBEM Project Figure Modeling Supplemental : VC"
        fig.title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.05}
        return fig
