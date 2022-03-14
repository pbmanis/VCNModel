import datetime
import importlib
import pickle
import string
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import matplotlib
import matplotlib.pyplot as mpl
import numpy as np
import pandas as pd
import scipy.spatial
import seaborn as sns
import toml
import vcnmodel.util.fixpicklemodule as FPM
import vcnmodel.util.readmodel as readmodel
from matplotlib import image as mpimg
from matplotlib.lines import Line2D
from pylibrary.plotting import plothelpers as PH
from pylibrary.tools import cprint as CP
from rich.console import Console
from rich.text import Text
from vcnmodel.analyzers import analyze_data
from vcnmodel.analyzers import isi_cv as ISI
from vcnmodel.analyzers import sac as SAC
from vcnmodel.plotters import SAM_VS_vplots
from vcnmodel.plotters import efficacy_plot as EF
from vcnmodel.plotters import (
    figure_data as FD,
)  # table of simulation runs used for plotting figures
from vcnmodel.plotters import plot_z as PZ

from . import plot_functions as PF

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
spike_colors = {"passive": "r", "normal": "r", "active": "r"}

print(FD)


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

    gradeA: list = field(default_factory=grAList)
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
    filename: Union[str, Path] = ""
    title: dict = field(default_factory=title_data)
    title2: dict = field(default_factory=title_data)


class Figures(object):
    """
    This class generates final figures for the SBEM manuscript by reaching back to the
    original simulation data, including, in some cases, refitting.
    Both primary "examplar" figures, and supplemental figures, are generated.
    The figures are made consistent by using both sns.set_style and
    mpl.style for figures.mplstyle, which overrides some defaults in mpl.
    The resulting figures are editable in Illustrator without any missing fonts.

    This is overall ugly code, but it gets the job done.
    Note that some plots are modifications of what is present in plot_sims,
    and also that plot_sims is used here, accessed through self.parent.PLT.
    Some plotting routines have been moved into plot_functions.

    This is part of the DataTablesVCN interactive analysis tool for the simulations
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
        self.axis_offset = -0.02
        config = toml.load(open("wheres_my_data.toml", "r"))
        self.ReadModel = readmodel.ReadModel()
        self.ReadModel.set_parent(
            parent=self.parent, my_parent=self
        )  # pass our parent and US to the reader

    def newPData(self):
        """
        Return Pdata with the paths set from self.config

        Returns:
            PData: dataclass
        """
        return PData(
            basepath=self.config["baseDataDirectory"],
            renderpath=str(Path(self.config["codeDirectory"], "Renderings")),
            revcorrpath=self.config["revcorrDataDirectory"],
        )

    def reset_style(self):
        sns.set_style(rc={"pdf.fonttype": 42})
        mpl.style.use("~/.matplotlib/figures.mplstyle")

    def make_figure(self, figure_name: Union[str, None] = None):
        self.reset_style()
        print("Making_figure:", figure_name)
        # dispatch
        dispatch_table = {
            # "Fig 3 Ephys-1 main": self.Figure3Main,
            "Figure3-Ephys_1_Main": self.Figure3_Main,
            "Figure3-Supplemental1_VC-KLTCalibration": self.plot_VC_gKLT,
            "Figure3-Supplemental1_VC_Rin_Taum": self.Figure3_Supplemental1,
            "Figure3-Supplemental2_CC": self.Figure3_Supplemental2_CC,
            "Figure3-Supplemental3_Zin": self.Figure3_Supplemental3_Zin,
            "Figure3-Supplemental4_PSTH": self.Figure3_Supplemental4_PSTH,
            "Figure4-Ephys_2_Main": self.Figure4_Main,
            "Figure4-Ephys_2_Supplemental1": self.Figure4_Supplmental1,
            "Figure7-Ephys-3 Main": self.Figure7_Main,
            # "Figure7-Ephys_3_Main": self.Figure7_Main,
            "Figure7-Ephys-3_Supplemental1": self.Figure7_Supplemental1,
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
            "VS-SAM Tone (no figure file - analysis only)": self.plot_VS_SAM,
        }
        if figure_name in list(dispatch_table.keys()):
            fig = dispatch_table[figure_name]()
            if fig is not None:
                mpl.show()
                self.save_figure(fig)
        else:
            cprint("r", f"Figure name '{figure_name:s}' was not in dispatch table.")

    def save_figure(self, fig):
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
        ofile = Path(self.config["baseDataDirectory"], "Figures", fig.filename)
        ofile.parent.mkdir(exist_ok=True)
        cprint(
            "g",
            f"Saving to: {str(Path(self.config['baseDataDirectory'])):s}, Figures: {str(fig.filename):s}",
        )
        cprint("g", f"   Figure title: {fig.title['title']:s}")
        mpl.savefig(
            Path(self.config["baseDataDirectory"], "Figures", fig.filename),
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
        for cell in grAList():
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

        """
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
        print("ivaxis: ", ivaxis)
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
                ivaxis=ivaxis,  # accumulate IV's in bottom left plot
                ivcolor=colors[iax],
                calx=120.0,
                caly=-10.0,
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

        if not toponly:
            labels = ["Passive", "Half-active", "Active", "Spk Thresh"]
            self.P2.axarr[0, 0].legend(
                lines, labels, bbox_to_anchor=(0.32, 0.35), fontsize=7
            )
            res_label = r"$\mathregular{R_{in} (M\Omega)}$"
            tau_label = r"$\mathregular{\tau_{m} (ms)}$"
            phase_label = r"$\mathregular{\phi (radians)}$"

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
                PH.nice_plot(secax, position=self.axis_offset, direction="outward")

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

    def Figure3_Main(self):
        message = """
        This Figure was made by using Illustrator to combine parts of other figures/files as follows:
        Panel A came from a syglass rendition of the obj (mesh).
        Panel B (and the inset in C) came from a vispy rendition (neuronvis) of the HOC file.
        Panel C is a combination of Illustrator text, and cnmodel output (see the notebook folder,
        nb/Figure3_main_PanelC.ipynb, for the generation code; the matplotlib windows were then pulled
        into Illustrator to make the figure)
        
        Panels D and E are taken from Figure3_Supplemental2_CC.pdf, for BC 17.
        Panels F, G and H are taken from Figure3_Supplemental4_PSTH.pdf for BC 17.
        """
        cprint("y", message)

    def Figure3_Supplemental1(self):
        """
        Figure 3, Supplemental Figure 1
        Combined voltage clamp traces, IV with Rin, taum plots
        """
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
        fig.filename = f"Figure3/Figure3-Supplemental1_VC_Rin_Taum.pdf"
        fig.title["title"] = "SBEM Project Figure4 Supplemental Figure 1 VC_Rin_Taum"
        return fig

    def Figure3_Supplemental2_CC(self, parent_figure=None, show_pngs=False):
        """
        Plot all of the IVS, for a supplemental figure
        Passive, normal, active, plus the crossed IV
        Also put the PNG for the cell on the left.
        """
        nivs = len(FD.figure_AllIVs)
        cprint("c", "Plotting Figure3_Supplemental2_CC")
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

        for rax, iv in enumerate(FD.figure_AllIVs.keys()):
            # if iv not in [9, 10]:
            #      continue
            cprint(
                "c", f"    Doing Cell BC{iv:02d} -----------------------------------"
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
                if rax > 0:
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
                    iax=iax,
                    figure=self.P.figure_handle,
                    show_title=False,
                    ivaxis=self.P.axarr[rax, jrax],  # accumulate IV's in right side
                    ivcolor=colors[iax],
                    iv_spike_color=spike_colors[dendmode],
                    spike_marker_size=1.5,
                    spike_marker_color=spike_colors[dendmode],
                    calx=calx,
                    caly=-10.0,
                    axis_index=iax,
                )
                if rax == 0:
                    self.P.axarr[rax, jax].set_title(dendmode.title())
                if iax == 0 and not show_pngs:
                    self.P.axarr[rax, 0].text(
                        -35.0, 0.6, f"BC{iv:02d}", fontweight="bold", fontsize=12
                    )
        if parent_figure is None:
            fig = FigInfo()
            fig.P = self.P
            fig.filename = f"Figure3/Figure3_supp/Figure3-Supplemental2_CC.pdf"
            timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
            fig.title[
                "title"
            ] = f"SBEM Project Figure 3 Modeling (Supplemental 2) ({timestamp_str:s})"
            return fig
        else:
            return self.P

    def Figure3_Supplemental3_Zin(self):
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
        PD = self.newPData()
        calx = 800.0
        if parent_figure is None:
            parent_figure = self.make_eff_fig()
        EFP = EF.EfficacyPlots(parent_figure=parent_figure)
        EFP.plot_efficacy(
            datasetname="Full", ax=EFP.parent_figure.axdict["B"], loc=loc, clean=True
        )
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
        """
        What it says: plot traces in a stacked format

        Parameters
        ----------
        cells : a list of ints (if none, default is to use all of the cells in the
            gradeAList defined at the top of this file)
        figure: Figure handle for the plot passed to plot_traces in plot_sims.py)
        axes : a list
             of axes passed to plot_traces in plot_sims.py
        calxp : float
            the position for the calibration bar in x (time), in msec
        calv : float
            the size of the calibration bar in voltage (in mV)
        maxstack: int (default 9)
            determines the height of the stack, relative to the offset (see code)

        Returns
        -------
        Nothing
        """
        if cells is None:
            cells = grAList()
        trace_ht = 80  # mV

        simulation_experiment = "Full"
        print("Cells: ", cells)
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
            PD = self.newPData()

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
                if n == (len(fn) - 1) and ic == 0:  # (len(cells)-1):
                    calxv = calxp
                    calyv = -50.0
                    iax = n
                else:
                    iax = None
                    calxv = None
                    calyv = -50.0
                # cprint("c", f"iax: {str(iax):s}, calxv = {str(calxv):s}  calyv = {str(calyv):s}")
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
                    caly=calyv,
                    calt=50.0,
                    calv=20.0,
                )
                axes[ic].annotate(
                    text=f"{n+1:d} ",
                    xy=(400.0, y0 - 60.0),  # ms
                    fontsize=8,
                    horizontalalignment="right",
                    verticalalignment="center",
                    # transform=axes[ic].transAxes,
                )

                if n == 0:
                    axes[ic].set_title(
                        f"BC{cellN:02d}",
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
            self.config["baseDataDirectory"],
            self.config["revcorrDataDirectory"],
            f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl",
        )
        with open(rc_datafile, "rb") as fh:
            d = FPM.pickle_load(fh)
        # print(d.keys())
        return d

    def Figure4_Supplmental1(self):
        fig = self.Figure4_Main(supplemental1=True)
        return fig

    def Figure4_assign_panel(self, supplemental1: bool = False, index: int = 0):
        if not supplemental1:
            revcorr_panel = f"D{index:d}"
            vm_panel = f"E{index:d}"
        else:
            revcorr_panel = f"B{index:d}"
            vm_panel = f"C{index:d}"
        return revcorr_panel, vm_panel

    def Figure4_Main(self, supplemental1=False):
        """
        Generate Figure 4 for the paper. Combined bits from various other plots
        and revcorrs/singles
        """
        print(f"Figure 4 main: supplemental1={str(supplemental1):s}")
        if supplemental1:
            example_cells = [2, 6, 10, 11, 13, 18]
        else:
            example_cells = [5, 9, 17, 30]

        start_letter = "A"
        parent_figure = None

        xw = 2.2
        xl = 0.75
        if not supplemental1:
            sizer = {
                # "B": {"pos": [6.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                # "C": {"pos": [9.5, 2.2, 4.25, 2.5], "labelpos": (-0.15, 1.02),},
                # "F": {"pos": [6.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
                # "G": {"pos": [9.5, 2.2, 0.5, 2.5], "labelpos": (-0.15, 1.02),},
                "B": {"pos": [xl, xw, 0.5, 2.5], "labelpos": (-0.15, 1.02)},
                "C": {
                    "pos": [xl + 2.75, xw, 0.5, 2.5],
                    "labelpos": (-0.15, 1.02),
                },
                "F": {
                    "pos": [6.5, 2.2, 0.5, 2.5],
                    "labelpos": (-0.15, 1.02),
                },
                "G": {
                    "pos": [9.5, 2.2, 0.5, 2.5],
                    "labelpos": (-0.15, 1.02),
                },
            }
            figsize = (12, 8)
        else:
            sizer = {}
            figsize = (9, 8)

        xw = 1.1
        xw2 = 1.0
        trace_axes = []
        if supplemental1:
            yh2 = 1.1
            yb2 = 2.25
            yb3 = 0.5
            yb1 = 3.75
            yh1 = 3.75
        else:
            yh2 = 1.2
            yb2 = 3.5 + 2.5
            yb3 = 3.5 + 0.5
            yb1 = 3.25
            yh1 = 4.25
        for j in range(len(example_cells)):
            i = j + 1
            pan_rev, pan_vm = self.Figure4_assign_panel(supplemental1, i)
            if not supplemental1:
                xl1 = j * 1.25 + 0.75
                xl2 = j * 1.25 + 6.5  # set panels to the right
            else:
                xl1 = j * 1.25 + 0.75
                xl2 = j * 1.25 + 0.75  # set panels on the bottom matching columns
            axn = f"A{i:d}"
            trace_axes.append(axn)
            sizer[axn] = {
                "pos": [xl1, xw, yb1, yh1],
                "labelpos": (-0.15, 1.03),
                "noaxes": True,
            }
            sizer[pan_rev] = {  # reverse correlatoin
                "pos": [xl2, xw2, yb2, yh2],
                "labelpos": (-0.15, 1.05),
                # "noaxes": True,
            }
            sizer[pan_vm] = {
                "pos": [xl2, xw2, yb3, yh2],
                "labelpos": (-0.15, 1.05),
                "noaxes": True,
            }
        # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            units="in",
            figsize=figsize,
            label=True,
            fontsize={"tick": 8, "label": 10, "panel": 13},
            showgrid=False,
            parent_figure=parent_figure,
        )
        # Efficacy plot
        if not supplemental1:
            EFP = EF.EfficacyPlots(parent_figure=P)
            EFP.plot_efficacy(
                "Full", datasetname_added="Added", ax=P.axdict["C"], clean=True
            )

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
                label=f"BC{n:02d}",
                clip_on=False,
                s=12,
            )
            ax.set_xlabel(r"Input ASA (${\mu m^2}$)")
            ax.set_xlim(0, 350)
            ax.set_ylim(0, 3)
            ax.set_ylabel(f"Participation ratio {dB:2d}/{0:2d} dBSPL")
            PH.nice_plot(ax, position=self.axis_offset, direction="outward")
            PH.talbotTicks(
                ax,
                density=(1.0, 1.5),
                insideMargin=0,
                tickPlacesAdd={"x": 0, "y": 1},
                floatAdd={"x": 0, "y": 1},
                axrange={"x": (0, 300.0), "y": (0, 3)},
                pointSize=None,
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

        # Traces
        axl = [P.axdict[axi] for axi in trace_axes]
        self.plot_stacked_traces(
            cells=example_cells, figure=P.figure_handle, axes=axl, maxstack=10
        )

        dB = 30
        if not supplemental1:
            # clusters
            plot_clustering(P.axdict["B"])

            # participation
            ds = self._load_rcdata("Spont")
            drc = self._load_rcdata(f"{dB:2d}dB")
            sns.set_palette(palette="tab10", n_colors=10)
            palette = sns.color_palette(palette="tab10", n_colors=len(ds.keys()))
            for i, c in enumerate(ds.keys()):
                # plot_participation(P.axdictax[0], c, ds, drc, dB=dB, color=palette[i])
                plot_diff_participation(
                    P.axdict["F"], c, ds, drc, dB=dB, color=palette[i], legend=False
                )

            # Cumulative plots
            self.plot_revcorr_compare(
                parent_figure=P,
                axlist=[P.axdict["G"]],  # P.axdict["G"]],
                dBSPLs=["Spont", "30dB"],
                legend=False,
            )
            PH.nice_plot(P.axdict["G"], position=self.axis_offset, direction="outward")
            synlabel_num = 5
        else:
            synlabel_num = 2

        # revcorrs and traces
        self.plot_revcorr_supplement(
            cells=example_cells,
            parent_figure=P,
            supplemental1=supplemental1,
            dBSPL="30dB",
            synlabel_num=synlabel_num,
        )
        # self.plot_efficacy_supplement(cells=example_cells, parent_figure=P, traces=False)

        # Revcorr axes cleanup
        for j in range(len(example_cells)):
            pan_rev, pan_vm = self.Figure4_assign_panel(supplemental1, j + 1)
            ax = P.axdict[pan_rev]
            ax.set_ylim(0, 0.8)
            ax.set_xlim(-5.0, 2.5)
            ax2 = P.axdict[pan_vm]
            ax2.set_xlim(-5.0, 2.5)
            ax2.set_ylim(-70, 0)

            if j == 0:
                ax.set_ylabel("Presynaptic\nCoinc. Rate (Hz)", ha="center", fontsize=10)
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
            fig.filename = "Figure4_Ephys2_main_v8.pdf"
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 (main) Modeling: singles inputs, efficacy and revcorr, revised version 8"
        else:
            fig.filename = "Figure4/Figure4_supp/Figure4-Supplemental1_Revcorr_V3.pdf"
            fig.title[
                "title"
            ] = "SBEM Project Figure 4 Modeling: Supplemental 1: other cells single inputs and revcorr"

        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        # fig.title2 = title2
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
        recompute=False,  # if you need to recompute the revcorrs for all the grade A cells, just set this True
    ) -> tuple:
        """
        Get the revcorr data associated with the cell number
        and the stimulus level (dbSPL, as a string)
        """
        cprint("c", "Calling _get_revcorr")
        cell_revcorr = FD.figure_revcorr[cell_number]
        run_calcs = False
        rc_datafile = Path(
            self.newPData().revcorrpath, f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl"
        )
        cprint("c", f"Rcdatafile: {str(rc_datafile):s}")
        if rc_datafile.is_file() and not recompute:
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = FPM.pickle_load(fh)
            RCD = all_RCD_RCP[cell_number][0]
            RCP = all_RCD_RCP[cell_number][1]
            PD = self.newPData()
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
            P, PD, RCP, RCD = self.parent.PLT.compute_revcorr(
                P=None,  # no plotting (and return P is None)
                gbc=str(cell_number),
                fn=fn[0],
                PD=self.newPData(),
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
        supplemental1=False,
        synlabel_num: int = 0,
        colormap: str = "magma",  # default color map
        save_calcs: bool = False,  # set to True if need to update.
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
        rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
        # if rc_datafile.is_file() and not run_calcs:
        #     with open(rc_datafile, "rb") as fh:
        #         all_RCD_RCP = FPM.pickle_load(fh)
        # else:
        #     print(
        #         "Must run revcorr_supplement plot first, then the file we need will be present"
        #     )
        #     run_calcs = False
        all_RCD_RCP = {}  # all revcorr data

        i_plot = 0
        if cells is None:
            cells = grAList()

        for i, cell_number in enumerate(cells):

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
            fn = sorted(list(sfi.glob("*")))
            print("Revcorr supplement file: ", fn)

            if i == 0:
                tcal = True  # point font for cal bar
            else:
                tcal = False  # no cal bar
            if synlabel_num == 0:  # labels all of them
                synlabel = True
            else:
                if cell_number == synlabel_num:
                    synlabel = True
                else:
                    synlabel = False

            all_RCD_RCP[cell_number] = [RCD, RCP]
            ax_top_row_name, ax_bot_row_name = self.Figure4_assign_panel(
                supplemental1, index=i_plot + 1
            )
            ax_top_row = P.axdict[ax_top_row_name]
            ax_bot_row = P.axdict[ax_bot_row_name]
            if cell_number in cells:
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
                    axarray=axlist,
                    calbar_show=tcal,
                    calbar_fontsize=7,
                    yaxis_label=False,
                    synlabel=synlabel,
                    colormap=colormap,
                )

                print(
                    f"  Mean pre: {np.mean(RCD.mean_pre_intervals):5.3f} ms  ({1./np.mean(RCD.mean_pre_intervals):7.3f} Hz)"
                )
                print(
                    f"  Mean Post: {RCD.mean_post_intervals:5.3f} ms  ({1./RCD.mean_post_intervals:7.3f} Hz)"
                )
                # if parent_figure is None:
                ax_top_row.text(
                    0.5,
                    1.0,
                    f"BC{cell_number:02d}",
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

        title = (
            "SBEM Project Supplemental Figure 3 Modeling : Reverse Correlation Summary",
        )
        save_file = f"Fig_M3_supplemental_Full_{dBSPL:s}.pdf"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        # fig.title["title"] = title
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
            print(len(r_keys))
            print("palette in G: ", palette)
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
                    label=f"BC{cell_number:02d}",
                )
                if i == 0:
                    hpoints.extend(
                        [[xnspike[k], RCD.ynspike[k]] for k in range(len(xnspike))]
                    )
            ax.set_xlim(0, 12)
            ax.set_ylim(0, 1.0)
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

            # ax.set_title(dBSPL)

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

    def Figure7_Main(self, parent_figure=None):
        sizer = OrderedDict(  # define figure layout
            [
                ("A", {"pos": [0.08, 0.4, 0.87, 0.09]}),
                ("B", {"pos": [0.08, 0.4, 0.73, 0.09]}),
                ("C", {"pos": [0.08, 0.4, 0.59, 0.09]}),
                ("D", {"pos": [0.08, 0.4, 0.45, 0.09]}),
                ("E", {"pos": [0.08, 0.4, 0.31, 0.09]}),
                ("F", {"pos": [0.08, 0.4, 0.17, 0.09]}),
                ("G", {"pos": [0.08, 0.4, 0.05, 0.07]}),
                ("H", {"pos": [0.57, 0.4, 0.87, 0.09]}),
                ("I", {"pos": [0.57, 0.4, 0.73, 0.09]}),
                ("J", {"pos": [0.57, 0.4, 0.59, 0.09]}),
                ("K", {"pos": [0.57, 0.4, 0.45, 0.09]}),
                ("L", {"pos": [0.57, 0.4, 0.31, 0.09]}),
                ("M", {"pos": [0.57, 0.4, 0.17, 0.09]}),
                ("N", {"pos": [0.57, 0.4, 0.05, 0.07]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.

        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            label=True,
            figsize=(6.0, 8.0),
            # labelposition=(-0.05, 1.02),
        )
        # Column 1: AM
        P.axdict["A"].set_ylabel("mV", fontsize=8)

        P.axdict["B"].set_title("Bushy Spike Raster", fontsize=9)
        P.axdict["B"].set_ylabel("Trial")

        P.axdict["C"].set_title(
            "Bushy PSTH", fontsize=9, verticalalignment="top", y=0.95
        )
        P.axdict["C"].set_ylabel("Spikes/second", fontsize=9)

        P.axdict["D"].set_title("ANF Spike Raster", fontsize=9)
        P.axdict["D"].set_ylabel("Trial")

        P.axdict["E"].set_title("ANF PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["E"].set_ylabel("Spikes/second", fontsize=9)

        P.axdict["F"].set_title("Phase", fontsize=9)
        P.axdict["F"].set_ylabel("Spike Count", fontsize=9)
        P.axdict["F"].set_title("Angle (radians)", fontsize=9)

        P.axdict["G"].set_title("Stimulus", fontsize=9)
        P.axdict["G"].set_ylabel("Amplitude (Pa)", fontsize=8)
        P.axdict["G"].set_xlabel("T (s)", fontsize=9)

        # column 2 Click train
        P.axdict["H"].set_ylabel("mV", fontsize=8)

        P.axdict["I"].set_title("Bushy Spike Raster", fontsize=9)
        P.axdict["I"].set_ylabel("Trial")

        P.axdict["J"].set_title(
            "Bushy PSTH", fontsize=9, verticalalignment="top", y=0.95
        )
        P.axdict["J"].set_ylabel("Spikes/second", fontsize=9)

        P.axdict["K"].set_title("ANF Spike Raster", fontsize=9)
        P.axdict["K"].set_ylabel("Trial")

        P.axdict["L"].set_title("ANF PSTH", fontsize=9, verticalalignment="top", y=0.95)
        P.axdict["L"].set_ylabel("Spikes/second", fontsize=9)

        P.axdict["M"].set_title("SAC", fontsize=9)
        P.axdict["M"].set_ylabel("CI", fontsize=9)
        P.axdict["M"].set_title("Time (ms)", fontsize=9)

        P.axdict["N"].set_title("Stimulus", fontsize=9)
        P.axdict["N"].set_ylabel("Amplitude (Pa)", fontsize=8)
        P.axdict["N"].set_xlabel("T (s)", fontsize=9)

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

        cell_number = 11
        self.Figure7_one_column(
            "SAM",
            cell_number,
            FD.figure_SAM_SAC,
            P,
            ["A", "B", "C", "D", "E", "F", "G"],
        )
        self.Figure7_one_column(
            "SAC",
            cell_number,
            FD.figure_SAM_SAC,
            P,
            ["H", "I", "J", "K", "L", "M", "N"],
        )

        fig = FigInfo()
        if parent_figure is not None:
            fig.P = parent_figure
        else:
            fig.P = P
        fig.filename = "Figure7/Figure7_Ephys3_main_v1_left.pdf"
        fig.title["title"] = "SBEM Project Figure 7 Modeling: SAM, SAC"
        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        fig.title2 = title2
        print("returnin fig: ", fig)
        return fig

    def Figure7_one_column(self, mode, cell_number, dataset, P, pan):

        PD = self.newPData()
        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )

        sfi = Path(cellpath, Path(dataset[mode][0]).name)
        if not sfi.is_dir():
            print("file not found: ", str(sfi))
            return

        fn = sorted(list(sfi.glob("*")))[0]
        changetimestamp = get_changetimestamp()
        X = self.ReadModel.get_data_file(fn, changetimestamp, PD)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if X is None:
            print("No simulation found that matched conditions")
            print(fn)
            return
        # unpack x
        par, stitle, ivdatafile, filemode, d = X
        axon_name = par["axonExpt"]
        protocol = "runANPSTH"
        AR, SP, RMA = analyze_data.analyze_data(ivdatafile, filemode, protocol)
        i = 0
        plot_win = (0.4, 0.550)
        psth_binw = 0.0005
        i = 0  # Voltage for first trial
        self.plot_voltage(ax=P.axdict[pan[0]], ntrace=i, d=d, AR=AR, time_win=plot_win)
        self.plot_stim_waveform(
            ax=P.axdict[pan[6]], ntrace=i, d=d, AR=AR, stim_win=plot_win
        )
        all_bu_st = self.get_bu_spikearray(AR, d)
        self.plot_spiketrain_raster(all_bu_st, ax=P.axdict[pan[1]], plot_win=plot_win)
        self.plot_psth_psth(
            ax=P.axdict[pan[2]],
            data=all_bu_st,
            ri=d["runInfo"],
            psth_binw=psth_binw,
            ntr=len(all_bu_st),
            psth_win=plot_win,
        )

        an_st_by_input, all_an_st, an_st_grand = self.get_an_spikearray(AR, d)
        ninputs = len(an_st_by_input)

        self.plot_stacked_spiketrain_rasters(
            an_st_by_input, ax=P.axdict[pan[3]], si=d["Params"], plot_win=plot_win
        )
        P.axdict[pan[3]].set_xlabel("Time (s)")

        self.plot_psth_psth(
            ax=P.axdict[pan[4]],
            data=all_an_st,
            ri=d["runInfo"],
            psth_binw=psth_binw,
            psth_win=plot_win,
            ntr=len(all_an_st),
            ninputs=ninputs,
        )
        P.axdict[pan[4]].set_xlabel("Time (sec)")
        P.axdict[pan[4]].set_title("AN")
        ri = d["runInfo"]
        si = d["Params"]
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
            P.axdict[pan[5]].plot(
                bu_sacbins,
                bu_sac,
                "k-",
                # label=sac_label,
            )
            an_sac, an_sacbins = S.SAC_with_histo(
                an_st_grand,
                pars=spars,
                engine=sac_engine,
                dither=1e-3 * si.dtIC / 2.0,
            )
            P.axdict[pan[5]].plot(
                an_sacbins,
                an_sac,
                "r-",
                # label=sac_label,
            )
        else:
            phasewin = [
                pip_start + 0.25 * pip_duration,
                pip_start + pip_duration,
            ]
            y = []
            for z in all_bu_st:
                y.extend(z)
            x = np.array(y)
            v = np.where((phasewin[0] <= x) & (x < phasewin[1]))[0]
            bu_spikesinwin = x[v]
            vs_bu = self.parent.PLT.VS.vector_strength(bu_spikesinwin, fmod)
            vs_an = self.parent.PLT.VS.vector_strength(all_an_st, fmod)
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
            self.parent.PLT.plot_psth(
                vs_bu.circ_phase,
                run_info=ri,
                max_time=2 * np.pi,
                bin_width=est_binw,
                ax=P.axdict[pan[5]],
            )
            self.parent.PLT.plot_psth(
                vs_an.circ_phase,
                run_info=ri,
                max_time=2 * np.pi,
                bin_width=est_binw,
                ax=P.axdict[pan[5]],
                bin_fill=False,
                edge_color="r",
                alpha=0.5,
            )
            self.parent.PLT.plot_psth(
                vs_bu.circ_phase,
                run_info=ri,
                max_time=2 * np.pi,
                bin_width=est_binw,
                ax=P.axdict[pan[5]],
                alpha=0.5,
            )
            # P.axdict["E"].hist(
            #     vs["ph"],
            #     bins=2 * np.pi * np.arange(30) / 30.0,
            #     facecolor="k",
            #     edgecolor="k",
            # )
            P.axdict[pan[5]].set_xlim((0.0, 2 * np.pi))
            P.axdict[pan[5]].set_title(
                f"VS: AN = {vs_an.vs:.3f} BU = {vs_bu.vs:.3f}",
                fontsize=8,
                horizontalalignment="center",
            )

    def Figure7_Supplemental1(self):
        V = SAM_VS_vplots.VS_Plots()
        V.make_figure()

    def plot_psth_psth(
        self,
        ax: object,
        data: object,
        ri: dict,
        psth_binw: float = 0.0005,
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
            0.20,
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
        for cell in grAList():
            fig = self.plot_PSTH(cellN=cell)
        return fig

    def plot_PSTH(self, cellN=None):
        dBSPL = "30dB"
        if cellN is None:
            cell_number = 17
        else:
            cell_number = cellN
        print(f"Plotting PSTH for BC{str(cell_number):s}")
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
        title = ("SBEM Project Figure 3 Modeling : PSTH Summary",)
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
        ax.plot((time_base - time_win[0]), vtrial, "k-", linewidth=0.5)
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
    ):
        # stimulus waveform
        trd = d["Results"][ntrace]
        waveform = trd["stimWaveform"].tolist()
        stb = trd["stimTimebase"]
        stimdt = np.mean(np.diff(stb))
        sttb_beg = int(stim_win[0] / stimdt)
        sttb_end = int(stim_win[1] / stimdt)

        ax.plot(
            (stb[sttb_beg:sttb_end] - stim_win[0]),
            np.array(waveform)[sttb_beg:sttb_end] * 1e3,
            "k-",
            linewidth=0.5,
        )  # stimulus underneath
        ax.set_xlim(0, stim_win[1] - stim_win[0])
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 0.5),
            insideMargin=0.02,
            # pointSize=ticklabelsize,
            tickPlacesAdd={"x": 2, "y": 1},
            floatAdd={"x": 2, "y": 1},
        )
        ax.set_ylabel("mPa")

    def plot_spiketrain_raster(self, spike_times, ax=None, plot_win: tuple = (0, 0.25)):
        # print(len(spike_times))
        # print(plot_win)
        for i in range(len(spike_times)):
            ispt = [
                j
                for j in range(len(spike_times[i]))
                if spike_times[i][j] >= plot_win[0] and spike_times[i][j] < plot_win[1]
            ]
            ax.plot(
                np.array(spike_times[i])[ispt] - plot_win[0],
                i * np.ones(len(ispt)),
                "|",
                markersize=1.5,
                color="b",
            )

    def plot_stacked_spiketrain_rasters(
        self,
        spike_times_by_input,
        ax=None,
        si=None,
        plot_win: tuple = (0, 0.25),
        max_trials=5,
        use_colors: bool = True,
        colormap="Set3",
        cbar_vmax: float = 300.0,
    ):
        """
        Spike trains are plotted as a raster for all inputs in the AN data

        Parameters
        ----------
        spike_times_by_input : list
            list by input of spike times by trial

        ax : matplotlib axis to place the plat

        si : Params

        plot_win : tuple
            time window of data to show

        max_trials : int
            number of trials to show (the first max_trials are plotted)

        use_colors : bool (default True)
            Plot the raster ticks with colors scaled by input surface area

        colormap : str
            color map to use for plotting when use_colors is True

        cbar_vmax : float
            maximum for the colorbar scale

        """
        n_inputs = len(spike_times_by_input)
        n_trials = len(spike_times_by_input[0])
        trial_spc = 1.0 / (n_inputs * 2)
        input_spc = 1

        if use_colors:
            cmx = sns.color_palette(colormap, as_cmap=True)
            cell_n = int(si.cellID[-2:])
            SC, syninfo = self.parent.PLT.get_synaptic_info(cell_n)
            syn_ASA = np.array([syninfo[1][isite][0] for isite in range(n_inputs)])
            max_ASA = np.max(syn_ASA)
        print("syn asa: ", syn_ASA)
        for k in range(n_inputs):  # raster of input spikes by input
            for i in range(n_trials):  # and by trial
                if i > max_trials - 1:
                    continue
                if np.max(spike_times_by_input[k][i]) > 2.0:  # probably in msec...
                    tk = np.array(spike_times_by_input[k][i]) * 1e-3
                else:
                    tk = np.array(spike_times_by_input[k][i])

                in_spt = [
                    i
                    for i in range(tk.shape[0])
                    if (tk[i] >= plot_win[0]) and (tk[i] < plot_win[1])
                ]
                y = (i * input_spc + ((n_trials - k - 1) * trial_spc)) * np.ones(
                    len(in_spt)
                )
                if use_colors:
                    color = cmx.colors[int(cmx.N * syn_ASA[k] / cbar_vmax) - 1]
                else:
                    color = "k"
                ax.plot(
                    tk[in_spt] - plot_win[0],
                    y,
                    "|",
                    markersize=2.5,
                    color=color,
                    linewidth=0.5,
                )
        ax.set_xlim(plot_win)

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
        label_x_axis=True,
    ):
        dataset = FD.figure_psth[cell_number]
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
            return

        fn = sorted(list(sfi.glob("*")))[0]
        changetimestamp = get_changetimestamp()
        X = self.ReadModel.get_data_file(fn, changetimestamp, PD)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if X is None:
            print("No simulation found that matched conditions")
            print(fn)
            return
        # unpack x
        par, stitle, ivdatafile, filemode, d = X
        axon_name = par["axonExpt"]
        protocol = "runANPSTH"
        AR, SP, RMA = analyze_data.analyze_data(
            ivdatafile=ivdatafile, filemode=filemode, protocol=protocol
        )
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

        all_bu_st = self.get_bu_spikearray(AR, d)
        self.all_bu_st = all_bu_st
        psth_binw = 0.5e-3
        ninputs = 1

        if bupsth_ax is not None:
            self.plot_psth_psth(
                ax=bupsth_ax,
                data=all_bu_st,
                ri=ri,
                psth_binw=psth_binw,
                psth_win=psth_win,
                ntr=ntr,
                ninputs=1,
            )
            PH.nice_plot(bupsth_ax, direction="outward", ticklength=3.0)

            if label_x_axis:
                bupsth_ax.set_xlabel("Time (s)")

        an_st_by_input, all_an_st, an_st_grand = self.get_an_spikearray(AR, d)

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

        if bufsl_ax is not None:
            PF.plot_fsl_ssl(
                all_bu_st,
                run_info=ri,
                max_time=25.0,
                bin_width=0.25,
                fsl_win=bu_fsl_win,
                ax=bufsl_ax,
                zero_time=psth_win[0] * 1e-3,
                cellID=cell_number,
            )
            PH.nice_plot(bufsl_ax, direction="outward", ticklength=3.0)

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

    def plot_one_CV(
        self,
        cell_number: int,
        dbSPL: float,
        cv_win: Union[list, tuple] = (0.220, 0.3),
        reftime: float = 0.0,
        t_grace: float = 0.0,
        cv_binw: float = 0.001,
        cv_ax: object = None,
        label_x_axis: bool = False,
    ):
        # use the data from PSTH. If not present, skip the CV plot
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
        cv_ax.plot((cvt - reftime) * 1e3, cvs / cvm, "k-")

        PH.nice_plot(cv_ax, direction="outward", ticklength=3.0)
        cv_ax.set_ylim(0, 1.2)
        cv_ax.set_xlim(0, 1e3 * (cv_win[1] - cv_win[0])- 20.0)
        PH.talbotTicks(
            cv_ax,
            axes="xy",
            density=(3, 1.5),
            insideMargin=0.02,
            # pointSize=ticklabelsize,
            tickPlacesAdd={"x": 2, "y": 1},
            floatAdd={"x": 2, "y": 1},
        )
        # compute the mean CV from 10-60 msec for display on the plot
        itwin = np.where(((cvt - reftime) > 0.01) &( (cvt - reftime) < 0.06))[0]
        CVp = np.nanmean(cvs[itwin] / cvm[itwin])
        cvlabel = r"$CV\prime$"
        cv_ax.text(
            1.0,
            1.15,
            f"{cvlabel:s}: {CVp:4.2f}",
            transform=cv_ax.transAxes,
            fontsize=7,
            horizontalalignment='right',
            verticalalignment='top',
        )
        cv_ax.set_ylabel("CV")
        if label_x_axis:
            cv_ax.set_xlabel("Latency (ms)")

    def Figure3_Supplemental4_PSTH(self):
        print("Plotting Figure 3 Supplement 4 PSTH")
        dBSPL = "30dB"
        lmar = 0.125
        rmar = 0.1
        hspc = 0.08
        P = PH.regular_grid(
            rows=len(grAList()),
            cols=4,
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

        for i, cell_number in enumerate(grAList()):
            print("Plotting psth for cell: ", cell_number)
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
            self.plot_one_CV(
                cell_number=cell_number,
                dbSPL=dBSPL,
                cv_win=(0.0, 0.100),
                cv_ax=P.axarr[i, 3],
                label_x_axis=show_label,
                t_grace=0.0,
            )

            P.axarr[i, 0].text(
                -0.40,
                0.5,
                f"BC{cell_number:02d}",
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

        save_file = f"Figure3/Figure3_supp/Figure3-Supplemental4_PSTH.pdf"
        title = "SBEM Project Figure 3 Modeling Supplemental : PSTH and FSL, all grade A cells"
        fig = FigInfo()
        fig.P = P
        fig.filename = save_file
        fig.title["title"] = title
        return fig

    def plot_VS_SAM(self):
        self.generate_VS_data_file()
        pass

    def analyze_VS_data(
        self, VS_data, cell_number, fout, firstline=False, sac_flag=False
    ):
        """
        Generate tables of vs measures for all cells
        across the frequencies listed
        """
        self.parent.PLT.textclear()  # just at start
        PD = self.newPData()
        P = None
        self.parent.cellID = cell_number
        for i, filename in enumerate(VS_data.samdata[cell_number]):
            print(f"Cell: {cell_number:d}  Filename: {filename:s}")
            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, filename + ".pkl")
            print(f"Opening pkl file: {str(sfi):s}")
            with open(sfi, "rb") as fh:
                d = FPM.pickle_load(fh)

            self.parent.PLT.plot_AN_response(
                P, d.files[0], PD, "runANPSTH", sac_flag=sac_flag
            )
            with open(fout, "a") as fth:
                if firstline:
                    fth.write(self.parent.PLT.VS_colnames, sac_flag=sac_flag)
                    fth.write("\n")
                    firstline = False
                fth.write(self.parent.PLT.VS_line)
                fth.write("\n")

        print("Analyze VS data completed for cell: ", cell_number)

    #
    def generate_VS_data_file(self, testmode=False, dB=30, append=True):
        """
        Write the VS_data file from the selected datasets.
        The datasets are in VS_datasets_nndB.py
        """

        if dB == 30:
            if f"VS_datasets_{dB:d}dB" not in list(dir()):
                import VS_datasets_30dB as VS_datasets
        if dB == 15:
            if f"VS_datasets_{dB:d}dB" not in list(dir()):
                from vcnmodel import VS_datasets_15dB as VS_datasets
        importlib.reload(VS_datasets)
        print("VS_datasets: ", VS_datasets)
        print(f"Data set keys found: {str(list(VS_datasets.samdata.keys())):s}")

        config = toml.load(open("wheres_my_data.toml", "r"))

        """
        Generate the table in VS_data.py by analyzing the data from 
        VS_datasets.py
        """
        cprint("g", f"Generate VS Data for {dB:d} dB")

        fout = f"VS_data_{dB:d}dB.py"  # we will generate this automatically
        if not append:
            with open(fout, "w") as fh:
                fh.write(f'"""\n')
                fh.write(
                    "    Vector strength for models with SAM tones, different input configurations.\n"
                )
                fh.write("    17 Aug 2021 version.\n")
                fh.write(
                    "    Results are a printout from DataTablesVCN after selecting the data runs.\n"
                )
                fh.write(
                    "NOTE: This table is automatically written by figures.py and should not be\n"
                )
                fh.write("      directly edited.")
                fh.write(
                    f"To Regenerate:\n  After running the simulations, use 'Print File Info' for each cell, "
                )
                fh.write(
                    f"  selecting the relevant simulations, and copy the text in the 'Reports' box in DataTablesVCN"
                )
                fh.write(
                    f"  into a 'VS_datasets_xxdB.py' file, where xx is the sound pressure level.\n"
                )
                fh.write(
                    f"Then select 'VS-SAMTone-no figure' in DataTables, and 'Create Figure."
                )
                fh.write(
                    f"  No figure will be generated, but the VS_data_xxdB.py file will be created.\n"
                )
                fh.write(
                    f"The VS_data file holds all of the vector-strength information, in a text format, "
                )
                fh.write(f"  and is read by the plotting programs.\n")
                fh.write(f'--pbm 2014-2021\n"""\n')  # end of the comment
                fh.write('\ndata = """')  # start of the data

        fl = True
        for i, celln in enumerate(grAList()):
            if testmode and i != 1:
                continue
            self.analyze_VS_data(VS_datasets, celln, fout, firstline=fl, sac_flag=True)
            fl = False
        with open(fout, "a") as fh:
            fh.write(f'"""\n')  # close the data text.
        cprint("g", "The VS_data file {str(fout):s} has been generated.")

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
