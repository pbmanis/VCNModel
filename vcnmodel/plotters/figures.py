import datetime
import importlib
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import matplotlib
import numpy as np
import pandas as pd
import rich as RI
import seaborn as sns
from cycler import cycler
from matplotlib import image as mpimg
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH
from pylibrary.tools import cprint as CP
from rich.console import Console
from rich.text import Text

import toml
from vcnmodel.plotters import efficacy_plot as EF
from vcnmodel.plotters import plot_z as PZ

config = toml.load(open("wheres_my_data.toml", "r"))
cprint = CP.cprint


def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """
    return [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]


syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x"]
colors = ["c", "k", "m", "r"]

title_text = {"passive": "Passive", "normal": "Half-active", "active": "Active"}
font_colors = {"passive": "c", "normal": "k", "active": "m"}
spike_colors = font_colors


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
    thiscell: str = ""


def title_data():
    return {'title': "", "x":None, "y":None}


@dataclass
class FigInfo:
    """
    data class to hold figure information for savefig
    """
    P: object=None
    filename: Union[str, Path]=""
    title: dict = field(default_factory = title_data)
    title2: dict = field(default_factory = title_data)


figure_IV = {
    "Cell": 17,
    "normal": "runIV-all-2020-07-30.18-29-53",
    "passive": "runIV-all-2020-07-30.18-43-29",
    "active": "runIV-all-2020-07-30.18-56-13",
    "Z_normal": "VCN_c17_Full_normal_Z.pkl",
    "Z_passive": "VCN_c17_Full_pasdend_Z.pkl",
    "Z_active": "VCN_c17_Full_actdend_Z.pkl",
}
figure_AllIVs = {
    2: {
        "normal": "runIV-all-2020-07-30.18-20-16",
        "passive": "runIV-all-2020-07-30.18-33-07",
        "active": "runIV-all-2020-07-30.18-46-46",
    },
    5: {
        "normal": "runIV-all-2020-07-30.18-21-32",
        "passive": "runIV-all-2020-07-30.18-34-20",
        "active": "runIV-all-2020-07-30.18-47-59",
    },
    6: {
        "normal": "runIV-all-2020-07-30.18-23-12",
        "passive": "runIV-all-2020-07-30.18-36-04",
        "active": "runIV-all-2020-07-30.18-49-38",
    },
    9: {
        "normal": "runIV-all-2020-07-30.18-24-19",
        "passive": "runIV-all-2020-07-30.18-37-10",
        "active": "runIV-all-2020-07-30.18-50-44",
    },
    10: {
        "normal": "runIV-all-2020-07-30.18-25-42",
        "passive": "runIV-all-2020-07-30.18-38-31",
        "active": "runIV-all-2020-07-30.18-52-07",
    },
    11: {
        "normal": "runIV-all-2020-07-30.18-27-24",
        "passive": "runIV-all-2020-07-30.18-40-45",
        "active": "runIV-all-2020-07-30.18-53-50",
    },
    13: {
        "normal": "runIV-all-2020-07-30.18-28-30",
        "passive": "runIV-all-2020-07-30.18-42-00",
        "active": "runIV-all-2020-07-30.18-54-51",
    },
    17: {
        "normal": "runIV-all-2020-07-30.18-29-53",
        "passive": "runIV-all-2020-07-30.18-43-29",
        "active": "runIV-all-2020-07-30.18-56-13",
    },
    18: {
        "normal": "runIV-all-2020-11-16.10-54-53",
        "active": "runIV-all-2020-11-16.10-55-23",
        "passive": "runIV-all-2020-11-16.10-55-08",
    },
    30: {
        "normal": "runIV-all-2020-07-30.18-31-35",
        "passive": "runIV-all-2020-07-30.18-45-12",
        "active": "runIV-all-2020-07-30.18-57-54",
    },
}

"""
The efficacy data is taken from runs using the latest measure
of synapse density, 0.7686 syn/um2  11/15/2020
"""
figure_efficacy_supplement = {
    2: {
        "NoDend": "runANSingles-all-2020-11-16.17-04-23",
        "Full": "runANSingles-all-2020-11-16.17-08-55",
    },
    5: {
        "NoDend": "runANSingles-all-2020-11-16.17-19-30",
        "Full": "runANSingles-all-2020-11-16.17-25-11",
    },
    6: {
        "NoDend": "runANSingles-all-2020-11-16.17-40-50",
        "Full": "runANSingles-all-2020-11-16.17-46-05",
    },
    9: {
        "NoDend": "runANSingles-all-2020-11-16.17-56-43",
        "Full": "runANSingles-all-2020-11-16.18-04-06",
    },
    10: {
        "NoDend": "runANSingles-all-2020-11-16.18-20-31",
        "Full": "runANSingles-all-2020-11-16.18-28-43",
    },
    11: {
        "NoDend": "runANSingles-all-2020-11-16.18-51-40",
        "Full": "runANSingles-all-2020-11-16.18-57-43",
    },
    13: {
        "NoDend": "runANSingles-all-2020-11-16.19-09-30",
        "Full": "runANSingles-all-2020-11-16.19-14-35",
    },
    17: {
        "NoDend": "runANSingles-all-2020-11-16.19-27-02",
        "Full": "runANSingles-all-2020-11-16.19-33-47",
    },
    18: {
        "NoDend": "runANSingles-all-2020-11-16.19-50-22",
        "Full": "runANSingles-all-2020-11-16.19-57-52",
    },
    30: {
        "NoDend": "runANSingles-all-2020-11-16.20-09-36",
        "Full": "runANSingles-all-2020-11-16.20-16-56",
    },
}


figure_revcorr_example = {
    17: {
        "Spont": "runANPSTH-all-2020-11-24.09-12-49",
        "40dB": "runANPSTH-all-2020-08-20.14-54-39",
    }
}

figure_revcorr = {
    2: {
        "Spont": "runANPSTH-all-2020-11-24.06-59-18",
        "40dB": "runANPSTH-all-2020-11-24.11-13-17",
    },
    5: {
        "Spont": "runANPSTH-all-2020-11-24.07-15-41",
        "40dB": "runANPSTH-all-2020-11-24.11-19-05",
    },
    6: {
        "Spont": "runANPSTH-all-2020-11-24.07-38-28",
        "40dB": "runANPSTH-all-2020-11-24.11-27-20",
    },
    9: {
        "Spont": "runANPSTH-all-2020-11-24.07-55-04",
        "40dB": "runANPSTH-all-2020-11-24.11-33-12",
    },
    10: {
        "Spont": "runANPSTH-all-2020-11-24.08-15-28",
        "40dB": "runANPSTH-all-2020-11-24.11-40-19",
    },
    11: {
        "Spont": "runANPSTH-all-2020-11-24.08-39-48",
        "40dB": "runANPSTH-all-2020-11-24.11-47-36",
    },
    13: {
        "Spont": "runANPSTH-all-2020-11-24.08-55-43",
        "40dB": "runANPSTH-all-2020-11-24.11-52-44",
    },
    17: {
        "Spont": "runANPSTH-all-2020-11-24.09-12-49",
        "40dB": "runANPSTH-all-2020-11-24.11-58-11",
        "50dB": "runANPSTH-all-2020-08-20.14-54-39",
    },
    18: {
        "Spont": "runANPSTH-all-2020-11-24.09-37-01",
        "40dB": "runANPSTH-all-2020-11-24.12-05-53",
    },
    30: {
        "Spont": "runANPSTH-all-2020-11-24.09-51-06",
        "40dB": "runANPSTH-all-2020-11-24.12-10-28",
    },
}


figure_VClamp = {
    17: [
        "runVC-all-2020-07-29.10-36-59",
        "runVC-all-2020-07-29.10-30-30",
        "runVC-all-2020-07-29.12-17-13",
    ],
}

figure_psth = {
    2: {"40dB": "runANPSTH-all-2020-11-24.15-39-05",},
    5: {"40dB": "runANPSTH-all-2020-11-24.15-46-45",},
    6: {"40dB": "runANPSTH-all-2020-11-24.15-57-32",},
    9: {"40dB": "runANPSTH-all-2020-11-24.16-05-10",},
    10: {"40dB": "runANPSTH-all-2020-11-24.16-14-34",},
    11: {"40dB": "runANPSTH-all-2020-11-24.16-25-29",},
    13: {"40dB": "runANPSTH-all-2020-11-24.16-32-59",},
    17: {"40dB": "runANPSTH-all-2020-11-24.16-40-56",},
    18: {"40dB": "runANPSTH-all-2020-11-24.16-51-50",},
    30: {"40dB": "runANPSTH-all-2020-11-24.16-58-24",},
}

all_figures = [
    figure_psth,
    figure_VClamp,
    figure_revcorr,
    figure_revcorr_example,
    figure_efficacy_supplement,
    figure_IV,
    figure_AllIVs,
]


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
            "IV Figure": self.plotIV,
            "All IVs": self.allIVs,
            "IV Supplement": self.plot_IVS,
            "Zin Supplement": self.plot_ZinS,
            "Efficacy": self.plot_efficacy,
            "Efficacy Supplement": self.plot_efficacy_supplement,
            "Revcorr Ex": self.plot_revcorr,
            "All Revcorr": self.plot_all_revcorr,
            "Revcorr Supplement": self.plot_revcorr_supplement,
            "Revcorr Compare": self.plot_revcorr_compare,
            "PSTH-FSL": self.plot_PSTH,
            "All PSTHs": self.plot_All_PSTH,
            "PSTH-FSL Supplement": self.plot_PSTH_supplement,
            "VS-SAM Tone": self.plot_VS_SAM,
            "VC-KLTCalibration": self.plot_VC_gKLT,
            "CombinedVCIV": self.plot_combined_VC_IV,
        }
        print(figure_name, dispatch_table.keys())
        if figure_name in list(dispatch_table.keys()):
            print("disp")
            fig = dispatch_table[figure_name]()
            if fig is not None:
                self.save_figure(fig)

        
    def save_figure(self, fig):
        #P, save_file, title, title2={"title": None, 'x': 0.0, 'y': 0.0}):

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
        print('saving to: ', str(Path(config["baseDataDirectory"], "Figures", fig.filename)))
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", fig.filename),
            metadata={"Creator": "Paul Manis", "Author": "Paul Manis", "Title": fig.title},
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

    def plotIV(self, cell=None, parent_figure=None, loc:Union[None, tuple]=None):
        # style = STY.styler(journal="JNeurosci", figuresize='full', font='stixsans')  # two colukn...
        if cell is None:
            cellN = figure_IV["Cell"]
            d1 = figure_IV["passive"]
        else:
            cellN = cell
            d1 = figure_AllIVs[cellN]["passive"]
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
            tmar = fsize[1]-loc[3]
        else:
            bmar = 0.
            tmar = 0.

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
                "bottommargin": 2.25+bmar,
                "leftmargin": 0.5,
                "rightmargin": 0.5,
                "topmargin": 0.5+tmar,
            },
            labelposition=(-0.05, 1.05),
            parent_figure=parent_figure,
            panel_labels=p0panels,
        )

        self.P2 = PH.regular_grid( # lower row of 4
            1,
            4,
            order="rowsfirst",
            units="in",
            showgrid=False,
            figsize=fsize,
            verticalspacing=0.25,
            horizontalspacing=0.65,
            margins={
                "bottommargin": 0.5+bmar,
                "leftmargin": 0.5,
                "rightmargin": 0.5,
                "topmargin": 2.0+tmar,
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
                "bottommargin": 0.5+bmar,
                "leftmargin": 2.35,
                "rightmargin": 4.15,
                "topmargin": 2.0+tmar,
            },
            labelposition=(-0.25, 1.05),
            parent_figure=self.P,
            panel_labels=p3panels,
        )

        for iax, iv in enumerate(["passive", "normal", "active"]):
            if cell is None:
                dfile = figure_IV[iv]
            else:
                dfile = figure_AllIVs[cellN][iv]
            sfi = Path(cellpath, Path(dfile).name)
            if not sfi.is_dir():
                 print("Missing file: sfi 1: ", sfi)
                 return None
            fn = list(sfi.glob("*"))

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
            self.P.axarr[0, iax].set_title(title_text[iv], color="k", fontweight="normal")

        from matplotlib.lines import Line2D
        ls = ['-', '-', '-', '']
        mc = ['w', 'w', 'w', 'r']
        lines = [Line2D([0], [0], color=c, linewidth=0.75, linestyle=ls[ic], marker='o', mfc=mc[ic], mew=0.5) for ic, c in enumerate(colors)]
        labels = ["Passive", "Half-active", "Active", "Spk Thresh" ]
        self.P2.axarr[0,0].legend(lines, labels, bbox_to_anchor=(0.32, 0.35), fontsize=7)
            # self.P2.axarr[0, 0].line.set_label(iv)
        # self.P2.axarr[0,0].legend(["passive", "normal", "active", "spike thresholds"])
        
        res_label = r"$\mathregular{R_{in} (M\Omega)}$"
        tau_label = r"$\mathregular{\tau_{m} (ms)}$"
        phase_label = r"$\mathregular{\phi (radians)}$"
        
        # plot overlays of all cell z/phase
        for iax, mode in enumerate(["Z_passive", "Z_normal", "Z_active"]):
            if mode not in figure_IV.keys():
                continue
            # cprint('r', f"doing iv: {str(mode):s}")
            sfi = Path(config["codeDirectory"], figure_IV[mode])
            if not sfi.is_file():
                cprint('r', f"File not found!!!!!!\n->>> {str(sfi):s}")
                return None
                continue
            # cprint('c', f"sfi Z: {str(sfi):s}")
            # ax = self.P2.axarr[0, 1]
            ax = self.P3.axdict[phaseplot1]
            label = sfi.name  # .replace("_", "\_")
            with open(sfi, "rb") as fh:
                d = pickle.load(fh)
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
        for rax, iv in enumerate(figure_AllIVs.keys()):
            # cellpath = Path(self.config['cellDataDirectory'], f"VCN_c{iv:02d}", 'Simulations', 'IV')
            for iax, dendmode in enumerate(["passive", "normal", "active"]):
                sfi = Path(
                    self.config["cellDataDirectory"],
                    f"VCN_c{iv:02d}",
                    "Simulations",
                    "IV",
                    figure_AllIVs[iv][dendmode],
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
            data=df_rin, x="dendrites", y="Rin", ax=self.P2.axdict[rinplot], color="0.6",
            size=4.0, edgecolor='k',
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
            data=df_tau, x="dendrites", y="taum", ax=self.P2.axdict[tauplot], color="0.6",
            size=4.0, edgecolor='k',
        )
        self.P2.axdict[tauplot].set_ylim(0, 2.0)
        self.P2.axdict[tauplot].set_ylabel(tau_label)
        self.P2.axdict[tauplot].set_xlabel("Dendrite Decoration")

        fig = FigInfo()
        fig.P = self.P
        if cellN == None:
            fig.filename = "Fig_M1.pdf"
        else:
             fig.filename  = f"Fig_IV/IV_cell_VCN_c{cellN:02d}.png"
        fig.title["title"] = "SBEM Project Figure 1 Modeling (Main)"
        # self.save_figure(self.P, save_file, title)
        return fig

    def plot_combined_VC_IV(self):
        P0 = PH.regular_grid(  # dummy figure space
            1,
            1,
            figsize=(8., 9.),
            order="rowsfirst",
            units="in",
            showgrid=True,
            parent_figure=None,
        )
        
        figp1 = self.plot_VC_gKLT(parent_figure=P0, loc=(0, 8, 4., 9.))
        P2 = self.plotIV(parent_figure=figp1.P, loc=(0, 8, 0., 4.))
        mpl.show()
        


    def plot_IVS(self, parent_figure=None):
        """
        All of the IVS, for a supplemental figure
        Passive, normal, active, plus the crossed IV
        Also put the PNG for the cell on the left.
        """
        nivs = len(figure_AllIVs)
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
        for rax, iv in enumerate(figure_AllIVs.keys()):
            celln = Path(png_path, f"VCN_c{iv:02d}.png")
            if celln.is_file():  # add images from png files
                img = mpimg.imread(str(celln))
                self.P.axarr[rax, 0].imshow(img, aspect="equal")
                ylim = self.P.axarr[rax, 0].get_ylim()
                self.P.axarr[rax, 0].set_xlim(900, 1500)
                PH.noaxes(self.P.axarr[rax, 0])
            # plot 3 dendrite decorations
            for iax, dendmode in enumerate(["passive", "normal", "active"]):
                sfi = Path(
                    cellpath,
                    f"VCN_c{iv:02d}",
                    "Simulations",
                    "IV",
                    figure_AllIVs[iv][dendmode],
                )
                if not sfi.is_dir():
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
        save_file = "Fig_M1A_Supplemental.pdf"
        title = "SBEM Project Figure 1 Modeling (Supplemental A)"
        if parent_figure is None:
            self.save_figure(self.P)

    def plot_ZinS(self):
        """
        All of the Zins for a supplemental figure
        """
        PZ.PlotZ()

    def plot_efficacy(self):
        cell_number = 17
        example = figure_efficacy_supplement[cell_number]

        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "AN",
        )
        sfi = Path(cellpath, Path(example["Full"]).name)
        if not sfi.is_dir():
            return
        fn = sorted(list(sfi.glob("*")))
        print("fn: ", fn)
        PD = PData()
        calx = 800.0
        figsize = (8.0, 3.0)

        EFP = EF.EfficacyPlots(None, hold=True, cols=4)
        EFP.P.figure_handle.set_size_inches(figsize[0], figsize[1])

        self.P_Eff1 = PH.regular_grid(
            len(fn),
            1,
            order="rowsfirst",
            figsize=(8.0, 4.5),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.1,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.05,
                "width": 0.17,
                "height": 0.8,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=EFP.P,
            # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        self.P_Eff2 = PH.regular_grid(
            len(fn),
            1,
            order="rowsfirst",
            figsize=(8.0, 4.5),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.1,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.54,
                "width": 0.17,
                "height": 0.8,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=EFP.P,
            # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        calx = 800.0
        nfiles = len(fn)
        for n in range(nfiles):
            sfile = Path(sfi, fn[n])
            print(n, sfile)
            if n == nfiles - 1:
                calxv = calx
            else:
                calxv = None
            print("Calxv: ", calxv)
            self.parent.PLT.plot_traces(
                ax=self.P_Eff1.axarr[n, 0],
                fn=sfile,
                PD=PD,
                protocol="runANSingles",
                ymin=-90.0,
                ymax=20.0,
                xmin=400.0,
                xmax=900.0,
                iax=n,
                nax=len(fn),
                rep=0,
                figure=self.P_Eff1.figure_handle,
                longtitle=True,
                ivaxis=None,
                ivcolor="k",
                iv_spike_color="r",
                spike_marker_size=1.5,
                spike_marker_color="c",
                calx=calxv,
                caly=-20.0,
                calt=50.0,
                calv=20.0,
            )

            self.P_Eff1.axarr[n, 0].text(
                -0.02,
                0.5,
                f"Syn {n:d}",
                fontsize=8,
                horizontalalignment="right",
                verticalalignment="center",
                transform=self.P_Eff1.axarr[n, 0].transAxes,
            )

        sfi = Path(cellpath, Path(example["NoDend"]).name)
        if not sfi.is_dir():
            return
        fn = sorted(list(sfi.glob("*")))
        print("fn: ", fn)
        PD = PData()
        calx = 800.0
        nfiles = len(fn)
        for n in range(nfiles):
            sfile = Path(sfi, fn[n])
            print("N, sfile: ", n, sfile)
            if n == (nfiles - 1):
                calxv = calx
            else:
                calxv = None
            print("Calxv: ", calxv)
            self.parent.PLT.plot_traces(
                ax=self.P_Eff2.axarr[n, 0],
                fn=sfile,
                PD=PD,
                protocol="runANSingles",
                ymin=-90.0,
                ymax=20.0,
                xmin=400.0,
                xmax=900.0,
                iax=n,
                nax=len(fn),
                rep=0,
                figure=self.P_Eff2.figure_handle,
                longtitle=True,
                ivaxis=None,
                ivcolor="k",
                iv_spike_color="r",
                spike_marker_size=1.5,
                spike_marker_color="c",
                calx=calxv,
                caly=-20.0,
                calt=50.0,
                calv=20.0,
            )
            calx = None
            self.P_Eff2.axarr[n, 0].text(
                -0.02,
                0.5,
                f"Syn {n:d}",
                fontsize=8,
                horizontalalignment="right",
                verticalalignment="center",
                transform=self.P_Eff2.axarr[n, 0].transAxes,
            )

        EFP.P.figure_handle.text(
            0.99,
            0.01,
            f"Cell {cell_number:d}",
            fontsize=10,
            horizontalalignment="right",
            verticalalignment="bottom",
        )
        x = [1, 3]
        for i, data in enumerate(EFP.datasets):
            EFP.plot_dataset(
                data, plotno=i, ax=EFP.P.axarr[0, x[i]]
            )  # , title=EFP.titles[x[i]])
            EFP.P.figure_handle.text(
                0.25 + 0.5 * i,
                0.99,
                EFP.titles[x[i]],  # .replace('_', '\_'),
                transform=EFP.P.figure_handle.transFigure,
                fontsize=11,
                horizontalalignment="center",
                verticalalignment="top",
            )
        label = ["A", "B", "C", "D"]
        for i, pl in enumerate(label):
            EFP.P.figure_handle.text(
                0.05 + 0.25 * i,
                0.9,
                label[i],  # .replace('_', '\_'),
                transform=EFP.P.figure_handle.transFigure,
                horizontalalignment="right",
                verticalalignment="bottom",
                fontsize=11,
            )

        save_file = "Fig_M2.pdf"
        title = "SBEM Project Figure 2 Modeling : Efficacy"
        self.save_figure(EFP.P, save_file, title)

    def plot_efficacy_supplement(self):
        cells = grAList()
        EFP = EF.EfficacyPlots(None, hold=True, cols=1, rows=1)
        EFP.P.figure_handle.set_size_inches((8.0, 8.0))
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
        # print(xp)
        #        print(yp)

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
                Path(Figure_efficacy_supplement[cellN][simulation_experiment]).name,
            )
            if not sfiles.is_dir():
                return
            fn = sorted(list(sfiles.glob("*")))
            print("fn: ", fn)
            PD = PData()

            row_n = 1 - int(ic / 5)
            col_n = ic % 5
            calxp = 800.0
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
                parent_figure=EFP.P,
                # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
            )
            effp.append(P_Eff)
            nfn = len(fn)
            ymax = 20.0
            if simulation_experiment != "Full":
                ymax = 40.0
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
                    ymin=-90.0,
                    ymax=ymax,
                    xmin=400.0,
                    xmax=900.0,
                    iax=n,
                    nax=len(fn),
                    rep=0,
                    figure=P_Eff.figure_handle,
                    longtitle=True,
                    ivaxis=None,
                    ivcolor="k",
                    iv_spike_color="r",
                    spike_marker_size=1.5,
                    spike_marker_color="c",
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
        self.save_figure(EFP, save_file, title)

    def plot_all_revcorr(self):
        for cell in grAList():
            self.plot_revcorr(cell)
    
    def plot_revcorr(self, cellN=None):
        if cellN == None:
            cell_number = 17
            example = figure_revcorr_example[cell_number]
        else:
            cell_number = cellN
            example = figure_revcorr[cell_number]
        dBSPL = "Spont"

        run_calcs = False
        rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
        if rc_datafile.is_file():
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = pickle.load(fh)
            RCD = all_RCD_RCP[cell_number][0]
            RCP = all_RCD_RCP[cell_number][1]
            PD = PData()
        else:
            print(
                "Must run revcorr_supplement plot first, then the file we need will be present"
            )
            run_calcs = True
            cellpath = Path(
                self.config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, Path(example[dBSPL]).name)
            if not sfi.is_dir():
                return
            fn = sorted(list(sfi.glob("*")))
            (P, PD, RCP, RCD) = self.parent.PLT.compute_revcorr(
                P=None,
                gbc="17",
                fn=fn[0],
                PD=PData(),
                protocol="runANPSTH",
                revcorrtype="RevcorrSPKS",
                thr=-20.0,
                width=4.0,
            )
        box_size = 0.32
        sizer = {
            "A": {
                "pos": [0.1, box_size, 0.56, box_size],
                "labelpos": (-0.15, 1.02),
                "noaxes": True,
            },
            "B": {"pos": [0.1, box_size, 0.12, box_size], "labelpos": (-0.15, 1.02)},
            # "C": {"pos": [0.52, 0.20, 0.52, 0.28], "labelpos": (-0.15, 1.00)},
            # "D": {"pos": [0.52, 0.20, 0.08, 0.28], "labelpos": (-0.15, 1.00)},
            "C": {"pos": [0.57, box_size, 0.56, box_size], "labelpos": (-0.15, 1.02)},
            "D": {"pos": [0.57, box_size, 0.12, box_size], "labelpos": (-0.15, 1.02)},
        }  # dict pos elements are [left, width, bottom, height] for the axes in the plot. gr = [(a, a+1, 0, 1) for a in range(0, 8)] # just generate subplots - shape do not matter axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.arbitrary_grid(
            sizer,
            order="columnsfirst",
            figsize=(6, 6),
            label=True,
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

        ax = P.axdict["A"]
        summarySiteTC = self.parent.PLT.plot_revcorr2(P, PD, RCP, RCD)

        maxp = np.max(RCD.pairwise)
        psh = RCD.pairwise.shape
        pos = np.zeros((psh[0], psh[1], 2))
        for i in range(RCP.ninputs):
            for j in range(RCP.ninputs):
                # print(f"{pairwise[i,j]:.3f}", end='  ')
                pos[i, j, 0] = i + 1
                pos[i, j, 1] = j + 1

        sax = P.axdict

        sax["C"].plot(RCD.sites, RCD.participation / RCD.nspikes, "kx", markersize=5)
        sax["C"].set_ylim((0, 1.0))
        sax["C"].set_xlim(left=0)
        sax["C"].set_xlabel("# Release Sites")
        sax["C"].set_ylabel("Participation")

        sax["D"].plot(np.arange(len(RCD.ynspike)) + 1, RCD.ynspike, "k^-", markersize=5)
        sax["D"].set_ylim(0, 1.05)
        sax["D"].set_xlabel(
            f"# Inputs in [{RCD.pre_w[0]:.1f} to {RCD.pre_w[1]:.1f}] before spike"
        )
        sax["D"].set_ylabel("Cumulative Bushy Spikes witn N AN inputs")

        PH.cleanAxes(P.axarr.ravel())
        PH.noaxes(P.axdict["A"])

        PH.talbotTicks(
            sax["B"],
            tickPlacesAdd={"x": 1, "y": 2},
            floatAdd={"x": 1, "y": 2},
            # pointSize=7,
        )
        PH.talbotTicks(
            sax["C"],
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            # pointSize=7,
        )
        PH.talbotTicks(
            sax["D"],
            tickPlacesAdd={"x": 0, "y": 1},
            floatAdd={"x": 0, "y": 2},
            # pointSize=7,
        )

        if dBSPL == "Spont":
            if cellN is None:
                save_file = "Fig_M3.pdf"
            else:
                save_file = f"Fig_Revcorr/Revcorr_VCN_c{cell_number:02d}.png"
        else:
            save_file = f"Fig_M3_{dBSPL:s}.pdf"
        title = 'Reverse Correlation'
        title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.01}
        self.save_figure(P, save_file, title=title, title2=title2)
        return (summarySiteTC, RCD.sites)

    def plot_revcorr_supplement(self):

        dBSPL = "Spont"
        P = PH.regular_grid(
            rows=len(grAList()),
            cols=4,
            order="rowsfirst",
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
        PH.cleanAxes(P.axarr.ravel())
        ncells = len(grAList())
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

        run_calcs = False
        rc_datafile = Path(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl")
        if rc_datafile.is_file():
            with open(rc_datafile, "rb") as fh:
                all_RCD_RCP = pickle.load(fh)
        else:
            print(
                "Must run revcorr_supplement plot first, then the file we need will be present"
            )
            run_calcs = True
            all_RCD_RCP = {}  # all revcorr data

        for i, cell_number in enumerate(grAList()):
            # if i > 1:
            #                 continue
            dataset = figure_revcorr[cell_number]

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
            if run_calcs:
                (P0, PD, RCP, RCD) = self.parent.PLT.compute_revcorr(
                    P=None,
                    gbc=str(cell_number),
                    fn=fn[0],
                    PD=PData(),
                    protocol="runANPSTH",
                    revcorrtype="RevcorrSPKS",
                    thr=-20.0,
                    width=4.0,
                )
            else:
                PD = PData()
                RCD = all_RCD_RCP[cell_number][0]
                RCP = all_RCD_RCP[cell_number][1]

            if i == 0:
                tcal = True  # point font for cal bar
            else:
                tcal = False  # no cal bar

            all_RCD_RCP[cell_number] = [RCD, RCP]
            summarySiteTC = self.parent.PLT.plot_revcorr2(
                P,
                PD,
                RCP,
                RCD,
                axarray=P.axarr[i, 0:2],
                calbar_show=tcal,
                calbar_fontsize=7,
                yaxis_label=False,
            )

            P.axarr[i, 0].text(
                -0.25,
                0.5,
                f"VCN_c{cell_number:02d}",
                fontsize=9,
                color="k",
                transform=P.axarr[i, 0].transAxes,
                horizontalalignment="right",
            )

            maxp = np.max(RCD.pairwise)
            psh = RCD.pairwise.shape
            pos = np.zeros((psh[0], psh[1], 2))
            for j in range(RCP.ninputs):
                for k in range(RCP.ninputs):
                    # print(f"{pairwise[i,j]:.3f}", end='  ')
                    pos[j, k, 0] = j + 1
                    pos[j, k, 1] = k + 1

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

            PH.noaxes(P.axarr[i, 1])

            sax[2].plot(
                RCD.sites,
                RCD.participation / RCD.nspikes,
                "kx",
                markersize=5,
                clip_on=False,
            )
            sax[2].set_ylim((0, 1.0))
            sax[2].set_xlim(0, 225)

            if i == ncells - 1:
                sax[2].set_xlabel("# Release Sites")
            # sax[2].set_ylabel("Participation")
            PH.talbotTicks(  # Participation
                sax[2],
                tickPlacesAdd={"x": 0, "y": 2},
                floatAdd={"x": 0, "y": 2},
                # pointSize=7,
            )
            sax[2].tick_params(direction="in", length=3.0, width=1.0)

            sax[3].plot(
                np.arange(len(RCD.ynspike)) + 1,
                RCD.ynspike,
                "k^-",
                markersize=5,
                clip_on=False,
            )
            sax[3].set_ylim(0, 1.05)
            if i == ncells - 1:
                sax[3].set_xlabel(
                    f"# Inputs in [{RCD.pre_w[0]:.1f} to {RCD.pre_w[1]:.1f}] before spike"
                )
            PH.talbotTicks(  # cumulative inputs before spike
                sax[3],
                tickPlacesAdd={"x": 0, "y": 1},
                floatAdd={"x": 0, "y": 2},
                # pointSize=7,
            )
            sax[3].tick_params(direction="in", length=3.0, width=1.0)

        # save the accumulated RCD data
        with open(rc_datafile, "wb") as fh:
            pickle.dump(all_RCD_RCP, fh)
        title = ("SBEM Project Figure 3 Modeling : Reverse Correlation Summary",)
        save_file = f"Fig_M3_supplemental_Full_{dBSPL:s}.pdf"
        self.save_figure(P, save_file, title=title)

    def plot_revcorr_compare(self):
        dBSPLs = ["Spont", "40dB"]
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
        for i, dBSPL in enumerate(dBSPLs):
            with open(f"GradeA_RCD_RCP_all_revcorrs_{dBSPL:s}.pkl", "rb") as fh:
                try:
                    R = pickle.load(fh)
                except:
                    print(
                        "Must run revcorr_supplement plot first, then the file we need will be present"
                    )
                    return

            for cell_number in R.keys():
                RCD = R[cell_number][0]
                RCP = R[cell_number][1]
                P.axarr[0, i].plot(
                    np.arange(len(RCD.ynspike)) + 1,
                    RCD.ynspike,
                    "^-",
                    markersize=4,
                    clip_on=False,
                    label=f"VCN_c{cell_number:02d}",
                )
            P.axarr[0, i].set_xlim(0, 12)
            P.axarr[0, i].set_ylim(0, 1.0)
            P.axarr[0, i].set_xlabel("Number of inputs prior to spike")
            P.axarr[0, i].set_ylabel("Cumulative Fraction of Spikes")
            if i == 0:
                P.axarr[0, i].legend()
            P.axarr[0, i].set_title(dBSPL)

        save_file = f"Fig_M4.pdf"
        title = "SBEM Project Figure 4 Modeling : Reverse Correlation Comparison"
        self.save_figure(P, save_file, title=title)

    def plot_psth_psth(
        self,
        ax: object,
        data: object,
        ri: dict,
        time_base: np.ndarray,
        psth_binw: float = 0.5,
        wstart: float = 0.0,
        wdur: float = 1.0,
        ntr: int = 1,
        ninputs: int = 1,
    ):

        self.parent.PLT.plot_psth(
            data,
            run_info=ri,
            max_time=np.max(time_base),
            bin_width=psth_binw,
            ax=ax,
            zero_time=wstart,
            scale=1.0 / ntr / psth_binw / ninputs,
        )
        ax.set_xlim(0, wdur)
        PH.talbotTicks(
            ax,
            axes="xy",
            density=(1.0, 1.0),
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
        label_x_axis=True,
    ):

        grand_fsl, grand_ssl = self.parent.PLT.plot_fsl_ssl(
            an_st_grand,
            run_info=ri,
            max_time=25.0,
            min_fsl=1.2e-3,
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
        dBSPL = "40dB"
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

        self.plot_one_PSTH(
            cell_number=cell_number,
            dBSPL=dBSPL,
            tr_ax=tr_ax,
            st_ax=st_ax,
            bupsth_ax=bupsth_ax,
            anpsth_ax=anpsth_ax,
            bufsl_ax=bufsl_ax,
            anfsl_ax=anfsl_ax,
        )

        if cellN is None:
            save_file = f"Fig_M5.pdf"
        else:
            save_file = f"All_PSTH/PSTH_VCN_c{cell_number:02d}.png"
        title = "SBEM Project Figure 5Modeling : PSTH and FSL, Example"
        title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.01}
        self.save_figure(P, save_file, title=title, title2=title2)

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
        label_x_axis=True,
    ):
        dataset = figure_psth[cell_number]
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
        protocol = "runANPSTH"
        AR, SP, RMA = self.parent.PLT.analyze_data(ivdatafile, filemode, protocol)
        ntr = len(AR.traces)  # number of trials
        v0 = -160.0
        trstep = 25.0 / ntr
        inpstep = 5.0 / ntr
        sz = 50.0 / ntr
        # print(dir(AR))
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

        ntr = len(AR.traces)
        all_bu_st = []
        wstart = 0.15
        wdur = 0.3
        for i in range(ntr):  # for all trials in the measure.
            vtrial = AR.traces[i] * 1e3
            time_base = AR.time_base / 1000.0  # convert to seconds
            dt = si.dtIC / 1000.0  # convert from msec to seconds
            trd = d["Results"][i]
            ninputs = len(trd["inputSpikeTimes"])
            if i == 0:
                all_an_st = [[] for x in range(ninputs)]  # by input, then trial
                an_st = [[] for x in range(ntr)]  # by trial, collapsed across inputs

            waveform = trd["stimWaveform"].tolist()
            stb = trd["stimTimebase"]
            if not isinstance(trd["spikeTimes"], list):
                cprint("r", "spiketimes is not a list")
                return
            all_bu_st.append(trd["spikeTimes"])
            if i == 0:
                tr_ax.plot((time_base - wstart) * 1e3, vtrial, "k-", linewidth=0.5)
                spiketimes = np.array(trd["spikeTimes"])
                # trim spike mark array so we don't get spots at the edge of the plot
                spikeindex = [
                    int(t / dt)
                    for t in spiketimes
                    if (t >= wstart and t < ((wstart + wdur)))
                ]
                tr_ax.plot(
                    (time_base[spikeindex] - wstart) * 1e3,
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
                tr_ax.set_xlim(0.0, wdur * 1e3)
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
                st_ax.plot(
                    stb - wstart, np.array(waveform) * 1e3, "k-", linewidth=0.5
                )  # stimulus underneath
                st_ax.set_xlim(0, wdur)
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
                wstart=wstart,
                wdur=wdur,
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
            self.plot_psth_psth(
                ax=anpsth_ax,
                data=all_an_st,
                ri=ri,
                time_base=time_base,
                psth_binw=psth_binw,
                wstart=wstart,
                wdur=wdur,
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
                min_fsl=2e-3,
                ax=bufsl_ax,
                zero_time=wstart * 1e-3,
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
            self.plot_psth_ANFSL(
                anfsl_ax, wstart, cell_number, ri, an_st_grand, label_x_axis
            )

    def plot_PSTH_supplement(self):
        print("PSTH supplement")
        dBSPL = "40dB"
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
            self.plot_one_PSTH(
                cell_number=cell_number,
                dBSPL=dBSPL,
                tr_ax=P.axarr[i, 0],
                st_ax=None,
                bupsth_ax=P.axarr[i, 1],
                anpsth_ax=None,
                bufsl_ax=P.axarr[i, 2],
                anfsl_ax=None,
                label_x_axis=show_label,
            )

            P.axarr[i, 0].text(
                -0.25,
                0.5,
                f"VCN_c{cell_number:02d}",
                fontsize=9,
                color="k",
                transform=P.axarr[i, 0].transAxes,
                horizontalalignment="right",
            )

        save_file = f"Fig_M5_Supplmental.pdf"
        title = "SBEM Project Figure 5 Modeling Supplemental : PSTH and FSL, All cells"
        self.save_file(P, save_file, title)

    def plot_VS_SAM(self):
        self.generate_VS_data()
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
            print(filename)
            cellpath = Path(
                config["cellDataDirectory"],
                f"VCN_c{cell_number:02d}",
                "Simulations",
                "AN",
            )
            sfi = Path(cellpath, filename + ".pkl")
            with open(sfi, "rb") as fh:
                d = pickle.load(fh)

            self.parent.PLT.plot_AN_response(P, d.files[0], PD, "runANPSTH")
            with open(fout, "a") as fth:
                if firstline:
                    fth.write(self.parent.PLT.VS_colnames)
                    fth.write("\n")
                    firstline = False
                fth.write(self.parent.PLT.VS_line)
                fth.write("\n")

        print("anvs done for cell: ", cell_number)

    #
    def generate_VS_data(self):
        if "vcnmodel.VS_datasets" not in list(dir()):
            from vcnmodel import VS_datasets as VS_datasets
        print(dir())
        importlib.reload(VS_datasets)
        config = toml.load(open("wheres_my_data.toml", "r"))

        """
        Generate the table in VS_data.py by analyzing the data from 
        VS_datasets.py
        """

        fout = "VS_data.py"  # we will generate this automatically
        with open(fout, "w") as fh:
            fh.write(f'"""\n')
            fh.write(
                "    Vector strength for models with SAM tones, different input configurations.\n"
            )
            fh.write("    9 Dec 2020 version.\n")
            fh.write(
                "    Results are printout from DataTablesVCN after selecting the data runs.\n"
            )
            fh.write(f'    pbm\n"""\n')
            fh.write('\ndata = """')

        for i, celln in enumerate(grAList()):
            self.analyze_VS_data(VS_datasets, celln, fout, firstline=True)
        with open(fout, "a") as fh:
            fh.write(f'"""\n')

    def plot_VC_gKLT(self, parent_figure=None, loc:Union[None, tuple]=None):
        cell_number = 17
        dataset = figure_VClamp[cell_number]

        cellpath = Path(
            self.config["cellDataDirectory"],
            f"VCN_c{cell_number:02d}",
            "Simulations",
            "VC",
        )
        sfi = []
        for i, ds in enumerate(dataset):
            print(i, ds)
            sfd = Path(cellpath, Path(ds).name)
            print(sfd)
            if not sfd.is_dir():
                print("Directory not found: ", str(sfd))
                return
            fn = sorted(list(sfd.glob("*")))[0]
            sfi.append(fn)
        P = self.parent.PLT.plot_VC(sfi=sfi, show=False, parent_figure=parent_figure,
            loc=loc)
        save_file = f"Fig_M0_VC_Adjustment.pdf"
        fig = FigInfo()
        fig.P = P
        if cell_number == None:
            fig.filename = "Fig_M1.pdf"
        else:
             fig.filename  = save_file
        fig.title = "SBEM Project Figure Modeling Supplemental : VC"
        fig.title2 = {"title": f"Cell {cell_number:d}", "x": 0.99, "y": 0.05}
        return fig