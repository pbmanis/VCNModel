import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import matplotlib
import matplotlib.pyplot as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import image as mpimg
from pylibrary.plotting import plothelpers as PH

import toml

import seaborn as sns
import vcnmodel.plotters.plot_z as PZ
import vcnmodel.plotters.efficacy_plot as EF

config = toml.load(open("wheres_my_data.toml", "r"))


def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """
    return [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]


syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x"]
colors = ["c", "k", "m"]

title_text = {"passive": "Passive", "normal": "Half-active", "active": "Active"}
font_colors = {"passive": "c", "normal": "k", "active": "m"}
spike_colors = font_colors


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


Figure_IV = {
    "Cell": 17,
    "normal": "runIV-all-2020-07-30.18-29-53",
    "passive": "runIV-all-2020-07-30.18-43-29",
    "active": "runIV-all-2020-07-30.18-56-13",
    "Z_normal": "VCN_c17_Full_normal_Z.pkl",
    "Z_passive": "VCN_c17_Full_pasdend_Z.pkl",
    "Z_active": "VCN_c17_Full_actdend_Z.pkl",
}
Figure_AllIVs = {
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
    30: {
        "normal": "runIV-all-2020-07-30.18-31-35",
        "passive": "runIV-all-2020-07-30.18-45-12",
        "active": "runIV-all-2020-07-30.18-57-54",
    },
}


class Figures(object):
    """
    Entry point.
    This generates final figures for a  manuscript by reaching back to the 
    original simulation data, including, in some cases, refitting.
    """

    def __init__(self, parent):
        self.parent = parent  # point back to caller's space
        self.config = toml.load(
            open("wheres_my_data.toml", "r")
        )  # sorry, have to reload it here.

    def reset_style(self):
        sns.set_style(rc={"pdf.fonttype": 42})
        mpl.style.use('~/.matplotlib/figures.mplstyle')

    def make_figure(self, figure_name: Union[str, None] = None):
        self.reset_style()
        print("make_figure:", figure_name)
        # dispatch
        dispatch_table = {
            "IV Figure": self.plotIV,
            "IV Supplement": self.plot_IVS,
            "Zin Supplement": self.plot_ZinS,
            "Efficacy": self.plot_efficacy,
            "Efficacy Supplement": self.plot_efficacy_supplement
            #
            # "Revcorr Ex", "Revcorr Supplement", "Revcorr Compare",
            # "PSTH/VS", "FSL", "VS-SAM Tone"}
        }
        print(figure_name, dispatch_table.keys())
        if figure_name in list(dispatch_table.keys()):
            print("disp")
            dispatch_table[figure_name]()

    def force_log_ticks(self, ax):
        locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
        ax.xaxis.set_major_locator(locmaj)

        locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                              numticks=100)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        

    def plotIV(self):
        # style = STY.styler(journal="JNeurosci", figuresize='full', font='stixsans')  # two colukn...
        cellN = Figure_IV["Cell"]
        cellpath = Path(
            self.config["cellDataDirectory"], f"VCN_c{cellN:02d}", "Simulations", "IV"
        )
        d1 = Figure_IV["passive"]
        rows = 2
        cols = 3
        height = 4.0
        width = 8.0
        PD = PData()
        ymin = -125.0
        ymax = 40.0
        self.P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            figsize=(width, height),
            showgrid=False,
            verticalspacing=0.06,
            horizontalspacing=0.06,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.08,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=None,
            panel_labels=["A", "B", "C", "D", "E", "F"],
        )

        self.P2 = PH.regular_grid(
            1,
            4,
            order="rowsfirst",
            showgrid=False,
            verticalspacing=0.06,
            horizontalspacing=0.075,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.52,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=self.P,
            panel_labels=["D", "E", "F", "G"],
        )
        # x=-0.05
        # y =1.05
        # self.P3 = PH.Plotter(
        #     (2, 1),
        #     arrangement={
        #         "E1": {"pos": [0.27, 0.3, 0.25, 0.18], "labelpos": (x,y), "noaxes": False},
        #         "E2": {"pos": [0.27, 0.1, 0.25, 0.18], "labelpos": (x,y), "noaxes": False},
        #     },
        #     parent_figure=self.P2,
            # order="rowsfirst",
            # showgrid=True,
            # verticalspacing=0.02,
            # horizontalspacing=0.075,
            # margins={
            #     "bottommargin": 0.1,
            #     "leftmargin": 0.3,
            #     "rightmargin": 0.58,
            #     "topmargin": 0.52,
            # },
            # labelposition=(-0.1, 1.06),
            # parent_figure=self.P,
            # panel_labels=["E1", "E2"],
        # )

        self.P3 = PH.regular_grid(
            2,
            1,
            order="rowsfirst",
            showgrid=False,
            verticalspacing=0.02,
            horizontalspacing=0.075,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.3,
                "rightmargin": 0.50,
                "topmargin": 0.52,
            },
            labelposition=(-0.25, 0.95),
            parent_figure=self.P,
            panel_labels=["E1", "E2"],
        )

        for iax, iv in enumerate(["passive", "normal", "active"]):
            sfi = Path(cellpath, Path(Figure_IV[iv]).name)
            print("sfi: ", sfi)
            print("is it a dir? ", sfi.is_dir())
            if not sfi.is_dir():
                continue
            fn = list(sfi.glob("*"))
            print("fn: ", fn)
            sfi = Path(sfi, fn[0])
            print("sfi: ", sfi)
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
                calx = 120.,
                caly = -10.,
            )
            self.P.axarr[0, iax].set_title(
                title_text[iv], color=font_colors[iv]
            )

        res_label = r"$\mathregular{R_{in} (M\Omega)}$"
        tau_label = r"$\mathregular{\tau_{m} (ms)}$"
        phase_label = r"$\mathregular{\phi (radians)}$"
        for iax, iv in enumerate(["Z_passive", "Z_normal", "Z_active"]):
            if iv not in Figure_IV.keys():
                continue
            sfi = Path(config["codeDirectory"], Figure_IV[iv])
            if not sfi.is_file():
                continue
            print("sif: ", sfi)
            # ax = self.P2.axarr[0, 1]
            ax = self.P3.axdict['E1']
            label = sfi.name # .replace("_", "\_")
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
            secax= self.P3.axdict['E2']
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
        for rax, iv in enumerate(Figure_AllIVs.keys()):
            # cellpath = Path(self.config['cellDataDirectory'], f"VCN_c{iv:02d}", 'Simulations', 'IV')
            for iax, dendmode in enumerate(["passive", "normal", "active"]):
                sfi = Path(
                    self.config["cellDataDirectory"],
                    f"VCN_c{iv:02d}",
                    "Simulations",
                    "IV",
                    Figure_AllIVs[iv][dendmode],
                )
                print(sfi)
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
            ax=self.P2.axdict["F"],
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_rin, x="dendrites", y="Rin", ax=self.P2.axdict["F"], color="0.6"
        )
        self.P2.axdict["F"].set_ylim(0, 30.0)
        self.P2.axdict["F"].set_ylabel(res_label)

        sns.boxplot(
            data=df_tau,
            x="dendrites",
            y="taum",
            ax=self.P2.axdict["G"],
            saturation=0.3,
            palette=colors,
        )
        sns.stripplot(
            data=df_tau, x="dendrites", y="taum", ax=self.P2.axdict["G"], color="0.6"
        )
        self.P2.axdict["G"].set_ylim(0, 2.0)
        self.P2.axdict["G"].set_ylabel(tau_label)


        save_file = "Fig_M1.pdf"
        self.P.figure_handle.text(
            0.98,
            0.98,
            save_file, # .replace('_', '\_'),
            transform=self.P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure 1 Modeling (Main)",
            },
        )
        self.P.figure_handle.show()

    def plot_IVS(self):
        """
        All of the IVS, for a supplemental figure
        Passive, normal, active, plus the crossed IV
        Also put the PNG for the cell on the left.
        """
        nivs = len(Figure_AllIVs)
        rows = nivs
        cols = 5
        height = 1.5 * nivs
        width = 8.5
        PD = PData()
        ymin = -125.0
        ymax = 40.0
        calx = 120.
        
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
            parent_figure=None,
            # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        cellpath = config["cellDataDirectory"]
        png_path = Path(config["baseDataDirectory"], config["pngDirectory"])
        for rax, iv in enumerate(Figure_AllIVs.keys()):
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
                    Figure_AllIVs[iv][dendmode],
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
                    caly=-10.,
                )
                if rax == 0:
                    self.P.axarr[rax, iax + 1].set_title(dendmode)
                if iax == 0:
                    self.P.axarr[rax, 0].text(-0.1, 0.5, str(iv))
        save_file = "Fig_M1A_Supplemental.pdf"
        self.P.figure_handle.text(
            0.98,
            0.98,
            save_file, # .replace('_', '\_'),
            transform=self.P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure 1 Modeling (Supplemental A)",
            },
        )
        self.P.figure_handle.show()

    def plot_ZinS(self):
        """
        All of the Zins for a supplemental figure
        """
        PZ.PlotZ()
        
    def plot_efficacy(self):
        example = {'Cell': 17,
                   'normal': 'runANSingles-all-2020-11-04.14-18-08',
                   'NoDend': 'runANSingles-all-2020-11-04.09-44-00',
                    }
        cellN = example['Cell']
        cellpath = Path(
            self.config["cellDataDirectory"], f"VCN_c{cellN:02d}", "Simulations", "AN"
        )        
        sfi = Path(cellpath, Path(example['normal']).name )
        if not sfi.is_dir():
            return
        fn = sorted(list(sfi.glob("*")))
        print("fn: ", fn)
        PD = PData()
        calx = 20.
        
        EFP = EF.EfficacyPlots(None, hold=True, cols=4)

        self.P_Eff1 = PH.regular_grid(
            len(fn),
            1,
            order="rowsfirst",
            figsize=(8., 4.5),
            showgrid=True,
            verticalspacing=0.01,
            horizontalspacing=0.1,
            margins={
             "bottommargin": 0.1,
             "leftmargin": 0.05,
             "rightmargin": 0.75,
             "topmargin": 0.08,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=EFP.P,
         # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        self.P_Eff2 = PH.regular_grid(
            len(fn),
            1,
            order="rowsfirst",
            figsize=(8., 4.5),
            showgrid=True,
            verticalspacing=0.01,
            horizontalspacing=0.1,
            margins={
             "bottommargin": 0.1,
             "leftmargin": 0.5,
             "rightmargin": 0.3,
             "topmargin": 0.08,
            },
            labelposition=(-0.05, 1.06),
            parent_figure=EFP.P,
         # panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        )
        for n in range(len(fn)):
            sfile = Path(sfi, fn[n])
            print(n, sfile)
            self.parent.PLT.plot_traces(
                ax = self.P_Eff1.axarr[n,0],
                fn = sfile,
                PD = PD,
                protocol= 'runANSingles',
                ymin=-90.0,
                ymax=20.0,
                xmin = 400., 
                xmax = 900.,
                iax=n,
                nax=len(fn),
                rep=0,
                figure=self.P_Eff1.figure_handle,
                longtitle=True,
                ivaxis=None,
                ivcolor='k',
                iv_spike_color='r',
                spike_marker_size=2.5,
                spike_marker_color='g',
                calx = calx,
                caly = 0.,
            )
            calx = None
            self.P_Eff1.axarr[n, 0].text(-0.02, 0.5, f"Input {n:d}",
                fontsize=8,
                horizontalalignment='right',
                verticalalignment="center",
                transform = self.P_Eff1.axarr[n, 0].transAxes)

        sfi = Path(cellpath, Path(example['NoDend']).name )
        if not sfi.is_dir():
            return
        fn = sorted(list(sfi.glob("*")))
        print("fn: ", fn)
        PD = PData()
        calx = 20.
        
        for n in range(len(fn)):
            sfile = Path(sfi, fn[n])
            print(n, sfile)
            self.parent.PLT.plot_traces(
                ax = self.P_Eff2.axarr[n,0],
                fn = sfile,
                PD = PD,
                protocol= 'runANSingles',
                ymin=-90.0,
                ymax=20.0,
                xmin = 400., 
                xmax = 900.,
                iax=n,
                nax=len(fn),
                rep=0,
                figure=self.P_Eff2.figure_handle,
                longtitle=True,
                ivaxis=None,
                ivcolor='k',
                iv_spike_color='r',
                spike_marker_size=2.5,
                spike_marker_color='g',
                calx = calx,
                caly = 0.,
            )
            calx = None
            self.P_Eff2.axarr[n, 0].text(-0.02, 0.5, f"Inputx {n:d}",
                fontsize=8,
                horizontalalignment='right',
                verticalalignment="center",
                transform = self.P_Eff2.axarr[n, 0].transAxes)
        x = [1, 3]
        for i, data in enumerate(EFP.datasets):
            EFP.plot_dataset(data, plotno=i, ax=EFP.P.axarr[0,x[i]], title=EFP.titles[x[i]])

        save_file = "Fig_M2.pdf"
        EFP.P.figure_handle.text(
            0.99,
            0.99,
            save_file, # .replace('_', '\_'),
            transform=EFP.P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure 2 Modeling : Efficacy",
            },
        )
        EFP.P.figure_handle.show()

    
    def plot_efficacy_supplement(self):
        pass
        
