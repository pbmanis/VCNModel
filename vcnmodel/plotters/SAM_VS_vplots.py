"""
Plot responses to SAM stimuli and the vector strengths
Assumes that the analysis has already been done, and that the
location of that data is given in the 'wheres_my_data.toml' file

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""
import importlib
import io
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
import seaborn as sns
from vcnmodel.util.set_figure_path import set_figure_path

sys.path.insert(0, os.path.abspath("nb"))
# import matplotlib.collections
import plotnine as PN
import toml
import VS_data_15dB as VS_data_15dB
import VS_data_15dB_BC09 as VS_data_15dB_BC09
import VS_data_30dB as VS_data_30dB
from matplotlib import pyplot as mpl
from matplotlib import ticker
from matplotlib.lines import Line2D
from pylibrary.plotting import plothelpers as PH

# seaborn default palette, first 10 colors
colors = [
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

colors_swarm = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
]

def reset_style():
    #     sns.set_style(rc={"pdf.fonttype": 42})
    mpl.style.use("~/.matplotlib/figures.mplstyle")


def defemptylist():
    return []


@dataclass
class PlotInfo:
    plotwhat: str = "VectorStrength"
    labels: list = field(default_factory=defemptylist)
    anpt: str = ""
    labels: list = field(default_factory=defemptylist)
    yscale: list = field(default_factory=defemptylist)
    freqs: list = field(default_factory=defemptylist)
    fr: list = field(default_factory=defemptylist)
    vsscale: list = field(default_factory=defemptylist)


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


class VS_Plots:
    def __init__(
        self, output="P9", sels=[2, 5, 6, 9, 10, 11, 13, 17, 18, 30], dBSPL=15, 
        dends:str=""
    ):
        """Vector strength plotting class

        Parameters
        ----------
        output : str, optional
            which plotter to use, by default "P9"
        sels : list, optional
            cell selection to plot, by default [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
        dBSPL : int, optional
            The base intensity for the VS stimuli, by default 15
        """
        self.dBSPL = dBSPL
        if dends == "9I9U":
            importlib.reload(VS_data_15dB_BC09)
            self.datas = VS_data_15dB_BC09.data
        elif dBSPL == 15:
            importlib.reload(VS_data_15dB)
            self.datas = VS_data_15dB.data
        elif dBSPL == 30:
            importlib.reload(VS_data_30dB)
            self.datas = VS_data_30dB.data
        self.output = output
        self.cell_inputs = {
            "2": "subthreshold",
            "5": "subthreshold",
            "6": "subthreshold",
            "10": "subthreshold",
            "30": "subthreshold",
            "9": "one",
            "13": "one",
            "18": "one",
            "11": "two",
            "17": "two",
        }
        self.cell_list = []
        with open("wheres_my_data.toml", "r") as fh:
            self.config = toml.load(fh)

        df = self.prepare_data(self.datas)

        df["strength"] = df["Cell"].apply(
            lambda x: "subthreshold"
            if x in [5, 10, 2, 6, 30]
            else "one"
            if x in [9, 13, 18]
            else "two"
        )
        df = df.sort_values(
            ["Cell", "frequency", "Configuration"], ascending=(True, True, True)
        )
        df = df.reset_index()

        self.PI = PlotInfo()

        reset_style()

        # plotwhat = "phasesdtime"
        # plotwhat = 'phasesd'
        plotwhat = "VectorStrength"
        # plotwhat = 'SpikeCount'
        # plotwhat = 'logSpikeCount'
        self.PI.plotwhat = plotwhat

        if self.PI.plotwhat == "phasesdtime":
            df["phasesdtime"] = 1e3 * df["phasesd"] / (2 * np.pi * df["frequency"])
            df["AN_phasesdtime"] = (
                1e3 * df["AN_phasesd"] / (2 * np.pi * df["frequency"])
            )

        if self.PI.plotwhat == "logSpikeCount":
            df["logSpikeCount"] = np.log10(df["SpikeCount"])


        labels = set(df["frequency"])
        labelnames = [float(label) for label in labels]

        self.PI.labels = [str(int(l)) for l in sorted(labelnames)]

        self.PI.anpt = "AN_phasesd"
        if plotwhat == "phasesdtime":
            self.PI.yscale = [0.0, 1.5]
            self.PI.anpt = "AN_phasesdtime"
        elif plotwhat == "SpikeCount":
            self.PI.yscale = [0, 20000]
        elif plotwhat == "logSpikeCount":
            self.PI.yscale = [0, 5]
        elif plotwhat == "phasesd":
            self.PI.yscale = [0, 2 * np.pi]
            self.PI.anpt = "AN_phasesd"
        else:
            self.PI.yscale = [-0.1, 1.0]
            self.PI.anpt = "AN_VS"

        self.PI.freqs = [50, 100, 200, 300, 400, 500, 750, 1000]
        self.PI.fr = self.fscale(self.PI.freqs)
        self.PI.vsscale = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        self.df = df[df["Cell"].isin(sels)]

    def make_figure(self, figure=None, axs=None):
        if self.output == "P9":
            fig, P = self.plot_with_p9(self.df, self.PI, figure=figure, axs=axs)
            return fig, P

    def prepare_data(self, datas):
        sio = io.StringIO(datas)
        df = pd.read_table(sio, sep=",")
        return df

    def fscale(self, x):
        return [str(int(y)) for y in x]

    def vscale(self, x):
        return [f"{y:.1f}" for y in x]

    def getCell(self, cell):
        self.cell_list.append(cell)
        return f"BC{int(cell):02d}"

    def scalefun(self, x):
        return [f"{y:.2f}" for y in x]

    def plot_with_p9(self, df, PI, figure=None, axs=None):
        if figure is None:
            PN.options.figure_size = (10, 10)
        dodge = PN.positions.position_dodge(width=0.6)
        themer = PN.themes.theme
        gg = (
            PN.ggplot(
                df,
                PN.aes(
                    x="factor(frequency)",
                    y=PI.plotwhat,
                    group="factor(frequency)",
                    color="factor(Configuration)",
                ),
            )
            + PN.scale_color_manual(values=colors)
            + PN.scale_x_discrete(breaks=PI.freqs, labels=PI.fr)
            + PN.scale_y_continuous(breaks=PI.vsscale, labels=self.vscale(PI.vsscale))
            # + PN.scales.scale_y_log10()
            + PN.geom_point(
                PN.aes(
                    x="factor(frequency)",
                    y=PI.plotwhat,
                    inherit=True,
                    #                         size=1.5,
                ),
                stroke=0,
                position=PN.positions.position_dodge2(width=1),
                # position=PN.positions.position_jitterdodge(jitter_width=1, dodge_width=1),
                alpha=1,
                size=2.5,
            )
            + PN.labs(
                x="Frequency (Hz)", y="Vector Strength", color="Synapse Configuration"
            )
            + PN.geom_line(
                PN.aes(
                    x="factor(frequency)",
                    y=PI.plotwhat,
                ),
                position=PN.positions.position_dodge2(width=1),
                alpha=1,
                color="grey",
            )
            + PN.geom_line(
                PN.aes(
                    x="factor(frequency)",
                    y=PI.anpt,
                ),
                position=PN.positions.position_dodge2(width=1),
                alpha=1,
                color="red",
            )
            + PN.scales.scale_colour_brewer(type="qual", palette="Dark2")
            + PN.facet_wrap("Cell", labeller=self.getCell, nrow=5, ncol=2, dir="v")
            + PN.theme_minimal()  # (style='whitegrid')
            + PN.theme(
                text=PN.element_text(
                    family="Arial",
                )
            )
            + PN.themes.theme(panel_grid_minor=PN.element_blank())
            + PN.themes.theme(panel_grid_major=PN.element_blank())
            + PN.themes.theme(
                panel_spacing_x=0.5,
            )
            + PN.themes.theme(
                panel_spacing_y=0.25,
            )
            + PN.themes.theme(axis_line=PN.element_rect())
            + PN.themes.theme(axis_line=PN.element_rect())
            + PN.themes.theme(
                axis_ticks_length=20,
            )
            + PN.themes.theme(
                axis_ticks_length_minor=0,
            )
            + PN.themes.theme(
                axis_ticks_length_major=3,
            )
            + PN.themes.theme(axis_ticks_direction="in")
            + PN.themes.theme(
                axis_ticks=PN.element_line(size=1),
            )
            + PN.themes.theme(
                axis_ticks_pad=20,
            )
            + PN.themes.theme(
                axis_ticks_major=PN.element_line(color="k"),
            )
            + PN.themes.theme(
                axis_ticks_major_x=PN.element_line(size=1),
            )
            + PN.themes.theme(
                axis_ticks_major_y=PN.element_line(color="k"),
            )
            + PN.themes.theme(
                axis_ticks_minor=PN.element_line(color="k"),
            )
            + PN.themes.theme(
                axis_ticks_minor_x=PN.element_line(size=0),
            )
            + PN.themes.theme(
                axis_ticks_minor_y=PN.element_line(color="k"),
            )
            + PN.themes.theme(axis_text_x=PN.element_text(angle=45))
            #       + PN.scale_x_discrete(labels=labels)
            + PN.scales.ylim(PI.yscale)
            + PN.themes.theme(legend_position=(0.5, 0.93), legend_entry_spacing=0.001)
            + PN.theme(figure_size=(8, 10))
            # + PN.scales.scale_y_continuous(breaks=[f"{f:.1f}" for f in [0, 0.2, 0.4 0.6, 0.8, 1.0]]) # scalefun) # PN.scales.number_format(accuracy = 0.01))
        )
 
        fig, P = gg.draw(return_ggplot=True)
 
        # save_file = Path(f"Figure7/Figure7_supp/Figure7_Main_RHS_top.pdf")
        save_file = Path(f"Figure7/Figure7_supp/Figure7_Supplemental2_V3.pdf")
        save_file_png = save_file.with_suffix(".png")

        title = "SBEM Project Supplemental Figure 7 Modeling : Vector Strength Summary"

        mpl.savefig(
            Path(self.config["baseDataDirectory"], self.config["figureDirectory"], save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )
        mpl.savefig(
            Path(self.config["baseDataDirectory"], self.config["figureDirectory"], save_file_png),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )

        fig.filename = save_file
        if hasattr(fig, "title"):
            fig.title["title"] = title
        return fig, P

    def plot_VS_Data(self, axin=None, legendflag=True):
        """Plot the VS response to SAM stimuli for each cells and frequency in a swarmplot

        Parameters
        ----------
        axin : matplotlib axis instance, optional
            The axis to plot the figure into, by default None

        legend: bool
            Set true to include legend in panel
        Returns
        -------
        matplotlib axis instance
            the axis into which the data was plotted; if none was
            specified, then this is the axis that was generated.
        """
        label_font = {"fontsize": 8, "fontweight": "normal"}
        title_font = {"fontsize": 9, "fontweight": "normal"}
        if axin is None:
            fig, ax = mpl.subplots(1, 1, figsize=(10, 4))
        else:
            ax = axin
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3.0)
        ax.set_ylim(0, 1)
        ax.set_xlim(-1, 10)
        sax = sns.swarmplot(
            x="frequency",
            y="VectorStrength",
            hue="Configuration",
            palette=colors_swarm,
            data=self.df,
            dodge=True,
            size=4,
            ax=ax,
        )
        nfreq = len(self.df["frequency"].unique())
        nconfig = len(self.df["Configuration"].unique())
        x = np.arange(nfreq)
        if nconfig == 3:
            anline = 320
        else:
            anline = 75
        # print(self.df["AN_VS"][:24:3])
        y = self.df["AN_VS"][::nconfig]
        # print(x)
        # print(y)
        # ax.scatter(
        #     x,
        #     y,
        #     s=anline,
        #     color="firebrick",
        #     marker="_",
        #     linewidth=2,
        #     legend="AN"
        # )
        xl = ax.get_xticklabels()
        newxl = [f"{int(float(x.get_text())):d}" for x in xl]
        ax.set_xticklabels(newxl)
        ax.set_ylabel("Vector Strength)", fontdict=label_font)
        ax.set_xlabel("Modulation Frequency (Hz)", fontdict=label_font)
        ax.set_title(
            f"{int(self.dBSPL):2d} dBSPL, 100% modulation", fontdict=title_font
        )
        if not legendflag:
            ax.legend(handles=[], labels= [])
            ax.get_legend().remove()
        else:
            ax.legend(markerscale=0.5)
        if axin is None:
            mpl.show()
        else:
            ax = axin
        return ax

    def plot_VS_summary(self, cell, axin=None, legendflag=True, legend_type="inputs",
                        barwidth=240, show_cell=True,
                        mode:str="line"):
        """Summarize the VS for one cell for all conditions/frequencies.
        This is similar to the plotnine routine above, but just
        uses matplotlib to acheive the same plot

        Parameters
        ----------
        cell : Union[int, str]
            cell number to pull data from pandas data frame to plot
            If str, probably 9U or 9I, then use those 2 datasets insteead of "all" etc.
        axin : matplotlib axis object, optional
            The axis to plot into, by default None
            if None, we just make our own new figure and plot to that.
        legendflag : bool, optional
            Defines whether to include a legend on the plot, by default True
        legendtypes : str
            Inputs: use the input experiment description
            Pruning:  use the pruning experiment description
        """
        label_font = {"fontsize": 8, "fontweight": "normal"}
        title_font = {"fontsize": 9, "fontweight": "normal"}

        if axin is None:
            fig, ax = mpl.subplots(1, 1, figsize=(5, 5))
        else:
            ax = axin

        dfl = self.df[self.df["Cell"] == cell]  # just get the cell's data

        # make categorical swarm plot for each configuration
        offs = 0.25
        nconfigs = len(set(dfl["Configuration"]))
        cfg_color = {'all': colors[0], "largestonly": colors[1], "removelargest": colors[2], "removetwolargest": colors_swarm[3],
                    'Intact': colors[0], "Pruned": colors[1]}
        legs = {'Intact': None, 'Pruned': None, "ANF": None, "all": None, "largestonly": None, "removelargest": None, "removetwolargest": None}
        
        if mode == 'V':
            for ifr, fr in enumerate(sorted(set(dfl['frequency']))):
                xc = ifr
                x = np.zeros(nconfigs)
                y = np.zeros(nconfigs)
                ye = np.zeros(nconfigs)
                yc = [None]*nconfigs
                mfc = [None]*nconfigs
                x_an = np.zeros(nconfigs)
                y_an = np.zeros(nconfigs)
                for jcfg, cfg in enumerate(sorted(set(dfl['Configuration']))):
                    x[jcfg] = xc+(jcfg-1)*offs
                    run = dfl[(dfl['frequency'] == fr) & (dfl['Configuration'] == cfg)]
                    y[jcfg] = run['VS_mean']
                    ye[jcfg] = run['VS_SD']
                    yc[jcfg] = cfg
                    mfc[jcfg] = cfg_color[cfg]
                    x_an[jcfg] = x[jcfg]
                    y_an[jcfg] = run['AN_VS']
                # print(cell, fr, x, y, ye, yc, mfc)
                ax.errorbar(x, y, yerr=ye, marker='o', mfc='none', ms=1, mec='none', mew=0, color='grey')
                ax.scatter(x, y, marker='o', c=mfc, s=12)
                ax.plot(x_an, y_an, '-', color="firebrick", lw=1.5, zorder=-100)  # in back of data
        elif mode == "line":
            for icfg, cfg in enumerate(sorted(set(dfl['Configuration']))):
                run = dfl[dfl['Configuration'] == cfg]
                x = np.arange(8) # run['frequency'].values
                y = run['VS_mean'].values
                ye = run['VS_SD'].values
                x_an = np.arange(8)
                y_an = run['AN_VS'].values
                mfc = cfg_color[cfg]
                ax.errorbar(x, y, yerr=ye, marker='o', mfc=mfc, ms=4, mec='none', mew=0, color='grey', label=cfg)
                if icfg == 0:
                    label = "ANF"
                else:
                    label = None
                ax.plot(x_an, y_an, '-', color="firebrick", lw=0.5, zorder=-100, label=label)  # in back of data
        else:
            raise ValueError("Mode must be one of 'V' or 'line' for VS plot")
            # print(x_an, y_an)
        freqs = sorted(set(dfl['frequency'].values))
        if mode == "V":
            freq_list = range(len(freqs))
            xticks_str = [f"{int(f):d}" for f in freqs]
            minor_list = []
        else:
            freq_list = np.arange(8) # [0, 50, 250, 500,  750, 1000]
            xticks_str=['50', '100', '200', '300', '400', '500', '750', '1000']
            minor_list = [] # [100, 200, 400, 500]
        configs = sorted(set(dfl['Configuration']))

        PH.set_axes_ticks(
            ax=ax,
            xticks=freq_list,
            xticks_str=xticks_str,
            xticks_pad=None,
            x_minor = minor_list,
            major_length=3.0,
            minor_length=1.5,
            yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            yticks_str= ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], 
            yticks_pad=None,
            y_minor=[0.1, 0.3, 0.5, 0.7, 0.9],
            y_rotation=0., 
            fontsize=8,
        )
        
        # sns.swarmplot(
        #     x="frequency",
        #     y = "VS_mean", # y="VectorStrength",
        #     hue="Configuration",
        #     data=dfl,
        #     dodge=True,  # offset configurations from each other
        #     size=3.5,
        #     ax=ax,
        # )
        # # find out how many configurations and frequencies are present
        # n_configurations = len(dfl["Configuration"].unique())
        # n_frequencies = len(dfl["frequency"].unique())

        # # get the first set of children plotted - the rest are legend etc
        # locs = ax.get_children()[: n_configurations * n_frequencies]
        # axis_items = []
        # for child in ax.get_children():
        #     if isinstance(child, matplotlib.collections.PathCollection):
        #         axis_items.append(child)
        # x = np.arange(n_frequencies)
        # for i in range(0, len(locs), n_configurations):
        #     xl = []  # hold locations of points that were plotted in this frequency
        #     yl = []
        #     yerr = []
        #     for j in range(n_configurations):
        #         # axis_item = ax.get_children()[i + j]
        #         # if not isinstance(axis_item, matplotlib.collections.PathCollection):
        #         #     continue
        #         offs = axis_items[i + j].get_offsets()
        #         xl.append(offs[0][0])  # build up the line
        #         yl.append(offs[0][1])
        #     ax.plot(xl, yl, color="slategrey", linewidth=1.25, linestyle="-")
        #     ax.scatter(  # plot the AN vector strength bars for reference
        #         x,
        #         dfl["AN_VS"][: n_configurations * n_frequencies : n_configurations],
        #         s=barwidth,
        #         color="firebrick",
        #         marker="_",
        #         linewidth=2,
        #         alpha=0.3,
        #     )
        ax.set_ylim(0, 1)
        xl = ax.get_xticklabels()  # do this before niceplot...
        xl = ax.get_xticklabels()
        newxl = [f"{int(float(x.get_text())):d}" for x in xl if x.get_text() != ""]
        ax.set_xticklabels(newxl)
        ax.set_xlabel("Modulation Frequency (Hz)", fontdict=label_font)
        ax.set_ylabel("Vector Strength", fontdict=label_font)
        PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)

        if legendflag:
            if legend_type == 'Dendrites':
                custom_legend = [Line2D([0], [0], marker="_", markersize=2, color="firebrick", lw=2, label='AN'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[0], markersize=5, label="Intact"),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[1], markersize=5, label="Pruned"),
                    ]
                ax.legend(handles=custom_legend, handlelength=1, loc="lower left", fontsize=7, labelspacing=0.33, markerscale=0.5)
            else:
                ax.legend()
        if isinstance(cell, str):
            cellname = f"BC{int(cell[0]):02d}"
        else:
            cellname = f"BC{cell:02d}"
        if show_cell:
            ax.text(
                x=0.5,
                y=1.0,
                s=cellname,
                fontdict=title_font,
                #    "horizontalalignment":"center","verticalalignment":"top"},
                transform=ax.transAxes,
            )
        if axin is None:
            mpl.show()
        return ax

    def Figure7_Supplemental2(self, mode:str="line"):
        V1 = VS_Plots(dBSPL=15)

        cells = [5, 6, 10, 11, 13, 18]
        P = PH.regular_grid(
            2,
            3,
            figsize=(8, 6),
            panel_labels=["A", "B", "C", "D", "E", "F"],
            margins={"leftmargin": 0.08, "rightmargin": 0.08, "bottommargin": 0.1, "topmargin":0.05},
            verticalspacing=0.15,
            horizontalspacing=0.07,
            labelposition=(-0.12, 1.05),
        )
        axl = P.axarr.ravel()

        for i, cell in enumerate(cells):
            if i == 0:
                legend = True
            else:
                legend = False
            self.plot_VS_summary(cell, axin=axl[i], legendflag=legend, show_cell=True,
                barwidth=180, mode=mode)

        fig = FigInfo()
        fig.P = P

        fig.filename = set_figure_path(
            fignum=7, filedescriptor="VS_Summary_V3", suppnum=2
        )
        title = "SBEM Project Supplemental Figure 7 Modeling : Vector Strength Summary"
        fig.title["title"] = title
        mpl.savefig(
            fig.filename,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )
        fig.filename = fig.filename.with_suffix(".png")
        mpl.savefig(
            fig.filename,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )
        return fig, P

    def Figure8_M(self, ax=None):
        """Plot Figure 8 panel M: Vector strength to SAM tones at
        different modulation frequencies.

        Args:
            ax (matplotlib axis, optional): axis. Defaults to None.
        """
        # V1 = VS_Plots(dBSPL=15, dends="9I9U")

        cells = [9]
        if ax is None:
            P = PH.regular_grid(
                1,
                1,
                figsize=(6, 6),
                panel_labels=["M"],
                margins={"leftmargin": 0.08, "rightmargin": 0.08, "bottommargin": 0.1, "topmargin":0.05},
                verticalspacing=0.15,
                horizontalspacing=0.10,
                labelposition=(-0.12, 1.05),
            )
            axl = P.axarr.ravel()
        else:
            axl = [ax]

        for i, cell in enumerate(cells):
            legend = False
            if i == 0:
                legend = True
            self.plot_VS_summary(cell, axin=axl[i], legendflag=legend, show_cell=True,
                barwidth=180)

        if ax is not None:
            return
        
        fig = FigInfo()
        fig.P = P

        fig.filename = set_figure_path(
            fignum=8, filedescriptor="BC09_Uninnervated_V1"
        )
        title = "SBEM Project Supplemental Figure 7 Modeling : Vector Strength Summary"
        fig.title["title"] = title
        mpl.show()
        mpl.savefig(
            fig.filename,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )
        fig.filename = fig.filename.with_suffix(".png")
        mpl.savefig(
            fig.filename,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )

    def Summarize(self, axin = None, legendflag = ''):
        label_font = {"fontsize": 8, "fontweight": "normal"}
        title_font = {"fontsize": 9, "fontweight": "normal"}

        if axin is None:
            fig, axs = mpl.subplots(3, 1, figsize=(5, 10))
        else:
            axs = axin
        superthreshset = [9, 11, 17, 18]
        for icell, cell in enumerate([2, 5, 6, 9, 10, 11, 13, 17, 18, 30]):
            dfl = self.df[self.df["Cell"] == cell]  # just get the cell's data

            cfg_color = {'all': colors[0], "largestonly": colors[1], "removelargest": colors[2], "removetwolargest": colors_swarm[3],
                        'Intact': colors[0], "Pruned": colors[1]}
            legs = {'Intact': None, 'Pruned': None, "ANF": None, "all": None, "largestonly": None, "removelargest": None, "removetwolargest": None}
            
            if cell in superthreshset:
                marker = 's'
            else:
                marker = 'o'
            mfc = colors[icell]
            for icfg, cfg in enumerate(sorted(set(dfl['Configuration']))):
                if icfg > 2:
                    continue
                ax = axs[icfg]
                # if cfg != "all":
                #     continue
                run = dfl[dfl['Configuration'] == cfg]
                x = np.log10([50, 100, 200, 300, 400, 500, 750, 1000]) # np.arange(8) # run['frequency'].values
                y = np.log10(run['VS_mean'].values)
                ye = run['VS_SD'].values
                x_an = x # np.arange(8)
                y_an = np.log10(run['AN_VS'].values)
                ax.errorbar(x, y, yerr=ye, marker=marker, mfc=mfc, ms=5, mec='none', mew=0, color=mfc, label=f"BC{cell:02d}")
                if icell == 0:
                    label = "ANF"
                else:
                    label = None
                ax.plot(x_an, y_an, '-', color="firebrick", lw=2, alpha=0.5, zorder=-100, label=label)  # in back of data
                ax.set_title(cfg)
        freqs = sorted(set(dfl['frequency'].values))
        freq_list = np.log10([50, 100, 200, 300, 400, 500, 750, 1000]) # np.arange(8) # [0, 50, 250, 500,  750, 1000]
        xticks_str=['50', '100', '200', '300', '400', '500', '750', '1000']
        minor_list = [] # [100, 200, 400, 500]
        for ax in axs:

            PH.set_axes_ticks(
                ax=ax,
                xticks=freq_list,
                xticks_str=xticks_str,
                xticks_pad=None,
                x_minor = minor_list,
                major_length=3.0,
                minor_length=1.5,
                yticks = np.log10(np.arange(0.1, 1.01, 0.1)), # [0, 0.2, 0.4, 0.6, 0.8, 1.0],
                yticks_str= [f"{x:.1f}" for x in np.arange(0.1, 1.01, 0.1)], # ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], 
                yticks_pad=None,
                y_minor=[0.1, 0.3, 0.5, 0.7, 0.9],
                y_rotation=0., 
                fontsize=8,
            )

            ax.set_ylim(-1, 0)
            xl = ax.get_xticklabels()  # do this before niceplot...
            xl = ax.get_xticklabels()
            newxl = [f"{int(float(x.get_text())):d}" for x in xl if x.get_text() != ""]
            ax.set_xticklabels(newxl)
            ax.set_xlabel("Modulation Frequency (Hz)", fontdict=label_font)
            ax.set_ylabel("Vector Strength", fontdict=label_font)
            PH.nice_plot(ax, position=-0.03, direction="outward", ticklength=3)

            ax.legend()
        fig.tight_layout()
        if axin is None:
            mpl.show()
        return ax        

if __name__ == "__main__":

    # V1 = VS_Plots(dBSPL=15)
    # V1.plot_VS_Data()
    # V1.plot_VS_summary(17)
    # exit()

    # V = VS_Plots()
    # fig, P = V.make_figure()

    # V1 = VS_Plots(sels=[9], dBSPL=15, dends="9I9U")
    # V1.Figure8_M()

    V1 = VS_Plots(dBSPL=15)
    V1.Summarize()

    # fm, new_ax = mpl.subplots(2,2)
    # new_ax = np.ravel(new_ax)

    # for i, ax in enumerate(P.axs):
    #     ax.tick_params(left=True, bottom=True)
    #     # print(ax.get_title())
    #     # print(self.cell_list[i], "\n")
    #     # print([a for a in ax.artists])
    #     # print([t for t in ax.texts])
    #     # print(dir(ax.texts))
    #     cn = ax.texts[0].get_text()
    #     # gg_ax = ax.figure
    #     pos = new_ax[i].get_position()
    #     buf = io.BytesIO()
    #     pickle.dump(fig, buf)
    #     buf.seek(0)
    #     # fig.delaxes(ax)

    #     fig2 = pickle.load(buf)
    #     fm.axes[i] = fig2.axes[i]
    #     # ax.remove()
    #     # ax.figure = fm

    #     # print(fm.axes)
    #     # fm.axes[j] = ax
    #     # fm.axes[j] = ax
    #     fm.axes[i].set_position(pos)

    # ax.text(0.1, 0.9, s=f"{self.cell_list[i]:s}", transform=ax.transAxes)
