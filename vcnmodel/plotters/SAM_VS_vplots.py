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
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, os.path.abspath("nb"))
import plotnine as PN
import toml
import VS_data_30dB as VS_data
from matplotlib import pyplot as mpl
from matplotlib import ticker
from pylibrary.plotting import plothelpers as PH


def reset_style():
    #     sns.set_style(rc={"pdf.fonttype": 42})
    mpl.style.use("~/.matplotlib/figures.mplstyle")


class VS_Plots:
    def __init__(self):
        importlib.reload(VS_data)
        self.datas = VS_data.data
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
        self.config = toml.load(open("wheres_my_data.toml", "r"))

    def prepare_data(self, datas):
        sio = io.StringIO(datas)
        df = pd.read_table(sio, sep=",")
        return df

    def make_figure(self):
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

        reset_style()

        # plotwhat = "phasesdtime"
        # plotwhat = 'phasesd'
        plotwhat = "VectorStrength"
        # plotwhat = 'SpikeCount'
        # plotwhat = 'logSpikeCount'

        if plotwhat == "phasesdtime":
            df["phasesdtime"] = 1e3 * df["phasesd"] / (2 * np.pi * df["frequency"])
            df["AN_phasesdtime"] = (
                1e3 * df["AN_phasesd"] / (2 * np.pi * df["frequency"])
            )

        if plotwhat == "logSpikeCount":
            df["logSpikeCount"] = np.log10(df["SpikeCount"])

        """
        Plot all the cells in two columns using ggplot
        """

        def getCell(cell):
            return f"BC{int(cell):02d}"

        labels = set(df["frequency"])
        labels = [str(int(l)) for l in sorted(labels)]

        def scalefun(x):
            return [f"{y:.2f}" for y in x]

        print(plotwhat)
        anpt = "AN_phasesd"
        if plotwhat == "phasesdtime":
            yscale = [0.0, 1.5]
            anpt = "AN_phasesdtime"
        elif plotwhat == "SpikeCount":
            yscale = [0, 20000]
        elif plotwhat == "logSpikeCount":
            yscale = [0, 5]
        elif plotwhat == "phasesd":
            yscale = [0, 2 * np.pi]
            anpt = "AN_phasesd"
        else:
            yscale = [-0.1, 1.0]
            anpt = "AN_VS"

        freqs = [50, 100, 200, 300, 400, 500, 750, 1000]

        def fscale(x):
            return [str(int(y)) for y in x]

        fr = fscale(freqs)
        vsscale = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

        def vscale(x):
            return [f"{y:.1f}" for y in x]

        dodge = PN.positions.position_dodge(width=0.6)
        themer = PN.themes.theme
        gg = (
            PN.ggplot(
                df,
                PN.aes(
                    x="factor(frequency)",
                    y=plotwhat,
                    group="factor(frequency)",
                    color="factor(Configuration)",
                ),
            )
            + PN.scale_x_discrete(breaks=freqs, labels=fr)
            + PN.scale_y_continuous(breaks=vsscale, labels=vscale(vsscale))
            # + PN.scales.scale_y_log10()
            + PN.geom_point(
                PN.aes(
                    x="factor(frequency)",
                    y=plotwhat,
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
                    y=plotwhat,
                ),
                position=PN.positions.position_dodge2(width=1),
                alpha=1,
                color="grey",
            )
            + PN.geom_line(
                PN.aes(
                    x="factor(frequency)",
                    y=anpt,
                ),
                position=PN.positions.position_dodge2(width=1),
                alpha=1,
                color="red",
            )
            + PN.scales.scale_colour_brewer(type="qual", palette="Dark2")
            + PN.facet_wrap("Cell", labeller=getCell, nrow=5, ncol=2, dir="v")
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
            + PN.scales.ylim(yscale)
            + PN.themes.theme(legend_position=(0.5, 0.93), legend_entry_spacing=0.001)
            + PN.theme(figure_size=(6, 8))
            # + PN.scales.scale_y_continuous(breaks=[f"{f:.1f}" for f in [0, 0.2, 0.4 0.6, 0.8, 1.0]]) # scalefun) # PN.scales.number_format(accuracy = 0.01))
        )
        # print(dir(PN.scales))

        fig, P = gg.draw(return_ggplot=True)
        PN.options.figure_size = (10, 10)

        for i, ax in enumerate(P.axs):
            ax.tick_params(left=True, bottom=True)

        save_file = f"Figure7/Figure7_supp/Figure7_Supplemental2_V2.pdf"

        # fig.text(
        #     0.99,
        #     0.99,
        #     save_file,  # .replace('_', '\_'),
        #     # transform=fig.transFigure,
        #     horizontalalignment="right",
        #     verticalalignment="top",
        # )
        title = (
            "SBEM Project Supplemental Figure 7 Modeling : Vector Strength Summary",
        )
        mpl.savefig(
            Path(self.config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": title,
            },
        )

        fig.filename = save_file
        if hasattr(fig, "title"):
            fig.title["title"] = title
        return fig


if __name__ == "__main__":
    V = VS_Plots()
    V.make_figure()
