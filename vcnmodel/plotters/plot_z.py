"""
Plot the impedances as seen from the soma for VCN cells
for different dendritic configurations
Assumes that the analysis has already been done, and that the
location of that data is given in the 'wheres_my_data.toml' file
Generates: Figure3_Supplemental3_Zin.pdf

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""
import pickle
from pathlib import Path

import matplotlib
import numpy as np
import pyqtgraph as pg
import seaborn as sns
import toml
from matplotlib import pyplot as mpl
from matplotlib import rc
from pylibrary.plotting import plothelpers as PH

rc("text", usetex=True)
rc("text.latex", preamble=r"\usepackage{xcolor}")
rc("mathtext", fontset="stixsans")

sns.set_style(rc={"pdf.fonttype": 42})
mpl.style.use("~/.matplotlib/figures.mplstyle")


class PlotZ:
    def __init__(self, pg=False):
        self.pg = pg
        self.config = toml.load(open("wheres_my_data.toml", "r"))
        # f = ['VCN_c09_Full_Z.pkl', 'VCN_c09_NoUninnervated_Z.pkl', 'VCN_c09_NoDend_Z.pkl']
        self.fi = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
        # fi = [9, 11, 13, 30]
        # fi = [2, 5, 6, 10, 17]
        self.filenames = []
        for fin in self.fi:
            fin_r = self.checkfi(fin)
            self.filenames.append(fin_r)
        # cols = ['w', 'm', 'c', 'y', 'g', 'r', 'b', pg.mkPen()]
        self.syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x", "s"]

        if self.pg:
            self.pg_plot()
        else:
            self.mpl_plot()

    def checkfi(self, fin, cellmode:str="Full", decorate:str="normal"):
        """
        Map and generate filename - fixing cells 2 and 5
        cellmode is either "Full", "NoDend", or "AxonOnly"
        decorate is "normal, actdend, or pasdend
        """
        if fin not in [2, 5]:
            fname = f"VCN_c{fin:02d}_{cellmode:s}_Meshinflate_{decorate:s}_Z.pkl"
        else:
            fname = f"VCN_c{fin:02d}_{cellmode:s}_Meshinflate_standardized_axon_{decorate:s}_Z.pkl"
        return fname
        
    def pg_plot(self):
        """
        Pyqtgraph plotting version
        """
        pg.mkQApp()
        win = pg.GraphicsWindow()
        win.setGeometry(100, 100, 500, 1200)
        p1 = win.addPlot()
        win.nextRow()
        p2 = win.addPlot()
        win.nextRow()
        p3 = win.addPlot()
        p1.setLogMode(x=True, y=False)
        p2.setLogMode(x=True, y=False)
        p1.setLabels(title="Zin", left="Zin (Mohm)", bottom="Frequency(Hz)")
        p2.setLabels(title="Phase", left="Phase (pi radians)", bottom="Frequency(Hz)")
        p3.setLabels(title="Impedance Locus", left="-Zi (MOhm)", bottom="Zr (MOhm)")
        # legend = pg.LegendItem()
        # legend.setParentItem(p1)
        # legend = p2.addLegend((80, 50), offset=(-1, 1))
        # print(dir(legend))
        # legend.setProperty({textsize:'8pt'})

        for i, filename in enumerate(self.filenames):
            with open(
                Path(
                    self.config["cellDataDirectory"], self.config["impedanceDirectory"], filename
                ),
                "rb",
            ) as fh:
                d = pickle.load(fh)
            col = pg.intColor(i, hues=len(self.filenames))
            p1.plot(
                d["f"], d["zin"], pen=col, symbol=self.syms[i], symbolSize=3, name=filename
            )
            p2.plot(
                d["f"], d["phase"], pen=col, symbol=self.syms[i], symbolSize=3, name=filename
            )
            zr = d["zin"] * np.sqrt(1.0 / (1 + (np.tan(d["phase"]) ** 2.0)))
            zi = np.tan(d["phase"]) * zr

            p3.plot(zr, -zi, pen=col, symbol=self.syms[i], symbolSize=3)
        # legendLabelStyle = {'color': '#FFF', 'size': '12pt', 'bold': True, 'italic': False}
        # for item in legend.items:
        #     for single_item in item:
        #         if isinstance(single_item, pg.graphicsItems.LabelItem.LabelItem):
        #             single_item.setText(single_item.text, **legendLabelStyle)

        pg.QtGui.QApplication.instance().exec_()

    def mpl_plot(self):
        """
        matplotlib version (for production)
        """
        # style = STY.styler('JNeurosci', "single")
        P = PH.regular_grid(
            rows=3,
            cols=3,
            order="rowsfirst",
            figsize=(8, 8),
            verticalspacing=0.1,
            horizontalspacing=0.1,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.08,
                "rightmargin": 0.08,
                "topmargin": 0.12,
            },
            panel_labels=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"],
            labelposition=(-0.08, 1.1),
            fontsize={"tick": 8, "label": 10, "panel": 13},
            fontweight={"tick": "normal", "label": 'normal', "panel": "bold"},
        )
        axonly_scale = 5500.
        rho_scale=150.
        # style.style_apply()
        for i, cond in enumerate(["Full", "NoDend", "AxonOnly"]):
            f = []
            for fin in self.fi:
                f.append(self.checkfi(fin, cellmode=cond, decorate="normal"))
            self.plot_col(col=i, f=f, P=P, axonly_scale=axonly_scale)

        for i, fin in enumerate(self.fi):
            with open(
                Path(
                    self.config["cellDataDirectory"], self.config["impedanceDirectory"], self.checkfi(fin, cellmode="Full", decorate="normal")
                ),
                "rb",
            ) as fh:
                d1 = pickle.load(fh)
            with open(
                Path(
                    self.config["cellDataDirectory"], self.config["impedanceDirectory"], self.checkfi(fin, cellmode="AxonOnly", decorate="normal")
                ),
                "rb",
            ) as fh:
                d2 = pickle.load(fh)
            g_ds = (1./d1['zin'])  # g_dend _ g soma / g-axon
            g_axon = (1./d2['zin'])
            rho_axon = g_ds/g_axon

            P.axarr[1, 2].plot(d1['f'], rho_axon, marker=self.syms[i], markersize=1.5, label=f"BC{fin:02d}")
        P.axarr[1, 2].set_xlabel("Frequency (Hz)")
        rholabel = r"$\rho_{axon}$"
        P.axarr[1, 2].set_ylabel(rholabel)
        P.axarr[1, 2].set_ylim((0, rho_scale))

        save_file = "Figure3/Figure3_supp/Figure3_Supplemental3_Zin_undecorated.pdf"
        P.figure_handle.text(
            0.99,
            0.99,
            save_file,
            transform=P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )
        mpl.savefig(
            Path(self.config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure3, Supplemental Figure 3 Modeling, Zin",
            },
        )
        mpl.show()

    def plot_col(self, col, f, P, axonly_scale=3000.):
        import matplotlib.scale as MPLS
        import matplotlib.ticker
        
        ax = P.axarr

        for i, filename in enumerate(f):
            # print("col: ", col, i, "file: ", filename)
            label = f"BC{int(filename[5:7]):d}"
            with open(
                Path(
                    self.config["cellDataDirectory"], self.config["impedanceDirectory"], filename
                ),
                "rb",
            ) as fh:
                d = pickle.load(fh)
            print(f"File read: {str(filename):s}")
            # print("dkeys: ", d.keys())  # col = pg.intColor(i, hues=len(f))
            ax[0, col].plot(
                d["f"], d["zin"], marker=self.syms[i], markersize=1., label=label[:7]
            )

            ax[0, col].set_ylim(0, 100.0)
            ax[0, col].set_ylabel(r"R (M$\Omega$)")
            ax[0, col].set_xlabel("Frequency (Hz)")
            # ax[0, col].set_xscale("log")
            if col == 0:
                ax[0, col].set_title("Full Dendrite", {"y": 1.2, "fontweight": "bold"})
            if col == 1:
                ax[0, col].set_title("No Dendrite", {"y": 1.2, "fontweight": "bold"})
            if col == 2:
                ax[0, col].set_title("Axon Only", {"y": 1.2, "fontweight": "bold"})
                ax[0, col].set_ylim(0, axonly_scale)
            if col != 2:  # we plot rho instead here.
                ax[1, col].plot(
                    d["f"],
                    d["phase"],
                    marker=self.syms[i],
                    markersize=1.0,
                    # label=filename
                )
                ax[1, col].set_ylim(-1.1*np.pi/2., 1.1*np.pi/2)
                ax[1, col].set_ylabel(r"$\phi$ (radians)")
                ax[1, col].set_xlabel("Frequency (Hz)")

            zr = d["zin"] * np.sqrt(1.0 / (1 + (np.tan(d["phase"]) ** 2.0)))
            zi = np.tan(d["phase"]) * zr

            ax[2, col].plot(
                zr, -zi, marker=self.syms[i], markersize=1,
            )
            if col < 2:
                ax[2, col].set_ylim(-10.0, 60.0)
                ax[2, col].set_xlim(5, 100)
            else:
                ax[2, col].set_ylim(-1000.0, axonly_scale)
                ax[2, col].set_xlim(-1000.0, axonly_scale)
            
            ax[2, col].set_ylabel(r"-Im(Z) (M$\Omega$)")
            ax[2, col].set_xlabel(r"Re(Z) (M$\Omega$)")

        if col == 0:
            axbox = ax[0, col].get_position()
            ax[0, col].legend(
                fontsize=7,
                frameon=False,
                fancybox=False,
                bbox_to_anchor=[
                    axbox.x0 + 0.8 * axbox.width,
                    axbox.y0 + 0.27 * axbox.height,
                ],
                bbox_transform=P.figure_handle.transFigure,
            )

        for row in [0,1]:
            ax[row, col].tick_params('x', which='minor', direction='in', length=2.5, width=0.5)
            
            ax[row, col].set_xscale("log")
            ax[row, col].set_xlim(1, 10000.)
            xmajor = matplotlib.ticker.LogLocator(base=10, numticks=5)
            ticks = matplotlib.ticker.FixedLocator([1, 10, 100, 1000, 10000])
            ax[row, col].xaxis.set_major_locator(ticks)
            xminor = matplotlib.ticker.LogLocator(base=10, subs = np.arange(0, 10) * 0.1, numticks = 10)
            ax[row, col].xaxis.set_minor_locator(xminor)
            ticks = matplotlib.ticker.FixedLocator([1, 10, 100, 1000, 10000])
            ax[row, col].set_xticklabels(["1", "10", "100", "1000", "10000"])
        # MPLS.LogScale(ax[0, col], subs=[2,3,4,5,6,7,8,9])
 #        MPLS.LogScale(ax[1, col], subs=[2,3,4,5,6,7,8,9])
        for row in [0, 1, 2]:
            # for tick in ax[row, col].get_xticklabels():
#                 tick.set_fontname("Arial")
#             for tick in ax[row, col].get_yticklabels():
#                 tick.set_fontname("Arial")
            ax[row, col].tick_params('both', which='major', labelsize=8)
    
if __name__ == "__main__":
    Z = PlotZ()
