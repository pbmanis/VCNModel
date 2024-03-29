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

import numpy as np
import pyqtgraph as pg
import seaborn as sns
from matplotlib import pyplot as mpl
from matplotlib import rc
from pylibrary.plotting import plothelpers as PH
from vcnmodel.util.set_figure_path import set_figure_path
from vcnmodel.util.get_data_paths import get_data_paths
from vcnmodel.plotters import \
    figure_data as FD  # table of simulation runs used for plotting figures

rc("text", usetex=True)
rc("text.latex", preamble=r"\usepackage{xcolor}")
rc("mathtext", fontset="stixsans")

sns.set_style(rc={"pdf.fonttype": 42})
mpl.style.use("~/.matplotlib/figures.mplstyle")

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

class PlotZ:
    def __init__(self, pg=False):
        self.pg = pg
        self.config = get_data_paths()
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
            if i == 0:
                annotate=True
            else:
                annotate=False
            self.plot_col(col=i, f=f, P=P, axonly_scale=axonly_scale, annotate=annotate)

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

            P.axarr[1, 2].plot(d1['f'], rho_axon, color=sns_colors[i], marker=self.syms[i], markersize=1.5, label=f"{FD.BC_name:s}{fin:02d}")
        P.axarr[1, 2].set_xlabel("Frequency (Hz)")
        rholabel = r"$\rho_{axon}$"
        P.axarr[1, 2].set_ylabel(rholabel)
        P.axarr[1, 2].set_ylim((0, rho_scale))

        save_file = set_figure_path(fignum=3, filedescriptor="Zin_decorated", suppnum=3)
#        "Figure3/Figure3_supp/Figure3_Supplemental3_Zin_undecorated.pdf"
        # P.figure_handle.text(
        #     0.99,
        #     0.99,
        #     save_file.name,
        #     transform=P.figure_handle.transFigure,
        #     horizontalalignment="right",
        #     verticalalignment="top",
        # )
        savedir = Path(save_file).parent
        savedir.mkdir(parents=True, exist_ok=True)
        mpl.savefig(
            save_file,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure3, Supplemental Figure 3 Modeling, Zin",
            },
        )
        f2 = save_file.with_suffix(".png")
        mpl.savefig(
            f2,
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure3, Supplemental Figure 3 Modeling, Zin",
            },
        )
        mpl.show()

    def plot_col(self, col, f, P, axonly_scale=3000., annotate=True):
        import matplotlib.scale as MPLS
        import matplotlib.ticker
        
        ax = P.axarr

        for i, filename in enumerate(f):
            # print("col: ", col, i, "file: ", filename)
            label = f"{FD.BC_name:s}{int(filename[5:7]):d}"
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
                d["f"], d["zin"], color=sns_colors[i], marker=self.syms[i], markersize=1., label=label[:7]
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
                    color=sns_colors[i], 
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
                zr, -zi, color=sns_colors[i], marker=self.syms[i], markersize=1,
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
 #        MPLS.LogScale(ax[1, col], subs=[2,3,4,5,6,7,8,9])s
        for row in [0, 1, 2]:
            # for tick in ax[row, col].get_xticklabels():
#                 tick.set_fontname("Arial")
#             for tick in ax[row, col].get_yticklabels():
#                 tick.set_fontname("Arial")
            ax[row, col].tick_params('both', which='major', labelsize=8)

        # Annotations on panel H
        if annotate:
            axh = P.axdict["H"]
            axh.set_clip_on(False)
            # axh.text(x=20, y=0, s="10 kHz", 
            #     fontdict={"fontsize": 8, "fontweight": "normal",
            #         "ha": "center", "va": "top"})
            # axh.text(x=38, y=0, s="1 Hz", 
            #     fontdict={"fontsize": 8, "fontweight": "normal",
            #         "ha": "left", "va": "bottom"})
            # axh.text(x=100, y=30, s="95 Hz", 
            #     fontdict={"fontsize": 8, "fontweight": "normal",
            #         "ha": "center", "va": "center"})
            # arc around the plot to indicate directions
            npts = 100
            ang = np.linspace(-np.pi/8, 0.8*np.pi, npts)
            x_center = np.linspace(60, 45, npts)
            y_center = np.linspace(18, 18, npts)
            minrad = 38
            import scipy.signal.windows
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Rectangle
            radius = minrad + 12*(scipy.signal.windows.gaussian(200, 80)[85:185]-0.5)
            # np.linspace(40, 40, 100)
            x = x_center + radius*np.cos(ang)
            y = y_center + radius*np.sin(ang)
            axh.plot(x[1:], y[1:], 'k-', clip_on=False)
            # axh.plot(x_center, y_center, 'r-')

            axh.arrow(x[0], y[0], x[1]-x[0], y[1]-y[0], shape="full",
                overhang=0.3,
                lw = 1, width=0, head_length=6, head_width=5, color='k', clip_on=False)
            axh.arrow(x[-1], y[-1], x[-1]-x[-2], y[-1]-y[-2], shape="full",
                overhang=0.3,
                lw = 1, width=0, head_length=6, head_width=5, color='k', clip_on=False,)
            axh.text(x=20, y=0, s="10 kHz", 
                fontdict={"fontsize": 8, "fontweight": "normal",
                    "ha": "center", "va": "top"})
            axh.text(x=36, y=-0.5, s="1 Hz", 
                fontdict={"fontsize": 8, "fontweight": "normal",
                    "ha": "left", "va": "bottom"},
                    )
            # rect = [Rectangle((90, 27), 30, 8, facecolor='w', fill=True, edgecolor='k', clip_on=False)]
            # pc = PatchCollection(rect)
            # axh.add_collection(pc)
            axh.text(x=105, y=28, s="95 Hz", 
                fontdict={"fontsize": 8, "fontweight": "normal",
                    "ha": "center", "va": "center"},
                    bbox=dict(facecolor='white', alpha=1, edgecolor='none'))

if __name__ == "__main__":
    Z = PlotZ()
