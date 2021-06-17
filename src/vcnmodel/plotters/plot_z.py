import pickle
from pathlib import Path

import matplotlib
import numpy as np
import pyqtgraph as pg
import seaborn as sns
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH

import toml

matplotlib.rcParams["mathtext.fontset"] = "stixsans"
# matplotlib.rcParams['font.family'] = 'sans-serif'
# matplotlib.rcParams['font.sans-serif'] = ['stixsans'] #, 'Tahoma', 'DejaVu Sans',
#                                #'Lucida Grande', 'Verdana']
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["text.usetex"] = False
sns.set_style(rc={"pdf.fonttype": 42})

config = toml.load(open("wheres_my_data.toml", "r"))

"""
Plot the impedances as seen from the soma for VCN cells
for different dendritic configurations
Assumes that the analysis has already been done, and that the
location of that data is given in the 'wheres_my_data.toml' file
"""


class PlotZ:
    def __init__(self, pg=False):
        self.pg = pg

        # f = ['VCN_c09_Full_Z.pkl', 'VCN_c09_NoUninnervated_Z.pkl', 'VCN_c09_NoDend_Z.pkl']
        self.fi = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
        # fi = [9, 11, 13, 30]
        # fi = [2, 5, 6, 10, 17]
        self.filenames = []
        for fin in self.fi:
            self.filenames.append(f"VCN_c{fin:02d}_Full_normal_Z.pkl")

        # cols = ['w', 'm', 'c', 'y', 'g', 'r', 'b', pg.mkPen()]
        self.syms = ["s", "o", "x", "s", "o", "x", "s", "o", "x", "s"]

        if self.pg:
            self.pg_plot()
        else:
            self.mpl_plot()

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
                    config["cellDataDirectory"], config["impedanceDirectory"], filename
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
            3,
            2,
            order="rowsfirst",
            figsize=(6, 7),
            verticalspacing=0.08,
            horizontalspacing=0.12,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.08,
                "rightmargin": 0.08,
                "topmargin": 0.12,
            },
            panel_labels=["A", "B", "C", "D", "E", "F"],
            labelposition=(-0.05, 1.05),
        )

        # style.style_apply()
        for i, cond in enumerate(["Full_normal", "NoDend"]):
            f = []
            for fin in self.fi:
                f.append(f"VCN_c{fin:02d}_{cond:s}_Z.pkl")
            self.plot_col(col=i, f=f, P=P)
        save_file = "Fig_M1B_Supplemental_Zin.pdf"
        P.figure_handle.text(
            0.99,
            0.99,
            save_file,
            transform=P.figure_handle.transFigure,
            horizontalalignment="right",
            verticalalignment="top",
        )
        mpl.savefig(
            Path(config["baseDataDirectory"], "Figures", save_file),
            metadata={
                "Creator": "Paul Manis",
                "Author": "Paul Manis",
                "Title": "SBEM Project Figure1B Modeling (Supplemental)",
            },
        )
        mpl.show()

    def plot_col(self, col, f, P):
        ax = P.axarr
        for i, a in enumerate(ax.ravel()):
            if i < 4:
                a.set_xscale("log")
        for i, filename in enumerate(f):
            # print("col: ", col, i, "file: ", filename)
            label = filename[:]  # .replace('_', '\_')
            with open(
                Path(
                    config["cellDataDirectory"], config["impedanceDirectory"], filename
                ),
                "rb",
            ) as fh:
                d = pickle.load(fh)
            print(f"File read: {str(filename):s}")
            # print("dkeys: ", d.keys())  # col = pg.intColor(i, hues=len(f))
            ax[0, col].plot(
                d["f"], d["zin"], marker=self.syms[i], markersize=3, label=label[:7]
            )
            ax[0, col].set_ylim(0, 100.0)
            ax[0, col].set_ylabel(r"R (M$\Omega$)")
            ax[0, col].set_xlabel("Frequency (Hz)")

            if col == 0:
                ax[0, col].set_title("Full Dendrite", {"y": 1.15, "fontweight": "bold"})
            if col == 1:
                ax[0, col].set_title("No Dendrite", {"y": 1.15, "fontweight": "bold"})

            ax[1, col].plot(
                d["f"],
                d["phase"],
                marker=self.syms[i],
                markersize=3,
                # label=filename
            )
            ax[1, col].set_ylim(-1.5, 0.25)
            ax[1, col].set_ylabel(r"$\phi$ (radians)")
            ax[1, col].set_xlabel("Frequency (Hz)")

            zr = d["zin"] * np.sqrt(1.0 / (1 + (np.tan(d["phase"]) ** 2.0)))
            zi = np.tan(d["phase"]) * zr

            ax[2, col].plot(
                zr, -zi, marker=self.syms[i], markersize=3,
            )
            ax[2, col].set_ylim(-10.0, 60.0)
            ax[2, col].set_xlim(5, 100)
            ax[2, col].set_ylabel(r"-Im(Z) (M$\Omega$)")
            ax[2, col].set_xlabel(r"Re(Z) (M$\Omega$)")

        if col == 0:
            axbox = ax[0, col].get_position()
            ax[0, col].legend(
                fontsize=7,
                frameon=False,
                fancybox=False,
                bbox_to_anchor=[
                    axbox.x0 + 0.7 * axbox.width,
                    axbox.y0 + 0.27 * axbox.height,
                ],
                bbox_transform=P.figure_handle.transFigure,
            )


if __name__ == "__main__":
    Z = PlotZ()
