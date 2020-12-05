import json
from typing import Union
from collections import OrderedDict
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy.stats
import toml
from matplotlib import pyplot as mpl
from pylibrary.plotting import plothelpers as PH

"""
Data sets are from George Spirou & Co. (Michael Morehead, 
    Nathan Spencer, Matthew Kersting)

This version **imports** data from an excel spreadsheet in VCN_Cells.
The location of the spreadsheet is defined in wheres_my_data.toml

When spreadsheets are updated with new data, we can just read them to build a new table.
last modified: 11 Oct 2019  pbm
"""

"""cell_config generates the configuration of synaptic inputs and devines a cell

The synaptic input structure consists of a list of tuples.
Each element in the list corresponds to one terminal and all of it's active zones
each tuple consists of N sites (calculated from area * average synapses/um2)
delay, and SR
the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
the 6th entry is E or I (guess)
the 7th entry is the location to insert the synapse, defaulting to {'soma': [0, 0.5, 1.0]}
and the diameter. The brances are not included. distances are in microns

The make_dict routine converts these into the lists needed by model_run to population the cell
connections. Note that we also use this to define the cell type, which determines the ion channels
that population the cell.

Measurements:
distances are in microns
size is measured as radii (NOT diameter)
   [(ASA), nsyn(calculated), delay, SRgroup, delay2, axonlength, branch l
        ength, syntype, postlocationdescriptor]

syntype is "AN" for auditory nerve, "AMPA", "AMPA+NMDA", "glycine", "GABA",
    or can be more specific as to actual mechanism

location descriptor is as follows:
{{'soma': [0, 0.5, 1.0]}: [0, 0.5, 1.0]}
The synapse might be split across sections as follows:
{'nameofsection': {nrnsection1#: [location, fractionofgmax], nrnsection2#: [location, fractionofgmax]}}
if fraction of gmax is -1, then it is computed as the residual of the remaining gmax.
(this allows you to split up a big ending that crosses anatomically defined boundaries)
"""


matplotlib.use("Qt5Agg")
rcParams = matplotlib.rcParams
rcParams["text.usetex"] = False
rcParams["svg.fonttype"] = "none"  # No text as paths. Assume font installed.
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

# datafile_default = Path('MorphologyData', 'Dendrite Quality and Surface Areas_comparisons_pbm_15Mar2019_v2.xlsx')
# soma_area_data = 'Mesh Surface Area'
config = toml.load(open("wheres_my_data.toml", "r"))
dendqual = Path(config["baseMorphologyDirectory"], config["dendriteQualityFile"])

# specify the name of the columns in the excel sheet
# these keep getting edited and changing

dataset = "new-3-2020"  # or 'old'
if dataset == "old":
    hocsomarecons = "HOC10"
    hoddendrecons = "HOC Dendrite Area"
    soma_area_data = "Mesh Soma Area Smoothed"  # old
    dend_area_data = "Mesh Dendrite Area"
elif dataset == "new-3-2020":  # with split spreadsheet
    hocsomarecons = (
        "Hoc_Soma_Full_4-03-2020"  # "Hoc_Soma_3-15-2020" #"New 10-section recons"
    )
    hocdendrecons = (
        "Hoc_Dend_Full_4-03-2020"  # "Hoc_Dend_3-15-2020" # "Original reconstructions"
    )
    soma_area_data = "Somatic Surface Area"
    dend_area_data = "Dendritic Surface Area"
else:
    raise ValueError("cell_config: Data set must be either 'new-3-2020', 'old'")

cellsintable = [
    2,
    5,
    6,
    8,
    9,
    10,
    11,
    13,
    14,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    24,
    27,
    29,
    30,
    0,  # test cell with soma only
]
datafile = dendqual
inputs = [f"Input {i+1:d}" for i in range(12)]  # input column labels

# synperum2 = 0.65 # average density of synapses, synapses per micron squared
# Original value, from Spriou measurements in MNTB.
# value confirmed in 5 VCN gBCs (0.68, leaving largest out)
# 0.729 with all 5.

# synperum2 = 0.799  # new numbers from synapse counting: value is average of 15
# terminals from 9 cells in the VCN (Matthew Kersting)
# std is 0.109 4/2020

synperum2 = 0.7686  # new numbers form synapse counting: 23 terminals from 10 cells.
# cells: 2, 5, 6, 9, 10, 11, 13, 17, 18, 30
# SD is 0.1052


class CellConfig:
    def __init__(self, datafile=None, verbose=False):

        self.synperum2 = synperum2
        if datafile is None:
            datafile = dendqual
        self.verbose = verbose
        self.datafile = datafile
        if self.verbose:
            print("CellConfig: configuration dict: ")
        with open(datafile, "rb") as fh:
            self.SDSummary = pd.read_excel(
                fh, config["SomaAndDendriteData"], skiprows=3
            )

        with open(datafile, "rb") as fh:
            self.ASA = pd.read_excel(
                fh, config["asaData"], skiprows=config["asaHeaderSkip"]
            )

        self.VCN_Inputs = OrderedDict()

        for cellnum in cellsintable:
            self.build_cell(cellnum)
            r, ct = self.make_dict(f"VCN_c{cellnum:02d}", synperum2)
        if self.verbose:
            self.print_cell_inputs(self.VCN_Inputs)

    def build_cell(self, cellnum):
        dcell = self.ASA[self.ASA["Cell-Inputs"] == cellnum]
        celln = f"VCN_c{cellnum:02d}"
        self.VCN_Inputs[celln] = ["bushy", []]
        for i in inputs:
            inasa = dcell[i].values
            if pd.isnull(inasa):  # skip empty entries
                continue
            if len(inasa) == 0:
                continue
            self.VCN_Inputs[celln][1].append(
                [
                    (float(inasa)),
                    0.0,
                    2.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    "AN",
                    {"soma": [0, 0.5, 1.0]},
                ]
            )
        # print('dcell: ', dcell)

    def make_dict(self, cell, synperum2:Union[float, None]= None, velocity:float=3.0, areainflate:float=1.0):
        assert cell in self.VCN_Inputs
        indata = self.VCN_Inputs[cell][1]
        celltype = self.VCN_Inputs[cell][0]
        r = [None] * len(indata)
        if synperum2 is None:
            synperum2 = self.synperum2  # use value at top of table.
        for j, i in enumerate(indata):
            asa = (
                i[0] * areainflate
            )  # here is where we can control area by inflations from reference value
            r[j] = OrderedDict(
                [
                    ("input", j + 1),
                    ("asa", asa),
                    ("synperum2", synperum2),
                    ("nSyn", int(np.around(asa * synperum2))),
                    ("delay", i[1]),
                    ("SR", i[2]),
                    ("delay2", np.nansum([i[3], i[5]]) * 0.001 / velocity),
                    ("axonLen", i[3]),
                    ("axonR", i[4]),
                    ("branchLen", i[5]),
                    ("branchR", i[6]),
                    ("type", i[7]),
                    ("postlocations", i[8]),
                ]
            )
        return r, celltype

    def print_cell_inputs_json(self, r):
        if self.verbose:
            print("\nCell Inputs:")
            print(json.dumps(r, indent=4))

    def print_cell_inputs(self, r):
        chs = "Cell ID, ASA, nsyn(calculated), delay, SRgroup, delay2, axonlength, branch length, syntype, postlocation"
        cht = chs.split(", ")
        slen = max([len(c) for c in cht]) + 2
        ch = ""
        for i in range(len(cht)):
            sfmt = "{0:>{1}s}".format(cht[i], slen)
            ch += sfmt

        print(f"{ch:s}")  # header
        for v in list(r.keys()):
            t = f"{v:<{slen}s}"
            # print(r[v])
            for inp in r[v][1]:
                # print('inp: ', inp)
                if isinstance(inp, tuple) or isinstance(inp, list):
                    # print('tu, li')
                    t += f"{inp[0]:{slen}.1f}"
                elif isinstance(inp, float):
                    t += f"{inp:{slen}.1f}"
                else:
                    t += f"{str(inp):^{slen}s}"
            print(t)

    def get_soma_ratio(self, cellID):
        if isinstance(cellID, str):
            cellnum = int(cellID[-2:])
        dcell = self.SDSummary[self.SDSummary["Cell Number"] == cellnum]
        # print(dcell.columns)
        # print(hocsomarecons)
        mesh_area = dcell[soma_area_data].values[0]
        hoc_soma_area = dcell[hocsomarecons].values[0]
        inflateratio = mesh_area / hoc_soma_area
        print(
            f"    Cell: {cellnum:02d}: Soma mesh area: {mesh_area:.2f}  Soma hoc area: {hoc_soma_area:.2f}  ",
            end="",
        )
        print(f"          Soma Inflation ratio: {inflateratio:.3f}")
        return inflateratio

    def get_dendrite_ratio(self, cellID):
        if isinstance(cellID, str):
            cellnum = int(cellID[-2:])
        dcell = self.SDSummary[self.SDSummary["Cell Number"] == cellnum]
        mesh_area = dcell[dend_area_data].values[0]
        hoc_dend_area = dcell[hocdendrecons].values[0]
        if hoc_dend_area > 0.0:  # assuming there are any dendrites...
            inflateratio = mesh_area / hoc_dend_area
        else:
            inflateratio = 1.0
        print(
            f"    Cell: {cellnum:02d}: Dendrite mesh area: {mesh_area:.2f}  HOC Dendrite area: {hoc_dend_area:.2f}  ",
            end="",
        )
        print(f"          Dendrite Inflation ratio: {inflateratio:.3f}")
        return inflateratio

    def _get_dict(self, cellnum, area=0.65):
        synperum2 = area
        r, celltype = self.make_dict(f"VCN_c{cellnum:02d}", synperum2=area)
        return r
        
    def summarize_release_sites(self):
        """
        Compare # of release sites under 3 different synaptic density assumptions
        0.65 syn/um2 (early simulations; value taken from MNTB calyx)
        0.799 syn/um2 (April 2020-10 Nov 2020 simulations; based on subset of endings)
        0.7686 syn/um2 (based on 23 endings, data from June 2020)
        """
        cells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]  # just grade A cells
        synsizes = [0.65, 0.799, 0.7686]
        ns = []
        for j, cell in enumerate(cells):
            print(F"Cell: VCN_c{cell:02d}")
            ns = []
            for i, synsize in enumerate(synsizes):
                nsA = self._get_dict(cell, area=synsize)
                if i == 0:
                    nsT = np.zeros((3, len(nsA)))
                for k, nsyn in enumerate(nsA):
                    nsT[i, k] = nsyn['nSyn']
            for i in range(nsT.shape[0]):
                print(f"{synsizes[i]:.4f}  ", end="")
                for x in nsT[i]:
                    print(f"{int(x):5d}", end="")
                print()

                    # for k, synn in enumerate(nsyn):
#                         print(k, synn['nSyn'])
        
    def summarize_inputs(self):
        P = PH.regular_grid(
            2,
            3,
            order="columnsfirst",
            figsize=(6, 4),
            showgrid=False,
            verticalspacing=0.18,
            horizontalspacing=0.08,
            margins={
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.06,
                "bottommargin": 0.1,
            },
            labelposition=(-0.22, 1.0),
            parent_figure=None,
            panel_labels=["A", "B", "C", "D", "E", "F"],
        )
        ax = [P.axdict[x] for x in P.axdict.keys()]
        #   PH.nice_plot(ax)
        allendings = []
        cellendings = {}
        for cell, v in self.VCN_Inputs.items():
            cellendings[cell] = []
            for s in v[1]:
                if isinstance(s, list) or isinstance(s, tuple):
                    allendings.append(s[0])
                    cellendings[cell].append(s[0])
                else:
                    continue
                    # print('not list or tuple: ', cell, s)
                    # allendings.append(s)
                    # cellendings[cell].append(s)
        tsize = 9
        ax[0].set_title("All ending areas", fontsize=tsize)
        ax[0].hist(allendings, bins=20)
        ax[0].set_ylim((0, 25))
        ax[0].set_ylabel("N")
        ax[0].set_xlim((0, 350))
        ax[0].set_xlabel("Area ($um^2$)")

        normd = []
        ratio1 = []
        ratio2 = []
        meansize = []
        maxsize = []
        convergence = []

        for cell in cellendings.keys():
            if len(cellendings[cell]) == 0:  # missing data in table
                continue
            normd.extend(cellendings[cell] / np.max(cellendings[cell]))
            ratio1.append(cellendings[cell][1] / cellendings[cell][0])
            ratio2.append(np.mean(cellendings[cell]) / cellendings[cell][0])
            meansize.append(np.mean(cellendings[cell]))
            maxsize.append(np.max(cellendings[cell]))
            convergence.append(len(cellendings[cell]))
        print("Convergence: ", convergence)

        ax[1].set_title("Normalized by largest", fontsize=tsize)
        ax[1].hist(normd, bins=20, range=(0, 1.0), align="mid")
        ax[1].set_xlabel("Area Ratio")
        ax[1].set_ylabel("N")
        PH.talbotTicks(
            ax[1],
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            pointSize=10,
            tickPlacesAdd={"x": 2, "y": 0},
            floatAdd={"x": 2, "y": 0},
            axrange={"x": None, "y": None},
        )

        ax[2].set_title("Ratio largest to next largest", fontsize=tsize)
        ax[2].hist(ratio1, bins=10, range=(0, 1.0), align="mid")
        ax[2].set_xlabel("Area Ratio")
        PH.talbotTicks(
            ax[2],
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            pointSize=10,
            tickPlacesAdd={"x": 1, "y": 0},
            floatAdd={"x": 1, "y": 0},
            axrange={"x": None, "y": None},
        )

        ax[3].set_title("Ratio of mean to largest", fontsize=tsize)
        ax[3].hist(ratio2, bins=10, range=(0, 1.0), align="mid")
        ax[3].set_xlabel("Area Ratio")
        PH.talbotTicks(
            ax[3],
            axes="xy",
            density=(1.0, 1.0),
            insideMargin=0.05,
            pointSize=10,
            tickPlacesAdd={"x": 2, "y": 0},
            floatAdd={"x": 2, "y": 0},
            axrange={"x": None, "y": None},
        )

        ax[4].set_title("Convergence vs. mean size", fontsize=tsize)
        ax[4].set_xlim((0.0, 200.0))
        ax[4].set_xlabel("Area ($um^2$)")
        ax[4].set_ylim((0.0, 15.0))
        ax[4].set_ylabel("Convergence")
        fit = np.polyfit(meansize, convergence, 1)
        fit_fn = np.poly1d(fit)
        r, p = scipy.stats.pearsonr(meansize, convergence)
        ax[4].text(x=0.05, y=0.95, s=f"r: {r:.3f}, p={p:.3e}", fontsize=8)
        ax[4].scatter(meansize, convergence)
        ax[4].plot(meansize, fit_fn(meansize), "--k")

        ax[5].set_title("Convergence vs max size", fontsize=tsize)
        fit = np.polyfit(maxsize, convergence, 1)
        fit_fn = np.poly1d(fit)
        r, p = scipy.stats.pearsonr(maxsize, convergence)
        ax[5].text(x=0.05, y=0.95, s=f"r: {r:.3f}, p={p:.3f}", fontsize=8)
        ax[5].scatter(maxsize, convergence)
        ax[5].plot(maxsize, fit_fn(maxsize), "--k")
        ax[5].set_xlim((0.0, 350.0))
        ax[5].set_xlabel("Area ($um^2$)")
        ax[5].set_ylim((0.0, 15.0))
        ax[5].set_xlabel("Convergence")

        mpl.show()


if __name__ == "__main__":
    # Check the formatting and display the results
    cc = CellConfig(datafile)
   # cc.summarize_inputs()
    cc.summarize_release_sites()
