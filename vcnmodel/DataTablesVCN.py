# -*- coding: utf-8 -*-
import functools
import importlib
import sys
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import pyqtgraph as pg
import toml
from pylibrary.plotting import plothelpers as PH
from pyqtgraph.parametertree import Parameter, ParameterTree
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

from vcnmodel import table_manager as table_manager
from vcnmodel.plotters import plot_sims as PS
from pylibrary.tools import cprint as CP

cprint = CP.cprint
"""
Use pyqtgraph tablewidget to build a table showing simulation
files/runs and enabling analysis via a GUI

"""

all_modules = [
    table_manager,
    PS,
]  

cellvalues = [
    2,
    5,
    6,
    9,
    10,
    11,
    13,
    17,
    24,
    29,
    30,
]
modeltypes = [
    "XM13_nacncoop",
    "mGBC",
    "XM13",
    "RM03",
]
runtypes = ["AN", "IV", "IO", "gifnoise"]
experimenttypes = [
    "default",
    "delays",
    "largestonly",
    "removelargest",
    "mean",
    "allmean",
    "twolargest",
    "threelargest",
    "fourlargest",
]
modetypes = ["find", "singles", "IO", "multi"]
analysistypes = ["traces", "PSTH", "revcorr", "SAC", "tuning"]
dendriteChoices = [
    "normal",
    "passive",
    "active",
]

class TableModel(QtGui.QStandardItemModel):
    _sort_order = QtCore.Qt.AscendingOrder

    def sortOrder(self):
        return self._sort_order

    def sort(self, column, order):
        if column == 0:
            self._sort_order = order
            QtGui.QStandardItemModel.sort(self, column, order)

class DataTables:
    def __init__(self):
        self.modeltypes = {
            "b": "Bushy",
            "t": "TStellate",
            "d": "DStellate",
            "v": "Tuberculoventral",
            "p": "Pyramidal",
        }
        self.basepath = None
        self.setPaths()
        self.app = pg.mkQApp()
        self.app.setStyle("fusion")
        # Enable High DPI display with PyQt5
        # self.app.setAttribute(Qt.AA_EnableHighDpiScaling)
        #         if hasattr(QStyleFactory, 'AA_UseHighDpiPixmaps'):
        #             self.app.setAttribute(Qt.AA_UseHighDpiPixmaps)

        self.win = pg.QtGui.QWidget()
        layout = pg.QtGui.QGridLayout()
        self.win.setLayout(layout)
        self.win.setWindowTitle("Model DataTables/FileSelector")
        self.win.resize(1280, 1024)
        self.table = pg.TableWidget(sortable=True)
        style = "::section {background-color: darkblue; }"
        self.selected_index_row = None
        self.table.horizontalHeader().setStyleSheet(style)
        self.model = None
        # self.table.sortingEnabled(True)
        self.voltage = False
        self.runtype = runtypes[0]
        self.cellID = 2
        self.modeltype = modeltypes[0]
        self.modetype = modetypes[0]
        self.experimenttype = experimenttypes[0]
        self.analysistype = analysistypes[0]
        self.dendriteChoices = dendriteChoices[0]
        self.selvals = {
            "ModelType": [modeltypes, self.modeltype],
            "Run Type": [runtypes, self.runtype],
            "Cells": [cellvalues, self.cellID],
            "Mode": [modetypes, self.modetype],
            "Experiment": [experimenttypes, self.experimenttype],
            "Analysis": [analysistypes, self.analysistype],
            "Dendrites": [dendriteChoices, self.dendriteChoices],
        }
        params = [
            # {"name": "Pick Cell", "type": "list", "values": cellvalues, "value": cellvalues[0]},
            {"name": "Scan Runs", "type": "action"},
            {"name": "Update Runs", "type": "action"},
            {
                "name": "Selections",
                "type": "group",
                "children": [
                    {
                        "name": "Run Type",
                        "type": "list",
                        "values": runtypes,
                        "value": runtypes[0],
                    },
                    {
                        "name": "Cells",
                        "type": "list",
                        "values": cellvalues,
                        "value": 2,
                    },
                    {
                        "name": "ModelType",
                        "type": "list",
                        "values": modeltypes,
                        "value": modeltypes[0],
                    },
                    {
                        "name": "Mode",
                        "type": "list",
                        "values": modetypes,
                        "value": modetypes[0],
                    },
                    {
                        "name": "Experiment",
                        "type": "list",
                        "values": experimenttypes,
                        "value": experimenttypes[0],
                    },
                    {
                        "name": "Analysis",
                        "type": "list",
                        "values": analysistypes,
                        "value": analysistypes[0],
                    },
                    {
                        "name": "Dendrites",
                        "type": "list",
                        "values": dendriteChoices,
                        "value": dendriteChoices[0],
                    },
                ],
            },
            {
                "name": "Analysis",
                "type": "group",
                "children": [
                    {
                        "name": "Traces",
                        "type": "action",
                        # "default": False,
                    },
                    {
                        "name": "Singles",
                        "type": "action",
                        # "default": False,
                    },
                    {"name": "Revcorr", "type": "action"},
                    {"name": "PSTH", "type": "action"},
                ],
            },
            {
                "name": "Tools",
                "type": "group",
                "children": [{"name": "Reload", "type": "action"}],
            },
        ]
        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name="Models", type="group", children=params)
        self.ptree.setParameters(self.ptreedata)
        self.ptree.setMaximumWidth(250)
        # build layout for plots and parameters
        layout.addWidget(self.ptree, 0, 0, 8, 1)  # Parameter Tree on left

        # add space for the graphs
        view = pg.GraphicsView()
        layout2 = pg.GraphicsLayout(border=(0, 0, 0))
        view.setCentralItem(layout2)
        layout.addWidget(view, 0, 1, 8, 8)
        layout.addWidget(self.table, 0, 1, 8, 8)  # data plots on right
        self.win.show()
        self.table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)

        self.table.doubleClicked.connect(functools.partial(self.on_Click, self.table))
        self.table.clicked.connect(functools.partial(self.on_Single_Click, self.table))
        self.ptreedata.sigTreeStateChanged.connect(self.command_dispatcher)

        self.table_manager = table_manager.TableManager(
            self.table, self.basepath, self.selvals, self.altColors
        )

    def setPaths(self, stimtype="AN", cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells",)}
        self.basepath = self.datapaths["baseDirectory"]

    def on_Click(self, w):
        index = w.selectionModel().currentIndex()
        self.selected_index_row = index.row()
        self.analyze_from_table(index.row())
        # print("row: ", self.selected_index_row)
     
    def on_Single_Click(self, w):
        index = w.selectionModel().currentIndex()
        self.selected_index_row = index.row()
        self.analyze_from_table(index.row())
        # print("row: ", self.selected_index_row)

    def handleSortIndicatorChanged(self, index, order):
        if index != 0:
            self.table.horizontalHeader().setSortIndicator(
                0, self.table.model().sortOrder())

    def command_dispatcher(self, param, changes):

        for param, change, data in changes:
            path = self.ptreedata.childPath(param)
            # print("Path: ", path[0])
            # if path[0] == "Pick Cell":
            #     self.selvals["Cell"] = data
            #     self.setPaths('AN', cell=data)
            #     self.table_manager.build_table(mode="scan")
            if path[0] == "Scan Runs":
                self.table_manager.build_table(mode="scan")
            if path[0] == "Update Runs":
                self.table_manager.build_table(mode="scan")
            # if path[0] == "Build Index":
            #     self.setPaths(self.runtype, cell=data)
            #     dpath = Path(self.datapaths['baseDirectory'],
            #         f"VCN_{self.cellID:02d}",
            #         "Simulations",
            #         self.runtype)
            #     # print(self.datapaths)
            #     # print(self.runtype)
            #     self.table_manager.find_build_indexfiles(dpath,
            #         force=False)
            # if path[0] == "Update Index":
            #     self.setPaths(self.runtype, cell=data)
            #     self.table_manager.find_build_indexfiles(
            #           Path(self.datapaths['baseDirectory'],
            #               f"VCN_{self.cellID:02d}",
            #               "Simulations",
            #               self.runtype),
            #          force=True)
            if path[0] == "Selections":
                self.selvals[path[1]][1] = str(data)
                self.cellID = self.selvals['Cells'][1]
                # if path[1] == "Run Type":
                #     self.runtype = str(data)
                # elif path[1] == "Cells":
                #     self.cellID = str(data)
                # elif path[1] == "ModelType":
                #     self.modeltype = str(data)
                # elif path[1] == ""
                self.setPaths("AN", cell=data)
                self.table_manager.build_table(mode="scan")
            
            if path[0] == "Analysis":
                if path[1] == "Traces":
                    self.analyze_traces()
                elif path[1] == "PSTH":
                    self.analyze_PSTH()
                elif path[1] == 'Revcorr':
                    self.analyze_revcorr()
                elif path[1] == "Singles":
                    self.analyze_singles()
                elif path[1] == "CMMR Summary":
                    self.analyze_cmmr_summary()
            
            if path[0] == "Tools":
                if path[1] == "Reload":
                    print("reloading...")
                    for module in all_modules:
                        print("reloading: ", module)
                        importlib.reload(module)
                    self.table_manager = table_manager.TableManager(
                        self.table, self.basepath, self.selvals, self.altColors
                    )
                    print("   reload ok")
                    print("-"*80)

                self.table.setSortingEnabled(True)
                self.table.horizontalHeader().sortIndicatorChanged.connect(
                    self.handleSortIndicatorChanged)

    def setColortoRow(self, rowIndex, color):
        for j in range(self.table.columnCount()):
            self.table.item(rowIndex, j).setBackground(color)

    def altColors(self, colors= [QtGui.QColor(0x00, 0x00, 0x00), QtGui.QColor(0x22, 0x22, 0x22)]):
        """
        Paint alternating table rows with different colors

        Parameters
        ----------
        colors : list of 2 elements
            colors[0] is for odd rows (RGB, Hex)
            colors[1] is for even rows
        """
        for j in range(self.table.rowCount()):
            if j % 2:
                self.setColortoRow(j, colors[0])
            else:
                self.setColortoRow(j, colors[1])

    def force_suffix(self, filename, suffix='.pkl'):
        fn = Path(filename)
        if fn.suffix != suffix:
            fn = str(fn)
            fn = fn + suffix
            fn = Path(fn)
        return fn
        
    def analyze_singles(self):
        if self.selected_index_row is None:
            return
        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        nfiles = len(selected.files)

        P = PH.regular_grid(
            nfiles,
            1,
            order="rowsfirst",
            figsize=(6.0, 10.0),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.01,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.03,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )


        PD = PS.PData()
        sfi = sorted(selected.files)

        df = pd.DataFrame(index=np.arange(0, nfiles), columns=['cell', 'syn#', 'nout', "nin", "efficacy"])
        for i in range(nfiles):
            synno, nout, nin = PS.plot_traces(
                P.axarr[i, 0], sfi[i], PD, selected.runProtocol
            )
            eff = float(nout)/nin
            df.iloc[i] = [self.cellID, synno, nout, nin, eff]
        u = df.head(n=nfiles)
        print(df.to_csv(sep='\t'))
        # print(u)

        
        P.figure_handle.show()

    def analyze_traces(self):
        if self.selected_index_row is None:
            return
        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()

        P = PH.regular_grid(
            1,
            1,
            order="rowsfirst",
            figsize=(6.0, 6.0),
            showgrid=False,
            verticalspacing=0.01,
            horizontalspacing=0.01,
            margins={
                "bottommargin": 0.1,
                "leftmargin": 0.07,
                "rightmargin": 0.05,
                "topmargin": 0.03,
            },
            labelposition=(0.0, 0.0),
            parent_figure=None,
            panel_labels=None,
        )

        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()

        sfi = Path(selected.simulation_path, selected.files[0])
        PS.plot_traces(P.axarr[0, 0], sfi, PD, selected.runProtocol)
        P.figure_handle.show()

    def analyze_PSTH(self):
        if self.selected_index_row is None:
            return
        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()

        sizer = OrderedDict(
            [
                ("A", {"pos": [0.08, 0.4, 0.71, 0.22]}),
                ("B", {"pos": [0.55, 0.4, 0.71, 0.22]}),
                ("C", {"pos": [0.08, 0.4, 0.39, 0.22]}),
                ("D", {"pos": [0.55, 0.4, 0.39, 0.22]}),
                ("E", {"pos": [0.08, 0.4, 0.07, 0.22]}),
                ("F", {"pos": [0.55, 0.4, 0.07, 0.22]}),
            ]
        )  # dict elements are [left, width, bottom, height] for the axes in the plot.
        n_panels = len(sizer.keys())
        gr = [
            (a, a + 1, 0, 1) for a in range(0, n_panels)
        ]  # just generate subplots - shape does not matter
        axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.Plotter((n_panels, 1), axmap=axmap, label=True, figsize=(8.0, 6.0))
        P.resize(sizer)  # perform positioning magic
        P.axdict["A"].set_ylabel("mV", fontsize=8)
        P.axdict["D"].set_xlabel("Phase", fontsize=8)
        P.axdict["C"].set_ylabel("Trial", fontsize=8)
        P.axdict["E"].set_ylabel("Trial, ANF", fontsize=8)
        P.axdict["B"].set_title("Stimulus", fontsize=9)
        P.axdict["E"].set_title("ANF Spike Raster", fontsize=9)
        P.axdict["C"].set_title("Bushy Spike Raster", fontsize=9)
        P.axdict["F"].set_title("PSTH", fontsize=9)

        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()
        sfi = Path(selected.simulation_path, selected.files[0])
        PS.plot_AN_response(P, sfi, PD, selected.runProtocol)
        P.figure_handle.show()

    def analyze_revcorr(self):
        if self.selected_index_row is None:
            return
        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()
        rows = 1
        cols = 1  # just a single selection
        sizex = 6.
        sizey = 6.

        plabels = [f"VCN_c{int(self.cellID):02d}"]
        pgbc = plabels[0]
        P = PH.regular_grid(
            rows,
            cols,
            order="rowsfirst",
            figsize=(sizex, sizey),
            panel_labels=plabels,
            labelposition=(0.05, 0.95),
            margins={
                "leftmargin": 0.1,
                "rightmargin": 0.01,
                "topmargin": 0.15,
                "bottommargin": 0.15,
            },
        )

        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]
        PD = PS.PData()
        sfi = Path(selected.simulation_path, selected.files[0])
        res = PS.compute_revcorr(
            P.axarr[0,0], pgbc, sfi, PD, selected.runProtocol)
        
        P.figure_handle.show()

    def analyze_from_table(self, i):
        if self.selected_index_row is None:
            return
        index_row = self.selected_index_row
        selected = self.table_manager.table_data[index_row]

        # map it:
        if selected.runProtocol == "runANSingles":  # subdirectory
            self.analyze_singles() 
            

def main():
    DataTables()
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QApplication.instance().exec_()


if __name__ == "__main__":
    main()
