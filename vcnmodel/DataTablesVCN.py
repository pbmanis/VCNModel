# -*- coding: utf-8 -*-
"""
Use pyqtgraph tablewidget to build a table showing simulation 
files/runs and enabling analysis via a GUI

"""

import sys
import numpy as np
import functools
from pathlib import Path
import pickle
import time
import importlib
import toml

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph import Qt
from pyqtgraph.parametertree import Parameter, ParameterTree

from pylibrary.tools import fileselector
import vcnmodel.table_manager as table_manager

import pylibrary.plotting.plothelpers as PH
import plotters.plot_sims as PS
            
# from . import read_physiology3 as read_physiology
# # thse are for the reload
# from . import show_network
# from . import analysis_reader_tools
# from . import analyze_psth
# from . import analysis_functions


all_modules = [
    table_manager,
    PS
]  # [read_physiology, show_network, analyze_psth, analysis_reader_tools, analysis_functions]

cellvalues = [
    2,
    5,
    6,
    9,
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
        self.table = pg.TableWidget(sortable=False)
        style = "::section {background-color: lightblue; }"
        self.table.horizontalHeader().setStyleSheet(style)
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
            {"name": "Rescan Runs", "type": "action"},
            {"name": "Build index", "type": "action"},
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
                        "name": "Show Sims",
                        "type": "action",
                        # "default": False,
                    },
                    # {
                    #     "name": "CMMR Summary",
                    #     "type": "action"
                    # },
                ],
            },
            {
                "name": "Tools",
                "type": "group",
                "children": [{"name": "Reload", "type": "action",}],
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
        l = pg.GraphicsLayout(border=(0, 0, 0))
        view.setCentralItem(l)
        layout.addWidget(view, 0, 1, 8, 8)
        layout.addWidget(self.table, 0, 1, 8, 8)  # data plots on right
        self.win.show()
        self.table.doubleClicked.connect(functools.partial(self.on_Click, self.table))
        self.ptreedata.sigTreeStateChanged.connect(self.command_dispatcher)

        self.table_manager = table_manager.TableManager(
            self.table, self.basepath, self.selvals
        )

        # print(dir(w.horizontalHeader()))
        #
        # w.horizontalHeader().sectionPressed.connect(functools.partial(on_doubleClick, w))

        # data = np.array([
        #     (1,   1.6,   'x'),
        #     (3,   5.4,   'y'),
        #     (8,   12.5,  'z'),
        #     (443, 1e-12, 'w'),
        #     ], dtype=[('Column 1', int), ('Column 2', float), ('Column 3', object)])
        #
        # self.table.setData(data)

    def setPaths(self, stimtype="AN", cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells",)}
        self.basepath = self.datapaths["baseDirectory"]

    def on_Click(self, w):
        index = w.selectionModel().currentIndex()
        self.analyze_from_table(index.row())
        # print("row: ", index.row())
        # for c in range(w.columnCount()):
        #     v = w.item(index.row(), c).value
        #     print('type: ', type(v))
        #     print('str: ', isinstance(v, str))
        #     print('float: ', isinstance(v, float))
        #     print('int: ', isinstance(v, int))

    def on_doubleClick(self, w):
        print(w.currentItem())

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
            if path[0] == "Rescan Runs":
                self.table_manager.build_table(mode="rescan")

            if path[0] == "Selections":
                self.selvals[path[1]][1] = str(data)
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
                if path[1] == "Show Voltages":
                    self.voltage = data

                print("Voltage Flag: ", self.voltage)
                if path[1] == "CMMR Summary":
                    self.analyze_cmmr_summary()
            if path[0] == "Tools":
                if path[1] == "Reload":
                    print("reloading...")
                    for module in all_modules:
                        print("reloading: ", module)
                        importlib.reload(module)
                    print("   reload ok")

    def setColortoRow(self, rowIndex, color):
        for j in range(self.table.columnCount()):
            self.table.item(rowIndex, j).setBackground(color)

    def altColors(self, colors):
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

    def analyze_from_table(self, i):
        selected = self.table_manager.table_data[i]

        if selected.runProtocol == "runANSingles":
            import pylibrary.plotting.plothelpers as PH
            import plotters.plot_sims as PS
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

            sfi = sorted(selected.files)
            import plotters.plot_sims as PS

            PD = PS.PData()
            changetimestamp = PS.get_changetimestamp()
            for i in range(nfiles):
                PS.plot_traces(
                    P.axarr[i, 0], sfi[i], PD, changetimestamp, selected.runProtocol
                )
            # for i in range(nfiles):
            #     with open(sfi[i], "rb") as fh:
            #         d = pickle.load(fh, encoding="latin1")
            #         for j in range(len(d["Results"]['somaVoltage'])):
            #             P.axarr[i,0].plot(d["Results"]['time'],d["Results"]['somaVoltage'][j], linewidth=0.6)
            #             P.axarr[i,0].set_ylim((-80., 10.))
            P.figure_handle.show()

    #     args = d["command_line"]
    #     protocol = args.protocol
    #     try:
    #         threshold = d["runInfo"]["threshold"]
    #     except:
    #         threshold = -20.0  # (default in read_physiology)
    #     try:
    #         cmmrs2n = d["runInfo"]["cmmrs2n"]
    #     except:
    #         cmmrs2n = None  # not defined?
    #     try:
    #         parallel = not args.noparallel
    #     except:
    #         parallel = True
    #     try:
    #         dsts = args.gly_dsts
    #     except:
    #         dsts = None
    #     print(f"Processing {len(d['files']):d} files")
    #     RP = read_physiology.ReadPhysiology(
    #         mode=protocol,
    #         model=args.model,
    #         stim=args.stim,
    #         files=d["files"],  # =args.filename,
    #         nreps=args.nreps,
    #         # dt=d['runInfo']['dt'],
    #         #     dbspl=d['runInfo']['dbspl'],
    #         voltage=self.voltage,  # show target cell voltages
    #         targetcelltype=args.model,
    #         threshold=threshold,  # args.threshold,
    #         parallel=parallel,
    #         cmmrs2n=cmmrs2n,  # args.cmmrs2n,
    #         gly_dsts=d["command_line"].gly_dsts,
    #         index=d,  # pass the index database.
    #     )
    #
    # def analyze_cmmr_summary(self):
    #     print(self.table.selectedIndexes())


def main():
    DT = DataTables()
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QApplication.instance().exec_()


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == "__main__":
    main()
