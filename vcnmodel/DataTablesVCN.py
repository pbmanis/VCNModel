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
import pyqtgraph.dockarea as PGD
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
from vcnmodel import table_manager as table_manager
from vcnmodel import cell_config
from vcnmodel.plotters import plot_sims
from vcnmodel.plotters import figures
from vcnmodel.plotters import plot_z
from vcnmodel.plotters import efficacy_plot
from pylibrary.tools import cprint as CP
# import vcnmodel.correlation_calcs
import vcnmodel.spikestatistics
import ephys


"""
This program provides a graphical interface for model results for the SBEM reconstruction project.
The display appears as 3 panels: One on the left with controls, one on the top right that is tabbed,
    showing either the current table, or the traces, and one on the bottom for text output.
The left panel provides a set of organized controls:
    Selections:
        Run Type (AN, IV) corresponding to auditory nerve input or current injection protocols
        Cells (from list of cells that are available to model)
            Selecting one of these will population the simulation table on the right with valid simulations.
        ModelType Mode, Experiment, Analysis, Dendrites: inactive.
    Analysis:
        This provides different fixed kinds of analysis for the model data.
        Traces: just plot the traces, stacked, for reference.
        IV : plot current-voltage relationships and calculate Rin, Taum, find spikes.
        VC : plot current voltage relationships in voltage clamp.
        Singles: For the "single" AN protocol, where only one input at a time is active, creates stacked plot
        Trace Viewer : dynamic plot of APs and preceding times for AN inputs in the "Traces" tab
         RevcorrSPKS : reverse correlation against postsynaptic spikes for each input. 
        RevcorrSimple : same as revcorrSPKS, but simpler algorithm. 
        RevcorrSTTC : not implemented.
        PSTH : Plot PSTH, raster, for bu cell and AN input; also compute phase locking to AM if needed.
    Filters:
        This provides data selection in the table. Most entries provide a drop-down list. The values that are not 
            None are applied with "and" logic. The filters are not applied until the Apply button is pressed;
            The filters are cleared by the Clear button.
    Options:
        These are options for the TraceViewer mode.
        Nubmer of traces
        Plot Vm or dVm/dt
        Movie button generates a movie through time.
        Frame interval sets the time between frames in the movie, in msec.
    Figures:
        Interface to figure generation. Figures are generated from the model data directly as much as possible.
        Some figures are generated from analysis data that is either compiled manually, or using a script.
    Tools:
        Reload: for all modules under DataTables, reload the code. Mostly used during development.
        View IndexFile: Print the index file in the text window.
        Print File Info: Prints the file info for the selected entries into the text window.
        Delete Selected Sim : for deleting simulations that are broken (or stopped early). 
    Quit:
        Exit the program.




Uses pyqtgraph tablewidget to build a table showing simulation
files/runs and enabling analysis via a GUI

    Supported primarily by R01DC015901 (Spirou, Manis, Ellisman),
    Early development: R01 DC004551 (Manis, until 2019)
    Later development: R01 DC019053 (Manis, 2020-2025)

"""


cprint = CP.cprint

# List reloadable modules
all_modules = [
    table_manager,
    plot_sims,
    figures,
    plot_z,
    efficacy_plot,
    cell_config,
    # vcnmodel.correlation_calcs,
    vcnmodel.spikestatistics,
    vcnmodel.analysis,
    ephys.ephysanalysis.SpikeAnalysis,
    ephys.ephysanalysis.Utility,
    ephys.ephysanalysis.MakeClamps,
    PH,
]

# Define the cell values in the Cells dropdown

cellvalues = [
    2,
    5,
    6,
    9,
    10,
    11,
    13,
    17,
    18,
    # 24,
    # 29,
    30,
]

# model types known to use - These define the decoration patterns.
# The data in the paper uses XM13A_nacncoop
modeltypes = [
    "None",
    "XM13_nacncoop",
    "XM13A_nacncoop",
    "mGBC",
    "XM13",
    "RM03",
]
# Types of runs - not all of these are used
runtypes = ["AN", "IV", "VC", "IO", "gifnoise"]

# types of experiments known to model_run2
experimenttypes = [
    "None",
    "all",
    "delays",
    "largestonly",
    "removelargest",
    "mean",
    "max=mean",
    "all=mean",
    "twolargest",
    "threelargest",
    "fourlargest",
]
# Modes for synapse model runs - not all are used.
modetypes = ["find", "singles", "IO", "multi"]
# For analysis - but not used.
analysistypes = ["traces", "PSTH", "revcorr", "SAC", "tuning", "traceviewer"]
# Known revers correlation types
revcorrtypes = [
    "RevcorrSPKS",
    "RevcorrSimple",
    "RevcorrSTTC",
]
# For Filter dropdown: 10 dB steps
dbValues = list(range(0, 91, 10))
dbValues.insert(0, "None")

# dendrite experiments.
dendriteChoices = [
    "None",
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
    """
    Main entry point for building the table and operating on the data
    """
    def __init__(self):
        self.PLT = plot_sims.PlotSims(parent=self)
        self.FIGS = figures.Figures(parent=self)
        self.QColor = QtGui.QColor  # for access in plotsims (pass so we can reload)
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

        # Define the table style for various parts
        # dark scheme
        from pyqtgraph.Qt import QtCore

        dark_palette = QtGui.QPalette()
        white = self.QColor(255, 255, 255)
        black = self.QColor(0, 0, 0)
        red = self.QColor(255, 0, 0)
        dark_palette.setColor(QtGui.QPalette.Window, self.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.WindowText, white)
        dark_palette.setColor(QtGui.QPalette.Base, self.QColor(25, 25, 25))
        dark_palette.setColor(QtGui.QPalette.AlternateBase, self.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.ToolTipBase, white)
        dark_palette.setColor(QtGui.QPalette.ToolTipText, white)
        dark_palette.setColor(QtGui.QPalette.Text, white)
        dark_palette.setColor(QtGui.QPalette.Button, self.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.ButtonText, white)
        dark_palette.setColor(QtGui.QPalette.BrightText, red)
        dark_palette.setColor(QtGui.QPalette.Link, self.QColor(42, 130, 218))
        dark_palette.setColor(QtGui.QPalette.Highlight, self.QColor(42, 130, 218))
        dark_palette.setColor(QtGui.QPalette.HighlightedText, black)

        self.app.setPalette(dark_palette)
        self.app.setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }")

        self.win = pg.QtGui.QMainWindow()
        # use dock system instead of layout.
        self.dockArea = PGD.DockArea()
        self.win.setCentralWidget(self.dockArea)
        self.win.setWindowTitle("Model DataTables/FileSelector")
        self.win.resize(1600, 1024)
        # Initial Dock Arrangment
        self.Dock_Params = PGD.Dock("Params", size=(250, 1024))
        self.Dock_Table = PGD.Dock("Simulation Table", size=(1000, 800))
        self.Dock_Report = PGD.Dock("Reporting", size=(1000, 200))
        self.Dock_Traces = PGD.Dock("Traces", size=(1000, 700))

        self.dockArea.addDock(self.Dock_Params, "left")
        self.dockArea.addDock(self.Dock_Table, "right", self.Dock_Params)
        self.dockArea.addDock(self.Dock_Traces, "below", self.Dock_Table)
        self.dockArea.addDock(self.Dock_Report, "bottom", self.Dock_Table)
        # self.dockArea.addDock(self.Dock_Traces_Slider, 'below', self.Dock_Traces)

        # self.Dock_Traces.addContainer(type=pg.QtGui.QGridLayout, obj=self.trace_layout)
        self.table = pg.TableWidget(sortable=True)
        self.Dock_Table.addWidget(self.table)
        self.Dock_Table.raiseDock()

        style = "::section {background-color: darkblue; }"
        self.selected_index_row = None  # for single selection mode
        self.selected_index_rows = None  # for multirow selection mode
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
        self.filters = {
            "Use Filter": False,
            "dBspl": None,
            "nReps": None,
            "pipDur": None,
            "fmod": None,
            "Protocol": None,
            "Experiment": None,
            "modelName": None,
            "dendMode": None,
            "DataTable": None,
            
        }
        self.trvalues = [1, 2, 4, 8, 16, 32]
        self.n_trace_sel = self.trvalues[2]
        self.V_disp = ["Vm", "dV/dt"]
        self.V_disp_sel = self.V_disp[0]
        self.movie_state = False
        self.frame_intervals = [0.033, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
        self.frame_interval = self.frame_intervals[3]
        self.target_figure = "IV Figure"
        
        # We use pyqtgraph's ParameterTree to set up the menus/buttons.
        # This defines the layout.
        self.params = [
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
                        "name": "IV",
                        "type": "action",
                        # "default": False,
                    },
                    {
                        "name": "VC",
                        "type": "action",
                        # "default": False,
                    },
                    {
                        "name": "Singles",
                        "type": "action",
                        # "default": False,
                    },
                    {"name": "Trace Viewer", "type": "action",},
                    {"name": "RevcorrSPKS", "type": "action"},
                    {"name": "RevcorrSimple", "type": "action"},
                    {"name": "RevcorrSTTC", "type": "action"},
                    {"name": "PSTH", "type": "action"},
                ],
            },
            {
                "name": "Filters",
                "type": "group",
                "children": [
                    # {"name": "Use Filter", "type": "bool", "value": False},
                    {
                        "name": "dBspl",
                        "type": "list",
                        "values": dbValues,
                        "value": "None",
                    },
                    {
                        "name": "nReps",
                        "type": "list",
                        "values": ["None", 1, 5, 10, 20, 25, 50, 100],
                        "value": "None",
                    },
                    {
                      "name": "pipDur",
                      "type": "list",
                      "values": ["None", 0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 10.0],
                      "value": "None",  
                    },
                    {
                        "name": "Protocol",
                        "type": "list",
                        "values": ["None", "runIV", "runANPSTH", "runANSingles"],
                        "value": "None",
                    },
                    {
                        "name": "Experiment",
                        "type": "list",
                        "values": experimenttypes,
                        "value": "None",
                    },
                    {
                        "name": "modelName",
                        "type": "list",
                        "values": modeltypes,
                        "value": "None",
                    },
                    {"name": "dataTable", "type": "str", "value": "None"},
                    {
                        "name": "dendMode",
                        "type": "list",
                        "values": ["None", "actdend", "normal", "pasdend"],
                        "value": "None",
                    },
                    {
                        "name": "soundType",
                        "type": "list",
                        "values": [
                            "None",
                            "tonepip",
                            "noise",
                            "SAM",
                            "CMMR",
                            "stationaryNoise",
                        ],
                        "value": "None",
                    },
                    {
                        "name": "fmod",
                        "type": "list",
                        "values": ["None", 10, 20, 50, 100, 200, 300, 400, 500, 750, 1000],
                        "value": "None",
                    },
                    {
                        "name": "Filter Actions",
                        "type": "group",
                        "children": [
                            {"name": "Apply", "type": "action"},
                            {"name": "Clear", "type": "action"},
                        ],
                    },
                ],
            },
            {
                "name": "Options",
                "type": "group",
                "children": [
                    {
                        "name": "Viewer #traces",
                        "type": "list",
                        "values": self.trvalues,
                        "value": self.n_trace_sel,
                    },
                    {
                        "name": "Vm or dV/dt",
                        "type": "list",
                        "values": self.V_disp,
                        "value": self.V_disp_sel,
                    },
                    {"name": "Movie", "type": "action"},
                    {
                        "name": "Frame Interval",
                        "type": "list",
                        "values": self.frame_intervals,
                        "value": self.frame_interval,
                    },
                ],
            },
            {
                "name": "Figures",
                "type": "group",
                "children": [
                    {"name": "Figures", "type": "list", 
                    "values": ["IV Figure", "All IVs", "IV Supplement", "Zin Supplement",
                               "Efficacy", "Efficacy Supplement",
                               "Revcorr Ex", "All Revcorr", "Revcorr Supplement", "Revcorr Compare", 
                               "PSTH-FSL", "All PSTHs", "PSTH-FSL Supplement",
                               "VS-SAM Tone", "VC-KLTCalibration"],
                    "value": "IV Figure",
                },
                    {"name": "Create Figure", "type": "action"},
                ],
            },
            {
                "name": "Tools",
                "type": "group",
                "children": [
                    {"name": "Reload", "type": "action"},
                    {"name": "View IndexFile", "type": "action"},
                    {"name": "Print File Info", "type": "action"},
                    {"name": "Delete Selected Sim", "type": "action"},
                ],
            },
            {"name": "Quit", "type": "action"},
        ]
        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(
            name="Models", type="group", children=self.params
        )
        self.ptree.setParameters(self.ptreedata)
        self.ptree.setMaximumWidth(300)
        self.ptree.setMinimumWidth(250)

        self.Dock_Params.addWidget(self.ptree)  # put the parameter three here

        self.trace_plots = pg.PlotWidget(title="Trace Plots")
        self.Dock_Traces.addWidget(self.trace_plots, rowspan=5, colspan=1)
        self.trace_plots.setXRange(-5.0, 2.5, padding=0.2)
        self.trace_plots.setContentsMargins(10, 10, 10, 10)
        # Build the trace selector
        self.trace_selector = pg.graphicsItems.InfiniteLine.InfiniteLine(
            0, movable=True, markers=[("^", 0), ("v", 1)]
        )
        self.trace_selector.setPen((255, 255, 0, 200))  # should be yellow
        self.trace_selector.setZValue(1)
        self.trace_selector_plot = pg.PlotWidget(title="Trace selector")
        self.trace_selector_plot.hideAxis("left")
        self.frameTicks = pg.graphicsItems.VTickGroup.VTickGroup(
            yrange=[0.8, 1], pen=0.4
        )
        self.trace_selector_plot.setXRange(0.0, 10.0, padding=0.2)
        self.trace_selector.setBounds((0, 10))
        self.trace_selector_plot.addItem(self.frameTicks, ignoreBounds=True)
        self.trace_selector_plot.addItem(self.trace_selector)
        self.trace_selector_plot.setMaximumHeight(100)
        self.trace_plots.setContentsMargins(10.0, 10.0, 10.0, 10.0)

        # print(dir(self.trace_selector_plot))
        self.Dock_Traces.addWidget(
            self.trace_selector_plot, row=5, col=0, rowspan=1, colspan=1
        )

        self.textbox = QtWidgets.QTextEdit()
        self.textbox.setReadOnly(True)
        self.textbox.setText("Text Edit box (RO)")
        self.Dock_Report.addWidget(self.textbox)

        self.win.show()
        self.table.setSelectionMode(pg.QtGui.QAbstractItemView.ExtendedSelection)
        self.table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.table_manager = table_manager.TableManager(
            parent=self,
            table=self.table,
            basepath=self.basepath,
            selvals=self.selvals,
            altcolormethod=self.altColors,
        )
        self.table.itemDoubleClicked.connect(
            functools.partial(self.on_double_Click, self.table)
        )
        self.table.clicked.connect(functools.partial(self.on_Single_Click, self.table))
        self.ptreedata.sigTreeStateChanged.connect(self.command_dispatcher)
        # Ok, we are in the loop - anything after this is menu-driven and handled
        # either as part of the TableWidget, the Traces widget, or through the
        # CommandDispatcher.

    def setPaths(self, stimtype="AN", cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"cellDataDirectory": Path("../VCN-SBEM-Data", "VCN_Cells",)}
        self.basepath = self.datapaths["cellDataDirectory"]

    def on_double_Click(self, w):
        index = w.selectionModel().currentIndex()
        self.selected_index_row = index.row()
        self.analyze_from_table(index.row())

    def on_Single_Click(self, w):
        selrows = w.selectionModel().selectedRows()
        self.selected_index_rows = selrows
        if len(selrows) == 0:
            self.selected_index_rows = None
        # for index in selrows:
        #     self.selected_index_row = index.row()
        #     self.analyze_from_table(index.row())

    def handleSortIndicatorChanged(self, index, order):
        # print(dir(self.table.model()))
        pass
        # if index != 0:

    #           self.table.horizontalHeader().setSortIndicator(
    #               0, self.table.model().sortOrder()
    #           )

    def command_dispatcher(self, param, changes):
        """
        Dispatcher for the commands from parametertree
        path[0] will be the command name
        path[1] will be the parameter (if there is one)
        path[2] will have the subcommand, if there is one
        data will be the field data (if there is any)
        """
        for param, change, data in changes:
            path = self.ptreedata.childPath(param)

            if path[0] == "Quit":
                exit()
            if path[0] == "Scan Runs":
                self.table_manager.build_table(mode="scan")
            if path[0] == "Update Runs":
                self.table_manager.build_table(mode="update")

            if path[0] == "Selections":
                self.selvals[path[1]][1] = str(data)
                self.cellID = self.selvals["Cells"][1]
                self.setPaths("AN", cell=data)
                self.table_manager.build_table(mode="scan")

            if path[0] == "Analysis":
                analyze_map = {
                    "Traces": self.analyze_traces,
                    "IV": self.analyze_IV,
                    "VC": self.analyze_VC,
                    "PSTH": self.analyze_PSTH,
                    "Singles": self.analyze_singles,
                    # 'CMMR_Summary': self.analyze_cmmr_summary,
                    "RevcorrSPKS": self.analyze_revcorr,
                    "RevcorrSimple": self.analyze_revcorr,
                    "RevcorrSTTC": self.analyze_revcorr,
                    "Trace Viewer": self.trace_viewer,
                }
                analyze_map[path[1]](path[1])  # call the analysis function

            if path[0] == "Filters":
                if path[1] == "Use Filter":  # currently not an option
                    # print(data)
                    self.filters["Use Filter"] = data
                elif path[1] in ["dBspl", "nReps", "fmod"]:
                    # print('dbspl/nreps: ', data)
                    if data != "None":
                        self.filters[path[1]] = int(data)
                    else:
                        self.filters[path[1]] = None
                elif path[1] in [
                    "Protocol",
                    "Experiment",
                    "modelName",
                    "dendMode",
                    "soundType",
                ]:
                    if data != "None":
                        self.filters[path[1]] = str(data)
                    else:
                        self.filters[path[1]] = None
                elif path[1] in ['pipDur']:
                    if data != "None":
                        self.filters[path[1]] = float(data)
                    else:
                        self.filters[path[1]] = None
                elif path[1] == "Filter Actions":
                    if path[2] in ["Apply"]:
                        self.filters["Use Filter"] = True
                        self.table_manager.apply_filter(QtCore=QtCore)
                    elif path[2] in ["Clear"]:
                        self.filters["Use Filter"] = False
                        self.table_manager.apply_filter(QtCore=QtCore)
                else:
                    if data != "None":
                        self.filters[path[1]] = float(data)
                    else:
                        self.fliters[path[1]] = None

                # print('Filters: ', self.filters)
            if path[0] == "Options":
                if path[1] == "Viewer #traces":
                    self.n_trace_sel = int(data)
                elif path[1] == "Vm or dV/dt":
                    self.V_disp_sel = str(data)
                elif path[1] == "Movie":
                    if not self.movie_state:
                        self.movie_state = True
                        self.trace_viewer()
                    else:
                        self.movie_state = False
                elif path[1] == "Frame Interval":
                    self.frame_interval = float(data)

            if path[0] == "Figures":
                if path[1] == "Figures":
                    self.target_figure = data
                elif path[1] == "Create Figure":
                    self.FIGS.make_figure(self.target_figure)

            if path[0] == "Tools":
                if path[1] == "Reload":
                    print("reloading...")
                    for module in all_modules:
                        print("reloading: ", module)
                        importlib.reload(module)
                    self.PLT = plot_sims.PlotSims(parent=self)
                    self.table_manager = table_manager.TableManager(
                        parent=self,
                        table=self.table,
                        basepath=self.basepath,
                        selvals=self.selvals,
                        altcolormethod=self.altColors,
                    )
                    

                    print("   reload ok")
                    print("-" * 80)
                    self.table_manager.build_table(mode="scan")
                    self.table.setSortingEnabled(True)
                    self.table.horizontalHeader().sortIndicatorChanged.connect(
                        self.handleSortIndicatorChanged
                    )
                    self.table_manager.apply_filter()
                    self.table.sortByColumn(1, QtCore.Qt.AscendingOrder)  # by date
                    self.altColors()  # reset the color list.
                    self.Dock_Table.raiseDock()
                    self.FIGS = figures.Figures(parent=self)
                    
                elif path[1] == "View IndexFile":
                    selected = self.table.selectionModel().selectedRows()
                    if selected is None:
                        return
                    index_row = selected[0]
                    self.table_manager.print_indexfile(index_row)
                elif path[1] == "Print File Info":
                    selected = self.table.selectionModel().selectedRows()
                    if selected is None:
                        return
                    self.PLT.print_file_info(selected)
                elif path[1] == "Delete Selected Sim":
                    if self.selected_index_rows is None:
                        return
                    self.selected_index_rows = self.table.selectionModel().selectedRows()
                    self.table_manager.remove_table_entry(self.selected_index_rows)
                    
                    
    def error_message(self, text):
        self.textbox.clear()
        color = "red"
        self.textbox.setTextColor(self.QColor(color))
        self.textbox.append(text)
        self.textbox.setTextColor(self.QColor("white"))
        
    def setColortoRow(self, rowIndex, color):
        for j in range(self.table.columnCount()):
            self.table.item(rowIndex, j).setBackground(color)

    def altColors(
        self, colors=[QtGui.QColor(0x00, 0x00, 0x00), QtGui.QColor(0x22, 0x22, 0x22)]
    ):
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

    def force_suffix(self, filename, suffix=".pkl"):
        fn = Path(filename)
        if fn.suffix != suffix:
            fn = str(fn)
            fn = fn + suffix
            fn = Path(fn)
        return fn

    # Here we provide a few management routines for
    # specific actions. These mostly call routins in 
    # plot_sims.py
    
    def analyze_singles(self, ana_name=None):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        index_row = self.selected_index_rows[0]
        selected = self.table_manager.get_table_data(index_row)  # table_data[index_row]
        if selected is None:
            return
        self.PLT.analyze_singles(index_row, selected)

    def trace_viewer(self, ana_name=None):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        index_row = self.selected_index_rows[0]
        selected = self.table_manager.get_table_data(index_row)  # table_data[index_row]
        if selected is None:
            return
        nfiles = len(selected.files)
        print(" nfiles: ", nfiles)
        print("selected files: ", selected.files)
        # if nfiles > 1:
        #     self.PLT.textappend('Please select only one file to view')
        # else:
        PD = plot_sims.PData()
        self.PLT.trace_viewer(selected.files[0], PD, selected.runProtocol)

    def analyze_traces(self, ana_name=None):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        self.plot_traces(
            rows=len(self.selected_index_rows),
            cols=1,
            height=len(self.selected_index_rows),
            width=6.0,
            stack=True,
            ymin=-80.0,
            ymax=20.0,
        )

    def analyze_IV(self, ana_name=None):
        if self.selected_index_rows is None:
            return
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        self.plot_traces(
            rows=1,
            cols=len(self.selected_index_rows),
            height=3.0,
            width=3 * len(self.selected_index_rows),
            stack=False,
            ymin=-120.0,
            ymax=20.0,
        )
        return


    def plot_traces(
        self, rows=1, cols=1, width=5.0, height=4.0, stack=True, ymin=-120.0, ymax=0.0
    ):
        """
        Just redirect to plotsims
        """
        self.PLT.simple_plot_traces(
            rows=rows, cols=cols, width=width, height=height,
                stack=stack, ymin=ymin, ymax=ymax
            )
    
    def analyze_VC(self, ana_name=None):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
            
        self.PLT.plot_VC(self.selected_index_rows)


    def analyze_PSTH(self, ana_name=None):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        P = self.PLT.setup_PSTH()
        PD = plot_sims.PData()
        index_row = self.selected_index_rows[0]
        selected = self.table_manager.get_table_data(index_row)  # table_data[index_row]
        if selected is None:
            return
        sfi = Path(selected.simulation_path, selected.files[0])
        self.PLT.plot_AN_response(P, sfi, PD, selected.runProtocol)
        P.figure_handle.show()

    def analyze_revcorr(self, ana_name):
        revcorrtype = ana_name
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        index_row = self.selected_index_rows[0]
        selected = self.table_manager.get_table_data(index_row)  # table_data[index_row]
        if selected is None:
            return
        self.PLT.plot_revcorr_figure(selected, revcorrtype)

    def analyze_from_table(self, i):
        self.selected_index_rows = self.table.selectionModel().selectedRows()
        if self.selected_index_rows is None:
            return
        index_row = self.selected_index_rows[0]
        selected = self.table_manager.get_table_data(index_row)  # table_data[index_row]
        if selected is None:
            return


def main():
    # Entry point.
    # Why do I do this ? 
    # It keeps sphinxdoc from running the code... in case we do that in the future.
    D = DataTables()  # must retain a pointer to the class, else we die!
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QApplication.instance().exec_()


if __name__ == "__main__":
    main()
