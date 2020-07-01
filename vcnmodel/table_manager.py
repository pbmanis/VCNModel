import sys
import subprocess

import numpy as np
import functools
from pathlib import Path
import pickle
import time
import importlib
import toml
import datetime
from dataclasses import dataclass, field
from typing import Union, Dict, List
import dataclasses

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph import Qt

import pylibrary.tools.cprint as CP
import cnmodel

"""
Fill in a table using pyqtgrph that summarizes the data in a directory.
Generates index files (pickled python dataclass) that summarize the results.

Set the force flag true to rebuild the pkl files.

pbm 7/1/2020

"""

cprint = CP.cprint

# Get the git has of he repositories so we know exactly what code was run (assumed the repo
# was updated and committed ahead of the run)

process = subprocess.Popen(
    ["git", "rev-parse", "HEAD"], shell=False, stdout=subprocess.PIPE
)
git_head_hash = process.communicate()[0].strip()
cnmodelpath = Path(cnmodel.__file__).parent
process = subprocess.Popen(
    ["git", "-C", str(cnmodelpath), "rev-parse", "HEAD"],
    shell=False,
    stdout=subprocess.PIPE,
)
cnmodel_git_hash = process.communicate()[0].strip()


def defemptylist():
    return []

#
# dataclass for the pkl file. Abstracted information from the main data files
# to speed up access.
#
@dataclass
class IndexData:
    datestr: str = ""
    cellType: str=""
    modelType: str=""
    modelName: str=""
    command_line: str = ""
    changed_data: Union[None, str] = None  # changed_data
    cnmodel_hash: str = cnmodel_git_hash  # save hash for the model code
    files: List = field(default_factory=defemptylist)
    vcnmodel_code_hash: str = git_head_hash  # this repository!
    simulationnpath: str = ""
    temperature: float = 0.0
    dt: Union[float, None] = None
    threshold: Union[float, None] = None
    dBspl: float=0.
    SRType: str=""
    elapsed: float = 0.0
    runProtocol: str = ""
    synapsetype: str = ""
    synapse_experiment: str = ""


class TableManager:
    def __init__(self, table, basepath, selvals):
        self.table = table
        self.basepath = basepath
        self.selvals = selvals

    def find_build_indexfiles(self, indexdir: Union[str, Path], force=False):
        """
        Given the indexdir, determine:
        a) whether an index file exists for this directory
            If it does, read it
            If it does not, then make one from the data in the 
            directory.
        
        Parameters
        ----------
        indexdir : (str or Path)
            The directory to check
        
        Returns
        ------
        Contents of the index file.
        """

        indexfile = Path(indexdir).with_suffix(".pkl")
        if indexfile.is_file() and not force:
            cprint("g", f"Found a pkl index file: {str(indexdir):s}")
            with open(indexfile, "rb") as fh:
                d = pickle.load(fh, encoding="latin1")
            return d
        else:
            cprint("c", f"Building a pkl index file for {str(indexdir):s}")
            runs = indexdir.glob("*.p")
            for r in runs:
                pars, runinfo = self.read_data_params(r)
            indexdata = self.write_indexfile(pars, runinfo, indexdir)
            return indexdata

    def read_data_params(self, datafile):
        """
        Reads the Params and runinfo entry from the simulation data file
        """
        with open(datafile, "rb") as fh:
            d = pickle.load(fh, encoding="latin1")
        # pars = dataclasses.asdict(d['Params'])
        return d["Params"], d["runInfo"]  # just return as dataclasses

    def write_indexfile(self, params, runinfo, indexdir):
        """
        Load up the index data class with selected information
        from the Param and runInfo data classes in the simulation
        run files, then write the index file.
        """
        Index_data = IndexData()
        Index_data.command_line = vars(params.commandline)
        Index_data.changed_data = None  # changed_data
        Index_data.cnmodel_hash = cnmodel_git_hash  # save hash for the model code
        Index_data.vcnmodel_code_hash = git_head_hash  # this repository!
        starttime = datetime.datetime.now()
        runtime = starttime.strftime("%Y-%m-%d.%H-%M-%S")
        Index_data.datestr = runtime
        Index_data.simulationnpath = str(
            Path(
                self.basepath,
                f"VCN_c{int(self.selvals['Cells'][1]):02d}",
                "Simulations",
                self.selvals["Run Type"][1],
            )
        )

        Index_data.cellType = params.cellType
        Index_data.modelType = params.modelType
        Index_data.modelName = params.modelName
        Index_data.temperature = params.celsius
        Index_data.dt = params.dt
        Index_data.threshold = None
        Index_data.elapsed = 0.0
        runs = indexdir.glob("*.p")
        for r in runs:
            Index_data.files.append(str(r))  # get the result file list here
        Index_data.runProtocol = runinfo.runProtocol
        Index_data.synapsetype = params.ANSynapseType
        Index_data.synapse_experiment = runinfo.Spirou
        Index_data.dBspl = runinfo.dB
        Index_data.SRType = params.SRType
        indexfile = Path(indexdir).with_suffix(".pkl")
        with open(indexfile, "wb") as fh:
            pickle.dump(Index_data, fh)
        return Index_data


    def read_indexfile(self, indexfilename):
        """
        Read the index file that we will use to populate the table, and to provide "hidden"
        information such as the file list, for analyses.
        
        """
        with open(indexfilename, "rb") as fh:
            indexdata = pickle.load(fh, encoding="latin1")
        return indexdata

        # then save off to indexfile.

    def build_table(self, mode="scan"):
        self.data = []
        self.table.setData(self.data)
        thispath = Path(
            self.basepath,
            f"VCN_c{int(self.selvals['Cells'][1]):02d}",
            "Simulations",
            self.selvals["Run Type"][1],
        )

        rundirs = thispath.glob("*")
        indexable_dirs = [r for r in rundirs if r.is_dir()]
        cprint('r', f"indexabls Directories: {str(indexable_dirs):s}")
        for indexdir in indexable_dirs:
            self.find_build_indexfiles(
                indexdir, force=False
            )  # get the indexed data and update as necessary.

        # load up the index files
        indexfiles = list(thispath.glob("*.pkl"))      
        indxs = []
        for i, f in enumerate(indexfiles):
            indxs.append(self.read_indexfile(f))
        self.table_data = indxs
        #transfer to the data array for the table
        self.data = np.array(
            [
                (
                    str(Path("/".join(indexfiles[i].parts[-4:]))),
                    indxs[i].datestr,
                    indxs[i].cellType,
                    indxs[i].modelType,
                    indxs[i].modelName,
                    indxs[i].runProtocol,
                    indxs[i].dBspl,
                    len(indxs[i].files),
                    
                )
                for i in range(len(indexfiles))
            ],
            dtype=[
                ("Filename", object),
                ("datestr", object),
                ("cellType", object),
                ("modelType", object),
                ("modelName", object),
                ("runProtocol", object),
                # ("inputtype", object),
                # ("modetype", object),
                # ("scaling", object),
                # ("Freq", object),
                ("dBspl", object),
                ("# Files", object),
                # ("SRType", object),
                # ("# files", object),
                # ("synapsetype", object),
            ],
        )
        style = "section:: {font-size: 4pt; color:black; font:TimesRoman;}"
        self.table.setData(self.data)
        self.table.setStyleSheet(style)
        self.altColors([QtGui.QColor(0x00, 0x00, 0x00), QtGui.QColor(0x22, 0x22, 0x22)])
        # self.table.setStyle(QtGui.QFont('Arial', 6))
        self.table.resizeRowsToContents()
        self.table.resizeColumnsToContents()

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
