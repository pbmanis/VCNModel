import dataclasses
from dataclasses import dataclass, field
import datetime
import pickle
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Union

import cnmodel
import numpy as np
from pylibrary.tools import cprint as CP


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
    filetype: str = "F"
    cellType: str = ""
    modelType: str = ""
    modelName: str = ""
    dendriteMode: str=""
    command_line: str = ""
    changed_data: Union[None, str] = None  # changed_data
    cnmodel_hash: str = cnmodel_git_hash  # save hash for the model code
    files: List = field(default_factory=defemptylist)
    vcnmodel_code_hash: str = git_head_hash  # this repository!
    simulation_path: str = ""
    temperature: float = 0.0
    dt: Union[float, None] = None
    threshold: Union[float, None] = None
    dBspl: float = 0.0
    nReps: int = 0
    SRType: str = ""
    elapsed: float = 0.0
    runProtocol: str = ""
    synapsetype: str = ""
    synapse_experiment: str = ""


class TableManager:
    def __init__(self, table, basepath, selvals, altcolormethod):
        self.table = table
        self.basepath = basepath
        self.selvals = selvals
        self.altColors = altcolormethod

    def force_suffix(self, filename, suffix='.pkl'):
        fn = Path(filename)
        if fn.suffix != suffix:
            fn = str(fn)
            fn = fn + suffix
            fn = Path(fn)
        return fn
        
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
            The directory to check.
        
        Returns
        ------
        Contents of the index file.
        """

        cprint('b', indexdir)
        indexfile = self.force_suffix(indexdir)
        cprint("c", f"Checking for index file: {str(indexfile):s}")
        if indexfile.is_file() and not force:
            cprint("g", f"    Found index file, reading")
            with open(indexfile, "rb") as fh:
                d = pickle.load(fh, encoding="latin1")
            return d
        if force or not indexfile.is_file():
            cprint("c", f"Building a new .pkl index file for {str(indexdir):s}")
            # print(indexfile)
  #           print(indexfile.is_file())
            dpath = Path(indexfile.parent, indexfile.stem)
            # cprint("c", dpath)
            runs = list(dpath.glob("*.p"))
            # print('runsa: ', runs)
            if len(runs) == 0:
                return None
            for r in runs:
                p = self.read_pfile_params(r)
                if p is None:
                    return None
                else:
                    pars, runinfo, indexfile = p
            cprint('m', f"to build indexdir: {str(indexdir):s}")
            indexdata = self.write_indexfile(pars, runinfo, indexdir)
            return indexdata

    def read_pfile_params(self, datafile) -> Union[tuple, None]:
        """
        Reads the Params and runinfo entry from the simulation data file
        """
        print("Reading pfile: ", str(datafile.name))
        try:
            with open(datafile, "rb") as fh:
                d = pickle.load(fh, encoding="latin1")
        except IOError:
            raise IOError('SKIPPING: File is too old; re-run for new structure')
            return None
        if 'runInfo' not in list(d.keys()):
            cprint('r', 'SKIPPING: File is too old (missing "runinfo"); re-run for new structure')
            return(None)
        if "Params" not in list(d.keys()):
            cprint('r', 'SKIPPING: File is too old (missing "Params"); re-run for new structure')
            return(None)
        return (d["Params"], d["runInfo"], str(datafile.name))  # just the dicts

    def make_indexdata(self, params, runinfo, fn=None, indexdir=None):
        """
        Load up the index data class with selected information
        from the Param and runInfo data classes in the simulation
        run files, then write the index file.
        """
        Index_data = IndexData()
        usedict = False
        if indexdir is not None:
            Index_data.command_line = vars(params.commandline)
            Index_data.filetype = "D"  # directory
        else:
            Index_data.filetype = "F"  # j ust one file
            if isinstance(params, OrderedDict):
                usedict = True
                Index_data.command_line = params["commandline"]
            else:
                Index_data.command_line = params.commandline

        Index_data.changed_data = None  # changed_data
        Index_data.cnmodel_hash = cnmodel_git_hash  # save hash for the model code
        Index_data.vcnmodel_code_hash = git_head_hash  # this repository!
        Index_data.simulation_path = str(
            Path(
                self.basepath,
                f"VCN_c{int(self.selvals['Cells'][1]):02d}",
                "Simulations",
                self.selvals["Run Type"][1],
            )
        )        
        if Index_data.filetype == "F":

            mtime = Path(Index_data.simulation_path, fn).stat().st_mtime
        elif Index_data.filetype == "D":
            mtime = Path(Index_data.simulation_path, indexdir).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d-%H:%M")
        Index_data.datestr = timestamp_str


        Index_data.threshold = None
        Index_data.elapsed = 0.0
        if indexdir is not None:
            runs = indexdir.glob("*.p")
            for r in runs:
                Index_data.files.append(str(r))  # get the result file list here
        else:
            Index_data.files = [fn]
        print('Params: ' , params)
        print('runinfo: ', runinfo)
        if usedict:  # old style files with dictionaries
            Index_data.cellType = params["cellType"]
            Index_data.modelType = params["modelType"]
            Index_data.modelName = params["modelName"]
            Index_data.temperature = 37.0  # params["Celsius"]  Not always in the file information we retrieve

            try:
                Index_data.dt = params["dt"]
            except:
                Index_data.dt = 20e-6
            Index_data.runProtocol = params["runProtocol"]
            Index_data.synapsetype = params["ANSynapseType"]
            Index_data.synapse_experiment = params["spirou"]
            Index_data.dBspl = params["dB"]
            Index_data.SRType = params["SRType"]
        else:  # new style files with dataclasses
            Index_data.cellType = params.cellType
            Index_data.modelType = params.modelType
            Index_data.modelName = params.modelName
            try:
                Index_data.dendriteMode = params.dendriteMode
            except:
                Index_data.dendriteMode = "normal"
            Index_data.temperature = params.celsius
            Index_data.dt = params.dt
            Index_data.runProtocol = runinfo.runProtocol
            Index_data.synapsetype = params.ANSynapseType
            Index_data.synapse_experiment = runinfo.Spirou
            Index_data.dBspl = runinfo.dB
            Index_data.nReps = runinfo.nReps
            Index_data.SRType = params.SRType

        return Index_data
        
    def write_indexfile(self, params, runinfo, indexdir):
        Index_data = self.make_indexdata(params, runinfo, indexdir=indexdir)
        indexfile = self.force_suffix(indexdir)
        with open(indexfile, "wb") as fh:
            pickle.dump(Index_data, fh)
        return Index_data

    def read_indexfile(self, indexfilename):
        """
        Read the index file that we will use to populate the table, and to provide "hidden"
        information such as the file list, for analyses.
        
        """
        indexfilename = self.force_suffix(indexfilename)
        with open(indexfilename, "rb") as fh:
            indexdata = pickle.load(fh, encoding="latin1")
        return indexdata

    def build_table(self, mode="scan"):
        if mode == 'scan':
            force = False
        if mode == 'update':
            force = True
        self.data = []
        self.table.setData(self.data)
        # cprint('m', self.basepath)
        # cprint('y', f"VCN_c{int(self.selvals['Cells'][1]):02d}")
        # cprint('g', self.selvals["Run Type"][1])
        thispath = Path(
            self.basepath,
            f"VCN_c{int(self.selvals['Cells'][1]):02d}",
            "Simulations",
            self.selvals["Run Type"][1],
        )
        # cprint('c', thispath)
        rundirs = list(thispath.glob("*"))
        indxs = []
        # first build elements with directories
        indexable_dirs = sorted([r for r in rundirs if r.is_dir()])
        indexfiles = []
        if len(indexable_dirs) > 0:
            cprint("c", "Indexable Directories: ")
            for d in indexable_dirs:
                cprint("c", f"     {str(d):s}")
                self.find_build_indexfiles(
                    d, force=force
                )  # get the indexed data and update as necessary.
            # load up the index files
            indexfiles = list(thispath.glob("*.pkl"))
            for i, f in enumerate(indexfiles):
                indxs.append(self.read_indexfile(f))
                # if indxs[-1] is None:
                # cprint('y', f"None in #{i:d} :{str(f):s}")

        runfiles = list(thispath.glob("*.p"))
        dfiles = sorted(list(runfiles))  # regular data files
        for i, f in enumerate(dfiles):
            p = self.read_pfile_params(f)
            if p is None:
                # indxs.append(None)
                cprint('y', f"None in #{i:d} :{str(f):s}")
            else:
                params, runinfo, fn = p  # ok to unpack
                indxs.append(self.make_indexdata(params, runinfo, fn))
        indexfiles = indexfiles + dfiles

        # print(indxs)
        # for i in range(len(indexfiles)):
        #     try:
        #         print(   indxs[i].command_line["nReps"],)
        #     except:
        #         print("......", indxs[i].command_line)
        # print('*'*80)
        # for ix in indxs:
        #            print(ix, end="\n\n")
        self.table_data = indxs
        # transfer to the data array for the table
        self.data = np.array(
            [
                (
                    indxs[i].filetype,
                    indxs[i].datestr,
                    indxs[i].cellType,
                    indxs[i].modelType,
                    indxs[i].modelName,
                    indxs[i].runProtocol,
                    indxs[i].dendriteMode,
                    indxs[i].synapse_experiment,
                    indxs[i].dBspl,
                    indxs[i].nReps, # command_line["nReps"],
                    len(indxs[i].files),
                    str(Path("/".join(indexfiles[i].parts[-4:]))),
                )
                for i in range(len(indxs))
            ],
            dtype=[
                ("type", object),
                ("date", object),
                ("cell", object),
                ("modelType", object),
                ("modelName", object),
                ("Protocol", object),
                ("DendMode", object),
                ("Experiment", object),
                # ("inputtype", object),
                # ("modetype", object),
                # ("scaling", object),
                # ("Freq", object),
                ("dBspl", object),
                ("nReps", object),
                ("# Files", object),
                ("Filename", object),
                # ("SRType", object),
                # ("# files", object),
                # ("synapsetype", object),
            ],
        )
        self.table.setData(self.data)
       
        style = "section:: {font-size: 4pt; color:black; font:TimesRoman;}"
        self.table.setStyleSheet(style)
        self.altColors()
        # self.table.setStyle(QtGui.QFont('Arial', 6))
        self.table.resizeRowsToContents()
        self.table.resizeColumnsToContents()

    def setColortoRow(self, rowIndex, color):
        for j in range(self.table.columnCount()):
            self.table.item(rowIndex, j).setBackground(color)

    def get_table_data(self, index_row):
        """
        Regardless of the sort, read the current index row and
        map it back to the data in the table.
        """
        index = self.table.currentIndex()
        value = index.sibling(index.row(), 1).data()
        for i, d in enumerate(self.data):
            if self.data[i][1] == value:
                return self.table_data[i]
        return None
