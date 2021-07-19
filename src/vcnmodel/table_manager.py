import dataclasses
from dataclasses import dataclass, field
import datetime
import functools
import pickle
import pprint
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Union

import cnmodel
import numpy as np
from pylibrary.tools import cprint as CP
import pylibrary.tools
import src.vcnmodel.util.fixpicklemodule as FPM


"""
Fill in a table using pyqtgrph that summarizes the data in a directory.
Generates index files (pickled python dataclass) that summarize the results.

Set the force flag true to rebuild the pkl files.

pbm 7/1/2020

"""

cprint = CP.cprint
PP = pprint.PrettyPrinter(indent=8, width=80)

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
    dendriteExpt: str=""
    dendriteMode: str=""
    axonExpt: str=""
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
    pipDur: float = 0.0
    soundType: str = ""
    fmod: Union[float, None]= None
    SRType: str = ""
    ANSynapticDepression: Union[int] = 0
    elapsed: float = 0.0
    runProtocol: str = ""
    synapsetype: str = ""
    synapseExperiment: str = ""
    dataTable: str = ""
    hocfile: str=""


def winprint(func):
    """
    Wrapper decorator for functions that print to the text area
    Clears the print area first,
    and puts a line of '*' when the function returns
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        self.textclear()
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        self.textappend("*" * 80)
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print


def winprint_continuous(func):
    """
    Wrapper decorator for functions that print to the text area
    DOES NOT clear the print area first,
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print


class TableManager:
    def __init__(self, parent=None, table:object=None, basepath:Union[str, Path]="", selvals:dict={}, altcolormethod:object=None):
        assert parent is not None
        self.parent = parent
        self.table = table
        self.basepath = basepath
        self.selvals = selvals
        self.altColors = altcolormethod

    def textclear(self):
        if self.parent is None:
            print("parent is None")
            raise ValueError()
        else:
            self.parent.textbox.clear()

    def textappend(self, text, color="white"):
        if self.parent is None:
            cprint(color, text)  # just go straight to the terminal
        else:
            self.parent.textbox.setTextColor(self.parent.QColor(color))
            self.parent.textbox.append(text)
            self.parent.textbox.setTextColor(self.parent.QColor("white"))

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
            print('          indexfile: ', indexfile)
            with open(indexfile, "rb") as fh:
                # d = FPM.pickle_load(fh, fix_imports=True) # , encoding="latin1")
                d = pickle.load(fh, fix_imports=True) # , encoding="latin1")
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

    @winprint_continuous
    def read_pfile_params(self, datafile) -> Union[tuple, None]:
        """
        Reads the Params and runinfo entry from the simulation data file
        """
        self.textappend(f"Reading pfile: {str(datafile.name):s}", color='white')
        try:
            with open(datafile, "rb") as fh:
                d = pickle.load(fh, fix_imports=True, encoding="latin1")
        except (ModuleNotFoundError, IOError, pickle.UnpicklingError):
            self.textappend('SKIPPING: File is too old; re-run for new structure', color="red")
            self.textappend(f"File: {str(fh):s}", color="red")
            return None
        # except ModuleNotFoundError:
        #     raise IOError('SKIPPING: File has old structure; re-run for new structure')
        #     return None
            
        if 'runInfo' not in list(d.keys()):
            self.textappend('SKIPPING: File is too old (missing "runinfo"); re-run for new structure', color="red")
            # print('  Avail keys: ', list(d.keys()))
            return None
        if "Params" not in list(d.keys()):
            self.textappend('SKIPPING: File is too old (missing "Params"); re-run for new structure', color="red")
            # print('  Avail keys: ', list(d.keys()))
            # print(d['Results'][0])
            return None
        # print(d["Params"])
        # print(d["runInfo"])
        # exit()
        return (d["Params"], d["runInfo"], str(datafile.name))  # just the dicts

    def get_sim_runtime(self, filename):
        """
        Switch the time stamp to different format
        Here the initial value is a string, which we convert to a datetime
      
        """
        with open(filename, "rb") as fh:
            d = pickle.load(fh, fix_imports=True) # d = FPM.pickle_load(fh)
        if d['runInfo'] is None:
            cprint('r', f"runinfo is None? file = {str(filename):s}")
            return None
        if isinstance(d['runInfo'], dict):
            ts = d['runInfo']['runTime']
        else:
            ts = d['runInfo'].runTime
        times = datetime.datetime.strptime(ts,"%a %b %d %H:%M:%S %Y") 
        return times

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
        cprint('r', f"simulation path: {str(Index_data.simulation_path):s}")
        if Index_data.filetype == "F":
            runtime = self.get_sim_runtime(Path(Index_data.simulation_path, fn))
        elif Index_data.filetype == "D":
            # get the first tile in the directory
            ddir = list(Path(Index_data.simulation_path, indexdir).glob('*.p'))
            if len(ddir) > 0:
                runtime = self.get_sim_runtime(ddir[0])
            else:
                timestamp_str = "No Valid Data Files Found"
        if runtime is not None:
            timestamp_str = datetime.datetime.strftime(runtime, "%F-%T") # "%Y-%m-%d-%H:%M:%S")
        else:
            timestamp_str = "No Valid Data Files Found"
        Index_data.datestr = timestamp_str


        Index_data.threshold = None
        Index_data.elapsed = 0.0
        if indexdir is not None:
            runs = indexdir.glob("*.p")
            for r in runs:
                Index_data.files.append(str(r))  # get the result file list here
        else:
            Index_data.files = [fn]
        # print('Params: ' , params)
        # print('runinfo: ', runinfo)
        if usedict:  # old style files with dictionaries
            Index_data.cellType = params["cellType"]
            Index_data.modelType = params["modelType"]
            Index_data.modelName = params["modelName"]
            Index_data.runProtocol = params["runProtocol"]
            Index_data.synapsetype = params["ANSynapseType"]
            Index_data.synapseExperiment = params["spirou"]
            Index_data.dBspl = params["dB"]
            Index_data.SRType = params["SRType"]
            Index_data.soundType = params["soundtype"]
            Index_data.fmod = params["fmod"]
            Index_data.hocfile = params['hocfile']

            try:
                Index_data.ANSynapticDepression = params["ANSynapticDepression"]
            except:
                Index_data.ANSynapticDepression = 0
                
            try:
                Index_data.dt = params["dt"]
            except:
                Index_data.dt = 20e-6
            try:
                Index_data.temperature = params["Celsius"] #  Not always in the file information we retrieve
            except:
                Index_data.temperature = 0. # mark if missing
            try:
                Index_data.dataTable = params['dataTable']
            except:
                Index_data.dataTable = params["modelName"]
                cprint('r', 'Inserted dataTable with usedict')
            try:
                Index_data.pipDur = params['pip_duration']
            except:

                cprint('r', 'cant find pipdur')
                print(params)
                return

        else:  # new style files with dataclasses
            # print('runinfo: ', runinfo)
            Index_data.cellType = params.cellType
            Index_data.modelType = params.modelType
            Index_data.modelName = params.modelName
            Index_data.hocfile = params.hocfile
            Index_data.runProtocol = runinfo.runProtocol
            Index_data.synapsetype = params.ANSynapseType
            Index_data.ANSynapticDepression = params.ANSynapticDepression
            Index_data.synapseExperiment = runinfo.Spirou
            Index_data.dBspl = runinfo.dB
            Index_data.nReps = runinfo.nReps
            Index_data.pipDur = runinfo.pip_duration
            Index_data.SRType = params.SRType
            
            Index_data.soundType = runinfo.soundtype
            Index_data.fmod = runinfo.fmod
            Index_data.temperature = params.celsius
            Index_data.simulation_path = str(
                    Path(
                        self.basepath,
                        # f"VCN_c{int(self.selvals['Cells'][1]):02d}",
   #                      "Simulations",
   #                      self.selvals["Run Type"][1],
                    )
                )    
            try:
                Index_data.stimVC = runinfo.stimVC
            except:
                Index_data.stimVC = None
            try:
                Index_data.dataTable = params.dataTable
            except:
                Index_data.dataTable = params.modelName
                cprint('r', 'Inserted dataTable with dataclasses')
            try:
                Index_data.dendriteExpt = params.dendriteExpt
            except:
                Index_data.dendriteExpt = "Full"
            try:
                Index_data.axonExpt = params.axonExpt
            except:
                Index_data.axonExpt = "default"

            try:
                Index_data.dendriteMode = params.dendriteMode
            except:
                Index_data.dendriteMode = "normal"
            try:
                Index_data.dt = params.dt
            except:
                if runinfo.postMode == 'CC':
                    Index_data.dt = params.dtIC
                elif runinfo.postMode == 'VC':
                    Index_data.dt = params.dtVC
                else:
                    raise ValueError('Rininfo postmode not known: ', runinfo.postMode)
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
            indexdata = pickle.load(fh, fix_imports=True) # FPM.pickle_load(fh) # , encoding="latin1")
        return indexdata

    def remove_table_entry(self, indexrow):
        if len(indexrow) != 1:
            self.parent.error_message("Selection Error: Can only delete one row at a time")
            return
        indexrow = indexrow[0]
        data = self.get_table_data(indexrow)
        for f in data.files:
            fn = Path(data.simulation_path, f)
            print(f"Would delete: {str(fn):s}")
            print(fn.is_file())
            Path(fn).unlink()
        indexfilename = self.force_suffix(data.datestr)
        print(Path(indexfilename).is_file())
        if Path(indexfilename).is_file():
            print(f" and index file: {str(indexfilename):s}")
        # now update the table
        # print(indexrow)
        ind = self.get_table_data_index(indexrow)
        if ind is None:
            return
        # if ind is not None:
            # print(ind)
            # print(self.table_data[ind])
            # print(type(self.table_data[ind]))
        self.table.removeRow(ind)
        # print(dir(self.table))
        # print(self.table.viewport)
        # print(dir(self.table.viewport()))
        self.table.viewport().update()
        
        

    def print_indexfile(self, indexrow):
        """
        Print the values in the index file
        """
        print('='*80)
        print('\nIndex file and data file params')
        cprint('c', f"Index row: {str(indexrow.row):s}")
        data = self.get_table_data(indexrow)
        print('Data: ', data)
        for f in data.files:
            with open(f, "rb") as fh:
                fdata = pickle.load(fh, fix_imports=True) # FPM.pickle_load(fh) # , encoding="latin1")
                print('*'*80)
                cprint('c', f)
                cprint('y', 'File Data Keys')
                print(fdata.keys())
                # print(dir(fdata['Params']))
                cprint('y', 'Parameters')
                PP.pprint(fdata['Params'].__dict__) # , sort_dicts=True)
                print('-'*80)
                cprint('y', 'RunInfo')
                PP.pprint(fdata['runInfo'].__dict__) # , sort_dicts=True
                print('-'*80)
                cprint('y', 'ModelPars')
                PP.pprint(fdata['modelPars']) # , sort_dicts=True
                
                print('*'*80)
    
    @winprint_continuous
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
                # print('\nfile: ', str(f))
                # print('     index spriou: ', indxs[i].synapseExperiment)
                try:
                    x = indxs[i].dataTable
                except:
                    indxs[i]

        runfiles = list(thispath.glob("*.p"))
        dfiles = sorted(list(runfiles))  # regular data files
        valid_dfiles = []
        for i, f in enumerate(dfiles):
            p = self.read_pfile_params(f)
            if p is None:
                # indxs.append(None)
                self.textappend(f"None in #{i:d} :{str(f):s}", "yellow")
            else:
                params, runinfo, fn = p  # ok to unpack
                indxs.append(self.make_indexdata(params, runinfo, fn))
                valid_dfiles.append(f)  # keep track of accepted files
        dfiles = valid_dfiles  # replace with shorter list with only the valid files
        indexfiles = indexfiles + dfiles
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
                    indxs[i].dendriteExpt,
                    indxs[i].dendriteMode,
                    indxs[i].axonExpt,
                    indxs[i].synapseExperiment,
                    indxs[i].SRType,
                    indxs[i].ANSynapticDepression,
                    indxs[i].dBspl,
                    indxs[i].nReps, # command_line["nReps"],
                    indxs[i].pipDur,
                    indxs[i].soundType,
                    indxs[i].fmod,
                    len(indxs[i].files),
                    indxs[i].dataTable,
                    indxs[i].hocfile,
                    indxs[i].simulation_path,
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
                ("DendExpt", object),
                ("DendMode", object),
                ("AxonExpt", object),
                ("Experiment", object),
                # ("inputtype", object),
                # ("modetype", object),
                # ("scaling", object),
                # ("Freq", object),
                ('SRType', object),
                ("Depr", object),
                ("dBspl", object),
                ("nReps", object),
                ("pipDur", object),
                ("soundType", object),
                ("fmod", object),
                ("# Files", object),
                ("DataTable", object),
                ("HocFile", object),
                ('simulation_path', object),
                ("Filename", object),
                # ("# files", object),
                # ("synapsetype", object),
            ],
        )
        self.update_table(self.data)
        
    def update_table(self, data, QtCore=None):
        self.table.setData(data)
        style = "section:: {font-size: 4pt; color:black; font:TimesRoman;}"
        self.table.setStyleSheet(style)
        if QtCore is not None:
            # print('sorting by a column')
            self.table.sortByColumn(1, QtCore.Qt.AscendingOrder)
        self.altColors()  # reset the coloring for alternate lines
        # self.table.setStyle(QtGui.QFont('Arial', 6))
        self.table.resizeRowsToContents()
        self.table.resizeColumnsToContents()
        self.current_table_data = data
        self.parent.Dock_Table.raiseDock()

    def apply_filter(self,QtCore=None):
        """
        self.filters = {'Use Filter': False, 'dBspl': None, 'nReps': None, 'Protocol': None,
                'Experiment': None, 'modelName': None, 'dendMode': None, "dataTable": None,}
        """
        if not self.parent.filters['Use Filter']:  # no filter, so update to show whole table
            self.update_table(self.data)
        
        else:
            self.filter_table(self.parent.filters, QtCore)
        
    def filter_table(self, filters, QtCore):
            
            coldict = {'modelName': 4, 'Protocol': 5, 'dendMode': 6, 
                        'dendExpt':7, 'Experiment': 8, 'SRType': 9,
                        'Depr':10,'dBspl': 11, 'nReps': 12,
                         "pipDur": 13, "soundType": 14, "fmod": 15, "DataTable": 17,}
            filtered_table = self.data.copy()
            matchsets = dict([(x, set()) for x in filters.keys() if x != 'Use Filter'])
            for k, d in enumerate(self.data):
                for f, v in filters.items():
                    if (v is not None) and (coldict.get(f, False)) and (self.data[k][coldict[f]] == v):
                        #and f in list(coldict.keys())) and
                        # if (self.data[k][coldict[f]] == v):
                        # print("f: ", f, "   v: ", v)
                        matchsets[f].add(k)

            baseset = set()
            for k, v in matchsets.items():
                if len(v) > 0:
                    baseset = v
                    break
            # and logic:
            finds = [v for k, v in matchsets.items() if len(v) > 0]
            keep_index = baseset.intersection(*finds)
            self.keep_index = keep_index  # we might want this later!
            # print('Filter index: ', keep_index)
            filtered_table = [filtered_table[ft] for ft in keep_index]
            self.update_table(filtered_table, QtCore)

    def setColortoRow(self, rowIndex, color):
        for j in range(self.table.columnCount()):
            self.table.item(rowIndex, j).setBackground(color)

    def get_table_data_index(self, index_row):
        value = index_row.sibling(index_row.row(), 1).data()
        for i, d in enumerate(self.data):
            if self.data[i][1] == value:
                return i
        return None
        
    def get_table_data(self, index_row):
        """
        Regardless of the sort, read the current index row and
        map it back to the data in the table.
        """
        ind = self.get_table_data_index(index_row)
        if ind is not None:
            return self.table_data[ind]
        else:
            return None
            
        # # print('gettabledata: indexrow: ', index_row, index_row.row())
        # value = index_row.sibling(index_row.row(), 1).data()
        # for i, d in enumerate(self.data):
        #     if self.data[i][1] == value:
        #         return self.table_data[i]
        # return None
