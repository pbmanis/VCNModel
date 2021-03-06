"""
Make a master file list for a group of runs
This list is created after the runs are done, and in a subdirectory.

"""

from pathlib import Path
import pickle
import toml
import datetime
from dataclasses import dataclass, field
from typing import Union, Dict, List
import dataclasses

# @dataclass
# class Rundata:
#
#     rundata["command_line"] = args
#     rundata["changed_data"] = changed_data
#     rundata["cnmodel_hash"] = cnmodel_git_hash  # save hash for the model code
#     rundata["network_code_hash"] = git_head_hash  # this repository!
#     starttime = datetime.datetime.now()
#     runtime = starttime.strftime("%Y-%m-%d.%H-%M-%S")
#     rundata["datestr"] = runtime
#     rundata["cachepath"] = fullcachepath
#     rundata["files"] = []  # list of the files generated by this run (under the path)
#     rundata["temp"] = prot.network.temp
#     rundata["dt"] = prot.network.dt
#     rundata["threshold"] = args.spikethreshold
#     pp = pprint.PrettyPrinter(indent=4, width=160)
#     rundata["tables"] = pp.pformat(changed_data)  # save table used for this run.
#
#
    
    
class MasterFiles():
    def __init__(self):
        self.set_paths()
        self.get_dirs()
        self.get_info(self.expt_files[0])
        
    def set_paths(self, stimtype='AN', cell=11):
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path(
            "../VCN-SBEM-Data", "VCN_Cells")
        } 
        self.baseDirectory = self.datapaths["baseDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"
        print('basedir: ', self.baseDirectory)
        self.datapath = Path(self.baseDirectory, f"VCN_c{cell:02d}", self.simDirectory, stimtype)
        print('datapath: ', self.datapath)
        fs = self.datapath.glob('*')
        for f in list(fs):
            if f.is_dir():
                print(f)
        rundata = {}
        rundata['datestr'] = datetime.datetime.now()
        self.masterfn = Path(self.baseDirectory, stimtype, f"{rundata['datestr']:s}_master.pkl")
    
    def get_dirs(self):
        self.expts = []
        self.expt_files = []
        fs = self.datapath.glob('*')
        for f in list(fs):
            if f.is_dir():
                self.expts.append(f)
                self.expt_files = list(f.glob("*.p"))
        print(f.name, 'with ', len(self.expt_files), 'files')

    def get_file(self, f):
        pass

    def get_info(self, f):
        with open(f, "rb") as fh:
            runinfo = pickle.load(fh)
        print(runinfo.keys())
        print(runinfo['Params'])
        print(runinfo['modelPars'])
        print(runinfo['mode'])
    # if len(args.script) > 0:
    #     script_data = Path(args.script).read_text()
    #     rundata["script"] = {"name": args.script, "data": script_data}
    # else:
    #     rundata["script"] = {"name": None, "data": None}

    # for tr in list(results.keys()):
    #     f = results[tr][2]
    #     rundata["files"].append(f)  # get the result file list here
    #
    # mfpath = Path(masterfn).parent
    # mfpath.mkdir(parents=True, exist_ok=True)
    # print("fullcachepath: ", str(fullcachepath))
    # print("fullcachepath parent: ", str(Path(fullcachepath).parent))
    # print("mfpath: ", str(mfpath))
    # with open(masterfn, "wb") as fh:
    #     pickle.dump(rundata, fh)



if __name__ == '__main__':
    MF = MasterFiles()
    