import copy
import dataclasses
import datetime
import importlib
import json
from math import floor, log10
import multiprocessing as MPROC
import pickle
import re
import sys
import time
import timeit
from typing import Union
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import ephys.tools.Utility
import matplotlib
import numpy as np
import toml
from cnmodel import cells
from cnmodel.decorator import Decorator
from cnmodel.util import sound
from matplotlib import pyplot as mpl

# from neuronvis import hoc_graphics as hoc_graphics
from pylibrary.tools import cprint as CP

# from pylibrary.tools import utility as pu  # access to a spike finder routine
from pyqtgraph import multiprocess as MP

import vcnmodel.model_params
from vcnmodel import cell_config as cell_config
from vcnmodel import cellInitialization as cellInit
from vcnmodel.generate_run import GenerateRun

EPU = ephys.tools.Utility.Utility()
matplotlib.use("Qt5Agg")


def powerise10(x):
    """
    Returns x as a*10**b with 0 <= a < 10
    """
    if x == 0:
        return 0, 0
    Neg = x < 0
    if Neg:
        x = -x
    a = 1.0 * x / 10 ** (floor(log10(x)))
    b = int(floor(log10(x)))
    if Neg:
        a = -a
    return a, b


def eng(x):
    """
    Return a string representing x in an engineer-friendly notation
    """
    a, b = powerise10(x)
    if -3 < b < 3:
        return "%.4g" % x
    a = a * 10 ** (b % 3)
    b = b - b % 3
    return "%8.4gE%s" % (a, b)


showCell = True

cprint = CP.cprint


class ModelRun:
    """
    model_run2.py

    Run a model based on a hoc cell structure, decorating the structure with ion channels and synapses.

    Requires:

    1. Python 3.7, 3.8 or 3.9 (anaconda distribution or a local environment - preferred).
    2. Neuron7.7 - 8.0 (neuron.yale.edu)
    3. pyqtgraph (Luke Campagnola, from github)
    4. neuronvis (Campagnola/Manis, from github)
    5. cnmodel (Campagnola/Manis, from github)
    6. ephys (Manis, from github)
    7. pylibrary (Manis, from github)
    8. cochlea : Rudniki and Hemmert python implementation of Zilany model. A very slightly modified version that works with Python 3, "cochlea-1", is in Manis' github repository.
    9. thorns  # required for cochlea

    See the requirements.txt file, or preferably build an environment using make_local_env.sh.

    Where should we look for data? The paths to the directories are defined in **wheres_the_data.toml**
    
    Supported primarily by R01DC015901 (Spirou, Manis, Ellisman). Other development was supported by 
    R01 DC004551 (Manis, until 2019) and R01 DC019053 (Manis, 2020-2025)
    
    This program expects the following directory structure to hold simulation results
    and the morphology files::

        VCN_SBEM_Data/ # top level for
           VCN_Cells/   # top level for data
                For each cell there is a separate directory:
                VCN_c18/    # Cell ID for example (first argument in call should be this directory name)
                    MorphologyFiles/  # location for swc and hoc files
                        VCN_c18.hoc (cell body scaled version)
                        VCN_c18_unscaled.hoc  (unscaled version)
                        VCN_c18.swc  (various swc files that were translated to hoc files)
                    InitializationFiles/  # location for Init files
                        IV/  # initialization for just the IV with no synaptic input
                            VCN_c18....   # different base structures
                        AN/  # initialization for synaptic inputs
                            (ditto)  # different base structures of input arrangements
                    Simulations/
                        IV/  # results from IV simulations
                        AN/  # results from AN simulations

    The VCN_Cells directory should be placed in a folder outside of the program repository.
    The additional directory structure
    is defined in 'wheres_my_data.toml', which is used by cell_config.py to 
    set directory locations and read the excel file with ASA information.

    Usage (updated 11/10/2020)::

         model_run2.py cell
                        [-h]
                        [--type {Bushy,TStellate,DStellate}]
                        [--model {XM13,XM13_nacncoop,XM13A_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}]
                        [--modeltype {II,II-I,I-II,I-c,I-t,II-o}]
                        [--dendritemode {normal,passive,active,allpassive}]
                        [--hocfile HOCFILE]
                        [-D {default,Full,NoDend,NoDistal,NoUninnervated}]
                        [--datatable DATATABLE] [--sgcmodel {Zilany,cochlea}]
                        [--protocol {initIV,Zin,runIV,initandrunIV,initVC,runVC,initAN,runANPSTH,runANIO,runANSingles,gifnoise}]
                        [--style {cylinders,graph,volume,surface}]
                        [--mechanism DISPLAYMECHANISM]
                        [--displaymode {None,vm,sec-type,mechanism}]
                        [--displayscale] [--inputpattern INPUTPATTERN]
                        [--stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}]
                        [--check] [--testsetup] [-C CONFIGFILE] [-d DB] [-f F0]
                        [-r NREPS] [--seed SEED] [-S {LS,MS,HS,fromcell}]
                        [--synapsetype {simple,multisite}] [--depression {0,1}]
                        [--pip_start PIP_START] [--pip_duration PIP_DURATION]
                        [--pip_offduration PIP_OFFDURATION] [--fmod FMOD]
                        [--dmod DMOD] [--S2M SIGNALTOMASKER]
                        [--cmmrmode {CM,CD,REF}]
                        [--Spirou {all,max=mean,all=mean,removelargest,removetwolargest,largestonly,twolargest,threelargest,fourlargest}]
                        [--soma_inflate SOMA_INFLATION] [--soma_autoinflate]
                        [--dendrite_inflate DENDRITE_INFLATION]
                        [--dendrite_autoinflate] [--dendrite_from_soma]
                        [--ASA_from_soma] [--hold VSTIMHOLDING]
                        [--tagstring TAGSTRING] [-a AMPASCALE] [--allmodes]
                        [--sequence SEQUENCE] [--plot] [--workers NWORKERS]
                        [--noparallel] [--auto] [--saveall] [--verbose]
                        [--gifi GIF_I0] [--gifsigma GIF_SIGMA]
                        [--giffmod GIF_FMOD] [--giftau GIF_TAU]
                        [--gifdur GIF_DUR] [--gifskew GIF_SKEW]

    Optional arguments::
    
        -h, --help            show this help message and exit
        -T, --type {Bushy,TStellate,DStellate}
                              Define the cell type (default: Bushy)
         -M --model] {XM13,XM13_nacncoop,XM13A_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr},
                              Define the model type (default: XM13)
        --modeltype {II,II-I,I-II,I-c,I-t,II-o}
                              Define the model type (default: XM13)
        --dendritemode {normal,passive,active,allpassive}
                              Choose dendrite table (normal, active, passive)
        --hocfile HOCFILE     hoc file to use for simulation (default is the
                              selected "cell".hoc)
        -D {default,Full,NoDend,NoDistal,NoUninnervated}, --dendriteexpt {default,Full,FullSA,NoDend,AxonOnly,NoDistal,NoUninnervated}
                              Choose dendrite experiment (default, Full, NoDend,
                              NoDistal, NoUninnervated)
        --datatable DATATABLE
                              Specify the data table for this run
        --sgcmodel {Zilany,cochlea}
                              Define the SGC model type (default: Zilany)
        -P --protocol {initIV,Zin,runIV,initandrunIV,initVC,runVC,initAN,runANPSTH,runANIO,runANSingles,gifnoise},
                              Protocol to use for simulation (default: IV)
        --style {cylinders,graph,volume,surface}
                              Render cell with neuronvis with style:
                              {str(CmdChoices.displayStyleChoices):s}
        --mechanism DISPLAYMECHANISM
                              Render channel mechanisms: ['klt', 'kht', 'ihvcn',
                              'nacncoop', 'nacn', 'najsr']
        --displaymode {None,vm,sec-type,mechanism}
                              Render cell with neuronvis : set mode to one of
                              ['None', 'vm', 'sec-type', 'mechanism']
        --displayscale        use display scale and orientation from table for
                              generating renderings
        --inputpattern INPUTPATTERN, -i INPUTPATTERN
                              cell input pattern to use (substitute) from
                              cell_config.py
        -s --stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}
                              Define the stimulus type (default: tonepip)
        --check, -/           Only check command line for valid input; do not run
                              model
        --testsetup           Test all setup, but do not run simulations
        -C CONFIGFILE, --configfile CONFIGFILE
                              Read a formatted configuration file (JSON, TOML) for
                              commands
        -d DB, --dB DB        Set sound intensity dB SPL (default 30)
        -f F0, --frequency F0
                              Set tone frequency, Hz (default 4000)
        -r NREPS, --reps NREPS
                              # repetitions
        --seed SEED           AN starting seed
        -S {LS,MS,HS,fromcell}, --SRType {LS,MS,HS,fromcell}
                              Specify SR type (from: ['LS', 'MS', 'HS', 'fromcell'])
        --synapsetype {simple,multisite}
                              Specify AN synapse type (from: ['simple',
                              'multisite'])
        --depression {0,1}    Specify AN depression flag for multisite synapses
                              (from: [0, 1])
        --pip_start PIP_START
                              Set delay to onset of acoustic stimulus
        --pip_duration PIP_DURATION
                              Set duration of acoustic stimulus
        --pip_offduration PIP_OFFDURATION
                              Time to continue simulation AFTER sound ends
        --fmod FMOD           Set SAM modulation frequency
        --dmod DMOD           Set SAM modulation depth (in percent)
        --S2M SIGNALTOMASKER  Signal to Masker ratio (dB)
        --cmmrmode {CM,CD,REF}
                              Specify mode (from: ['CM', 'CD', 'REF'])
        --Spirou {all,max=mean,all=mean,removelargest,removetwolargest,largestonly,twolargest,threelargest,fourlargest}
                              Specify Spirou experiment type....
        --soma_inflate SOMA_INFLATION
                              Specify factor by which to inflate soma AREA
        --soma_autoinflate    Automatically inflate soma based on table
        --dendrite_inflate DENDRITE_INFLATION
                              Specify factor by which to inflate total dendritic
                              AREA
        --dendrite_autoinflate
                              Automatically inflate dendrite area based on table
        --dendrite_from_soma  Automatically inflate dendrite area based on soma
                              inflation
        --ASA_from_soma       Automatically inflate dendrite area based on soma
                              inflation
        --hold VSTIMHOLDING   Holding voltage in VClamp (mV) (default: -80 mV)
        --tagstring TAGSTRING
                              Add a tag string to the output filename to distinguish
                              it
        -a AMPASCALE, --AMPAScale AMPASCALE
                              Set AMPAR conductance scale factor (default 1.0)
        --allmodes            Force run of all modes (CMR, CMD, REF) for stimulus
                              configuration.
        --sequence SEQUENCE   Specify a sequence for the primary run parameters
        --plot                Plot results as they are generated - requires user
                              intervention...
        --workers NWORKERS    Number of "workers" for parallel processing (default:
                              4)
        --noparallel          Use parallel or not (default: True)
        --auto                Force auto initialization if reading the state fails
                              in initialization
        --saveall             Save data from all sections in model
        --verbose             Print out extra stuff for debugging
        --gifi GIF_I0         Set Noise for GIF current level (default 0 nA)
        --gifsigma GIF_SIGMA  Set Noise for GIF variance (default 0.2 nA)
        --giffmod GIF_FMOD    Set Noise for GIF fmod (default 0.2 Hz)
        --giftau GIF_TAU      Set Noise for GIF tau (default 0.3 ms)
        --gifdur GIF_DUR      Set Noise for GIF duration (default 10 s)
        --gifskew GIF_SKEW    Set Noise for GIF to have skewed distribution (0 =
                              normal)
    
    Example::

        Set up initialization:
        python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol initIV --model XM13
        Then run model:
        python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol runIV --model XM13


    This module is part of *vcnmodel*.

    Support::

        NIH grants:
        DC R01DC015901 (Spirou, Manis, Ellisman),
        DC R01 DC004551 (Manis, 2013-2019, Early development)
        DC R01 DC019053 (Manis, 2020-2025, Later development)

    """
    
    def __init__(self, params: dataclass = None, runinfo: dataclass = None, args=None):
        """
        The two dataclass parameters are defined in model_params.py, and must be passed here
        
        Parameters
        ----------
        params: a dataclass
            The parameter class that holds variables that control simulation
                parameters. Defined in model_params.py
        runinfo: dataclass
            The run information dataclass that holds variables governing the
                runs. Defined in model_params.py
        args: argparser object
            The command line arguments

        Returns
        -------
        Nothing
        """
        self.Params = params
        self.RunInfo = runinfo

        if self.Params.checkcommand:
            self.Params.commandline = sys.argv[1:]
            print(
                "Parameters: ", json.dumps(dataclasses.asdict(self.Params), indent=4)
            )  # pprint doesn't work well with ordered dicts
            print("RunInfo: ", dataclasses.asdict(self.RunInfo))
            print("Command line: ", self.Params.commandline)
            exit()

        self.Params.commandline = args
        self.Params.commands = sys.argv[1:]

        if self.Params.verbose:
            self.print_modelsetup()
        self.cconfig = cell_config.CellConfig(
            verbose=self.Params.verbose,
            spont_mapping=self.Params.SRType,
            add_inputs=self.RunInfo.SpirouSubs,  # input substitutions/additions
        )

        # find out where our files live
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {
                "cellDataDirectory": Path("../VCN-SBEM-Data", "VCN_Cells")
            }
        self.baseDirectory = self.datapaths["cellDataDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"

    def print_modelsetup(self):
        """
        Print out all of the parameters in the model
        """
        print("Params:\n", self.Params)
        print("RunInfo:\n", self.RunInfo)
        print("-----------")

    def set_celltype(self, cellType: str):
        """
        Set the cell Type, as requested. The cell type must be in the cell choices

        Parameters
        ----------
        cellType : string
            The type of the cell that will be the basis for the model

        Returns
        -------
            Nothing

        """
        if cellType not in self.cellChoices:
            choices = str(self.cellChoices)
            print(f"Celltype must be one of: {choices:s}. Got: {cellType:s}")
            exit()
        self.Params.cellType = cellType

    def set_model_name(self, model_name: str):
        """
        Set the model Type, as requested. The model type must be in the model choices

        Parameters
        ----------
        model_name : string
            The type of the model that will be used (condutance settings)

        Returns
        -------
            Nothing

        """
        if model_name not in self.modfelNameChoices:
            print(
                f"Model type must be one of: {', '.join(self.modelNameChoices):s}. Got: {self.Params.modelType:s}"
            )
            exit()
        self.Params.modelName = model_name

    def _make_filenames(self):
        """
        Define program-wide file names (used also in cell_
        initialization and generate_run) one time for consistencye

        This routine generates two names:
            1. The name of the initizlization file. This name

                 includes the model name, the model type (for that name),
                 soma and dendrite inflation factors if relevenat.
                 Other parameters can be added if needed
            2. The name of the simulation data file. This name
                 is similar to the initizlizaiton name, but may include
                 information about the type of run, stimuli, etc.

        """
        if not self.Params.setup:
            raise ValueError("Model setup call must be run before generating filenames")

        self.Params.cellID = Path(
            self.Params.cell
        ).stem  # os.path.splitext(self.Params.cell)[0]
        # pick up parameters that should be in both init and run filenames:
        # Run / init type independent parameters:
        if self.Params.verbose:
            print(
                "inflations: ",
                self.Params.soma_inflation,
                self.Params.dendrite_inflation,
            )

        model_Pars = f"{str(self.Params.modelName):s}"
        model_Pars += f"_{str(self.Params.modelType):s}"
        model_Pars += f"_HF={str(Path(self.Params.hocfile).stem):s}"
        if self.Params.soma_inflation != 1.0:
            model_Pars += f"_soma={self.Params.soma_inflation:.3f}"
        if self.Params.dendrite_inflation != 1.0:
            model_Pars += f"_dend={self.Params.dendrite_inflation:.3f}"
        if self.Params.ASA_inflation != 1.0:
            model_Pars += f"_ASA={self.Params.ASA_inflation:.3f}"
        # if self.Params.dendriteMode != 'normal':
        model_Pars += f"_{self.Params.dendriteMode:s}"

        add_Pars = None
        if self.RunInfo.Spirou == "all":
            add_Pars = "all"
        elif self.RunInfo.Spirou == "max=mean":
            add_Pars = "mean"
        elif self.RunInfo.Spirou == "all=mean":
            add_Pars = "allmean"
        elif self.RunInfo.Spirou == "removelargest":
            add_Pars = "removelargest"
        elif self.RunInfo.Spirou == "removetwolargest":
            add_Pars = "removetwolargest"
        elif self.RunInfo.Spirou == "largestonly":
            add_Pars = "largestonly"
        elif self.RunInfo.Spirou == "twolargest":
            add_Pars = "twolargest"
        elif self.RunInfo.Spirou == "threelargest":
            add_Pars = "threelargest"
        elif self.RunInfo.Spirou == "fourlargest":
            add_Pars = "fourlargest"
        if add_Pars is None:
            raise ValueError("RunInfo has invalid value")
        if self.Params.inputPattern is None:
            inputs = "self"
        else:
            inputs = self.Params.inputPattern
        datestr = self.run_starttime.strftime(
            "%Y-%m-%d.%H-%M-%S"
        )  # get the actual start time for the top directory
        run_directory = f"{self.RunInfo.runProtocol:s}-{add_Pars:s}-{datestr:s}"

        print("\nRUNPROTOCOL: ", self.RunInfo.runProtocol)
        if self.RunInfo.runProtocol in [
            "initIV",
            "initandrunIV",
            "runIV",
            "Zin",
            "runIVSpikeThreshold",
        ]:
            simMode = "IV"
            self.RunInfo.postMode = "CC"
            initPath = Path(self.baseDirectory, self.Params.cellID, self.initDirectory)
            self.mkdir_p(initPath)  # confirm existence of the initialization directory
            if self.Params.initStateFile is None:
                fn0 = f"{simMode:s}_Init_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
                fn0 += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
                fn0 += f"_{self.Params.SRType:2s}"
                if self.Params.tagstring is not None:
                    fn0 += f"_{self.Params.tagstring:s}"
                fn0 += ".p"
                self.Params.initStateFile = Path(initPath, fn0)
                if (
                    self.RunInfo.runProtocol.startswith("init")
                    and self.Params.initStateFile.is_file()
                ):
                    self.Params.initStateFile.unlink()  # delete old initializaiton file first
            simPath = Path(
                self.baseDirectory,
                self.Params.cellID,
                self.simDirectory,
                simMode,
                run_directory,
            )
            self.mkdir_p(simPath)  # confirm that output path exists
            fn = f"{simMode:s}_Result_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
            fn += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
            fn += f"_{self.Params.SRType:2s}"
            if self.Params.tagstring is not None:
                fn += f"_{self.Params.tagstring:s}"
            fn += ".p"
            self.Params.simulationFilename = Path(simPath, fn)

        elif self.RunInfo.runProtocol in ["initVC", "runVC"]:
            simMode = "VC"
            self.RunInfo.postMode = simMode
            initPath = Path(self.baseDirectory, self.Params.cellID, self.initDirectory)
            self.mkdir_p(initPath)  # confirm existence of the initialization directory
            print("initpath: ", initPath)
            if self.Params.initStateFile is None:
                fn0 = f"{simMode:s}_Init_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
                fn0 += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
                fn0 += f"_{self.Params.SRType:2s}"
                if self.Params.tagstring is not None:
                    fn0 += f"_{self.Params.tagstring:s}"
                fn0 += ".p"
                self.Params.initStateFile = Path(initPath, fn0)
                if (
                    self.RunInfo.runProtocol.startswith("init")
                    and self.Params.initStateFile.is_file()
                ):
                    self.Params.initStateFile.unlink()  # delete old initializaiton file first
            simPath = Path(
                self.baseDirectory,
                self.Params.cellID,
                self.simDirectory,
                simMode,
                run_directory,
            )
            self.mkdir_p(simPath)  # confirm that output path exists
            fn = f"{simMode:s}_Result_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
            fn += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
            fn += f"_{self.Params.SRType:2s}"
            if self.Params.tagstring is not None:
                fn += f"_{self.Params.tagstring:s}"
            fn += ".p"
            self.Params.simulationFilename = Path(simPath, fn)

        elif (
            self.RunInfo.runProtocol.startswith("runAN")
            or self.RunInfo.runProtocol == "initAN"
            or self.RunInfo.runProtocol == "runANSingles"
            or self.RunInfo.runProtocol == "runANPSTH"
            or self.RunInfo.runProtocol == "runANThreshold"
        ):
            simMode = "AN"
            self.RunInfo.postMode = "CC"

            initPath = Path(self.baseDirectory, self.Params.cellID, self.initDirectory)
            self.mkdir_p(initPath)  # confirm existence of that file
            if self.Params.initStateFile is None:
                fn0 = f"{simMode:s}_Init_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
                fn0 += f"_{add_Pars:s}_{self.Params.ANSynapseType:s}"
                if self.RunInfo.soundtype in ["SAM", "sam"]:
                    fn0 += f"_{self.RunInfo.fmod:03.1f}_{int(self.RunInfo.dmod):03d}"
                fn0 += f"_{self.Params.SRType:2s}.p"
                self.Params.initStateFile = Path(initPath, fn0)
                print("initstatefile exists?: ", self.Params.initStateFile.is_file())
                if (
                    self.RunInfo.runProtocol.startswith("initAN")
                    and self.Params.initStateFile.is_file()
                ):
                    self.Params.initStateFile.unlink()  # delete old initializaiton file first
            simPath = Path(
                self.baseDirectory,
                self.Params.cellID,
                self.simDirectory,
                simMode,
                run_directory,
            )
            fn = f"{simMode:s}_Result_{self.Params.cellID:s}_inp={inputs:s}_{model_Pars:s}"
            fn += f"_{add_Pars:s}_{self.Params.ANSynapseType:s}"
            fn += f"_{self.RunInfo.nReps:03d}_{self.RunInfo.soundtype:s}"
            fn += f"_{int(self.RunInfo.dB):03d}dB_{self.RunInfo.F0:06.1f}"
            if self.RunInfo.soundtype in ["SAM", "sam"]:
                fn += f"_{self.RunInfo.fmod:03.1f}_{int(self.RunInfo.dmod):03d}"
            fn += f"_{self.Params.SRType:2s}"
            if self.Params.tagstring is not None:
                fn += f"_{self.Params.tagstring:s}"
            fn += ".p"
            self.Params.simulationFilename = Path(simPath, fn)
            self.mkdir_p(simPath)  # confirm that output path exists
        cprint(
            "cyan",
            f"{simMode:s} Initialization file :  {str(self.Params.initStateFile.name):s}",
        )
        cprint(
            "c",
            f"{simMode:s} Simulation data file:  {str(self.Params.simulationFilename.name):s}",
        )
        cprint("c", f"Full sim path: {str(simPath):s}")

    def set_spontaneousrate(self, spont_rate_type: int):
        """
        Set the SR, overriding SR in the cell_config file. The SR type must be in the SR choices

        Parameters
        ----------
        spont_rate_type : int
            The SR type that will be used for AN fibers.
            1 = Low, 2 = medium, 3  = high (following Zilany et al)
        Returns
        -------
            Nothing

        """
        if self.spont_rate_type not in self.spont_rate_choices:
            print(
                "SR type must be one of: {:s}. Got: {:s} ".format(
                    ", ".join(self.spont_rate_choices), spont_rate_type
                )
            )
            exit()
        self.Params.SRType = spont_rate_type

    def set_starttime(self, starttime: float):
        """
        store the simulation start time...
        """
        self.Params.StartTime = starttime

    def mkdir_p(self, path: [str, Path]):
        """
        Make a set of directories along the path if needed to ensure that the
        data has a home

        Parameters
        ----------
        path : string or pathlib object

        Returns
        -------
        Nothing
        """

        try:
            Path.mkdir(path, parents=True, exist_ok=True)
        except IOError:
            raise FileNotFoundError(f"Cannot create path: {str(path):s}")

    def fix_singlets(self, h: object, harsh: bool = False):
        """
        Sometimes the hoc file has sections that consist of a single point.
        Here we replace those sections with ones that include the last point
        of the left connecting, or parent, section( assuming you go left to right).

        """
        badsecs = []
        for i, sec in enumerate(h.allsec()):
            if sec.n3d() == 1:
                badsecs.append(sec)
                parent = sec.parentseg()
                if parent is None:  # until we can find a better way.
                    continue
                nparentseg = parent.sec.n3d() - 1
                x = parent.sec.x3d(nparentseg)
                y = parent.sec.y3d(nparentseg)
                z = parent.sec.z3d(nparentseg)
                d = parent.sec.diam3d(nparentseg)
                h.pt3dinsert(0, x, y, z, d, sec=sec)
        if harsh and len(badsecs) > 0:
            cprint("r", f"**** Fixed Singlet segments in {len(badsecs):d} sections")
            for b in badsecs:
                print(f"Bad sections: ', {b.name():s}")
            exit()
        else:
            cprint("g", "**** No singlets were found in this hoc file")
        badsecs = []
        for sec in h.allsec():
            if sec.n3d() == 1:
                badsecs.append(sec)
        if len(badsecs) == 0:
            cprint("g", "No bad sections")
        else:
            cprint(
                "r", f"Found {len(badsecs):d} bad sections with sec.n3d() = 1; repaired"
            )

    def setup_model(self, par_map: dict = None):
        """
        Main entry routine for running all models.
        Here we assig various parameters, deal
        with various morphological manipulations such as inflating areas, generate
        filenames, and read the relevant data tables into cnmodel.
        Probably too much for just one function, but hey, it all has to happen.

        NOTE: This is the ONLY place where cells are instantiated.
        
        Parameters
        ----------
        par_map : dict (default: None)
            A dictionary of parameters, passed to models that are run (not used).

        Returns
        -------
        Nothing

        """
        if self.Params.cell is None:
            raise ValueError("Cell must be defined before setup is called")
        # dendrite_names = [
        #     "Proximal_Dendrite",
        #     "Distal_Dendrite",
        #     "Dendritic_Hub",
        #     "Dendritic_Swelling",
        #     "dend",
        #     "proximaldendrite",
        #     "distaldendrite",
        # ]

        if self.Params.verbose:
            print("run_model entry")
        if par_map is not None and "id" in list(par_map.keys()):
            self.idnum = par_map["id"]
        else:
            self.idnum = 9999
        print(f"Cell ID: {self.Params.cell:s}")

        self.Params.cellID = Path(
            self.Params.cell
        ).stem  # os.path.splitext(self.Params.cell)[0]
        print("self.Params.cellID: ", self.Params.cellID)
        if self.Params.verbose:
            print(f"Morphology directory: {self.morphDirectory:s}")

        if (
            self.Params.hocfile is not None
        ):  # -H flag: a passed hoc file will be used - overrides everything.
            label = "specified hocfile"

        else:  # build filename based on flags
            label = ""
            self.Params.hocfile = self.Params.cell
            if self.Params.dendriteExpt == "default":  # -dendriteExpt flags
                self.Params.hocfile += "_Full"
                label = self.Params.dendriteExpt
            if self.Params.dendriteExpt == "NoUninnervated":
                self.Params.hocfile += "_NoUninnervated"
                label = self.Params.dendriteExpt
            else:
                self.Params.hocfile += f"_{self.Params.dendriteExpt:s}"
                label += self.Params.dendriteExpt
            if not self.Params.meshInflate:  # not using raw? specify mesh inflated file
                self.Params.hocfile += "_MeshInflate"

            if self.Params.axonExpt == "standardized":
                self.Params.hocfile += "_standardized_axon"
            if self.Params.axonExpt.startswith("AIS"):
                self.Params.hocfile += "_" + self.Params.axonExpt
                label += ", " + self.Params.axonExpt

            self.Params.hocfile += ".hoc"  # now we can add extension

        cprint("c", f"Using {label:s} hoc file: {self.Params.hocfile:s}")

        hoc_filename = Path(
            self.baseDirectory,
            self.Params.cellID,
            self.morphDirectory,
            self.Params.hocfile,
        )

        cprint("c", f"Using hoc file at path: {str(hoc_filename):s}")
        # confirm hoc file exists
        if not hoc_filename.is_file():
            cprint("r", f"No Hoc file matching: {str(hoc_filename):s}")
            exit()  # fatal error
        else:
            with(open(hoc_filename, 'r')) as fh:
                self.Params.hocstring = fh.read() 
        
        # instantiate cells
        if self.Params.cellType in ["Bushy", "bushy"]:
            print("Creating a bushy cell (run_model) ")
            # We will use dynamic imports to get the table associated with the
            # selected channel mapping onto dendrites and sodium channel type
            from cnmodel import data

            dmodes = {
                "normal": "",
                "passive": "_pasdend",
                "active": "_actdend",
                "allpassive": "_allpassive",
            }
            changes = None
            nach = None  # uses default
            if self.Params.dataTable == "":
                table_name = f"vcnmodel.model_data.data_{self.Params.modelName:s}"
                table_name += f"{dmodes[self.Params.dendriteMode]:s}"
            else:
                table_name = f"vcnmodel.model_data.{self.Params.dataTable:s}"
                cprint("c", f"**** USING SPECIFIED DATA TABLE: {str(table_name):s}")
                knownmodes = ["normal", "actdend", "pasdend"]
                self.Params.dendriteMode = "normal"
                for mode in knownmodes:
                    if table_name.find(mode) > 0:
                        self.Params.dendriteMode = mode

                cprint("c", f"Dendrite mode: {self.Params.dendriteMode:s}")
            name_parts = self.Params.modelName.split("_")
            if len(name_parts) > 1:
                nach = name_parts[1]
            else:
                nach = "nav11"
            CHAN = importlib.import_module(table_name)
            """
            save the channel data and compartment data tables, as pulled from the 
            file we just imported, directly into the simulations
            This makes sure the file has even more (if not all) of the data needed to specify the
            run parameters
            """
            self.RunInfo.tableData["ChannelData"] = CHAN.ChannelData
            self.RunInfo.tableData["ChannelCompartments"] = CHAN.ChannelCompartments

            channels = f"{name_parts[0]:s}_{nach:s}_channels"
            compartments = f"{name_parts[0]:s}_{nach:s}_channels_compartments"
            # print("table name: ", table_name)
            # print(self.Params.dataTable)
            # print('nach: ', nach)
            # print("Channels: ", channels)
            # print("Compartments: ", compartments)
            # print(CHAN.ChannelData)
            # print(CHAN.ChannelCompartments)
            changes = data.add_table_data(
                channels,
                row_key="field",
                col_key="model_type",
                species="mouse",
                data=CHAN.ChannelData,
            )
            cprint("r", f"Channel Changes: {str(changes):s}")
            changes_c = data.add_table_data(
                compartments,
                row_key="parameter",
                col_key="compartment",
                species="mouse",
                model_type="II",
                data=CHAN.ChannelCompartments,
            )
            cprint("r", f"Compartment Changes: {str(changes_c):s}")
            if changes is not None:
                data.report_changes(changes)
                data.report_changes(changes_c)
            
            self.post_cell = cells.Bushy.create(
                morphology=str(hoc_filename),
                decorator=Decorator,
                species=self.Params.species,
                modelName=self.Params.modelName,
                modelType=self.Params.modelType,
                nach=nach,
            )

        elif self.Params.cellType in ["tstellate", "TStellate"]:
            print("Creating a t-stellate cell (run_model) ")
            self.post_cell = cells.TStellate.create(
                morphology=str(self.Params.hocfile),
                decorator=Decorator,
                species=self.Params.species,
                modelType=self.Params.modelName,
            )

        elif self.Params.cellType in ["dstellate", "DStellate"]:
            print("Creating a D-stellate cell (run_model)")
            self.post_cell = cells.DStellate.create(
                morphology=str(self.Params.hocfile),
                decorator=Decorator,
                species=self.Params.species,
                modelType=self.Params.modelName,
            )
        else:
            raise ValueError(f"cell type {self.Params.cellType:s} not implemented")

        # Set up run parameters
        print(
            f"Requested temperature (deg C): {self.post_cell.status['temperature']:.2f}"
        )
        self.post_cell.hr.h.celsius = self.post_cell.status[
            "temperature"
        ]  # this is set by prepareRun in generateRun. Only place it should be changed
        self.post_cell.hr.h.Ra = self.Params.Ra
        print(f"Specified Temperature = {self.post_cell.hr.h.celsius:8.1f} degC ")
        print("Ra (ohm.cm) = {:8.1f}".format(self.post_cell.hr.h.Ra))
        
        ###############################################
        # remove all non-parented sections, except the first.
        ###############################################
        for i, section in enumerate(self.post_cell.hr.h.allsec()):
            sref = self.post_cell.hr.h.SectionRef(sec=section)
            psec = self.post_cell.hr.h.parent_section(0, sec=section)
            if sref.has_parent():
                pass
                # print("i: ", i, "  sec: ", section.name(), "  parent sec: ", psec)
            else:
                if i > 0:
                    self.post_cell.hr.h.delete_section(sec=section)
                    print("Section index i: ", i, "  sec: ", section.name(), "has no parent")
        
        ###############################################
        # Fix singlet section (one point)
        ###############################################
        self.fix_singlets(self.post_cell.hr.h)

        if self.Params.ASA_fromsoma:
            self.Params.ASA_inflation = self.Params.soma_inflation
            if self.RunInfo.runProtocol not in ["initVC", "testVC", "runVC"]:
                rtau = self.post_cell.compute_rmrintau(
                    auto_initialize=True, vrange=[-80.0, -60.0]
                )
                print(f"    New Rin after somatic inflation: {rtau['Rin']:.2f}", end="")
                print(f" , tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}")

        # adjustments to areas are done with adjust_areas.py, and stored as a hoc file
        # no adjustments to areas are done on the fly here.
        # Set Ra for all regions, uniformly
        for group in list(self.post_cell.hr.sec_groups.keys()):
            g = self.post_cell.hr.sec_groups[group]
            for section in list(g):
                self.post_cell.hr.get_section(section).Ra = self.post_cell.hr.h.Ra
                if self.Params.verbose:
                    print("Section: ", section)
                    print("Ra: ", self.post_cell.hr.get_section(section).Ra)

        if 'soma' in self.post_cell.hr.sec_groups.keys() and len(self.post_cell.hr.sec_groups["soma"]) > 0:
            electrode_section = list(self.post_cell.hr.sec_groups["soma"])[0]
            secname = 'soma'
            self.RunInfo.electrodeSectionName = "soma"
        else:
            electrode_section = list(self.post_cell.hr.sec_groups["Axon_Initial_Segment"])[0]
            self.RunInfo.electrodeSectionName = "Axon_Initial_Segment"
        self.RunInfo.electrodeSection = self.post_cell.hr.get_section(electrode_section)
        # self.hg = hoc_graphics
        self.get_hoc_file(self.post_cell.hr)

        if self.Params.verbose:
            if par_map is not None:
                print("Listing par_map (run_model): ", par_map)
            self.post_cell.hr.h.topology()
        
        self.post_cell.set_nseg(freq=self.Params.lambdaFreq)

        self.Params.setup = True
        self._make_filenames()  # make filenames AFTER all manipulations to the cell

    def configure_cell(self, thisCell: object, synapseConfig: dict, celltype: str):
        """
        Configure the cell. This routine builds the cell in NEURON, adds presynaptic inputs
        as described in the synapseConfig, and configures those according to parameters in
        self.Params and self.RunInfo.
        It is used only for runs with auditory nerve input; basic IV/VC runs
        do not need to call this.

        Parameters
        ----------
        thisCell : hoc_reader object
            Access to neuron and file information

        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.

        celltype : string
            A string describing the cell type. Determines how the channels are populated
            and how the cell is built, even if it comes from a .hoc structure.

        Returns
        -------
            preCell : list
                A list of the preCell hoc objects attached to the synapses
            this_cell : cells object
                The target cell
            synapse : list
                A list of the synapses that were created and connected to this cell
            electrode_site : Neuron hoc object
                The section and location where the electrode is attached to the cell.
        """

        debug = False
        if debug:
            print("hf.sec_groups : ", list(thisCell.sec_groups.keys()))
            print(thisCell.print_all_mechs())
        preCell = []
        synapse = []
        # reconfigure syanpses to set the spont rate group(s) on a per-synapses basis
        for i, syn in enumerate(synapseConfig):
            if self.Params.SRType == "fromcell":  # use the one in the table
                preCell.append(cells.DummySGC(cf=self.RunInfo.F0, sr=syn["SR"]))
                preCell[-1]._synapsetype = self.Params.ANSynapseType
                print(
                    f"(configure_cell: fromcell): Retrieved SR index {syn['SR']:d} with SR type {self.Params.SRType:s}",
                    end=" ",
                )
            elif self.Params.SRType == "mixed1":  # use the one in the table
                # print('syn: ', i, syn["SR"])
                preCell.append(cells.DummySGC(cf=self.RunInfo.F0, sr=syn["SR"]))
                preCell[-1]._synapsetype = self.Params.ANSynapseType
                print(
                    f"(configure_cell: Mixed): Retrieved SR index {syn['SR']:d} with SR type {self.Params.SRType:s}",
                    end=" ",
                )
            else:
                try:
                    srindex = self.Params.srnames.index(self.Params.SRType)
                    print(
                        f"(configure_cell: NOT mixed): Retrieved SR index {srindex:d}", end=" "
                    )
                    print(f" with SR type {self.Params.SRType:s}")
                except ValueError:
                    raise ValueError(
                        "\nSR type '%s' not found in SR type list" % self.Params.SRType
                    )

                preCell.append(
                    cells.DummySGC(cf=self.RunInfo.F0, sr=srindex)
                )  # override

            if self.Params.verbose:
                print("  SRtype, srindex: ", self.Params.SRType, srindex)
            # print(self.Params.ANSynapseType)

            # note that we provide the opportunity to split the number of zones
            # between multiple sites
            for pl in syn["postlocations"]:
                postsite = syn["postlocations"][pl]
                plsecs = list(thisCell.hr.sec_groups[pl])
                # print('plsecs: ', plsecs)
                firstplsec = plsecs[0]  # get the first one (for soma)
                plsecn = re.split(r"[\[\]+]", firstplsec)
                # note this will need to be mapped later to put synapses on the right
                # sections in dendrites. But for now, this will have to do.
                postsite[0] = int(plsecn[1])  # from ['sections[118]'] for example
                # print('postsite: ', postsite)

                if self.Params.ANSynapseType == "simple":
                    print("  Synapsetype: simple", end=" ")
                    synapse.append(
                        preCell[-1].connect(
                            thisCell,
                            type="simple",
                            post_opts={"postsec": pl, "postsite": postsite[0:2]},
                        )
                    )
                else:
                    print(f" Synapsetype: multisite, with {int(syn['nSyn']):4d} zones", end=" ")
                    synapse.append(
                        preCell[-1].connect(
                            thisCell,
                            type="multisite",
                            pre_opts={
                                "nzones": int(syn["nSyn"] * postsite[2]),
                                "delay": syn["delay2"],
                            },
                            post_opts={"postsec": pl, "postsite": postsite[0:2]},
                        )
                    )

        if self.Params.ANSynapseType == "multisite":
            for i, s in enumerate(synapse):
                s.terminal.relsite.Dep_Flag = (
                    self.Params.ANSynapticDepression
                )  # turn on or off depression computation

        electrodeSection = list(thisCell.hr.sec_groups["soma"])[0]
        electrode_site = thisCell.hr.get_section(electrodeSection)
        return (preCell, synapse, electrode_site)

    def view_model(self, par_map: dict = None):
        """
        This method provides a sub connection to neuronvis that lets you
        view the hoc structure that has been requested
        """
        if not self.Params.setup:
            self.setup_model(par_map=par_map)
        return  # until we figure out opengl on mac osx Big Sur
        # import neuronvis as NV
        #
        # if self.Params.cellID in list(vcnmodel.model_params.display_orient_cells):
        #     angles = vcnmodel.model_params.display_orient_cells[self.Params.cellID]
        # else:
        #     angles = [200.0, 0.0, 0.0]
        # print("Rendering via model_run2")
        # NV.hocRender.Render(
        #     hoc_file=self.post_cell.hr.h,  # hoc_file,
        #     display_style=self.Params.displayStyle,  # display_mode,
        #     display_renderer="pyqtgraph",  # args["renderer"],
        #     display_mode=self.Params.displayMode,  # 'sec-type', # args["display"],
        #     mechanism=self.Params.displayMechanism,  # args["mechanism"],
        #     initial_view=angles,
        #     alpha=1.0,  # args["alpha"],
        #     sim_data=None,  # sim_data,
        # )
        # exit()

    def run_model(self, par_map: dict = None):
        """
        Main routine for dispatching the control to various
        run routines.

        """
        self.run_starttime = datetime.datetime.now()
        if not self.Params.setup:
            self.setup_model(par_map=par_map)

        dispatcher = {
            "initIV": self.initIV,
            "runIV": self.iv_run,
            "Zin": self.Zin,
            "runIVSpikeThreshold": self.iv_run_spike_threshold,
            "initVC": self.initVC,
            "runVC": self.VC_run,
            "initAN": self.an_run_init,  # (self.post_cell, make_an_intial_conditions=True),
            "runANPSTH": self.an_run,
            "runANIO": self.an_run_IO_gSyn,
            "runANSingles": self.an_run_singles,
            "runANOmitOne": self.an_run_omit_one,  # .(self.post_cell, exclude=True),
            "runANThreshold": self.an_run_find_syn_thr,
            "gifnoise": self.noise_run,
        }
        if self.RunInfo.runProtocol in [
            "runANPSTH",
            "runANSingles",
            "runANThreshold",
            "runANOmitOne",
            "runANIO",
        ]:
            self.RunInfo.run_duration = (
                np.sum(self.RunInfo.pip_start)
                + self.RunInfo.pip_duration
                + self.RunInfo.pip_offduration
            )

        dispatcher[self.RunInfo.runProtocol]()

    # ===============================================================================
    # The following methods are all of the initialization and run routines for the
    # model.
    # ===============================================================================

    def initIV(self):
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        cellInit.get_initial_condition_state(
            self.post_cell,
            mode=self.RunInfo.postMode,
            tdur=self.Params.initialization_time,
            filename=self.Params.initStateFile,
            electrode_site=self.R.electrode_site,
        )
        print(f"Ran to get initial state for {self.post_cell.hr.h.t:.1f} msec")

    def Zin(self, measure_section:Union[str, None]=None):
        """
        Measure Zin from the soma or specific compartment for a compartmental cell.

        Parameters
        ----------
        measure_section: object, default: None
             is the section (neuron hoc) where the
            measurement will be made. if None, will make
            measurement from the soma.
            The section should be named by "string[#][position]"
            e.g., section[9][0.5]
        """
        if measure_section is None:
            measure_section = self.post_cell.soma

        self.Params.initialization_time = 500.0
        self.initIV()
        # print(dir(self.post_cell.hr.h))
        Z = self.post_cell.hr.h.Impedance(0.5, measure_section)
        freqs = np.logspace(0, 4, 50)
        Zin = np.zeros_like(freqs)
        Zphase = np.zeros_like(freqs)
        for i, freq in enumerate(freqs):
            Z.compute(freq, 1)
            Zin[i] = Z.input(0.5, measure_section)
            Zphase[i] = Z.input_phase(0.5, measure_section)
            cprint("y", f"f:{freq:.3f} Vm: { self.post_cell.soma.v:.3f}")
            # print(f"Fr: {freq:.1f}  Zin: {Zin[i]:8.2f}, Phase: {Zphase[i]:6.2f}")
        fo = Path(self.Params.hocfile).stem
        fo += f"_{self.Params.dendriteMode:s}"
        fo += "_Z.pkl"
        fo = Path(self.baseDirectory, "Impedance_Calculations", fo)
        d = {"f": freqs, "zin": Zin, "phase": Zphase, "Vm": measure_section.v}
        with open(fo, "wb") as fh:
            pickle.dump(d, fh)

        if self.Params.plotFlag:
            import pyqtgraph as pg

            pg.mkQApp()
            win = pg.GraphicsWindow()
            p1 = win.addPlot()
            p2 = win.addPlot()
            p1.setLogMode(x=True, y=False)
            p2.setLogMode(x=True, y=False)
            p1.plot(freqs, Zin, pen="k", symbol="o", symbolsize=1,)
            p2.plot(freqs, Zphase, pen="b", symbol="s", symbolsize=1)
            pg.QtGui.QApplication.instance().exec_()

        exit()

    def iv_run(self, par_map: dict = None):
        """
        Main entry routine for running all IV (current-voltage relationships with somatic electrode)

        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).

        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times

        """
        print("iv_run: starting")
        start_time = timeit.default_timer()

        if self.RunInfo.sequence != "":  # replace sequence?
            self.RunInfo.stimInj = {"pulse": eval(self.Params.sequence)}
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        self.R.RunInfo.folder = Path(
            self.baseDirectory, self.Params.cellID, self.simDirectory, "IV"
        )
        if self.Params.verbose:
            print("iv_run: calling do_run")
        nworkers = self.Params.nWorkers
        #        print(self.Params.Parallel)
        if self.Params.Parallel is False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(
            self.Params.hocfile,
            parMap=self.RunInfo.stimInj,
            save="monitor",
            restore_from_file=True,
            initfile=self.Params.initStateFile,
            workers=nworkers,
        )
        if self.Params.verbose:
            print("   iv_run: do_run completed")
        elapsed = timeit.default_timer() - start_time
        print(f"   iv_rin: Elapsed time: {elapsed:2f} seconds")
        isteps = self.R.IVResult["I"]
        if self.Params.verbose:
            for k, i in enumerate(self.R.IVResult["tauih"].keys()):
                print(
                    "   ih: %3d (%6.1fnA) tau: %f"
                    % (i, isteps[k], self.R.IVResult["tauih"][i]["tau"])
                )
                print("           dV : %f" % self.R.IVResult["tauih"][i]["a"])
            for k, i in enumerate(self.R.IVResult["taus"].keys()):
                print(
                    "   i: %3d (%6.1fnA) tau: %f"
                    % (i, isteps[k], self.R.IVResult["taus"][i]["tau"])
                )
                print("          dV : %f" % (self.R.IVResult["taus"][i]["a"]))

        # print('   Nspike, Ispike: ', self.R.IVResult['Nspike'], self.R.IVResult['Ispike'])
        print("   N spikes:   {0:d}".format(int(np.sum(self.R.IVResult["Nspike"]))))
        print("   Rinss:      {0:.1f} Mohm".format(self.R.IVResult["Rinss"]))
        print(
            "   Tau(mean):  {0:.3f} ms".format(
                np.mean(
                    [
                        self.R.IVResult["taus"][i]["tau"]
                        for i in range(len(self.R.IVResult["taus"]))
                    ]
                )
            )
        )
        print("   Vm:         {0:.1f} mV".format(np.mean(self.R.IVResult["Vm"])))
        if len(list(self.R.IVResult["taus"].keys())) == 0:
            taum_mean = 0.0
            tauih_mean = 0.0
        else:
            taum_mean = np.mean(
                [
                    self.R.IVResult["taus"][i]["tau"]
                    for k, i in enumerate(self.R.IVResult["taus"].keys())
                ]
            )
            tauih_mean = np.mean(
                [
                    self.R.IVResult["tauih"][i]["tau"]
                    for k, i in enumerate(self.R.IVResult["tauih"].keys())
                ]
            )
        # construct dictionary for return results:
        self.IVSummary = {
            "basefile": self.R.basename,
            "par_map": par_map,
            "ID": self.idnum,
            "sequence": self.RunInfo.stimInj,
            "Vm": np.mean(self.R.IVResult["Vm"]),
            "Rin": self.R.IVResult["Rinss"],
            "taum": taum_mean,
            "tauih": tauih_mean,
            "spikes": {"i": self.R.IVResult["Ispike"], "n": self.R.IVResult["Nspike"]},
        }
        return self.IVSummary

    def iv_run_spike_threshold(self, par_map: dict = None):
        """
        find spike threshold - lowest current level (within a delta I) that generates a
        spike.

        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).

        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times

        """
        print("iv_run_spike_threshold: starting")
        # start_time = timeit.default_timer()

        # if self.RunInfo.sequence != "":  # replace sequence?
        #      self.RunInfo.stimInj = {"pulse": np.linspace(0.01, 1.60, endpoint=True)}
        #  self.R = GenerateRun(
        #      self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        #  )
        self.RunInfo.folder = Path(
            self.baseDirectory, self.Params.cellID, self.simDirectory, "IV"
        )
        if self.Params.verbose:
            print("iv_run_spike_threshold: calling do_run")
        nworkers = self.Params.nWorkers
        #        print(self.Params.Parallel)
        if self.Params.Parallel is False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
        # coarse run first:
        self.RunInfo.stimDur = 20.0  # msec short pulses are sufficient

        current_I = 0.4  # start near middle of the range
        precision = 0.001  # measure to nearest 10 pA
        step = current_I / 2.0
        nworkers = 1
        niter = 0
        last_spike_thr = current_I
        while step > precision:
            niter += 1
            self.RunInfo.stimInj = {
                "pulse": [current_I]
            }
            self.R = GenerateRun(
                self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
            )
            self.R.doRun(
                self.Params.hocfile,
                parMap=self.RunInfo.stimInj,
                save="monitor",
                restore_from_file=True,
                initfile=self.Params.initStateFile,
                workers=nworkers,
            )
            fspk = np.where(self.R.IVResult["Nspike"])[0]
            if len(fspk) == 0: # no spikes, increase current
                current_I += step
            else:
                current_I -= step
                last_spike_thr = current_I
                step /= 2.0  # only decrease step size when spike detected

        thrstr = (
            f"{datetime.datetime.now().strftime('%m.%d.%Y %H:%M:%S'):s} Prec. {precision:.3f} pA, " +
            f" {last_spike_thr:.3f} niter={niter:3d} {str(self.Params.hocfile):s}"
        )
        cprint("m", thrstr)
        with (open("thrrun.txt", "a")) as fh:
            fh.write(thrstr + "\n")
        return

    def initVC(self):
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        cellInit.get_initial_condition_state(
            self.post_cell,
            mode=self.RunInfo.postMode,
            tdur=self.Params.initialization_time,
            filename=self.Params.initStateFile,
            electrode_site=self.R.electrode_site,
        )
        print(f"VC Ran to get initial state for {self.post_cell.hr.h.t:.1f} msec")

    def testVC(self):
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        cellInit.test_initial_conditions(
            self.post_cell,
            filename=self.Params.initStateFile,
            electrode_site=self.R.electrode_site,
        )
        self.R.testRun(initfile=self.Params.initStateFile)
        return  # that is ALL, never make testVC/init and then keep running.

    def VC_run(self, par_map: dict = None):
        """
        Main entry routine for running all IV (current-voltage relationships with somatic electrode)

        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).

        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times

        """
        print("VC_run: starting")
        start_time = timeit.default_timer()

        if self.RunInfo.sequence != "":  # replace sequence?
            self.RunInfo.stimInj = {"pulse": eval(self.Params.sequence)}
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        self.R.RunInfo.folder = Path(
            self.baseDirectory, self.Params.cellID, self.simDirectory, "VC"
        )
        if self.Params.verbose:
            print("VC_run: calling do_run")
        nworkers = self.Params.nWorkers
        #        print(self.Params.Parallel)
        if self.Params.Parallel is False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(
            self.Params.hocfile,
            parMap=self.RunInfo.stimInj,
            save="monitor",
            restore_from_file=True,
            initfile=self.Params.initStateFile,
            workers=nworkers,
        )
        if self.Params.verbose:
            print("   VC_run: do_run completed")
        elapsed = timeit.default_timer() - start_time
        print(f"   VC_rin: Elapsed time: {elapsed:2f} seconds")
        # isteps = self.R.VCResult["I"]
        # if self.Params.verbose:
        #     for k, i in enumerate(self.R.IVResult["tauih"].keys()):
        #         print(
        #             "   ih: %3d (%6.1fnA) tau: %f"
        #             % (i, isteps[k], self.R.IVResult["tauih"][i]["tau"])
        #         )
        #         print("           dV : %f" % self.R.IVResult["tauih"][i]["a"])
        #     for k, i in enumerate(self.R.IVResult["taus"].keys()):
        #         print(
        #             "   i: %3d (%6.1fnA) tau: %f"
        #             % (i, isteps[k], self.R.IVResult["taus"][i]["tau"])
        #         )
        #         print("          dV : %f" % (self.R.IVResult["taus"][i]["a"]))
        #
        # # print('   Nspike, Ispike: ', self.R.IVResult['Nspike'], self.R.IVResult['Ispike'])
        # print("   N spikes:   {0:d}".format(int(np.sum(self.R.IVResult["Nspike"]))))
        # print("   Rinss:      {0:.1f} Mohm".format(self.R.IVResult["Rinss"]))
        # print(
        #     "   Tau(mean):  {0:.3f} ms".format(
        #         np.mean(
        #             [
        #                 self.R.IVResult["taus"][i]["tau"]
        #                 for i in range(len(self.R.IVResult["taus"]))
        #             ]
        #         )
        #     )
        # )
        # print("   Vm:         {0:.1f} mV".format(np.mean(self.R.IVResult["Vm"])))
        # if len(list(self.R.IVResult["taus"].keys())) == 0:
        #      taum_mean = 0.0
        #      tauih_mean = 0.0
        #  else:
        #      taum_mean = np.mean(
        #          [
        #              self.R.IVResult["taus"][i]["tau"]
        #              for k, i in enumerate(self.R.IVResult["taus"].keys())
        #          ]
        #      )
        #      tauih_mean = np.mean(
        #          [
        #              self.R.IVResult["tauih"][i]["tau"]
        #              for k, i in enumerate(self.R.IVResult["tauih"].keys())
        #          ]
        #      )
        # construct dictionary for return results:
        self.VCSummary = {
            "basefile": self.R.basename,
            "par_map": par_map,
            "ID": self.idnum,
            "sequence": self.RunInfo.stimInj,
            "Vm": np.mean(self.R.VCResult["Vm"]),
            # "Rin": self.R.VCResult["Rinss"],
            # "taum": taum_mean,
            # "tauih": tauih_mean,
        }
        return self.VCSummary

    def noise_run(self, par_map: dict = {}):
        """
        Main entry routine for stimulating a cell with noise
        via current injection (may be useful for generating GIF models)

        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).

        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times

        """
        print("noise_run: starting")
        # parse i_test_range and pass it here
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        # ivinitfile = Path(self.baseDirectory, self.Params.cellID,
        #                         self.initDirectory, self.Params.initStateFile)
        self.R.runInfo.folder = Path(
            self.baseDirectory, self.Params.cellID, self.simDirectory, "Noise"
        )
        if self.Params.verbose:
            print("noise_run: calling do_run")
        nworkers = self.Params.nWorkers

        if self.Params.Parallel is False:
            nworkers = 1
        self.R.doRun(
            self.Params.hocfile,
            parMap=par_map,
            save="monitor",
            restore_from_file=True,
            initfile=self.Params.initStateFile,
            workers=nworkers,
        )
        if self.Params.verbose:
            print("   noise_run: do_run completed")

    def check_for_an_statefile(self):
        print("State file name: ", self.Params.initStateFile)
        print("Cell: ", self.Params.cell)
        print("Base Directory: ", self.baseDirectory)
        print("Initialization Directory: ", self.initDirectory)
        statefile = Path(
            self.baseDirectory,
            self.Params.cell,
            self.initDirectory,
            self.Params.initStateFile,
        )
        return statefile.is_file()

    def compute_seeds(self, nReps: int, synapseConfig: list):
        """
        Generate the seeds for the AN runs

        If Params['seed'] is None, we use the randomized seed method - every run is totally
        independent.
        If not, then we generate a consequetive sequence of seeds (so it is controlled and
        reusable) offset by the starting number specified as the seed
        
        Note: If you change the number of synapses, then the seeds for subsequent repetitions
        will be different. Probably should have made the seed array > # of synapses ever expected.
        
        Parameters
        ----------
        nReps: int (no default)
            Number of repetions of the stimulus that are needed

        synpaseConfig : list (no default)
            The length of the synpaseConfig list indicates how many synapses are on the
            cell, and each synapse is assigned a unique seed.

        Returns
        -------

        s : numpy array
            A 2-D numpy array (nReps, len(synapseConfig)) containing the seeds for each run and
            each synapse.
        """
        if self.Params.seed is None:  # every run is randomized
            seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        else:
            startingseed = self.Params.seed
            seeds = np.arange(0, nReps * len(synapseConfig)).reshape(
                (nReps, len(synapseConfig))
            )
            seeds = seeds + startingseed
        # print('AN Seeds: ', seeds)

        return seeds

    def print_timing(self):
        total_elapsed_time = time.time() - self.start_time
        #        total_run_time = time.time() - run_time
        print(
            f"{'Total Elapsed Time = ':>25s} {(total_elapsed_time / 60.0):8.2f} min ({total_elapsed_time:8.0f}s)"
        )
        print(
            f"{'Total Setup Time = ':>25s} {self.setup_time / 60.0:8.2f} min ({self.setup_time:8.0f}s)"
        )
        print(
            f"{'Total AN Calculation Time = ':>25s} {self.an_setup_time / 60.0:8.2f} min ({self.an_setup_time:8.0f}s)"
        )
        print(
            f"{'Total Neuron Run Time = ':>25s} {self.nrn_run_time / 60.0:8.2f} min ({self.nrn_run_time:8.0f}s)"
        )

    def get_synapses(self, synapses: list, printvalues=False):
        """
        Get values of settings in synapses
        """
        gMax = np.zeros(len(synapses))  # total g for each synapse
        nSyn = np.zeros(len(synapses))  # # sites each synapse
        ngMax = np.zeros(len(synapses))
        # print('# syn: ', len(synapses))
        if self.Params.ANSynapseType == "simple":
            for i, s in enumerate(synapses):
                gMax[i] = gMax[i] + float(s.psd.terminal.netcon.weight[0])
                nSyn[i] = nSyn[i] + 1

        else:
            for i, s in enumerate(synapses):
                nSyn[i] = len(s.psd.ampa_psd)
                for p in s.psd.ampa_psd:
                    gMax[i] = gMax[i] + p.gmax
                for p in s.psd.nmda_psd:
                    ngMax[i] = ngMax[i] + p.gmax
        if self.Params.verbose or printvalues:
            print("  getsyn")
            print("  Syn#    nsites    AMPA gmax    NMDA gmax   synperum2    SRType")
            for i, s in enumerate(synapses):
                # print(self.Params.SynapseConfig[i])
                print(
                    f"  {i:>4d}   {int(nSyn[i]):>5d}    {eng(gMax[i]):>9s}    {eng(ngMax[i]):>9s}",
                    end="",
                )
                print(
                    f"  {self.Params.SynapseConfig[i]['synperum2']}, {self.Params.SynapseConfig[i]['SR']:d}"
                )
        return (gMax, ngMax, nSyn)

    def set_synapse(self, synapse: list, gampa: float, gnmda: float):
        totalg = 0.0
        if self.Params.ANSynapseType == "simple":
            synapse.psd.terminal.netcon.weight[0] = gampa
            return
        # multisite:
        for p in synapse.psd.ampa_psd:
            p.gmax = gampa
            totalg = totalg + p.gmax
        totaln = 0.0
        for p in synapse.psd.nmda_psd:
            p.gmax = gnmda
            totaln = totaln + p.gmax
        if self.Params.verbose:
            print("setsyn: total ampa, nmda: ", totalg, totaln)

    def retrieve_data(self, tresults):
        """
        Do a bit of massaging from the data returned in tresults
        to build the result array and ancillary arrays for storage
        """
        celltime = tresults["time"]  # (self.time)
        spikeTimes = EPU.findspikes(
            x=tresults["time"] / 1000.0,  # in msec, convert to sec
            v=tresults["Vsoma"] * 1e-3,  # in mV
            thresh=self.RunInfo.threshold * 1e-3,  # in units of vsoma, so mV
            t0=0.0,  # / sec
            t1=self.RunInfo.run_duration,  # leave in sec
            dt=self.post_cell.hr.h.dt * 1e-3,  # convert to sec
            mode="peak",
            detector="Kalluri",
            refract=0.0007,  # min refractory period, msec
        )
        allDendriteVoltages = {}
        inputSpikeTimes = tresults[
            "ANSpikeTimes"
        ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
        somaVoltage = np.array(tresults["Vsoma"])
        dendriteVoltage = np.array(tresults["Vdend"])
        stimWaveform = np.array(tresults["stimWaveform"])
        stimTimebase = np.array(tresults["stimTimebase"])
        # print('input spikes: ', inputSpikeTimes)
        result = {
            "spikeTimes": spikeTimes,
            "inputSpikeTimes": inputSpikeTimes,
            "time": np.array(celltime),
            "somaVoltage": somaVoltage,
            "dendriteVoltage": dendriteVoltage,  # just for one select location
            "allDendriteVoltages": allDendriteVoltages,
            "stimWaveform": stimWaveform,
            "stimTimebase": stimTimebase,
        }
        if self.Params.save_all_sections:  # just save soma sections
            for sd in tresults["allsecVec"]:
                print(str(sd))
                result["allDendriteVoltages"][str(sd)] = np.array(tresults["allsecVec"][sd])
                # tresults["allsecVec"]
            # )  # convert out of neuron
        return celltime, result

    def an_run_init(self):
        self.an_run(make_an_intial_conditions=True)

    def an_run(
        self, verify: bool = False, make_an_intial_conditions: bool = False,
    ):
        """
        Establish AN inputs to soma, and run the model.
        Requires a synapseConfig list of dicts from cell_config.make_dict()
        each list element represents an AN fiber (SGC cell) with:
            * (N sites, delay (ms)
            * spont rate group [1=low, 2=high, 3=high],
            * type (e or i)
            * segment and segment location
        
        Note if segment is None, then the synapse is assigned to the soma.

        Parameters
        ----------
        verify : boolean (default: False)
            Flag to control printing of various intermediate results

        make_ANIntitialConditions : bool (default : False)
            Flag to control whether the initial conditions need to be recomputed and stored.

        Returns
        -------
            Nothing

        """
        if make_an_intial_conditions:
            cprint("g", f"\n*** Initializing for an_run, time={self.Params.initialization_time:.5f}\n")
        else:
            cprint("g", "\n*** an_run\n")
        if self.Params.inputPattern is not None:
            fromdict = self.Params.inputPattern
            print(f"Cell id: {self.Params.cellID:s} using input pattern: {fromdict:s}")
        else:
            print(f"Cell id: {self.Params.cellID:s} is using 'self' input pattern")
            fromdict = self.Params.cellID

        synapseConfig, celltype = self.cconfig.make_dict(
            fromdict, areainflate=self.Params.ASA_inflation,
        )

        self.Params.SynapseConfig = synapseConfig
        # print("SynapseConfig: ", self.Params.SynapseConfig)

        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        # for i, s in enumerate(synapseConfig):
        #     s[1] = s[3]*0.001/3.0
        #     print 'delay for input %d is %8.4f msec' % (i, s[1])

        nReps = self.RunInfo.nReps

        preCell, synapse, self.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig=self.Params.SynapseConfig, celltype=celltype
        )

        # see if we need to save the cell state now.
        if make_an_intial_conditions:
            cellInit.get_initial_condition_state(
                self.post_cell,
                mode=self.RunInfo.postMode,
                tdur=self.Params.initialization_time,
                filename=self.Params.initStateFile,
                electrode_site=self.electrode_site,
                reinit=self.Params.auto_initialize,
            )
            cprint(
                "c",
                f"Confirming Initial conditions with: {str(self.Params.initStateFile.name):s}",
            )
            cellInit.test_initial_conditions(
                self.post_cell,
                filename=self.Params.initStateFile,
                electrode_site=self.electrode_site,
            )
            return

        seeds = self.compute_seeds(nReps, synapseConfig)
        self.RunInfo.seeds = seeds  # keep the seed values too.

        celltime = [None] * nReps
        # allDendriteVoltages = (
        #     {}
        # )  # Saving of all dendrite voltages is controlled by --saveall flag
        result = {}
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0

        nWorkers = self.Params.nWorkers
        TASKS = [s for s in range(nReps)]
        tresults = [None] * len(TASKS)

        # manipulate syanptic strengths to test hypotheses here.
        allsf = self.Params.AMPAScale
        if self.Params.verbose:
            print("AMPA Scale: ", allsf)
        gMax, ngMax, nSyn = self.get_synapses(synapse)

        if allsf != 1.0:
            for i in range(len(nSyn)):
                gMax[i] = gMax[i] * allsf
                ngMax[i] = ngMax[i] * allsf
                self.set_synapse(synapse[i], gMax[i] / nSyn[i], ngMax[i] / nSyn[i])

        if self.RunInfo.Spirou != "all":
            print(f"   Spirou test type: {self.RunInfo.Spirou:s}")

        if self.RunInfo.Spirou in ["all=mean"]:
            cprint("c", "    Setting ALL to the mean of all inputs (no variance)")
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print("    AMPA: mean, gmax, imax: ", meangmax, gMax, imaxgmax)
            print("    NMDA: mean, gnmax: ", meangnmax, gnMax)
            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            for i in range(len(synapse)):
                self.set_synapse(synapse[i], meangmax / nSyn[i], meangnmax / nSyn[i])

        elif self.RunInfo.Spirou in ["max=mean"]:
            cprint("c", "    Setting largest to mean of all inputs")
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print("    AMPA: mean, gmax, imax: ", meangmax, gMax, imaxgmax)
            print("    NMDA: mean, gnmax: ", meangnmax, gnMax)
            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            self.set_synapse(
                synapse[imaxgmax], meangmax / nSyn[imaxgmax], meangnmax / nSyn[imaxgmax]
            )

        elif self.RunInfo.Spirou in ["removelargest"]:
            cprint("c", "    Setting largest terminal to 0")
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print("    AMPA:  gMax, imax: ", gMax, imaxgmax)
            print("    NMDA:  gnmax: ", gnMax)
            gMax[imaxgmax] = 0.0
            gnMax[imaxgmax] = 0.0
            self.set_synapse(
                synapse[imaxgmax], 0.0, 0.0,
            )
        elif self.RunInfo.Spirou in ["removetwolargest"]:
            cprint("c", "    Setting TWO largest terminals to 0")
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            igMax = np.argsort(gMax)
            print("igmax: ", igMax)
            for i in igMax[-2:]:
                gMax[i] = 0.0
                gnMax[i] = 0.0
                print("    removed : AMPA:  gMax, imax: ", gMax, igMax)
                print("                     NMDA:  gnmax: ", gnMax)
                self.set_synapse(
                    synapse[i], 0.0, 0.0,
                )
            gm, gnm, ns = self.get_synapses(synapse, printvalues=True)
            # exit()

        elif self.RunInfo.Spirou in ["largestonly"]:
            nlarge = 1
            cprint("c", "    Setting all terminal inputs to 0 except the one largest")
        elif self.RunInfo.Spirou in ["twolargest"]:
            nlarge = 2
            cprint("c", "    Setting all terminal inputs to 0 except the two largest")
        elif self.RunInfo.Spirou in ["threelargest"]:
            nlarge = 3
            cprint("c", "    Setting all terminal inputs to 0 except the three largest")
        elif self.RunInfo.Spirou in ["fourlargest"]:
            nlarge = 4
            cprint("c", "    Setting all terminal inputs to 0 except the four largest")
        if self.RunInfo.Spirou in [
            "largestonly",
            "twolargest",
            "threelargest",
            "fourlargest",
        ]:
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            i_gMaxSorted = np.argsort(gMax)  # synapse with largest conductance
            print("    AMPA:  gMax, imax: ", gMax, i_gMaxSorted)
            print("    NMDA:  gnmax: ", gnMax)
            # print(f">>>> experiment: {self.RunInfo.Spirou:s}")
            # print(f">>>> nlarge: {nlarge:d}")
            for i in range(len(gMax)):
                if i not in i_gMaxSorted[-nlarge:]:
                    self.set_synapse(
                        synapse[i], 0.0, 0.0,
                    )
        cprint("g", "    Revised synapse values:")
        self.get_synapses(synapse, printvalues=True)
        if self.Params.testsetup:
            return
        nWorkers = MPROC.cpu_count()
        # workers
        CP.cprint("m", f"anrun : Parallel with {nWorkers:d} processes")
        # run using pyqtgraph's parallel support
        if self.Params.Parallel:
            with MP.Parallelize(
                enumerate(TASKS), results=tresults, workers=nWorkers
            ) as tasker:
                for j, x in tasker:
                    tresults = self.single_an_run(
                        j, synapseConfig, seeds, preCell, self.an_setup_time,
                    )
                    tasker.results[j] = tresults

        else:
            # Non parallelized version (with --noparallel flag - useful for debugging):
            for j, N in enumerate(range(nReps)):
                print("Repetition %d" % N)

                tresults[j] = self.single_an_run(
                    j, synapseConfig, seeds, preCell, self.an_setup_time,
                )
        for j, N in enumerate(range(nReps)):
            celltime[N], result[N] = self.retrieve_data(tresults[j])

        self.print_timing()

        self.analysis_filewriter(self.Params.cell, result, tag="delays")
        if self.Params.plotFlag:
            print("plotting")
            self.plot_an(celltime, result)

    def an_run_singles(self, exclude=False, verify=False):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig is list of tuples gathered from the configuration file. 
        Each tuple represents an AN fiber (SGC cell) with:
        
            * N sites
            * delay (ms)
            * spont rate group [1=low, 2=high, 3=high])
        
        This routine is special - it runs nReps for each synapse, turning off all of the other synapses
        by setting the synaptic conductance to 0 (to avoid upsetting initialization)
        if the "Exclude" flag is set, it turns OFF each synapse, leaving the others running.
        The output file either says "Syn", or "ExclSyn"

        Parameters
        ----------

        exclude : boolean (default: False)
            Set false to do one synapse at a time.
            Set true to do all BUT one synapse at a time.

        verify : boolean (default: False)
            Flag to control printing of various intermediate results

        Returns
        -------
            Nothing

        """
        cprint("c", "\n****AN Singles****\n")
        self.start_time = time.time()

        if self.Params.inputPattern is not None:
            fromdict = self.Params.inputPattern
            print(
                "Cell id: %s  using input pattern: %s" % (self.Params.cellID, fromdict)
            )
        else:
            fromdict = self.Params.cellID
        synapseConfig, celltype = self.cconfig.make_dict(fromdict)
        self.Params.SynapseConfig = synapseConfig
        nReps = self.RunInfo.nReps
        preCell, synapse, self.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig, celltype
        )
        nSyns = len(preCell)
        seeds = self.compute_seeds(nReps, synapseConfig)
        self.RunInfo.seeds = seeds  # keep the seed values too.
        k = 0
        parallel = self.Params.Parallel
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0
        # get the gMax's
        # gMax = np.zeros(len(synapse))
        # gMaxNMDA = np.zeros(len(synapse))

        gMax, ngMax, nSyn = self.get_synapses(synapse)
        for k in range(nSyns):
            # only enable or disable gsyn on the selected input
            if exclude:  # disable the one
                tagname = "ExcludeSyn%03d" % k
                for i in range(nSyns):
                    if i != k:  # set values for all others
                        self.set_synapse(
                            synapse[k], gMax[k] / nSyn[k], ngMax[k] / nSyn[k]
                        )
                self.set_synapse(synapse[k], 0.0, 0.0)  # but turn this one off

            if not exclude:  # enable a single synapse
                tagname = "Syn%03d" % k
                self.set_synapse(
                    synapse[k], gMax[k] / nSyn[k], ngMax[k] / nSyn[k]
                )  # set value for selected one
                for i in range(nSyns):
                    if i != k:
                        self.set_synapse(
                            synapse[i], 0.0, 0.0
                        )  # but turn all others off

            self.get_synapses(synapse, printvalues=True)
            celltime = [None] * nReps
            result = {}

            tresults = [None] * nReps
            if self.Params.testsetup:
                continue

            if parallel and self.Params.nWorkers > 1:
                nWorkers = self.Params.nWorkers
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with MP.Parallelize(
                    enumerate(TASKS), results=tresults, workers=nWorkers
                ) as tasker:
                    for j, x in tasker:
                        tresults = self.single_an_run(
                            j, synapseConfig, seeds, preCell, self.an_setup_time,
                        )
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    tresults[j] = self.single_an_run(
                        # self.post_cell,
                        j,
                        synapseConfig,
                        seeds,
                        preCell,
                        self.an_setup_time,
                    )
            for j, N in enumerate(range(nReps)):
                celltime[N], result[N] = self.retrieve_data(tresults[j])

            self.analysis_filewriter(self.Params.cell, result, tag=tagname)
            self.print_timing()
            if self.Params.plotFlag:
                self.plot_an(celltime, result)
        if self.Params.testsetup:
            return

    def an_run_omit_one(self):
        self.an_run_singles(exclude=True)

    def _search(
        self,
        variable: float = 0.0,
        target: float = 0.0,
        current_value: float = 0.0,
        factors: list = [0.6, 0.3],
        last: int = 0,
    ) -> (float, int):

        if current_value > target:
            if last == 1:
                factor = factors[0]
            else:
                factor = factors[1]
            last = 1
            ndiff = int(
                variable * factor * (current_value - target)
            )  # decrease "some" fraction of distance
            variable -= ndiff
        else:
            if last == -1:
                factor = factors[0]
            else:
                factor = factors[1]
            last = -1
            ndiff = int(
                variable * factor * (target - current_value)
            )  # increase "some" fraction of distance
            variable += ndiff
        return variable, last, ndiff

    def an_run_find_syn_thr(self):
        """
        This routine finds the number of active zones needed to produce
        a target amount of (entrainment)
        Iterative process:
        1. Only one terminal is active; the number of active zones is varied.
        2. The number of zones is adjusted using a linear estimator until the
            the target_ratio of bushy spikes to AN spikes is attained.
        Note: this may get caught bouncing between two input sizes until it

        reaches the maxiter terminations. This is probably ok.

        """
        self.start_time = time.time()
        synapseConfig, celltype = self.cconfig.make_dict(self.Params.cellID)
        self.Params.SynapseConfig = synapseConfig
        max_nsyn = 0
        for i, sc in enumerate(synapseConfig):  # for each terminal
            if synapseConfig[i]["nSyn"] > max_nsyn:
                max_nsyn = synapseConfig[i]["nSyn"]
            synapseConfig[i]["nSyn"] = 0  # turn off all but one synapse
        # algorithm:
        # start with largest input. If spike count is 0,
        # increase by
        nsyn = 180
        done = False
        target_ratio = 0.7
        iteration = 0
        maxiter = 25
        last = 0
        ndiff = 0
        min_diff = 0
        datestr = datetime.datetime.now().strftime(
            "%Y-%m-%d.%H-%M-%S"
        )  # get the actual start time for the top directory
        f_log = f"thrtest_{str(self.Params.cellID):s}_{target_ratio:.2f}.txt"
        ostr = f"testing syn threshold, cell {str(self.Params.cellID):s}"
        ostr += f" target:{target_ratio:.1f} {datestr:s}\n"
        with (open(f_log, "w")) as fh:
            fh.write(ostr)
        while not done:
            result = self.an_run_syn_thr_finder(celltype, synapseConfig, nsyn=nsyn)
            iteration += 1
            bu_st = np.array(result[0]["spikeTimes"])
            bu_st = bu_st[np.where(bu_st > 0.01)[0]]
            an_st = result[0]["inputSpikeTimes"][0]
            an_st = an_st[an_st > 0.01]
            n_bu = len(bu_st)
            n_an = len(an_st)
            ratio = float(n_bu) / float(n_an)
            output_str = f"Niter: {iteration:3d}  nSyn: {nsyn:4d}  nbu: {n_bu:4d}  n_an: {n_an:4d}"
            output_str += f" ratio: {ratio:.2f} target: {target_ratio:.2f} ndiff: {ndiff:3d}"
            print(output_str)
            with (open(f_log, "a")) as fh:
                fh.write(output_str + "\n")

            cprint("y", output_str)
            nsyn, last, ndiff = self._search(
                variable=nsyn,
                target=target_ratio,
                current_value=ratio,
                factors=[0.33, 0.65],
                last=last,
            )
            # if ratio > target_ratio:
            #     if last == 1:
            #         factor = 0.6
            #     else:
            #         factor = 0.3
            #     last = 1
            #     ndiff = int(nsyn*factor*(ratio-target_ratio)) # decrease "some" fraction of distance
            #     nsyn -= ndiff
            # else:
            #     if last == -1:
            #         factor = 0.6
            #     else:
            #         factor = 0.3
            #     last = -1
            #     ndiff = int(factor*nsyn*(target_ratio-ratio))  # increase "some" fraction of distance
            #     nsyn += ndiff
            if iteration > maxiter or ndiff <= min_diff:
                print("done,iteration, ndiff: ", iteration, ndiff)
                done = True
        print("\nFinal: ")
        output_str = f"Iterations: {iteration:3d}  nSyn: {nsyn:4d}  nbu: {n_bu:4d}  n_an: {n_an:4d} "
        output_str += f"ratio: {ratio:.2f} target: {target_ratio:.2f} ndiff: {ndiff:3d}"
        print(output_str)
        with (open(f_log, "a")) as fh:
            fh.write(output_str + "\n")

    def an_run_syn_thr_finder(self, celltype, synapseConfig, nsyn=0):
        """
        Run a set of reps with one configuration of synapses

        """
        nReps = self.RunInfo.nReps
        synapseConfig[0]["nSyn"] = nsyn
        preCell, synapse, self.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig, celltype
        )
        if self.Params.ANSynapseType == "multisite":
            for i, s in enumerate(synapse):
                s.terminal.relsite.Dep_Flag = 0
                # print(dir(s.terminal.relsite))

        seeds = self.compute_seeds(nReps, synapseConfig)
        # special conditions:
        # no depression
        # no variability

        celltime = []
        parallel = self.Params.Parallel
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0

        # tagname = "SynIO%03d"
        tresults = [None] * nReps
        celltime = [None] * nReps
        result = {}

        if self.Params.testsetup:
            return None
        if parallel and self.Params.nWorkers > 1:
            nWorkers = self.Params.nWorkers
            TASKS = [s for s in range(nReps)]
            # run using pyqtgraph's parallel support
            with MP.Parallelize(
                enumerate(TASKS), results=tresults, workers=nWorkers
            ) as tasker:
                for j, x in tasker:
                    tresults = self.single_an_run(
                        j, synapseConfig, seeds, preCell, self.an_setup_time,
                    )
                    tasker.results[j] = tresults
        else:  # easier to debug witout parallelization
            for j, N in enumerate(range(nReps)):
                tresults[j] = self.single_an_run(
                    j, synapseConfig, seeds, preCell, self.an_setup_time,
                )
        for j, N in enumerate(range(nReps)):
            celltime[N], result[N] = self.retrieve_data(tresults[j])

        # self.analysis_filewriter(self.Params.cell, result, tag=tagname % k)
        if self.Params.testsetup:
            return
        if self.Params.plotFlag:
            self.plot_an(celltime, result)
        for rk in result:
            r = result[rk]
            print(len(r["spikeTimes"]), len(r["inputSpikeTimes"][0]))

            # mpl.figure()
            # mpl.eventplot([np.array(r['spikeTimes'])*1e3, r['inputSpikeTimes'][0]], ['r', 'k'])
            # mpl.show()
        return result

    def an_run_IO_gSyn(self):
        """
        Establish AN inputs to soma, and run the model adjusting gmax over the reps from 0.5 to 4x.

        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=medium, 3=high])
        This routine runs a series of nReps for each synapse, turning off all of the other synapses.
        and varying the synaptic conductance to 0 (to avoid upsetting initialization)
        The output file either says "SynIO" or

        Parameters
        ----------
        exclude : boolean (default: False)
            Set false to do one synapse at a time.
            Set true to do all BUT one synapse at a time.

        verify : boolean (default: False)
            Flag to control printing of various intermediate results

        Returns
        -------
            Nothing

        """
        print("\n*** an_run_IO\n")
        self.start_time = time.time()
        synapseConfig, celltype = self.cconfig.make_dict(self.Params.cellID)
        self.Params.SynapseConfig = synapseConfig
        nReps = self.RunInfo.nReps

        preCell, synapse, self.R.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig, celltype
        )
        seeds = self.compute_seeds()
        nSyns = len(preCell)
        k = 0
        celltime = []
        parallel = self.Params.Parallel
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0
        # get the gMax's
        gMax = np.zeros(len(synapse))
        for i, s in enumerate(synapse):
            for p in s.psd.ampa_psd:
                gMax[i] = p.gmax
        for k in range(nSyns):
            # only enable gsyn on the selected input
            tagname = "SynIO%03d"
            tresults = [None] * nReps
            if self.Params.testsetup:
                continue
            if parallel and self.Params.nWorkers > 1:
                nWorkers = self.Params.nWorkers
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with MP.Parallelize(
                    enumerate(TASKS), results=tresults, workers=nWorkers
                ) as tasker:
                    for j, x in tasker:
                        for i, s in enumerate(synapse):
                            for p in s.psd.ampa_psd:
                                if i != k:
                                    p.gmax = 0.0  # disable all others
                                else:
                                    p.gmax = (
                                        4.0 * (j + 1) * gMax[i] / float(nReps + 1)
                                    )  # except the shosen one
                        tresults = self.single_an_run(
                            j, synapseConfig, seeds, preCell, self.an_setup_time,
                        )
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    for i, s in enumerate(synapse):
                        for p in s.psd.ampa_psd:
                            if i != k:
                                p.gmax = 0.0  # disable all others
                            else:
                                p.gmax = (
                                    4.0 * float(j + 1) * gMax[i] / float(nReps + 1)
                                )  # except the chosen one
                    tresults[j] = self.single_an_run(
                        j, synapseConfig, seeds, preCell, self.an_setup_time,
                    )
            # gmaxs = [
            #     4.0 * float(j + 1) * gMax[0] / float(nReps + 1) for j in range(nReps)
            # ]
            result = [None] * nReps
            for j, N in enumerate(range(nReps)):
                celltime[N], result[N] = self.retrieve_data(tresults[j])

            self.analysis_filewriter(self.Params.cell, result, tag=tagname % k)
        if self.Params.testsetup:
            return
        if self.Params.plotFlag:
            self.plot_an(celltime, result)

    def set_dbspl(self, signal, dbspl):
        """Scale the level of `signal` to the given dB_SPL."""
        p0 = 20e-6
        rms = np.sqrt(np.sum(signal ** 2) / signal.size)
        scaled = signal * 10 ** (dbspl / 20.0) * p0 / rms
        return scaled

    def single_an_run(self, j, synapseConfig, seeds, preCell, an_setup_time):
        """
        Perform a single run with all AN input on the target cell turned off except for input j.

        Parameters
        ----------
        j : int
            The input that will be active in this run

        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.

        seeds : 2d numpy array
            An array listing the seeds for starting each auditory nerve input spike train

        preCell : list
            A list of the preCell hoc objects attached to the synapses

        an_setup_time : time object

        Returns
        -------
        anresult : dict
            A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'

        """
        print(f"\n*** single_an_run: j={j:4d}")

        # print("restoring state from : ", self.Params.initStateFile)
        # cellInit.restore_initial_conditions_state(
        #     cell = self.post_cell,
        #     filename = self.Params.initStateFile,
        #     electrode_site = self.electrode_site,
        # )
        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        if isinstance(self.RunInfo.pip_start, float):
            pips = [self.RunInfo.pip_start]
        if self.RunInfo.soundtype == "tonepip":
            stim = sound.TonePip(
                rate=self.RunInfo.Fs,
                duration=self.RunInfo.run_duration,
                f0=self.RunInfo.F0,
                dbspl=self.RunInfo.dB,
                ramp_duration=self.RunInfo.RF,
                pip_duration=self.RunInfo.pip_duration,
                pip_start=pips,
            )
        elif self.RunInfo.soundtype in ["stationaryNoise", "noise"]:
            if self.RunInfo.soundtype == "noise":   # non-stationary noise generator seed changes on per-run basis
                self.RunInfo.noise_seed = self.RunInfo.noise_seed + j
            print(f" **** Noise type: {self.RunInfo.soundtype:s}  seed={self.RunInfo.noise_seed}")
            stim = sound.NoisePip(
                rate=self.RunInfo.Fs,
                duration=self.RunInfo.run_duration,
                dbspl=self.RunInfo.dB,
                pip_duration=self.RunInfo.pip_duration,
                pip_start=pips,
                ramp_duration=self.RunInfo.RF,
                seed=self.RunInfo.noise_seed,
            )
        elif self.RunInfo.soundtype == "SAM":
            stim = sound.SAMTone(
                rate=self.RunInfo.Fs,
                duration=self.RunInfo.run_duration,
                f0=self.RunInfo.F0,
                dbspl=self.RunInfo.dB,
                ramp_duration=self.RunInfo.RF,
                fmod=self.RunInfo.fmod,
                dmod=self.RunInfo.dmod,
                pip_duration=self.RunInfo.pip_duration,
                pip_start=pips,
            )
        elif self.RunInfo.soundtype in ["regularClicks", "poissonClicks"]:
            if self.RunInfo.soundtype == "poissonClicks":
                eventintervals = np.random.exponential(
                    1.0 / self.RunInfo.clickRate, int(self.RunInfo.clickTrainDuration*self.RunInfo.clickRate)
                )
                events = np.cumsum(eventintervals)
                events = events[events < self.RunInfo.clickTrainDuration]
            else:
                events = np.linspace(self.RunInfo.clickStart, 
                                             self.RunInfo.clickTrainDuration, 
                                             int(self.RunInfo.clickTrainDuration*self.RunInfo.clickRate)
                                             )
            stim = sound.ClickTrain(
                rate=self.RunInfo.Fs,
                duration=self.RunInfo.clickTrainDuration,
                dbspl=self.RunInfo.dB,
                click_duration=self.RunInfo.clickDuration,
                click_starts=events,
            )
        else:
            raise ValueError(
                "RunInfo sound type %s not implemented" % self.RunInfo.soundtype
            )
        stimWaveform = stim.generate()
        stimTimebase = stim.time
        for i, syn in enumerate(synapseConfig):
            nseed = seeds[j, i]
            if self.Params.SGCmodelType in ["Zilany"]:
                preCell[i].set_sound_stim(
                    stim, seed=nseed, simulator="matlab"
                )  # generate spike train, connect to terminal
            elif self.Params.SGCmodelType in ["cochlea"]:
                wf = self.set_dbspl(stim.generate(), self.RunInfo.dB)
                stim._sound = wf
                preCell[i].set_sound_stim(
                    stim, seed=nseed, simulator="cochlea"
                )  # generate spike train, connect to terminal
            else:
                raise ValueError(
                    "SGC model type type %s not implemented" % self.Params.SGCmodelType
                )
            ANSpikeTimes.append(preCell[i]._spiketrain)

        an_setup_time += time.time() - an0_time
        nrn_start = time.time()
        Vsoma = self.post_cell.hr.h.Vector()
        Vdend = self.post_cell.hr.h.Vector()
        rectime = self.post_cell.hr.h.Vector()
        if (
            "dendrite" in self.post_cell.all_sections
            and len(self.post_cell.all_sections["dendrite"]) > 0
        ):
            dendsite = self.post_cell.all_sections["dendrite"][-1]
            Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        else:
            dendsite = None

        self.allsecVec = OrderedDict()
        # cprint('m', f"save all sections flag: {str(self.Params.save_all_sections):s}")
        if self.Params.save_all_sections:
            for group in list(
                self.post_cell.hr.sec_groups.keys()
            ):  # get morphological components
                g = self.post_cell.hr.sec_groups[group]
                for section in list(g):
                    sec = self.post_cell.hr.get_section(section)
                    self.allsecVec[sec.name()] = self.post_cell.hr.h.Vector()
                    self.allsecVec[sec.name()].record(
                        sec(0.5)._ref_v, sec=sec
                    )  # recording of voltage all set up here
            cprint('r', f"save all sections: n = {len(self.allsecVec):d}")
        Vsoma.record(self.post_cell.soma(0.5)._ref_v, sec=self.post_cell.soma)
        rectime.record(self.post_cell.hr.h._ref_t)
        # print("restoring state from : ", self.Params.initStateFile)
        # cellInit.restore_initial_conditions_state(
        #     cell = self.post_cell,
        #     filename = self.Params.initStateFile,
        #     electrode_site = self.electrode_site,
        # )
        self.post_cell.hr.h.finitialize()
        self.post_cell.hr.h.tstop = self.RunInfo.run_duration * 1000.0
        self.post_cell.hr.h.t = 0.0
        self.post_cell.hr.h.batch_save()  # save nothing
        self.post_cell.hr.h.dt = self.Params.dtIC
        cprint("c", "single_an_run: Running")

        self.post_cell.hr.h.batch_run(
            self.post_cell.hr.h.tstop, self.post_cell.hr.h.dt, "an.dat"
        )
        nrn_run_time += time.time() - nrn_start
        # if len(rectime) == 0:
        #     npts = np.array(Vsoma).shape[0]
        #     rectime = np.linspace(0., npts*self.post_cell.hr.h.dt, num=npts)
        if dendsite is None:
            Vdend = np.zeros_like(Vsoma)

        anresult = {
            "Vsoma": np.array(Vsoma),
            "Vdend": np.array(Vdend),
            "time": np.array(rectime),
            "ANSpikeTimes": ANSpikeTimes,
            "stim": stim,
            "stimWaveform": stimWaveform,
            "stimTimebase": stimTimebase,
        }
        # print("max time: ", max(anresult["time"]))
        # mpl.plot(anresult["time"], anresult["Vsoma"], 'k-')
        # mpl.show()
        if self.Params.save_all_sections:
            anresult["allsecVec"] = self.allsecVec
        return anresult

    def plot_an(self, celltime, result):
        """
        Plot the cell's voltage reponse to the AN inputs

        Parameters
        ----------
        celltime : array (no default)
            time array for the cell voltage data
        result : dict from simulation

        Returns
        -------
        Nothing

        """
        if not self.Params.plotFlag:
            return
        # if self.Params.Parallel:
        #     return  # no plots when doing parallel runs...
        # nReps = self.RunInfo.nReps
        threshold = self.RunInfo.threshold
        fig, ax = mpl.subplots(3, 1, sharex=True)
        fig.suptitle("AN Inputs")
        # win = pgh.figure(title='AN Inputs')
        # layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
        for j, N in enumerate(range(len(result))):
            dt = np.mean(np.diff(celltime[N]))
            ax[0].plot(celltime[N], result[N]["somaVoltage"], c="k", linewidth=0.75)
            st = [int(x * 1e3 / dt) for x in result[N]["spikeTimes"]]
            ax[0].plot(
                celltime[N][st], result[N]["somaVoltage"][st], "ro", markersize=3.0
            )
            sta = [int(x / dt) for x in result[N]["inputSpikeTimes"][0]]
            ax[0].plot(
                celltime[N][sta], result[N]["somaVoltage"][sta], "b|", markersize=6
            )
            ax[0].plot(
                [np.min(celltime[N]), np.max(celltime[N])],
                [threshold, threshold],
                c=(0.5, 0.5, 0.5),
                linewidth=0.5,
            )
            dvdt = np.diff(result[N]["somaVoltage"]) / np.diff(
                celltime[N]
            )  # show in mV/ms
            ax[1].plot(celltime[N][:-1], dvdt, c="k", linewidth=0.75)
            # layout.getPlot(0).setXLink(layout.getPlot(1))
        if (
            "dendVoltage" in list(result[N].keys())
            and result[N]["dendVoltage"] is not None
        ):
            for j, N in enumerate(range(len(result))):
                ax[0].plot(
                    celltime[N],
                    result[N]["dendVoltage"][0],
                    c="b",
                    linewidth=0.6,
                    linestyle="--",
                )

        st = result[0]["inputSpikeTimes"][0]
        ax[2].eventplot(
            [st, np.array(result[0]["spikeTimes"]) * 1e3], colors=["k", "r"]
        )
        mpl.show()

    def cleanNeuronObjs(self):
        if not isinstance(self.RunInfo.electrodeSection, str):
            self.RunInfo.electrodeSection = str(
                self.RunInfo.electrodeSection.name()
            )  # electrodeSection
        # self.RunInfo.electrodeSectionName = str(self.RunInfo.electrodeSection.name())
        if self.RunInfo.dendriticElectrodeSection is not None and not isinstance(
            self.RunInfo.dendriticElectrodeSection, str
        ):
            self.RunInfo.dendriticElectrodeSection = str(
                self.dendriticElectrodeSection.name()
            )  # dendriticElectrodeSection,
        # dendriticSectionDistance = 100.0  # microns.

    def analysis_filewriter(self, filebase:str, result:dict, tag:str=""):
        """
        Write the analysis information to a pickled file

        Parameters
        ----------
        filebase : string (no default)
            base filename - *not used* (value replaced bye cellID)
        result : dict (no default)
            dict hlding results. Must be pickleable
        tag : string (default: '')
            tag to insert in filename string
        
        Returns
        -------
        Nothing
        """
        k = list(result.keys())
        # result will be a dict; each key is a repetition/run.
        requiredKeys = [
            "spikeTimes",
            "inputSpikeTimes",
            "somaVoltage",
            "time",
            "stimWaveform",
            "stimTimebase",
            # "allDendriteVoltages", # may be an empty list, or various dendrite V's
        ]

        res_mode = "reps"
        for rk in requiredKeys:
            if rk not in k:
                res_mode = "syn"
                break

        results = {}
        # results with be a dict with params, runinfo, modelpars and trials as keys
        print("\n*** analysis_filewriter\n")

        results["basename"] = self.Params.simulationFilename
        results["Params"] = self.Params  # include all the parameters of the run too
        # clean up Parame to remove PosixPath from filenames
        results["Params"].initStateFile = str(self.Params.initStateFile)
        results["Params"].initStateFile = str(self.Params.initStateFile)
        results["Params"].simulationFilename = str(self.Params.simulationFilename)
        results["Params"].hocfile = str(self.Params.hocfile)
        self.cleanNeuronObjs()
        results["runInfo"] = self.RunInfo
        results["modelPars"] = copy.deepcopy(self.post_cell.status)
        del results["modelPars"]["decorator"]  # remove neuron section objects
        results["Results"] = result
        results["mode"] = res_mode

        fout = (
            self.Params.simulationFilename
        )  # base name created by make-filename - do not overwrite
        if len(tag) > 0:
            fout = Path(str(fout).replace("_all_", "_" + tag + "_"))
        if self.Params.tagstring is not None:
            fnb = fout.stem
            fp = fout.parent
            fnb = str(fnb) + "_" + self.Params.tagstring + ".p"
            fout = Path(fp, fnb)
        with open(fout, "wb") as fh:
            pickle.dump(results, fh)
        print(f"**** Model output written to: ****\n   {str(fout):s}")
        self.ANFilename = str(fout)  # save most rexent value

    def get_hoc_file(self, hf:object):
        """
        Gets a section list from a previously loaded hoc file (morphology)
        and generates a list of section name (str)/reference (NEURON object) pairs
        
        Parameters
        ----------
        hf : object
            a NEURON object containing the geometry of the cell
        
        Returns
        -------
        Nothing
        """
        if hf.file_loaded is False:
            raise ValueError("No hoc file has been loaded")
        self.section_list = hf.get_section_prefixes()
        list(hf.sec_groups.keys())
        if len(hf.sec_groups) > 1:  # multiple names, so assign colors to structure type
            self.section_colors = {}
            # if self.Params.verbose:
            #     print("gethocfile # colors: ", len(self.hg.colorMap))
            #     print("gethocfile # sexcgroups: ", len(list(hf.sec_groups.keys())))
            # for i, s in enumerate(hf.sec_groups.keys()):
            #     if self.Params.verbose:
            #         print("gethocfile group: ", s)
            #     self.section_colors[s] = self.hg.colorMap[i]
        #        else: # single section name, assign colors to SectionList types:
        #        self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
        #            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

        (v, e) = hf.get_geometry()
        self.clist = []

        for si in hf.sections:  # self.section_list[s]:
            hf.h("access %s" % si)
            sr = hf.h.SectionRef()
            n1 = hf.h.cas().name()
            if sr.has_parent() == 1:
                x = sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])


def main():

    (parsedargs, params, runinfo,) = vcnmodel.model_params.getCommands(
        toml_dir="toml"
    )  # get from command line
    model = ModelRun(
        params=params, runinfo=runinfo, args=parsedargs
    )  # create instance of the model
    if parsedargs.displayMode != "None":
        model.view_model()
    else:
        model.run_model()


if __name__ == "__main__":
    main()
