import copy
import dataclasses
import datetime
import json
import pickle
import re
import sys
import time
import timeit
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import toml
from cnmodel import cells
from cnmodel.decorator import Decorator
from cnmodel.util import sound
from matplotlib import pyplot as mpl
from neuronvis import hoc_graphics as hoc_graphics
from pylibrary.tools import cprint as CP
from pylibrary.tools import utility as pu  # access to a spike finder routine
from pyqtgraph import multiprocess as MP

import vcnmodel.model_params
from vcnmodel import cell_config as cell_config
from vcnmodel import cellInitialization as cellInit
from vcnmodel.adjust_areas import AdjustAreas
from vcnmodel.generate_run import GenerateRun

matplotlib.use("Qt5Agg")

"""
model_run.py

Run a model based on a hoc structure, decorating the structure with ion channels and synapses.

Requires:
Python 3.7 (anaconda distribution)
Neuron7.7 (neuron.yale.edu)
pyqtgraph (Campagnola, from github)
neuronvis (Campagnola/Manis, from github)
cnmodel (Campagnola/Manis, from github)
vcnmodel parts:
    generate_run
    cellInitialization
    analysis

cochlea # Rudniki python implementation of Zilany model
thorns  # required for cochlea


Expects the following directory structure:
(example)
VCN_Cells/   # top level for data
    cell_ID/    # VCN_c18, for example (first argument in call should be this directory name)
        MorphologyFiles/  # location for swc and hoc files
            VCN_c18_755V2.hoc (cell body scaled version)
            VCN_c18_755.hoc  (unscaled version)
            misc.swc  (various swc files that were translated to hoc files)
        InitializationFiles/  # location for Init files
            IV/  # initialization for just the IV with no synaptic input
                VCN_c18_755v2.ninit  # different base structures
                VCN_c18755.ninit
            AN/  # initialization for synaptic inputs
                (ditto)  # different base structures of input arrangements
        Simulations/
            IV/  # results from IV simulations
            AN/  # results from AN simulations


usage: model_run [-h] [--type {Bushy,TStellate,DStellate}]
                 [--model {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,
                        XM13PasDend,Calyx,MNTB,L23Pyr}]
                 [--modeltype {II,II-I,I-II,I-c,I-t,II-o}]
                 [--sgcmodel {Zilany,cochlea}]
                 [--protocol {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,
                        runANIO,runANSingles,runANOmitOne,gifnoise}]
                 [-H] [--hocfile HOCFILE] [-F] [--inputpattern INPUTPATTERN]
                 [--stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}]
                 [--check] [-C CONFIGFILE] [-d DB] [-f F0]
                 [--duration PIP_DURATION] [-r NREPS] [--seed SEED]
                 [-S {LS,MS,HS,fromcell}] [--synapsetype {simple,multisite}]
                 [--depression {0,1}] [--fmod FMOD] [--dmod DMOD]
                 [--S2M SIGNALTOMASKER] [--cmmrmode {CM,CD,REF}]
                 [--spirou {all,max=mean,all=mean}]
                 [--soma_inflate SOMA_INFLATION] [--soma_autoinflate]
                 [--dendrite_inflate DENDRITE_INFLATION]
                 [--dendrite_autoinflate] [--dendrite_from_soma]
                 [--ASA_from_soma] [--tagstring TAGSTRING] [-a AMPASCALE]
                 [--allmodes] [--sequence SEQUENCE] [--plot]
                 [--workers NWORKERS] [--noparallel] [--auto] [--saveall]
                 [--verbose] [--gifi GIF_I0] [--gifsigma GIF_SIGMA]
                 [--giffmod GIF_FMOD] [--giftau GIF_TAU] [--gifdur GIF_DUR]
                 [--gifskew GIF_SKEW]
                 cell

Simulate activity in a reconstructed model cell

positional arguments:
  cell                  Select the cell (no default)

optional arguments:
  -h, --help            show this help message and exit
  --type {Bushy,TStellate,DStellate}, -T {Bushy,TStellate,DStellate}
                        Define the cell type (default: Bushy)
  --model {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr},
                -M {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}
                        Define the model type (default: XM13)
  --modeltype {II,II-I,I-II,I-c,I-t,II-o}
                        Define the model type (default: XM13)
  --sgcmodel {Zilany,cochlea}
                        Define the SGC model type (default: Zilany)
  --protocol {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,runANIO,
                runANSingles,runANOmitOne,gifnoise}
                , -P {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,runANIO,
                runANSingles,runANOmitOne,gifnoise}
                        Protocol to use for simulation (default: IV)
  -H, --defaulthoc      Use default hoc file for this cell
  --hocfile HOCFILE     hoc file to use for simulation (default is the
                        selected "cell".hoc)
  -F, --full            Use "full" hoc file as in "VCN_c02_Full.hoc instead of
                        VCN_c02.hoc")
  --inputpattern INPUTPATTERN, -i INPUTPATTERN
                        cell input pattern to use (substitute) from
                        cell_config.py
  --stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}, -s {tonepip,noise,stationaryNoise,SAM,CMMR}
                        Define the stimulus type (default: tonepip)
  --check, -/           Only check command line for valid input; do not run
                        model
  -C CONFIGFILE, --configfile CONFIGFILE
                        Read a formatted configuration file (JSON, TOML) for
                        commands
  -d DB, --dB DB        Set sound intensity dB SPL (default 30)
  -f F0, --frequency F0
                        Set tone frequency, Hz (default 4000)
  --duration PIP_DURATION
                        Set sound stimulus duration (sec; default 0.1)
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
  --fmod FMOD           Set SAM modulation frequency
  --dmod DMOD           Set SAM modulation depth (in percent)
  --S2M SIGNALTOMASKER  Signal to Masker ratio (dB)
  --cmmrmode {CM,CD,REF}
                        Specify mode (from: ['CM', 'CD', 'REF'])
  --spirou {all,max=mean,all=mean}
                        Specify spirou experiment type....
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
Example:
Set up initialization:
python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol initIV --model XM13
Then run model:
python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol runIV --model XM13

Where should we look for data?

wheres_the_data.toml

"""


# import pylibrary.plotting.pyqtgraph_plothelpers as pgh


showCell = True

cprint = CP.cprint


class ModelRun:
    def __init__(self, params: dataclass = None, runinfo: dataclass = None, args=None):
        self.Params = params
        self.RunInfo = runinfo
        print(args)

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
        self.cconfig = cell_config.CellConfig(verbose=self.Params.verbose)

        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells")}
        self.baseDirectory = self.datapaths["baseDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"

    def print_modelsetup(self):
        """
        Print out all of the parameters in the model
        """

        # pp = pprint.PrettyPrinter(indent=4, width=120)
        print("Params:\n", self.Params)
        # pp.pprint(dataclasses.asdict(self.Params))
        print("RunInfo:\n", self.RunInfo)
        # pp.pprint(dataclasses.asdict(self.RunInfo))

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
        self.cellID = Path(
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

        model_Pars = f"{str(self.Params.modelName):s}_{str(self.Params.modelType):s}"
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
        if self.RunInfo.runProtocol in ["initIV", "initandrunIV", "runIV"]:
            simMode = "IV"
            initPath = Path(self.baseDirectory, self.cellID, self.initDirectory)
            self.mkdir_p(initPath)  # confirm existence of the initialization directory
            if self.Params.initStateFile is None:
                fn0 = f"{simMode:s}_Init_{self.cellID:s}_inp={inputs:s}_{model_Pars:s}"
                fn0 += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
                fn0 += f"_{self.Params.SRType:2s}"
                if self.Params.tagstring is not None:
                    fn0 += f"_{self.Params.tagstring:s}"
                fn0 += ".dat"
                self.Params.initStateFile = Path(initPath, fn0)
                if (
                    self.RunInfo.runProtocol.startswith("init")
                    and self.Params.initStateFile.is_file()
                ):
                    self.Params.initStateFile.unlink()  # delete old initializaiton file first
            simPath = Path(
                self.baseDirectory,
                self.cellID,
                self.simDirectory,
                simMode,
                run_directory,
            )
            self.mkdir_p(simPath)  # confirm that output path exists
            fn = f"{simMode:s}_Result_{self.cellID:s}_inp={inputs:s}_{model_Pars:s}"
            fn += f"_pulse_{add_Pars:s}_{self.Params.ANSynapseType:s}"
            fn += f"_{self.Params.SRType:2s}"
            if self.Params.tagstring is not None:
                fn += f"_{self.Params.tagstring:s}"
            fn += ".p"
            self.Params.simulationFilename = Path(simPath, fn)

        if (
            self.RunInfo.runProtocol.startswith("runAN")
            or self.RunInfo.runProtocol == "initAN"
            or self.RunInfo.runProtocol == "runANSingles"
            or self.RunInfo.runProtocol == "runANPSTH"
        ):
            simMode = "AN"

            initPath = Path(self.baseDirectory, self.cellID, self.initDirectory)
            self.mkdir_p(initPath)  # confirm existence of that file
            if self.Params.initStateFile is None:
                fn0 = f"{simMode:s}_Init_{self.cellID:s}_inp={inputs:s}_{model_Pars:s}"
                fn0 += f"_{add_Pars:s}_{self.Params.ANSynapseType:s}"
                if self.RunInfo.soundtype in ["SAM", "sam"]:
                    fn0 += f"_{self.RunInfo.fmod:03.1f}_{int(self.RunInfo.dmod):03d}"
                fn0 += f"_{self.Params.SRType:2s}.p"
                self.Params.initStateFile = Path(initPath, fn0)
                # print("AN Initialization file: ", self.Params.initStateFile.name)
                print("initstatefile exists?: ", self.Params.initStateFile.is_file())
                if (
                    self.RunInfo.runProtocol.startswith("initAN")
                    and self.Params.initStateFile.is_file()
                ):
                    self.Params.initStateFile.unlink()  # delete old initializaiton file first
            simPath = Path(
                self.baseDirectory,
                self.cellID,
                self.simDirectory,
                simMode,
                run_directory,
            )
            fn = f"{simMode:s}_Result_{self.cellID:s}_inp={inputs:s}_{model_Pars:s}"
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
        spont_rate_type : string
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
        store the start time...
        """
        self.Params.StartTime = starttime

    def mkdir_p(self, path: [str, Path]):
        try:
            Path.mkdir(path, parents=True, exist_ok=True)
        except IOError:
            raise FileNotFoundError(f"Cannot create path: {str(path):s}")

    def fix_singlets(self, h: object, harsh: bool = False):
        badsecs = []
        for i, sec in enumerate(h.allsec()):
            if sec.n3d() == 1:
                badsecs.append(sec)
                # print(f'Fixing Singlet Section: {str(sec):s}')
                # print(dir(sec))
                parent = sec.parentseg()
                # if (
                #     parent is None and i > 0
                # ):  # first section might not have a parent... obviously
                #     raise ValueError(
                #         f"Section: {str(sec):s} does not have a parent, and this must be fixed manually"
                #     )
                # if i == 0 and parent is None:
                if parent is None:  # until we can find a better way.
                    continue
                # print(f'    Section connected to {str(parent):s} ')
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
            cprint("g", "**** No singlets found in this hoc file")
        badsecs = []
        for sec in h.allsec():
            if sec.n3d() == 1:
                badsecs.append(sec)
        if len(badsecs) == 0:
            cprint("g", "No bad sections")
        else:
            cprint(
                "r", f"found {len(badsecs):d} bad sections with sec.n3d() = 1; repaired"
            )

    def setup_model(self, par_map: dict = None):
        """
        Main entry routine for running all models
        here we set up the model, assigning various parameters, and dealing
        with various morphological manipulations such as inflating areas.

        Parameters
        ----------
        par_map : dict (default: None)
            A dictionary of parameters, passed to models that are run (not used).

        Returns
        -------
            Nothing

        """
        dendrite_names = [
            "Proximal_Dendrite",
            "Distal_Dendrite",
            "Dendritic_Hub",
            "Dendritic_Swelling" "dend",
            "proximaldendrite",
            "distandendrite",
        ]

        if self.Params.verbose:
            print("run_model entry")
        if par_map is not None and "id" in list(par_map.keys()):
            self.idnum = par_map["id"]
        else:
            self.idnum = 9999
        print(f"Cell ID: {self.Params.cell:s}")
        if self.Params.cell is None:
            return
        self.cellID = Path(
            self.Params.cell
        ).stem  # os.path.splitext(self.Params.cell)[0]

        if self.Params.verbose:
            print(f"Morphology directory: {self.morphDirectory:s}")
        if self.Params.fullhocfile:
            self.Params.hocfile = self.Params.cell + "_Full.hoc"
            cprint("c", f"Using '_Full' hoc file: {self.Params.hocfile:s}")
        elif self.Params.usedefaulthoc or self.Params.hocfile is None:
            self.Params.hocfile = self.Params.cell + ".hoc"  # try to make a name
            cprint("c", f"Using default hoc file: {self.Params.hocfile:s}")
        else:  # a passed hoc file will be used
            cprint("c", f"Using passed hoc file: {self.Params.hocfile:s}")

        filename = Path(
            self.baseDirectory, self.cellID, self.morphDirectory, self.Params.hocfile
        )
        # if filename.with_suffix(".hocx").is_file():  # change to preferred if available
        #     filename = filename.with_suffix(".hocx")
        # cprint("yellow", f"Using hoc file: {str(filename):s}")

        # instantiate cells
        if self.Params.cellType in ["Bushy", "bushy"]:
            print("Creating a bushy cell (run_model) ")
            from cnmodel import data

            changes = None
            nach = None  # uses default
            if (
                self.Params.modelName == "XM13_nacncoop"
                and self.Params.dendriteMode == "normal"
            ):
                from vcnmodel.model_data import data_XM13_nacncoop as CHAN

                nach = "nacncoop"
                changes = data.add_table_data(
                    "XM13_nacncoop_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_nacncoop_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            if (
                self.Params.modelName == "XM13_nacncoop"
                and self.Params.dendriteMode == "passive"
            ):
                from vcnmodel.model_data import \
                    data_XM13_nacncoop_pasdend as CHAN

                nach = "nacncoop"
                changes = data.add_table_data(
                    "XM13_nacncoop_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_nacncoop_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            if (
                self.Params.modelName == "XM13_nacncoop"
                and self.Params.dendriteMode == "active"
            ):
                from vcnmodel.model_data import \
                    data_XM13_nacncoop_actdend as CHAN

                nach = "nacncoop"
                changes = data.add_table_data(
                    "XM13_nacncoop_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_nacncoop_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            elif self.Params.modelName == "XM13":
                from model_data import data_XM13 as CHAN

                nach = "nav11"
                changes = data.add_table_data(
                    "XM13_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            elif self.Params.modelName == "XM13_nacn":
                from model_data import data_XM13_nacn as CHAN

                nach = "nacn"
                changes = data.add_table_data(
                    "XM13_nacn_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_nacn_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            elif self.Params.modelName == "XM13_nabu":
                from model_data import data_XM13_nabu as CHAN

                nach = "nabu"
                changes = data.add_table_data(
                    "XM13_nabu_channels",
                    row_key="field",
                    col_key="model_type",
                    species="mouse",
                    data=CHAN.ChannelData,
                )
                changes_c = data.add_table_data(
                    "XM13_nabu_channels_compartments",
                    row_key="parameter",
                    col_key="compartment",
                    species="mouse",
                    model_type="II",
                    data=CHAN.ChannelCompartments,
                )

            if changes is not None:
                data.report_changes(changes)
                data.report_changes(changes_c)
            self.post_cell = cells.Bushy.create(
                morphology=str(filename),
                decorator=Decorator,
                species=self.Params.species,
                modelName=self.Params.modelName,
                modelType=self.Params.modelType,
                nach=nach,
            )

        elif self.Params.cellType in ["tstellate", "TStellate"]:
            print("Creating a t-stellate cell (run_model) ")
            self.post_cell = cells.TStellate.create(
                morphology=str(filename),
                decorator=Decorator,
                species=self.Params.species,
                modelType=self.Params.modelName,
            )

        elif self.Params.cellType in ["dstellate", "DStellate"]:
            print("Creating a D-stellate cell (run_model)")
            self.post_cell = cells.DStellate.create(
                morphology=str(filename),
                decorator=Decorator,
                species=self.Params.species,
                modelType=self.Params.modelName,
            )
        else:
            raise ValueError(f"cell type {self.Params.cellType:s} not implemented")

        AdjA = AdjustAreas()
        AdjA.sethoc_fromCNcell(self.post_cell)
        hoc_somaarea = AdjA.get_hoc_area(["soma"])
        hoc_dendritearea = AdjA.get_hoc_area(dendrite_names)
        cprint(
            "y",
            f"HOC: Original Soma area: {hoc_somaarea:.2f}  Dendrite Area: {hoc_dendritearea:.2f} um^2",
        )

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

        self.fix_singlets(self.post_cell.hr.h)
        if (
            self.Params.soma_autoinflate
        ):  # get values and inflate soma automatically to match mesh
            cprint("c", "Soma Autoinflation")

            inflateratio = self.cconfig.get_soma_ratio(self.cellID)
            if np.isnan(inflateratio):
                raise ValueError(
                    f"Soma Inflation Ratio is not defined for cell {self.cellID:s}"
                )
            self.Params.soma_inflation = inflateratio

        if self.Params.soma_inflation != 1.0:
            cprint("c", "    Inflating soma")
            self.post_cell.hr.h.finitialize()
            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -50.0]
            )

            print(
                f"    Original Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}"
            )
            # origdiam = {}

            AdjA = AdjustAreas(method="pt3d")
            AdjA.sethoc_fromCNcell(self.post_cell)
            AdjA.adjust_diameters(sectypes=["soma", "Soma"], inflateRatio=inflateratio)
            # AdjA.plot_areas(pt3d)

            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -60.0]
            )
            print(f"    New Rin after somatic inflation: {rtau['Rin']:.2f}", end="")
            print(f" , tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}")
        if self.Params.ASA_fromsoma:
            self.Params.ASA_inflation = self.Params.soma_inflation

        if self.Params.dendrite_fromsoma:
            self.Params.dendrite_inflation = self.Params.soma_inflation
            cprint("c", "Dendrite inflation uses soma inflation.")
        elif (
            self.Params.dendrite_autoinflate
        ):  # get values and inflate soma automatically to match mesh
            cprint("c", "Dendrite Autoinflation")
            inflateratio = self.cconfig.get_dendrite_ratio(self.cellID)
            if np.isnan(inflateratio):
                raise ValueError("Dendrite Inflation Ration is not defined!")
            self.Params.dendrite_inflation = inflateratio

        if self.Params.dendrite_inflation != 1.0:
            cprint("c", "    Inflating dendrite")
            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -55.0]
            )
            print(
                f"     Original Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}"
            )

            AdjA = AdjustAreas()
            AdjA.sethoc_fromCNcell(self.post_cell)
            AdjA.adjust_diameters(sectypes=dendrite_names, inflateRatio=inflateratio)
            # AdjA.plot_areas(pt3d)  # just to take a look at the adjustment

            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -55.0]
            )
            print(
                f"     New Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}"
            )

        for group in list(self.post_cell.hr.sec_groups.keys()):
            g = self.post_cell.hr.sec_groups[group]
            for section in list(g):
                self.post_cell.hr.get_section(section).Ra = self.post_cell.hr.h.Ra
                if self.Params.verbose:
                    print("Section: ", section)
                    print("Ra: ", self.post_cell.hr.get_section(section).Ra)

        electrode_section = list(self.post_cell.hr.sec_groups["soma"])[0]
        self.RunInfo.electrodeSection = self.post_cell.hr.get_section(electrode_section)
        self.RunInfo.electrodeSectionName = "soma"
        self.hg = hoc_graphics
        self.get_hoc_file(self.post_cell.hr)

        # if self.RunInfo.dendriticElectrodeSection is not None:
        #     dendritic_sites = list(
        #         self.post_cell.hr.sec_groups[self.Params.dendriteelectrodedendrite]
        #     )[0]
        #     self.dendriticelectrode_site = self.post_cell.hr.get_section(
        #         dendritic_electrode
        #     )
        #     self.dendriticElectrodeSection = self.Params.dendriteelectrodedendrite

        if self.Params.verbose:
            if par_map is not None:
                print("Listing par_map (run_model): ", par_map)
            self.post_cell.hr.h.topology()

        self.post_cell.set_nseg(freq=self.Params.lambdaFreq)

        # handle the following protocols:
        # ['initIV', 'initAN', 'runIV', 'run', 'runANSingles', 'gifnoise']

        self._make_filenames()  # make filenames AFTER all manipulations of the cell

        self.Params.setup = True

    def view_model(self, par_map: dict = None):
        if not self.Params.setup:
            self.setup_model(par_map=par_map)
        import neuronvis as NV

        if self.cellID in list(vcnmodel.model_params.display_orient_cells):
            angles = vcnmodel.model_params.display_orient_cells[self.cellID]
        else:
            angles = [200.0, 0.0, 0.0]
        print("Rendering via model_run2")
        NV.hocRender.Render(
            hoc_file=self.post_cell.hr.h,  # hoc_file,
            display_style=self.Params.displayStyle,  # display_mode,
            display_renderer="pyqtgraph",  # args["renderer"],
            display_mode=self.Params.displayMode,  # 'sec-type', # args["display"],
            mechanism=self.Params.displayMechanism,  # args["mechanism"],
            initial_view=angles,
            alpha=1.0,  # args["alpha"],
            sim_data=None,  # sim_data,
            # output_file=str(Path('Renderings',
            #     f'{self.cellID:s}_{self.Params.displayMechanism:s}_{self.Params.dendriteMode:s}.png')),
        )
        exit()

    # def view_cellparts(self, par_map: dict = None):
    #     cprint('r', 'cellparts')
    #     if not self.Params.setup:
    #         self.setup_model(par_map=par_map)
    #     import neuronvis as NV
    #     viewitem = "sec-type"
    #     mech = self.Params.displaychannel
    #     NV.hocRender.Render(
    #         hoc_file=self.post_cell.hr.h,# hoc_file,
    #         display_style='surface', # display_mode, (cylinders, graph, volume, surface)
    #         display_renderer='pyqtgraph', # args["renderer"],
    #         display_mode=viewitem, # 'sec-type', # args["display"],
    #         alpha=1, # args["alpha"],
    #         sim_data=None, # sim_data,
    #         output_file=str(Path('Renderings', f'{self.cellID:s}_parts.png')),
    #     )
    #     exit()

    def run_model(self, par_map: dict = None):
        self.run_starttime = datetime.datetime.now()
        if not self.Params.setup:
            self.setup_model(par_map=par_map)

        dispatcher = {
            "initIV": self.initIV,
            "runIV": self.iv_run,
            "testIV": self.testIV,
            "initAN": self.an_run_init,  # (self.post_cell, make_an_intial_conditions=True),
            "runANPSTH": self.an_run,
            "runANIO": self.an_run_IO,
            "runANSingles": self.an_run_singles,
            "runANOmitOne": self.an_run_omit_one,  # .(self.post_cell, exclude=True),
            "gifnoise": self.noise_run,
        }
        if self.RunInfo.runProtocol in [
            "runANPSTH",
            "runANSingles",
            "runANOmitOne",
            "runANIO",
        ]:
            self.RunInfo.run_duration = (
                np.sum(self.RunInfo.pip_start)
                + self.RunInfo.pip_duration
                + self.RunInfo.pip_offduration
            )

        dispatcher[self.RunInfo.runProtocol]()

    def initIV(self):
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        cellInit.get_initial_condition_state(
            self.post_cell,
            tdur=self.Params.initialization_time,
            filename=self.Params.initStateFile,
            electrode_site=self.R.electrode_site,
        )
        print(f"Ran to get initial state for {self.post_cell.hr.h.t:.1f} msec")

    def testIV(self):
        self.R = GenerateRun(
            self.Params, self.RunInfo, self.post_cell, idnum=self.idnum, starttime=None,
        )
        cellInit.test_initial_conditions(
            self.post_cell,
            filename=self.Params.initStateFile,
            electrode_site=self.R.electrode_site,
        )
        self.R.testRun(initfile=self.Params.initStateFile)
        return  # that is ALL, never make testIV/init and then keep running.

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
            self.baseDirectory, self.cellID, self.simDirectory, "IV"
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

    def noise_run(self, par_map: dict = {}):
        """
        Main entry routine for running noise into cell current injection (for generating GIF models)

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
        # ivinitfile = Path(self.baseDirectory, self.cellID,
        #                         self.initDirectory, self.Params.initStateFile)
        self.R.runInfo.folder = Path(
            self.baseDirectory, self.cellID, self.simDirectory, "Noise"
        )
        if self.Params.verbose:
            print("noise_run: calling do_run")
        nworkers = self.Params.nWorkers
        #        print(self.Params.Parallel)
        if self.Params.Parallel is False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
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
            for i, s in enumerate(synapses):
                print("getsyn: i nsyn gmax nmdagmax  : ", i, nSyn[i], gMax[i], ngMax[i])
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

    def an_run_init(self):
        self.an_run(make_an_intial_conditions=True)

    def an_run(
        self, verify: bool = False, make_an_intial_conditions: bool = False,
    ):
        """
        Establish AN inputs to soma, and run the model.
        Requires a synapseConfig list of dicts from cell_config.makeDict()
        each list element represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high],
                    type (e or i), segment and segment location)
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
            cprint("g", "\n*** Initializing for an_run\n")
        else:
            cprint("g", "\n*** an_run\n")

        if self.Params.inputPattern is not None:
            fromdict = self.Params.inputPattern
            print(f"Cell id: {self.cellID:s} using input pattern: {fromdict:s}")
        else:
            print(f"Cell id: {self.cellID:s} is using 'self' input pattern")
            fromdict = self.cellID

        synapseConfig, celltype = self.cconfig.makeDict(
            fromdict, self.Params.ASA_inflation
        )
        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        # for i, s in enumerate(synapseConfig):
        #     s[1] = s[3]*0.001/3.0
        #     print 'delay for input %d is %8.4f msec' % (i, s[1])

        nReps = self.RunInfo.nReps
        threshold = self.RunInfo.threshold  # spike threshold, mV

        preCell, synapse, self.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig=synapseConfig, celltype=celltype
        )

        # see if we need to save the cell state now.
        if make_an_intial_conditions:
            #            print('getting initial conditions for AN')
            # aninitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params.initStateFile)
            cellInit.get_initial_condition_state(
                self.post_cell,
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
        # spikeTimes = {}
        # inputSpikeTimes = {}
        # somaVoltage = {}
        # dendriteVoltage = {}
        celltime = []
        # stim = {}
        # stimWaveform = {}
        allDendriteVoltages = (
            {}
        )  # Saving of all dendrite voltages is controlled by --saveall flag
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0

        nWorkers = self.Params.nWorkers
        TASKS = [s for s in range(nReps)]
        tresults = [None] * len(TASKS)
        result = {}

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

        elif self.RunInfo.Spirou in ["all=mean"]:
            print("    setting ALL to the mean of all inputs (no variance)")
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print(self.RunInfo.Spirou)
            print("    AMPA: mean, gmax, imax: ", meangmax, gMax, imaxgmax)
            print("    NMDA: mean, gnmax: ", meangnmax, gnMax)
            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            for i in range(len(synapse)):
                self.set_synapse(synapse[i], meangmax / nSyn[i], meangnmax / nSyn[i])
            cprint("g", "    Revised synapse values:")
            self.get_synapses(synapse, printvalues=True)

        elif self.RunInfo.Spirou in ["max=mean"]:
            cprint("c", "Setting largest to mean of all inputs")
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
            cprint("g", "    Revised synapse values:")
            self.get_synapses(synapse, printvalues=True)

        elif self.RunInfo.Spirou in ["removelargest"]:
            cprint("c", "Setting largest terminal to 0")
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print("    AMPA:  gMax, imax: ", gMax, imaxgmax)
            print("    NMDA:  gnmax: ", gnMax)
            gMax[imaxgmax] = 0.0
            gnMax[imaxgmax] = 0.0
            self.set_synapse(
                synapse[imaxgmax], 0.0, 0.0,
            )
            cprint("g", "    Revised synapse values:")
            self.get_synapses(synapse, printvalues=True)

        elif self.RunInfo.Spirou in ["largestonly"]:
            nlarge = 1
            cprint("c", "Setting all terminal inputs to 0 except the one largest")
        elif self.RunInfo.Spirou in ["twolargest"]:
            nlarge = 2
            cprint("c", "Setting all terminal inputs to 0 except the two largest")
        elif self.RunInfo.Spirou in ["threelargest"]:
            nlarge = 3
            cprint("c", "Setting all terminal inputs to 0 except the three largest")
        elif self.RunInfo.Spirou in ["fourlargest"]:
            nlarge = 4
            cprint("c", "Setting all terminal inputs to 0 except the four largest")
        if self.RunInfo.Spirou.find("largest") >= 0:
            gMax, gnMax, nSyn = self.get_synapses(synapse)
            i_gMaxSorted = np.argsort(gMax)  # synapse with largest conductance
            print("    AMPA:  gMax, imax: ", gMax, i_gMaxSorted)
            print("    NMDA:  gnmax: ", gnMax)

            for i in range(len(gMax)):
                if i not in i_gMaxSorted[-nlarge:]:
                    self.set_synapse(
                        synapse[i], 0.0, 0.0,
                    )
            cprint("g", "    Revised synapse values:")
            self.get_synapses(synapse, printvalues=True)

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
            # retreive the data
            for j, N in enumerate(range(nReps)):
                celltime.append(tresults[j]["time"])  # (self.time)
                spikeTimes = pu.findspikes(
                    tresults[j]["time"],
                    tresults[j]["Vsoma"],
                    threshold,
                    t0=0.0,
                    t1=self.RunInfo.run_duration * 1000.0,
                    dt=1.0,
                    mode="peak",
                )
                spikeTimes = self.clean_spiketimes(spikeTimes)
                inputSpikeTimes = tresults[j][
                    "ANSpikeTimes"
                ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage = np.array(tresults[j]["Vsoma"])
                dendriteVoltage = np.array(tresults[j]["Vdend"])
                stimWaveform = np.array(tresults[j]["stimWaveform"])
                stimTimebase = np.array(tresults[j]["stimTimebase"])
                result[N] = {
                    "spikeTimes": spikeTimes,
                    "inputSpikeTimes": inputSpikeTimes,
                    "time": np.array(celltime[0]),
                    "somaVoltage": somaVoltage,
                    "dendriteVoltage": dendriteVoltage,
                    "allDendriteVoltages": allDendriteVoltages,
                    "stimWaveform": stimWaveform,
                    "stimTimebase": stimTimebase,
                }
                if self.Params.save_all_sections:  # just save soma sections
                    # for section in list(g):
                    allDendriteVoltages[N] = tresults[j]["allsecVec"]

        else:
            # Non parallelized version (with --noparallel flag - useful for debugging):
            for j, N in enumerate(range(nReps)):
                print("Repetition %d" % N)

                tresults[j] = self.single_an_run(
                    j, synapseConfig, seeds, preCell, self.an_setup_time,
                )

                celltime.append(tresults[j]["time"])  # (self.time)
                spikeTimes = pu.findspikes(
                    tresults[j]["time"],
                    tresults[j]["Vsoma"],
                    threshold,
                    t0=0.0,
                    t1=self.RunInfo.run_duration * 1000.0,
                    dt=1.0,
                    mode="peak",
                )
                spikeTimes = self.clean_spiketimes(spikeTimes)
                inputSpikeTimes = tresults[j][
                    "ANSpikeTimes"
                ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage = np.array(tresults[j]["Vsoma"])
                dendriteVoltage = np.array(tresults[j]["Vdend"])
                stimWaveform = np.array(tresults[j]["stimWaveform"])
                stimTimebase = np.array(tresults[j]["stimTimebase"])
                # stim = np.array(tresults[j]["stim"])  # save the stimulus
                result[N] = {
                    "spikeTimes": spikeTimes,
                    "inputSpikeTimes": inputSpikeTimes,
                    "time": np.array(celltime[0]),
                    "somaVoltage": somaVoltage,
                    "dendriteVoltage": dendriteVoltage,
                    "allDendriteVoltages": allDendriteVoltages,
                    "stimWaveform": stimWaveform,
                    "stimTimebase": stimTimebase,
                }
                if self.Params.save_all_sections:  # save data for all sections
                    allDendriteVoltages[N] = tresults[j]["allsecVec"]
        self.print_timing()

        if self.Params.save_all_sections:
            for n in range(len(allDendriteVoltages)):
                for s in list(allDendriteVoltages[n].keys()):
                    allDendriteVoltages[n][s] = np.array(allDendriteVoltages[n][s])
        self.analysis_filewriter(self.Params.cell, result, tag="delays")
        if self.Params.plotFlag:
            print("plotting")
            self.plot_an(celltime, result)

    def an_run_singles(self, exclude=False, verify=False):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
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
        print("\n****AN Singles****\n")
        self.start_time = time.time()

        if self.Params.inputPattern is not None:
            fromdict = self.Params.inputPattern
            print("Cell id: %s  using input pattern: %s" % (self.cellID, fromdict))
        else:
            fromdict = self.cellID
        synapseConfig, celltype = self.cconfig.makeDict(fromdict)
        nReps = self.RunInfo.nReps
        threshold = self.RunInfo.threshold  # spike threshold, mV

        preCell, synapse, self.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig, celltype
        )

        nSyns = len(preCell)
        seeds = self.compute_seeds(nReps, synapseConfig)

        self.RunInfo.seeds = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = self.Params.Parallel
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.0
        # get the gMax's
        gMax = np.zeros(len(synapse))
        gMaxNMDA = np.zeros(len(synapse))

        for i, s in enumerate(synapse):
            for p in s.psd.ampa_psd:
                gMax[i] = p.gmax
            for p in s.psd.nmda_psd:
                gMaxNMDA[i] = p.gmax
        print("gmaxAmpa: ", gMax)
        print("gmaxNMDA: ", gMaxNMDA)
        result = {}
        for k in range(nSyns):
            # only enable gsyn on the selected input
            if exclude:
                tagname = "ExcludeSyn%03d" % k
                for i, s in enumerate(synapse):
                    for p in s.psd.ampa_psd:
                        if i == k:
                            p.gmax = 0.0  # disable the one
                        else:
                            p.gmax = gMax[i]  # leaving the others set
            else:
                tagname = "Syn%03d" % k
                for i, s in enumerate(synapse):
                    for p in s.psd.ampa_psd:
                        if i != k:
                            p.gmax = 0.0  # disable all
                            p.gmax = 0.0
                        else:
                            p.gmax = gMax[i]  # except the chosen one
                    for p in s.psd.nmda_psd:
                        if i != k:
                            p.gmax = 0.0  # disable all
                            p.gmax = 0.0
                        else:
                            p.gmax = gMaxNMDA[i]  # except the chosen one
            print("Syn #%d   synapse gMax: %f " % (k, gMax[k]))
            # for i, s in enumerate(synapse):
            #     for p in s.psd.ampa_psd:
            #         print(f"Synapse: {i:d}  gmax={p.gmax:.6f}")
            # print("tag: %s" % tagname)

            tresults = [None] * nReps

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
            spikeTimes = {}
            inputSpikeTimes = {}
            somaVoltage = {}
            dendriteVoltage = {}
            celltime = []
            stim = {}
            for j, N in enumerate(range(nReps)):
                celltime.append(tresults[j]["time"])  # (self.time)
                spikeTimes[N] = pu.findspikes(
                    tresults[j]["time"],
                    tresults[j]["Vsoma"],
                    threshold,
                    t0=0.0,
                    t1=self.RunInfo.run_duration * 1000,
                    dt=1.0,
                    mode="peak",
                )
                spikeTimes[N] = self.clean_spiketimes(spikeTimes[N])
                inputSpikeTimes[N] = tresults[j][
                    "ANSpikeTimes"
                ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]["Vsoma"])
                dendriteVoltage[N] = np.array(tresults[j]["Vdend"])
                stim[N] = np.array(tresults[j]["stim"])
            stimWaveform = np.array(tresults[0]["stim"])
            stimTimebase = np.array(tresults[0]["stimTimebase"])

            result[k] = {
                "spikeTimes": spikeTimes,
                "inputSpikeTimes": inputSpikeTimes,
                "somaVoltage": somaVoltage,
                "dendriteVoltage": dendriteVoltage,
                "time": np.array(tresults[j]["time"]),
                "stimWaveform": stimWaveform,
                "stimTimebase": stimTimebase,
            }

            self.analysis_filewriter(self.Params.cell, result[k], tag=tagname)
            self.print_timing()
        if self.Params.plotFlag:
            self.plot_an(celltime, result)

    def an_run_omit_one(self):
        self.an_run_singles(exclude=True)

    def an_run_IO(self):
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
        synapseConfig, celltype = self.cconfig.makeDict(self.cellID)
        nReps = self.RunInfo.nReps
        threshold = self.RunInfo.threshold  # spike threshold, mV

        preCell, synapse, self.R.electrode_site = self.configure_cell(
            self.post_cell, synapseConfig, celltype
        )
        seeds = self.compute_seeds()
        nSyns = len(preCell)
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
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
            gmaxs = np.zeros(nReps)
            tresults = [None] * nReps

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
            gmaxs = [
                4.0 * float(j + 1) * gMax[0] / float(nReps + 1) for j in range(nReps)
            ]
            for j, N in enumerate(range(nReps)):
                celltime.append(tresults[j]["time"])  # (self.time)
                spikeTimes[N] = pu.findspikes(
                    tresults[j]["time"],
                    tresults[j]["Vsoma"],
                    threshold,
                    t0=0.0,
                    t1=self.RunInfo.run_duration * 1000,
                    dt=1.0,
                    mode="peak",
                )
                spikeTimes[N] = self.clean_spiketimes(spikeTimes[N])
                inputSpikeTimes[N] = tresults[j][
                    "ANSpikeTimes"
                ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]["Vsoma"])
                dendriteVoltage[N] = np.array(tresults[j]["Vdend"])
            self.RunInfo.gSyns = gmaxs
            self.print_timing()

            result = {
                "spikeTimes": spikeTimes,
                "inputSpikeTimes": inputSpikeTimes,
                "somaVoltage": somaVoltage,
                "dendriteVoltage": dendriteVoltage,
                "time": np.array(tresults[j]["time"]),
            }

            self.analysis_filewriter(self.Params.cell, result, tag=tagname % k)
        if self.Params.plotFlag:
            self.plot_an(celltime, result)

    def configure_cell(self, thisCell: object, synapseConfig: dict, celltype: str):
        """
        Configure the cell. This routine builds the cell in NEURON, adds presynaptic inputs
        as described in the synapseConfig, and configures those according to parameters in
        self.Params and self.RunInfo.

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
        # reconfigure syanpses to set the spont rate group
        for i, syn in enumerate(synapseConfig):
            if self.Params.SRType == "fromcell":  # use the one in the table
                preCell.append(cells.DummySGC(cf=self.RunInfo.F0, sr=syn["SR"]))
                preCell[-1]._synapsetype = self.Params.ANSynapseType
                self.Params.SRType = self.Params.srnames[
                    syn[2]
                ]  # use and report value from table
            else:
                try:
                    srindex = self.Params.srnames.index(self.Params.SRType)
                    print(
                        f"Retrieved index {srindex:d} with SR type {self.Params.SRType:s}",
                        end=" ",
                    )
                except ValueError:
                    raise ValueError(
                        "SR type '%s' not found in Sr type list" % self.Params.SRType
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
                    print("  Synapsetype: simple")
                    synapse.append(
                        preCell[-1].connect(
                            thisCell,
                            type="simple",
                            post_opts={"postsec": pl, "postsite": postsite[0:2]},
                        )
                    )
                else:
                    print(f" Synapsetype: multisite, with {int(syn['nSyn']):d} zones")
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

        # try:
        #     cellInit.restore_initial_conditions_state(
        #         self.post_cell,
        #         electrode_site=None,
        #         filename=self.Params.initStateFile,
        #         autoinit=self.Params.auto_initialize,
        #     )
        #     # print('### SUCCESS ###')
        # except:
        #     cprint('r',
        #         "single_an_run: could not restore initial conditions: will try to create again"
        #     )
        #     self.an_run(make_an_intial_conditions=True)
        #     cprint('g', "Success making initial conditions #2")
        #     try:
        #         cellInit.restore_initial_conditions_state(
        #             self.post_cell,
        #             electrode_site=None,
        #             filename=self.Params.initStateFile,
        #         )
        #     except:
        #         raise ValueError("Failed initialization for cell: ", self.cellID)

        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.0
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        if isinstance(self.RunInfo.pip_start, list):
            pips = self.RunInfo.pip_start
        else:
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
        rtime = self.post_cell.hr.h.Vector()
        if (
            "dendrite" in self.post_cell.all_sections
            and len(self.post_cell.all_sections["dendrite"]) > 0
        ):
            dendsite = self.post_cell.all_sections["dendrite"][-1]
            Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        else:
            dendsite = None

        self.allsecVec = OrderedDict()
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
        Vsoma.record(self.post_cell.soma(0.5)._ref_v, sec=self.post_cell.soma)
        rtime.record(self.post_cell.hr.h._ref_t)
        self.post_cell.hr.h.finitialize()
        self.post_cell.hr.h.tstop = self.RunInfo.run_duration * 1000.0
        self.post_cell.hr.h.t = 0.0
        self.post_cell.hr.h.batch_save()  # save nothing
        self.post_cell.hr.h.dt = self.Params.dt
        self.post_cell.hr.h.batch_run(
            self.post_cell.hr.h.tstop, self.post_cell.hr.h.dt, "an.dat"
        )
        nrn_run_time += time.time() - nrn_start
        if dendsite is None:
            Vdend = np.zeros_like(Vsoma)

        anresult = {
            "Vsoma": np.array(Vsoma),
            "Vdend": np.array(Vdend),
            "time": np.array(rtime),
            "ANSpikeTimes": ANSpikeTimes,
            "stim": stim,
            "stimWaveform": stimWaveform,
            "stimTimebase": stimTimebase,
        }
        if self.Params.save_all_sections:
            anresult["allsecVec"] = self.allsecVec
        return anresult

    def clean_spiketimes(self, spikeTimes, mindT=0.7):
        """
        Clean up spike time array, removing all less than mindT
        spikeTimes is a 1-D list or array
        mindT is difference in time, same units as spikeTimes
        If 1 or 0 spikes in array, just return the array
        """
        if len(spikeTimes) > 1:
            dst = np.diff(spikeTimes)
            st = np.array(spikeTimes[0])  # get first spike
            sok = np.where(dst > mindT)
            st = np.append(st, [spikeTimes[s + 1] for s in sok])
            spikeTimes = st
        return spikeTimes

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
        fig, ax = mpl.subplots(2, 1)
        fig.suptitle("AN Inputs")
        # win = pgh.figure(title='AN Inputs')
        # layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
        for j, N in enumerate(range(len(result))):
            ax[0].plot(celltime[N], result[N]["somaVoltage"], c="k", linewidth=0.75)
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

    def analysis_filewriter(self, filebase, result, tag=""):
        """
        Write the analysis information to a pickled file

        Parameters:
        filebase : string (no default)
            base filename - *not used* (value replaced bye cellID)
        result : dict (no default)
            dict hlding results. Must be pickleable
        tag : string (default: '')
            tag to insert in filename string
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
        # print(results)
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

    def get_hoc_file(self, hf):
        if hf.file_loaded is False:
            exit()
        self.section_list = hf.get_section_prefixes()
        list(hf.sec_groups.keys())
        if len(hf.sec_groups) > 1:  # multiple names, so assign colors to structure type
            self.section_colors = {}
            if self.Params.verbose:
                print("gethocfile # colors: ", len(self.hg.colorMap))
                print("gethocfile # sexcgroups: ", len(list(hf.sec_groups.keys())))
            for i, s in enumerate(hf.sec_groups.keys()):
                if self.Params.verbose:
                    print("gethocfile group: ", s)
                self.section_colors[s] = self.hg.colorMap[i]
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

    (
        parsedargs,
        params,
        runinfo,
    ) = vcnmodel.model_params.getCommands()  # get from command line
    model = ModelRun(
        params=params, runinfo=runinfo, args=parsedargs
    )  # create instance of the model
    if parsedargs.displayMode != "None":
        model.view_model()
    else:
        model.run_model()


if __name__ == "__main__":
    main()
