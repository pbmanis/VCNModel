#!/usr/bin/env python


__author__ = "pbmanis"
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
                 [--model {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}]
                 [--modeltype {II,II-I,I-II,I-c,I-t,II-o}]
                 [--sgcmodel {Zilany,cochlea}]
                 [--protocol {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}]
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
  --model {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}, -M {XM13,XM13_nacncoop,XM13_nacn,XM13_nabu,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}
                        Define the model type (default: XM13)
  --modeltype {II,II-I,I-II,I-c,I-t,II-o}
                        Define the model type (default: XM13)
  --sgcmodel {Zilany,cochlea}
                        Define the SGC model type (default: Zilany)
  --protocol {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}, -P {initIV,testIV,runIV,initandrunIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}
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

import sys
from pathlib import Path
import os
import errno
import copy
import pickle
import time
import argparse
from collections import OrderedDict
import pprint
import json
import toml
import numpy as np
import timeit
import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as mpl

# from neuronvis.hoc_viewer import HocViewer

import vcnmodel.cell_config as cell_config
import neuronvis.hoc_graphics as hoc_graphics
from vcnmodel.generate_run import GenerateRun
import vcnmodel.cellInitialization as cellInit
from vcnmodel.adjust_areas import AdjustAreas
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.decorator import Decorator
from cnmodel import data as DATA
import pylibrary.tools.utility as pu  # access to a spike finder routine
import pylibrary.tools.cprint as CP
import pyqtgraph as pg
import pyqtgraph.multiprocess as MP

# import pylibrary.plotting.pyqtgraph_plothelpers as pgh


showCell = True

cprint = CP.cprint


class ModelRun:
    def __init__(self, args:[object, None]=None):

        self.init_flags()
        self.init_parameters()

    def init_parameters(self):
        # The following are initial values. Some values of some of these are replaced after
        # parsing the command line.
        # See the main code section (at the end of this file).
        self.Params = OrderedDict()
        self.Params["cell"] = self.cellID
        self.Params["AMPAScale"] = 1.0  # Use the default scale for AMPAR conductances
        self.Params["ANSynapseType"] = "simple"  # or multisite
        self.Params[
            "ANSynapticDepression"
        ] = 0  # depression calculation is off by default
        self.Params["initIVStateFile"] = None  # 'IVneuronState_%s.dat'
        self.Params["initANStateFile"] = None  # 'ANneuronState_%s.dat'
        self.Params["simulationFilename"] = None
        self.Params["hocfile"] = None
        self.Params["fullhocfile"] = False
        self.Params["usedefaulthoc"] = False
        self.Params["cellType"] = self.cellChoices[0]
        self.Params["modelName"] = self.modelNameChoices[0]
        self.Params["modelType"] = self.modelTypeChoices[0]
        self.Params["SGCmodelType"] = self.SGCmodelChoices[0]
        self.Params["species"] = self.speciesChoices[0]
        self.Params["Ra"] = 150.0  # ohm.cm
        self.Params["dt"] = 0.025  # dt integration rate, ms
        self.Params["soma_inflation"] = 1.0  # factor to multiply soma section areas by
        self.Params["soma_autoinflate"] = False  #
        self.Params["dendrite_inflation"] = 1.0
        self.Params["dendrite_autoinflate"] = False
        self.Params["dendrite_fromsoma"] = False
        self.Params["ASA_inflation"] = 1.0
        self.Params["ASA_fromsoma"] = False
        self.Params["lambdaFreq"] = 2000.0  # Hz for segment number
        self.Params["dendriteelectrode"] = "Proximal_Dendrite"
        self.Params[
            "sequence"
        ] = ""  # sequence for run - may be string [start, stop, step]
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params["SRType"] = self.SRChoices[2]
        # self.Params['SR'] = self.Params['SRType']  # actually used SR: this might be cell-defined, rather than entirely specified from the command line
        self.Params[
            "inputPattern"
        ] = None  # ID of cellinput pattern (same as cellID): for substitute input patterns.
        self.Params["spirou"] = "all"
        self.Params["runProtocol"] = self.protocolChoices[
            2
        ]  # testIV is default because it is fast and should be run often
        self.Params["nReps"] = 1
        self.Params["seed"] = 100  # always the same start - avoids lots of recomutation
        # in production, vary this or set to None for random values
        self.Params[
            "initialization_time"
        ] = 50.0  # nominal time to let system settle, in msec
        self.Params["run_duration"] = 0.25  # in sec
        self.Params["soundtype"] = "SAM"  # or 'tonepip'
        self.Params["pip_duration"] = 0.1  # duration in seconds for tone pip
        self.Params["pip_start"] = [0.1]  # start (delay to start of pip)
        self.Params["pip_offduration"] = 0.05  # time after pip ends to keep running
        self.Params["Fs"] = 100e3  # cochlea/zilany model rate
        self.Params["F0"] = 16000.0  # stimulus frequency
        self.Params["dB"] = 30.0  # in SPL
        self.Params["RF"] = 2.5e-3  # rise-fall time
        self.Params["fmod"] = 20  # hz, modulation if SAM
        self.Params["dmod"] = 0  # percent if SAM
        self.Params["threshold"] = -35
        self.Params["signalToMasker"] = 0.0
        self.Params["CMMRmode"] = "CM"
        self.Params["all_modes"] = False
        # parameters for generating a noise signal to generating GIF model of cell
        self.Params["gif_i0"] = 0.0  # base current level
        self.Params["gif_sigma"] = 0.5  # std of noise
        self.Params["gif_fmod"] = 0.2  # mod in Hz
        self.Params["gif_tau"] = 3.0  # tau, msec
        self.Params["gif_dur"] = 10.0  # seconds
        self.Params["gif_skew"] = 0.0  # a in scipy.skew

        # general control parameters
        self.Params["plotFlag"] = False
        self.Params["auto_initialize"] = False
        self.Params["nWorkers"] = 8
        self.Params["Parallel"] = True
        self.Params["verbose"] = False
        self.Params["save_all_sections"] = False
        self.Params["commandline"] = ""  # store command line on run
        self.Params["checkcommand"] = False
        self.Params["configfile"] = None
        self.Params["tagstring"] = ""

    def init_flags(self):
        # use v2 files for model with rescaled soma
        self.cellChoices = ["Bushy", "TStellate", "DStellate"]
        self.modelNameChoices = [
            "XM13",
            "XM13_nacncoop",
            "XM13_nacn",
            "XM13_nabu",
            "RM03",
            "mGBC",
            "XM13PasDend",
            "Calyx",
            "MNTB",
            "L23Pyr",
        ]
        self.modelTypeChoices = ["II", "II-I", "I-II", "I-c", "I-t", "II-o"]
        self.SGCmodelChoices = [
            "Zilany",
            "cochlea",
        ]  # cochlea is python model of Zilany data, no matlab, JIT computation; Zilany model creates matlab instance for every run.
        self.cmmrModeChoices = ["CM", "CD", "REF"]  # comodulated, codeviant, reference
        self.SRChoices = [
            "LS",
            "MS",
            "HS",
            "fromcell",
        ]  # AN SR groups (assigned across all inputs)
        self.protocolChoices = [
            "initIV",
            "testIV",
            "runIV",
            "initandrunIV",
            "initAN",
            "runANPSTH",
            "runANIO",
            "runANSingles",
            "runANOmitOne",
            "gifnoise",
        ]
        self.soundChoices = ["tonepip", "noise", "stationaryNoise", "SAM", "CMMR"]
        self.speciesChoices = ["mouse", "guineapig"]
        self.spirouChoices = ["all", "max=mean", "all=mean"]
        self.ANSynapseChoices = ["simple", "multisite"]

        self.cellmap = {
            "bushy": cells.bushy,
            "tstellate": cells.tstellate,
            "dstellate": cells.dstellate,
        }
        self.srname = ["LS", "MS", "HS"]  # runs 0-2, not starting at 0
        self.cellID = (
            None  # ID of cell (string, corresponds to directory name under VCN_Cells)
        )
        self.setup = (
            False  # require setup ahead of run - but setup can be done separately
        )
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            datapaths = toml.load("wheres_my_data.toml")
        else:
            datapaths = {"baseDirectory": Path(
            "../VCN-SBEM-Data", "VCN_Cells")
        }  # "here"
        self.baseDirectory = datapaths["baseDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"

    def print_modelsetup(self):
        """
        Print out all of the parameters in the model
        """
        for p in list(self.Params.keys()):
            print("{:18s} = {:12}".format(p, self.Params[p]))
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
        self.Params["cellType"] = cellType

    def set_model_name(self, model_name: str):
        """
        Set the model Type, as requested. The model type must be in the model choices

        Parameters
        ----------
        model_type : string
            The type of the model that will be used (condutance settings)

        Returns
        -------
            Nothing

        """
        if model_name not in self.modfelNameChoices:
            print(
                "Model type must be one of: {:s}. Got: {:s}".format(
                    ", ".join(self.modelNameChoices), model_type
                )
            )
            exit()
        self.Params["modelName"] = model_name

    def _make_filenames(self):
        """
        Define program-wide file names (used also in cell_initialization and generate_run) one time for consistencye

        This routine generates two names:
            1. The name of the initizlization file. This name includes the model name, the model type (for that name),
                soma and dendrite inflation factors if relevenat. Other parameters can be added if needed
            2. The name of the simulation data file. This name is similar to the initizlizaiton name, but may include
                information about the type of run, stimuli, etc.


        """
        self.cellID = Path(
            self.Params["cell"]
        ).stem  # os.path.splitext(self.Params['cell'])[0]
        # pick up parameters that should be in both init and run filenames:
        # Run / init type independent parameters:
        namePars = (
            f"{str(self.Params['modelName']):s}_{str(self.Params['modelType']):s}"
        )
        print(
            "inflations: ",
            self.Params["soma_inflation"],
            self.Params["dendrite_inflation"],
        )

        if self.Params["soma_inflation"] != 1.0:
            namePars += f"_soma={self.Params['soma_inflation']:.3f}"
        if self.Params["dendrite_inflation"] != 1.0:
            namePars += f"_dend={self.Params['dendrite_inflation']:.3f}"
        if self.Params["ASA_inflation"] != 1.0:
            namePars += f"_ASA={self.Params['ASA_inflation']:.3f}"

        if self.Params["runProtocol"] in ["initIV", "initandrunIV", "runIV"]:
            if self.Params["initIVStateFile"] is None:
                fn = f"IVneuronState_{namePars:s}"
                ivinitdir = Path(self.baseDirectory, self.cellID, self.initDirectory)
                if self.Params["tagstring"] is not None:
                    self.Params["initIVStateFile"] = Path(
                        ivinitdir, f"{str(fn):s}_{self.Params['tagstring']:s}"
                    )
                else:
                    self.Params["initIVStateFile"] = Path(ivinitdir, fn)
            if not self.Params["initIVStateFile"].suffix == ".dat":
                self.Params["initIVStateFile"] = Path(
                    str(self.Params["initIVStateFile"]) + ".dat"
                )
            cprint(
                "cyan",
                "IV Initialization file:  {str(self.Params['initIVStateFile']):s}",
            )
            self.mkdir_p(ivinitdir)  # confirm existence of that file

        if self.Params["runProtocol"] in ["initandrunIV", "runIV"]:
            outPath = Path(self.baseDirectory, self.cellID, self.simDirectory, "IV")
            self.mkdir_p(outPath)  # confirm that output path exists
            fout = self.Params[
                "simulationFilename"
            ]  # base name created by make-filename - do not overwrite
            print("fout: ", fout)
            if self.Params["tagstring"] is not None:
                self.Params["simulationFilename"] = Path(
                    outPath,
                    f"{self.cellID:s}_pulse_{namePars:s}_monitor_{self.Params['tagstring']:s}.p",
                )
            else:
                self.Params["simulationFilename"] = Path(
                    outPath, f"{self.cellID:s}_pulse_{namePars:s}_monitor.p"
                )
            print("Simulation filename: ", self.Params["simulationFilename"])

        if (
            self.Params["runProtocol"].startswith("runAN")
            or self.Params["runProtocol"] == "initAN"
            or self.Params["runProtocol"] == "runANSingles"
        ):
            if self.Params["inputPattern"] is None:
                inputs = "self"
            else:
                inputs = self.Params["inputPattern"]
            if self.Params["initANStateFile"] is None:
                fn = f"ANneuronState_{namePars:s}_inp={inputs:s}_{self.Params['ANSynapseType']:s}.dat"
                aninitdir = Path(self.baseDirectory, self.cellID, self.initDirectory)
                self.Params["initANStateFile"] = Path(aninitdir, fn)
                self.mkdir_p(aninitdir)  # confirm existence of that file
                print("AN Initialization file: ", self.Params["initANStateFile"])

            outPath = Path(self.baseDirectory, self.cellID, self.simDirectory, "AN")
            self.mkdir_p(outPath)  # confirm that output path exists
            addarg = namePars
            if self.Params["spirou"] == "all":
                addarg += "_all"
            elif self.Params["spirou"] == "max=mean":
                addarg = "_mean"
            elif self.Params["spirou"] == "all=mean":
                addarg = "_allmean"
            print("soundtype: ", self.Params["soundtype"])
            fn = f"AN_Result_{self.cellID:s}_inp={inputs:s}_{addarg:s}_{self.Params['ANSynapseType']:s}"
            fn += f"_{self.Params['nReps']:03d}_{self.Params['soundtype']:s}"
            fn += f"_{int(self.Params['dB']):03d}dB_{self.Params['F0']:06.1f}"
            if self.Params["soundtype"] in ["SAM", "sam"]:
                fn += f"_{self.Params['fmod']:03.1f}_{int(self.Params['dmod']):03d}"
                fn += f"_{self.Params['SRType']:2s}.p"
                ofile = Path(outPath, fn)
            else:
                fn += f"_{self.Params['SRType']:2s}.p"
                ofile = Path(outPath, fn)
            self.Params["simulationFilename"] = ofile
            # if not self.Params["simulationFilename"].suffix == ".dat":
#                 self.Params["simulationFilename"] = Path(
#                     str(self.Params["simulationFilename"]) + ".dat"
#                 )

        cprint(
            "r", f"Output simulation file: {str(self.Params['simulationFilename']):s}"
        )
        # exit()

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
        self.Params["SRType"] = spont_rate_type

    def set_starttime(self, starttime: float):
        """
        store the start time...
        """
        self.Params["StartTime"] = starttime

    def mkdir_p(self, path: [str, Path]):
        try:
            Path.mkdir(path, parents=True, exist_ok=True)
        except:
            raise FileNotFoundError(f"Cannot create path: {str(path):s}")

    def fix_singlets(self, h: object):
        badsecs = []
        for i, sec in enumerate(h.allsec()):
            if sec.n3d() == 1:
                badsecs.append(sec)
                # print(f'Fixing Singlet Section: {str(sec):s}')
                # print(dir(sec))
                parent = sec.parentseg()
                if (
                    parent is None and i > 0
                ):  # first section might not have a parent... obviously
                    raise ValueError(
                        f"Section: {str(sec):s} does not have a parent, and this must be fixed manually"
                    )
                if i == 0 and parent is None:
                    continue
                # print(f'    Section connected to {str(parent):s} ')
                nparentseg = parent.sec.n3d() - 1
                x = parent.sec.x3d(nparentseg)
                y = parent.sec.y3d(nparentseg)
                z = parent.sec.z3d(nparentseg)
                d = parent.sec.diam3d(nparentseg)
                h.pt3dinsert(0, x, y, z, d, sec=sec)
        if len(badsecs) > 0:
            print(f"**** Fixed Singlet segments in {len(badsecs):d} sections")
        badsecs = []
        for sec in h.allsec():
            if sec.n3d() == 1:
                badsecs.append(sec)
        # if len(badsecs) == 0:
        #     print('no more bad sections')

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

        if self.Params["verbose"]:
            print(f"run_model entry")
        if par_map is not None and "id" in list(par_map.keys()):
            self.idnum = par_map["id"]
        else:
            self.idnum = 9999
        print(f"Cell ID: {self.Params['cell']:s}")
        if self.Params["cell"] is None:
            return
        self.cellID = Path(
            self.Params["cell"]
        ).stem  # os.path.splitext(self.Params['cell'])[0]

        print(f"Morphology directory: {self.morphDirectory:s}")
        if self.Params["usedefaulthoc"] or self.Params["hocfile"] == None:
            self.Params["hocfile"] = self.Params["cell"] + ".hoc"
            cprint("red", "Using default hoc file")
        if self.Params["usedefaulthoc"] and self.Params["fullhocfile"]:
            cprint("red", "Using _full hoc file")
            self.Params["hocfile"] = self.Params["cell"] + "_Full.hoc"

        filename = Path(
            self.baseDirectory, self.cellID, self.morphDirectory, self.Params["hocfile"]
        )
        print(f"Name of hoc file: {self.Params['hocfile']:s}")
        print(f"base hoc file: {str(filename):s}")
        print(f"hocx file: {str(filename.with_suffix('.hocx')):s}")
        if filename.with_suffix(".hocx").is_file():  # change to preferred if available
            filename = filename.with_suffix(".hocx")
        cprint("yellow", f"Using hoc file: {str(filename):s}")

        # instantiate cells
        if self.Params["cellType"] in ["Bushy", "bushy"]:
            print(f"Creating a bushy cell (run_model) ")
            from cnmodel import data

            changes = None
            nach = None  # uses default
            if self.Params["modelName"] == "XM13_nacncoop":
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
            elif self.Params["modelName"] == "XM13":
                import model_data.data_XM13 as CHAN

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
            elif self.Params["modelName"] == "XM13_nacn":
                import model_data.data_XM13_nacn as CHAN

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
            elif self.Params["modelName"] == "XM13_nabu":
                import model_data.data_XM13_nabu as CHAN

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
                species=self.Params["species"],
                modelName=self.Params["modelName"],
                modelType=self.Params["modelType"],
                nach=nach,
            )
        elif self.Params["cellType"] in ["tstellate", "TStellate"]:
            print(f"Creating a t-stellate cell (run_model) ")
            self.post_cell = cells.TStellate.create(
                morphology=str(filename),
                decorator=Decorator,
                species=self.Params["species"],
                modelType=self.Params["modelName"],
            )
        elif self.Params["cellType"] in ["dstellate", "DStellate"]:
            print(f"Creating a D-stellate cell (run_model)")
            self.post_cell = cells.DStellate.create(
                morphology=str(filename),
                decorator=Decorator,
                species=self.Params["species"],
                modelType=self.Params["modelName"],
            )
        else:
            raise ValueError(f"cell type {self.Params['cellType']:s} not implemented")

        AdjA = AdjustAreas()
        AdjA.sethoc_fromCNcell(self.post_cell)
        hoc_somaarea = AdjA.get_hoc_area(["soma"])
        hoc_dendritearea = AdjA.get_hoc_area(dendrite_names)
        cprint(
            "y",
            f"HOC: Soma area: {hoc_somaarea:.2f}  Dendrite Area: {hoc_dendritearea:.2f}",
        )

        # Set up run parameters
        print(
            f"Requested temperature (deg C): {self.post_cell.status['temperature']:.2f}"
        )
        self.post_cell.hr.h.celsius = self.post_cell.status[
            "temperature"
        ]  # this is set by prepareRun in generateRun. Only place it should be changed
        self.post_cell.hr.h.Ra = self.Params["Ra"]
        print("Ra (ohm.cm) = {:8.1f}".format(self.post_cell.hr.h.Ra))
        print(f"Specified Temperature = {self.post_cell.hr.h.celsius:8.1f} degC ")

        # for sec in self.post_cell.hr.h.allsec():
        #     print ('self.post_cell.hr.h: ', sec)
        self.fix_singlets(self.post_cell.hr.h)
        if self.Params[
            "soma_autoinflate"
        ]:  # get values and inflate soma automatically to match mesh
            cprint("c", "Soma Autoinflation")
            cconfig = cell_config.CellConfig()
            inflateratio = cconfig.get_soma_ratio(self.cellID)
            if np.isnan(inflateratio):
                raise ValueError(
                    f"Soma Inflation Ratio is not defined for cell {self.cellID:s}"
                )
            self.Params["soma_inflation"] = inflateratio

        if self.Params["soma_inflation"] != 1.0:
            print("!!!!!   Inflating soma")
            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -60.0]
            )
            print(
                f"     Original Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}"
            )
            # origdiam = {}
            AdjA = AdjustAreas(method="pt3d")
            AdjA.sethoc_fromCNcell(self.post_cell)
            pt3d = AdjA.adjust_diameters(
                sectypes=["soma", "Soma"], inflateRatio=inflateratio
            )
            # AdjA.plot_areas(pt3d)

            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -60.0]
            )
            print(
                f"     New Rin after somatic inflation: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}"
            )
        if self.Params["ASA_fromsoma"]:
            self.Params["ASA_inflation"] = self.Params["soma_inflation"]

        if self.Params["dendrite_fromsoma"]:
            self.Params["dendrite_inflation"] = self.Params["soma_inflation"]
        elif self.Params[
            "dendrite_autoinflate"
        ]:  # get values and inflate soma automatically to match mesh
            cprint("c", "Dendrite Autoinflation")
            cconfig = cell_config.CellConfig()
            inflateratio = cconfig.get_dendrite_ratio(self.cellID)
            if np.isnan(inflateratio):
                raise ValueError("Dendrite Inflation Ration is not defined!")
            self.Params["dendrite_inflation"] = inflateratio

        if self.Params["dendrite_inflation"] != 1.0:
            print("!!!!!   Inflating dendrite")
            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -55.0]
            )
            print(
                f"     Original Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e6:.2f}, RMP: {rtau['v']:.2f}"
            )

            AdjA = AdjustAreas()
            AdjA.sethoc_fromCNcell(self.post_cell)
            pt3d = AdjA.adjust_diameters(
                sectypes=dendrite_names, inflateRatio=inflateratio
            )
            # AdjA.plot_areas(pt3d)  # just to take a look at the adjustment

            rtau = self.post_cell.compute_rmrintau(
                auto_initialize=True, vrange=[-80.0, -55.0]
            )
            print(
                f"     New Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e6:.2f}, RMP: {rtau['v']:.2f}"
            )
            # print(f"     New Rin after dendrite inflation: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e6:.2f}, RMP: {rtau['v']:.2f}")
            # for section in list(dendrite_group):
            # print('      diam orig: ', origdiam[section])
            # print('      diam new:  ', ssection.diam)
            # for sec in badsecs:
        #      print("? bad diam: ", sec.diam, sec)
        #  print("# bad: ", len(badsecs))
        #  if len(badsecs) > 0:
        #      exit()

        for group in list(self.post_cell.hr.sec_groups.keys()):
            g = self.post_cell.hr.sec_groups[group]
            for section in list(g):
                self.post_cell.hr.get_section(section).Ra = self.post_cell.hr.h.Ra
                if self.Params["verbose"]:
                    print("Section: ", section)
                    print("Ra: ", self.post_cell.hr.get_section(section).Ra)

        electrode_section = list(self.post_cell.hr.sec_groups["soma"])[0]
        self.electrode_site = self.post_cell.hr.get_section(electrode_section)
        self.electrodeSection = "soma"
        self.hg = hoc_graphics
        self.get_hoc_file(self.post_cell.hr)

        try:
            dendritic_electrode = list(
                self.post_cell.hr.sec_groups[self.Params["dendriteelectrodedendrite"]]
            )[0]
            self.dendriticelectrode_site = self.post_cell.hr.get_section(
                dendritic_electrode
            )
            self.dendriticElectrodeSection = self.Params["dendriteelectrodedendrite"]
        except:
            dendritic_electrode = list(self.post_cell.hr.sec_groups["soma"])[0]
            self.dendriticelectrode_site = self.post_cell.hr.get_section(
                dendritic_electrode
            )
            self.dendriticElectrodeSection = "soma"

        if self.Params["verbose"]:
            if par_map is not None:
                print("Listing par_map (run_model): ", par_map)
            self.post_cell.hr.h.topology()

        self.post_cell.set_nseg(freq=self.Params["lambdaFreq"])

        # handle the following protocols:
        # ['initIV', 'initAN', 'runIV', 'run', 'runANSingles', 'gifnoise']

        self._make_filenames()  # make filenames AFTER all manipulations of the cell

        self.setup = True

    def run_model(self, par_map: dict = None):
        if not self.setup:
            self.setup_model(par_map=par_map)
        if self.Params["runProtocol"] in ["runANPSTH", "runANSingles"]:
            self.Params["run_duration"] = (
                self.Params["pip_duration"]
                + np.sum(self.Params["pip_start"])
                + self.Params["pip_offduration"]
            )

        if self.Params["runProtocol"] in ["initIV", "initandrunIV"]:
            if self.Params["verbose"]:
                print("run_model: protocol is initIV")
            # ivinitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(
                self.post_cell,
                idnum=self.idnum,
                celltype=self.Params["cellType"],
                starttime=None,
                electrodeSection=self.electrodeSection,
                dendriticElectrodeSection=self.dendriticElectrodeSection,
                iRange=self.post_cell.i_test_range,
                plotting=self.Params["plotFlag"],
                params=self.Params,
            )
            cellInit.get_initial_condition_state(
                self.post_cell,
                tdur=self.Params["initialization_time"],
                filename=self.Params["initIVStateFile"],
                electrode_site=self.electrode_site,
            )
            print(f"Ran to get initial state for {self.post_cell.hr.h.t:.1f} msec")

        if self.Params["runProtocol"] in ["runIV", "initandrunIV"]:
            if self.Params["verbose"]:
                print("Run IV")
            self.iv_run()

        if self.Params["runProtocol"] == "testIV":
            if self.Params["verbose"]:
                print("test_init")
            # ivinitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(
                self.post_cell,
                idnum=self.idnum,
                celltype=self.Params["cellType"],
                starttime=None,
                electrodeSection=self.electrodeSection,
                dendriticElectrodeSection=self.dendriticElectrodeSection,
                iRange=self.post_cell.irange,
                plotting=self.Params["plotFlag"],
                params=self.params,
            )
            cellInit.test_initial_conditions(
                self.post_cell,
                filename=self.Params["initIVStateFile"],
                electrode_site=self.electrode_site,
            )
            self.R.testRun(initfile=self.Params["initIVStateFile"])
            return  # that is ALL, never make testIV/init and then keep running.

        if self.Params["runProtocol"] == "runANPSTH":
            if self.Params["verbose"]:
                print("ANPSTH")
            self.an_run(self.post_cell)

        if self.Params["runProtocol"] == "initAN":
            if self.Params["verbose"]:
                print("Init AN")
            self.an_run(self.post_cell, make_an_intial_conditions=True)

        if self.Params["runProtocol"] == "runANIO":
            if self.Params["verbose"]:
                print("Run AN IO")
            self.an_run_IO(self.post_cell)

        if self.Params["runProtocol"] == "runANSingles":
            if self.Params["verbose"]:
                print("ANSingles")
            self.an_run_singles(self.post_cell)

        if self.Params["runProtocol"] == "runANOmitOne":
            if self.Params["verbose"]:
                print("ANOmitOne")
            self.an_run_singles(self.post_cell, exclude=True)

        if self.Params["runProtocol"] == "gifnoise":
            self.noise_run()

        # if showCell:
        #     print(dir(self.post_cell))
        #     self.render = HocViewer(self.post_cell.hr)
        #     cylinder=self.render.draw_cylinders()
        #     cylinder.set_group_colors(self.section_colors, alpha=0.8, mechanism=['nav11', 'gbar'])
        #     pg.show()

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
        if self.Params["verbose"]:
            print("iv_run: calling generateRun", self.post_cell.i_test_range)
        # parse i_test_range and pass it here
        if self.Params["sequence"] is "":
            if isinstance(
                self.post_cell.i_test_range, dict
            ):  # not defined, use default for the cell type
                iinjValues = self.post_cell.i_test_range
            else:  # if it was a list, try to put the list into the pulse range
                iinjValues = {"pulse": self.post_cell.i_test_range}
        else:
            iinjValues = {"pulse": eval(self.Params["sequence"])}
        self.R = GenerateRun(
            self.post_cell,
            idnum=self.idnum,
            celltype=self.Params["cellType"],
            starttime=None,
            electrodeSection=self.electrodeSection,
            dendriticElectrodeSection=self.dendriticElectrodeSection,
            iRange=iinjValues,
            plotting=self.Params["plotFlag"],
            params=self.Params,
        )
        self.R.runInfo.folder = Path(
            self.baseDirectory, self.cellID, self.simDirectory, "IV"
        )
        if self.Params["verbose"]:
            print("iv_run: calling do_run")
        nworkers = self.Params["nWorkers"]
        #        print(self.Params['Parallel'])
        if self.Params["Parallel"] == False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(
            self.Params["hocfile"],
            parMap=iinjValues,
            save="monitor",
            restore_from_file=True,
            initfile=self.Params["initIVStateFile"],
            workers=nworkers,
        )
        if self.Params["verbose"]:
            print("   iv_run: do_run completed")
        elapsed = timeit.default_timer() - start_time
        print(f"   iv_rin: Elapsed time: {elapsed:2f} seconds")
        isteps = self.R.IVResult["I"]
        if self.Params["verbose"]:
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
            "sequence": iinjValues,
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
            self.post_cell,
            idnum=self.idnum,
            celltype=self.Params["cellType"],
            starttime=None,
            electrodeSection=self.electrodeSection,
            dendriticElectrodeSection=self.dendriticElectrodeSection,
            stimtype="gifnoise",
            plotting=self.Params["plotFlag"],
            params=self.Params,
        )
        # ivinitfile = Path(self.baseDirectory, self.cellID,
        #                         self.initDirectory, self.Params['initIVStateFile'])
        self.R.runInfo.folder = Path(
            self.baseDirectory, self.cellID, self.simDirectory, "Noise"
        )
        if self.Params["verbose"]:
            print("noise_run: calling do_run")
        nworkers = self.Params["nWorkers"]
        #        print(self.Params['Parallel'])
        if self.Params["Parallel"] == False:
            nworkers = 1
        #        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(
            self.Params["hocfile"],
            parMap=par_map,
            save="monitor",
            restore_from_file=True,
            initfile=self.Params["initIVStateFile"],
            workers=nworkers,
        )
        if self.Params["verbose"]:
            print("   noise_run: do_run completed")

    def check_for_an_statefile(self):
        print("State file name: ", self.Params["initANStateFile"])
        print("Cell: ", self.Params["cell"])
        print("Base Directory: ", self.baseDirectory)
        print("Initialization Directory: ", self.initDirectory)
        statefile = Path(
            self.baseDirectory,
            self.Params["cell"],
            self.initDirectory,
            self.Params["initANStateFile"],
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
        if self.Params["seed"] is None:  # every run is randomized
            seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        else:
            startingseed = self.Params["seed"]
            seeds = np.arange(0, nReps * len(synapseConfig)).reshape(
                (nReps, len(synapseConfig))
            )
            seeds = seeds + startingseed
        # print('AN Seeds: ', seeds)

        return seeds

    def get_synapses(self, synapses: list):
        """
        Get values of settings in synapses
        """
        gMax = np.zeros(len(synapses))  # total g for each synapse
        nSyn = np.zeros(len(synapses))  # # sites each synapse
        ngMax = np.zeros(len(synapses))
        # print('# syn: ', len(synapses))
        if self.Params["ANSynapseType"] == "simple":
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
                # print('i nsyn gmax nmdagmax  : ', i, nSyn[i], gMax[i], ngMax[i])
        return (gMax, ngMax, nSyn)

    def set_synapse(self, synapse: list, gampa: float, gnmda: float):
        totalg = 0.0
        if self.Params["ANSynapseType"] == "simple":
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
        # print('total ampa, nmda: ', totalg, totaln)

    def an_run(
        self,
        post_cell: object,
        verify: bool = False,
        make_an_intial_conditions: bool = False,
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
        post_cell : Cell object
            Provides access to cell class, as well as neuron and hoc file information

        verify : boolean (default: False)
            Flag to control printing of various intermediate results

        make_ANIntitialConditions : bool (default : False)
            Flag to control whether the initial conditions need to be recomputed and stored.

        Returns
        -------
            Nothing

        """
        print("\n*** an_run\n")

        if self.Params["inputPattern"] is not None:
            fromdict = self.Params["inputPattern"]
            print("Cell id: %s  using input pattern: %s" % (self.cellID, fromdict))
        else:
            fromdict = self.cellID
        cconfig = cell_config.CellConfig()
        synapseConfig, celltype = cconfig.makeDict(
            fromdict, self.Params["ASA_inflation"]
        )
        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        # for i, s in enumerate(synapseConfig):
        #     s[1] = s[3]*0.001/3.0
        #     print 'delay for input %d is %8.4f msec' % (i, s[1])

        nReps = self.Params["nReps"]
        threshold = self.Params["threshold"]  # spike threshold, mV

        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(
            post_cell, synapseConfig, celltype, stimInfo
        )

        # see if we need to save the cell state now.
        if make_an_intial_conditions:
            #            print('getting initial conditions for AN')
            # aninitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initANStateFile'])
            cellInit.get_initial_condition_state(
                post_cell,
                tdur=self.Params["initialization_time"],
                filename=self.Params["initANStateFile"],
                electrode_site=self.electrode_site,
                reinit=self.Params["auto_initialize"],
            )
            print("TESTING; statefile = ", self.Params["initANStateFile"])
            cellInit.test_initial_conditions(
                post_cell,
                filename=self.Params["initANStateFile"],
                electrode_site=self.electrode_site,
            )
            return

        seeds = self.compute_seeds(nReps, synapseConfig)
        stimInfo["seeds"] = seeds  # keep the seed values too.
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

        nWorkers = self.Params["nWorkers"]
        TASKS = [s for s in range(nReps)]
        tresults = [None] * len(TASKS)
        result = {}

        # manipulate syanptic strengths to test hypotheses here.
        allsf = self.Params["AMPAScale"]
        gMax, ngMax, nSyn = self.get_synapses(synapse)
        for i in range(len(nSyn)):
            gMax[i] = gMax[i] * allsf
            ngMax[i] = ngMax[i] * allsf
            self.set_synapse(synapse[i], gMax[i] / nSyn[i], ngMax[i] / nSyn[i])

        if self.Params["spirou"] in ["max=mean"]:
            print("setting largest to mean of all inputs")
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print(self.Params["spirou"])
            print("AMPA: mean, gmax, imax: ", meangmax, gMax, imaxgmax)
            print("NMDA: mean, gnmax: ", meangnmax, gnMax)

            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            self.set_synapse(
                synapse[imaxgmax], meangmax / nSyn[imaxgmax], meangnmax / nSyn[imaxgmax]
            )
            print("revised values:")
            self.get_synapses(synapse)
            # for i, s in enumerate(synapse):
            #     if i == imaxgmax:
            #         p.gmax = gMax[i]/nSyn[i]  # except the chosen one
            # for p in s.psd.ampa_psd:
            #     gMax[i] = p.gmax
            #     print('revised gmax i : ', i, gMax[i])
        if self.Params["spirou"] in ["all=mean"]:
            print("setting ALL to the mean of all inputs (no variance)")
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print(self.Params["spirou"])
            print("AMPA: mean, gmax, imax: ", meangmax, gMax, imaxgmax)
            print("NMDA: mean, gnmax: ", meangnmax, gnMax)

            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            for i in range(len(synapse)):
                self.set_synapse(synapse[i], meangmax / nSyn[i], meangnmax / nSyn[i])
            print("revised values:")
            self.get_synapses(synapse)
            # for i, s in enumerate(synapse):
            #     if i == imaxgmax:
            #         p.gmax = gMax[i]/nSyn[i]  # except the chosen one
            # for p in s.psd.ampa_psd:
            #     gMax[i] = p.gmax
            #     print('revised gmax i : ', i, gMax[i])

        # run using pyqtgraph's parallel support
        if self.Params["Parallel"]:
            with MP.Parallelize(
                enumerate(TASKS), results=tresults, workers=nWorkers
            ) as tasker:
                for j, x in tasker:
                    tresults = self.single_an_run(
                        post_cell,
                        j,
                        synapseConfig,
                        stimInfo,
                        seeds,
                        preCell,
                        self.an_setup_time,
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
                    t1=stimInfo["run_duration"] * 1000.0,
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
                stim = np.array(tresults[j]["stim"])  # save the stimulus
                result[N] = {
                    "stimInfo": stimInfo,
                    "spikeTimes": spikeTimes,
                    "inputSpikeTimes": inputSpikeTimes,
                    "time": np.array(celltime[0]),
                    "somaVoltage": somaVoltage,
                    "dendriteVoltage": dendriteVoltage,
                    "allDendriteVoltages": allDendriteVoltages,
                    "stimWaveform": stimWaveform,
                    "stimTimebase": stimTimebase,
                }
                if self.Params[
                    "save_all_sections"
                ]:  # just save soma sections        for section in list(g):
                    allDendriteVoltages[N] = tresults[j]["allsecVec"]

        else:
            # Non parallelized version (with --noparallel flag - useful for debugging):
            for j, N in enumerate(range(nReps)):
                print("Repetition %d" % N)

                tresults[j] = self.single_an_run(
                    post_cell,
                    j,
                    synapseConfig,
                    stimInfo,
                    seeds,
                    preCell,
                    self.an_setup_time,
                )

                celltime.append(tresults[j]["time"])  # (self.time)
                spikeTimes = pu.findspikes(
                    tresults[j]["time"],
                    tresults[j]["Vsoma"],
                    threshold,
                    t0=0.0,
                    t1=stimInfo["run_duration"] * 1000.0,
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
                stim = np.array(tresults[j]["stim"])  # save the stimulus
                result[N] = {
                    "stimInfo": stimInfo,
                    "spikeTimes": spikeTimes,
                    "inputSpikeTimes": inputSpikeTimes,
                    "time": np.array(celltime[0]),
                    "somaVoltage": somaVoltage,
                    "dendriteVoltage": dendriteVoltage,
                    "allDendriteVoltages": allDendriteVoltages,
                    "stimWaveform": stimWaveform,
                    "stimTimebase": stimTimebase,
                }
                if self.Params["save_all_sections"]:  # save data for all sections
                    allDendriteVoltages[N] = tresults[j]["allsecVec"]

        total_elapsed_time = time.time() - self.start_time
        #        total_run_time = time.time() - run_time
        print(
            "Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(
                total_elapsed_time / 60.0, total_elapsed_time
            )
        )
        print(
            "Total Setup Time = {:8.2f} min ({:8.0f}s)".format(
                self.setup_time / 60.0, self.setup_time
            )
        )
        print(
            "Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(
                self.an_setup_time / 60.0, self.an_setup_time
            )
        )
        print(
            "Total Neuron Run Time = %{:8.2f} min ({:8.0f}s)".format(
                self.nrn_run_time / 60.0, self.nrn_run_time
            )
        )

        if self.Params["save_all_sections"]:
            for n in range(len(allDendriteVoltages)):
                for s in list(allDendriteVoltages[n].keys()):
                    allDendriteVoltages[n][s] = np.array(allDendriteVoltages[n][s])
        self.analysis_filewriter(self.Params["cell"], result, tag="delays")
        if self.Params["plotFlag"]:
            print("plotting")
            self.plot_an(celltime, result)

    def an_run_singles(self, post_cell: object, exclude=False, verify=False):
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
        post_cell : cnmodel cell object
            Access to neuron and file information

        exclude : boolean (default: False)
            Set false to do one synapse at a time.
            Set true to do all BUT one synapse at a time.

        verify : boolean (default: False)
            Flag to control printing of various intermediate results


        Returns
        -------
            Nothing

        """
        print("\n****singles****\n")
        self.start_time = time.time()
        # cconfig = cell_config.CellConfig()
        # synapseConfig, celltype = cconfig.makeDict(self.cellID)
        if self.Params["inputPattern"] is not None:
            fromdict = self.Params["inputPattern"]
            print("Cell id: %s  using input pattern: %s" % (self.cellID, fromdict))
        else:
            fromdict = self.cellID
        cconfig = cell_config.CellConfig()
        synapseConfig, celltype = cconfig.makeDict(fromdict)
        nReps = self.Params["nReps"]
        threshold = self.Params["threshold"]  # spike threshold, mV

        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(
            post_cell, synapseConfig, celltype, stimInfo
        )

        nSyns = len(preCell)
        seeds = self.compute_seeds(nReps, synapseConfig)

        stimInfo["seeds"] = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = self.Params["Parallel"]
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
            for i, s in enumerate(synapse):
                for p in s.psd.ampa_psd:
                    print(f"Synapse: {i:d}  gmax={p.gmax:.6f}")
            print("tag: %s" % tagname)

            tresults = [None] * nReps

            if parallel and self.Params["nWorkers"] > 1:
                nWorkers = self.Params["nWorkers"]
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with MP.Parallelize(
                    enumerate(TASKS), results=tresults, workers=nWorkers
                ) as tasker:
                    for j, x in tasker:
                        # post_cell, j, synapseConfig, stimInfo, seeds, preCell, an_setup_time
                        tresults = self.single_an_run(
                            post_cell,
                            j,
                            synapseConfig,
                            stimInfo,
                            seeds,
                            preCell,
                            self.an_setup_time,
                        )
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    tresults[j] = self.f(
                        post_cell,
                        j,
                        synapseConfig,
                        stimInfo,
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
                    t1=stimInfo["run_duration"] * 1000,
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
            total_elapsed_time = time.time() - self.start_time
            #        total_run_time = time.time() - run_time
            print(
                "Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(
                    total_elapsed_time / 60.0, total_elapsed_time
                )
            )
            print(
                "Total Setup Time = {:8.2f}min ({:8.0f}s)".format(
                    self.setup_time / 60.0, self.setup_time
                )
            )
            print(
                "Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(
                    self.an_setup_time / 60.0, self.an_setup_time
                )
            )
            print(
                "Total Neuron Run Time = {:8.2f} min ({:8.0f}s)".format(
                    self.nrn_run_time / 60.0, self.nrn_run_time
                )
            )

            result[k] = {
                "stimInfo": stimInfo,
                "spikeTimes": spikeTimes,
                "inputSpikeTimes": inputSpikeTimes,
                "somaVoltage": somaVoltage,
                "dendriteVoltage": dendriteVoltage,
                "time": np.array(tresults[j]["time"]),
                "stimWaveform": stimWaveform,
                "stimTimebase": stimTimebase,
            }

            self.analysis_filewriter(self.Params["cell"], result[k], tag=tagname)
        if self.Params["plotFlag"]:
            self.plot_an(celltime, result)

    def an_run_IO(self, post_cell):
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
        post_cell : cnmodel cell object
            Access to neuron and file information

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
        cconfig = cell_config.CellConfig()
        synapseConfig, celltype = cconfig.makeDict(self.cellID)
        nReps = self.Params["nReps"]
        threshold = self.Params["threshold"]  # spike threshold, mV

        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(
            post_cell, synapseConfig, celltype, stimInfo
        )

        nSyns = len(preCell)
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = self.Params["Parallel"]
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

            if parallel and self.Params["nWorkers"] > 1:
                nWorkers = self.Params["nWorkers"]
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
                            post_cell,
                            j,
                            synapseConfig,
                            stimInfo,
                            seeds,
                            preCell,
                            self.an_setup_time,
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
                        post_cell,
                        j,
                        synapseConfig,
                        stimInfo,
                        seeds,
                        preCell,
                        self.an_setup_time,
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
                    t1=stimInfo["run_duration"] * 1000,
                    dt=1.0,
                    mode="peak",
                )
                spikeTimes[N] = self.clean_spiketimes(spikeTimes[N])
                inputSpikeTimes[N] = tresults[j][
                    "ANSpikeTimes"
                ]  # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]["Vsoma"])
                dendriteVoltage[N] = np.array(tresults[j]["Vdend"])
            stimInfo["gSyns"] = gmaxs
            total_elapsed_time = time.time() - self.start_time
            #        total_run_time = time.time() - run_time
            print(
                "Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(
                    total_elapsed_time / 60.0, total_elapsed_time
                )
            )
            print(
                "Total Setup Time = {:8.2f}min ({:8.0f}s)".format(
                    self.setup_time / 60.0, self.setup_time
                )
            )
            print(
                "Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(
                    self.an_setup_time / 60.0, self.an_setup_time
                )
            )
            print(
                "Total Neuron Run Time = {:8.2f} min ({:8.0f}s)".format(
                    self.nrn_run_time / 60.0, self.nrn_run_time
                )
            )

            result = {
                "stimInfo": stimInfo,
                "spikeTimes": spikeTimes,
                "inputSpikeTimes": inputSpikeTimes,
                "somaVoltage": somaVoltage,
                "dendriteVoltage": dendriteVoltage,
                "time": np.array(tresults[j]["time"]),
            }

            self.analysis_filewriter(self.Params["cell"], result, tag=tagname % k)
        if self.Params["plotFlag"]:
            self.plot_an(celltime, result)

    # def single_an_run_fixed(self, post_cell, j, synapseConfig, stimInfo, preCell, an_setup_time):
    #     """
    #     Perform a single run with all AN input on the target cell turned off except for input j.
    #
    #     Parameters
    #     ----------
    #     hf : hoc_reader object
    #         Access to neuron and file information
    #
    #     j : int
    #         The input that will be active in this run
    #
    #     synapseConfig : dict
    #         A dictionary with information about the synapse configuration to use.
    #
    #     stimInfo : dict
    #         A dictionary whose elements include the stimulus delay
    #
    #     preCell : list
    #         A list of the preCell hoc objects attached to the synapses
    #
    #     post_cell : cells object
    #         The target cell
    #
    #     Returns
    #     -------
    #     anresult : dict
    #         A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'
    #
    #     """
    #     print('\n*** single_an_run_fixed\n')
    #     group = 'soma'
    #     hf = post_cell.hr
    #     try:
    #         cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'])
    #     except:
    #         self.an_run(post_cell, make_an_intial_conditions=True)
    #         print('return from an_run, initial conditions')
    #         try:
    #             cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'])
    #         except:
    #             raise ValueError('Failed initialization for cell: ', self.cellID)
    #
    #     # make independent inputs for each synapse
    #     ANSpikeTimes = []
    #     an0_time = time.time()
    #     nrn_run_time = 0.
    #     stim = {}
    #     #
    #     # Generate stimuli - they are always the same for every synaptic input, so just generate once
    #     #
    #     for i in range(len(preCell)):
    #         preCell[i].set_spiketrain([10., 20., 30., 40., 50.])
    #         ANSpikeTimes.append(preCell[i]._spiketrain)
    #
    #     an_setup_time += (time.time() - an0_time)
    #     nrn_start = time.time()
    #     self.allsecVec = OrderedDict()
    #
    #     if self.Params['save_all_sections']:
    #         for group in post_cell.hr.sec_groups.keys(): # get morphological components
    #             g = post_cell.hr.sec_groups[group]
    #             for section in list(g):
    #                 sec = post_cell.hr.get_section(section)
    #                 self.allsecVec[sec.name()] = post_cell.hr.h.Vector()
    #                 self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)  # recording of voltage all set up here
    #
    #
    #
    #     # if self.Params['save_all_sections']:  # just save soma sections        for section in list(g):
    #    #      print(hf.sec_groups.keys())
    #    #      g = hf.sec_groups[group]
    #    #      sec = hf.get_section(section)
    #    #      self.allsecVec[sec.name()] = hf.h.Vector()
    #    #      self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)  # recording of voltage all set up here
    #
    #     Vsoma = hf.h.Vector()
    #     Vdend = hf.h.Vector()
    #     rtime = hf.h.Vector()
    #     if 'dendrite' in post_cell.all_sections and len(post_cell.all_sections['dendrite']) > 0:
    #         dendsite = post_cell.all_sections['dendrite'][-1]
    #         Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
    #     else:
    #         dendsite = None
    #
    #     Vsoma.record(post_cell.soma(0.5)._ref_v, sec=post_cell.soma)
    #     rtime.record(hf.h._ref_t)
    #     hf.h.finitialize()
    #     hf.h.tstop = 50.
    #     hf.h.t = 0.
    #     hf.h.batch_save() # save nothing
    #     hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
    #     nrn_run_time += (time.time() - nrn_start)
    #     if dendsite == None:
    #         Vdend = np.zeros_like(Vsoma)
    #     anresult = {'Vsoma': np.array(Vsoma), 'Vdend': np.array(Vdend), 'time': np.array(rtime),
    #             'ANSpikeTimes': ANSpikeTimes, 'stim': stim, 'stimTimebase': np.array(rtime)}
    #
    #     # print('all sec vecs: ', self.allsecVec.keys())
    #     # print(' flag for all sections: ', self.Params['save_all_sections'])
    #     if self.Params['save_all_sections']:
    #         anresult['allsecVec'] = self.allsecVec
    #
    #       return anresult

    def configure_cell(self, thisCell, synapseConfig, celltype, stimInfo):
        """
        Configure the cell. This routine builds the cell in Neuron, adds presynaptic inputs
        as described in the synapseConfig, and configures those according to parameters in
        stiminfo.

        Parameters
        ----------
        thisCell : hoc_reader object
            Access to neuron and file information

        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.

        celltype : string
            A string describing the cell type. Determines how the channels are populated
            and how the cell is built, even if it comes from a .hoc structure.

        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)

        Returns
        -------
            preCell : list
                A list of the preCell hoc objects attached to the synapses
            post_cell : cells object
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
        stimInfo["SR"] = stimInfo["SRType"]
        for i, syn in enumerate(synapseConfig):
            if stimInfo["SRType"] is "fromcell":  # use the one in the table
                preCell.append(cells.DummySGC(cf=stimInfo["F0"], sr=syn["SR"]))
                preCell[-1]._synapsetype = self.Params["ANSynapseType"]
                stimInfo["SR"] = self.srname[syn[2]]  # use and report value from table
            else:
                try:
                    srindex = self.srname.index(stimInfo["SRType"])
                    print(
                        "Retrieved index {:d} with SR type {:s}".format(
                            srindex, stimInfo["SRType"]
                        )
                    )
                except:
                    raise ValueError(
                        'SR type "%s" not found in Sr type list' % stimInfo["SRType"]
                    )

                preCell.append(
                    cells.DummySGC(cf=stimInfo["F0"], sr=srindex)
                )  # override
            if self.Params["verbose"]:
                print("SRtype, srindex: ", stimInfo["SRType"], srindex)
            for pl in syn["postlocations"]:
                postsite = syn["postlocations"][pl]
                # note that we split the number of zones between multiple sites
                if self.Params["ANSynapseType"] == "simple":
                    print("*******Synapsetype is simple")
                    synapse.append(
                        preCell[-1].connect(
                            thisCell,
                            type="simple",
                            post_opts={"postsec": pl, "postsite": postsite[0:2]},
                        )
                    )
                else:
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
        if self.Params["ANSynapseType"] == "multisite":
            for i, s in enumerate(synapse):
                s.terminal.relsite.Dep_Flag = self.Params[
                    "ANSynapticDepression"
                ]  # turn on or off depression computation

        electrodeSection = list(thisCell.hr.sec_groups["soma"])[0]
        electrode_site = thisCell.hr.get_section(electrodeSection)
        return (preCell, synapse, electrode_site)

    def set_dbspl(self, signal, dbspl):
        """Scale the level of `signal` to the given dB_SPL."""
        p0 = 20e-6
        rms = np.sqrt(np.sum(signal ** 2) / signal.size)
        scaled = signal * 10 ** (dbspl / 20.0) * p0 / rms
        return scaled

    def single_an_run(
        self, post_cell, j, synapseConfig, stimInfo, seeds, preCell, an_setup_time
    ):
        """
        Perform a single run with all AN input on the target cell turned off except for input j.

        Parameters
        ----------
        hf : hoc_reader object
            Access to neuron and file information

        j : int
            The input that will be active in this run

        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.

        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)

        seeds : 2d numpy array
            An array listing the seeds for starting each auditory nerve input spike train

        preCell : list
            A list of the preCell hoc objects attached to the synapses

        post_cell : cells object
            The target cell

        an_setup_time : time object

        Returns
        -------
        anresult : dict
            A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'

        """
        print("\n*** single_an_run\n")
        # self.an_run(post_cell, make_an_intial_conditions=True)
        # cellInit.get_initial_condition_state(self.post_cell, tdur=self.Params['initialization_time'],
        #    filename=None, electrode_site=None)

        try:
            cellInit.restore_initial_conditions_state(
                post_cell,
                electrode_site=None,
                filename=self.Params["initANStateFile"],
                autoinit=self.Params["auto_initialize"],
            )
            # print('### SUCCESS ###')
        except:
            print(
                "single_an_run: could not restore initial conditions: will try to create again"
            )
            self.an_run(post_cell, make_an_intial_conditions=True)
            print("Return from inital run initial conditions #2")
            try:
                cellInit.restore_initial_conditions_state(
                    post_cell,
                    electrode_site=None,
                    filename=self.Params["initANStateFile"],
                )
            except:
                raise ValueError("Failed initialization for cell: ", self.cellID)

        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.0
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        if isinstance(stimInfo["pip_start"], list):
            pips = stimInfo["pip_start"]
        else:
            pips = [stimInfo["pip_start"]]
        if stimInfo["soundtype"] == "tonepip":
            stim = sound.TonePip(
                rate=stimInfo["Fs"],
                duration=stimInfo["run_duration"],
                f0=stimInfo["F0"],
                dbspl=stimInfo["dB"],
                ramp_duration=stimInfo["RF"],
                pip_duration=stimInfo["pip_duration"],
                pip_start=pips,
            )
        elif stimInfo["soundtype"] == "SAM":
            stim = sound.SAMTone(
                rate=stimInfo["Fs"],
                duration=stimInfo["run_duration"],
                f0=stimInfo["F0"],
                dbspl=stimInfo["dB"],
                ramp_duration=stimInfo["RF"],
                fmod=stimInfo["fmod"],
                dmod=stimInfo["dmod"],
                pip_duration=stimInfo["pip_duration"],
                pip_start=pips,
            )
        else:
            raise ValueError(
                "StimInfo sound type %s not implemented" % stimInfo["soundtype"]
            )
        stimWaveform = stim.generate()
        stimTimebase = stim.time
        for i, syn in enumerate(synapseConfig):
            nseed = seeds[j, i]
            if self.Params["SGCmodelType"] in ["Zilany"]:
                preCell[i].set_sound_stim(
                    stim, seed=nseed, simulator="matlab"
                )  # generate spike train, connect to terminal
            elif self.Params["SGCmodelType"] in ["cochlea"]:
                wf = self.set_dbspl(stim.generate(), stimInfo["dB"])
                stim._sound = wf
                preCell[i].set_sound_stim(
                    stim, seed=nseed, simulator="cochlea"
                )  # generate spike train, connect to terminal
            else:
                raise ValueError(
                    "SGC model type type %s not implemented"
                    % self.Params["SGCmodelType"]
                )
            ANSpikeTimes.append(preCell[i]._spiketrain)

        an_setup_time += time.time() - an0_time
        nrn_start = time.time()
        Vsoma = post_cell.hr.h.Vector()
        Vdend = post_cell.hr.h.Vector()
        rtime = post_cell.hr.h.Vector()
        if (
            "dendrite" in post_cell.all_sections
            and len(post_cell.all_sections["dendrite"]) > 0
        ):
            dendsite = post_cell.all_sections["dendrite"][-1]
            Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        else:
            dendsite = None

        self.allsecVec = OrderedDict()
        if self.Params["save_all_sections"]:
            for group in list(
                post_cell.hr.sec_groups.keys()
            ):  # get morphological components
                g = post_cell.hr.sec_groups[group]
                for section in list(g):
                    sec = post_cell.hr.get_section(section)
                    self.allsecVec[sec.name()] = post_cell.hr.h.Vector()
                    self.allsecVec[sec.name()].record(
                        sec(0.5)._ref_v, sec=sec
                    )  # recording of voltage all set up here
        Vsoma.record(post_cell.soma(0.5)._ref_v, sec=post_cell.soma)
        rtime.record(post_cell.hr.h._ref_t)
        post_cell.hr.h.finitialize()
        post_cell.hr.h.tstop = stimInfo["run_duration"] * 1000.0
        post_cell.hr.h.t = 0.0
        post_cell.hr.h.batch_save()  # save nothing
        post_cell.hr.h.dt = self.Params['dt']
        post_cell.hr.h.batch_run(post_cell.hr.h.tstop, post_cell.hr.h.dt, "an.dat")
        nrn_run_time += time.time() - nrn_start
        if dendsite == None:
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
        if self.Params["save_all_sections"]:
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
        somaVoltage : array (no default)
            Voltage recorded at the soma for this cell
        stimInfo : dict (no default)
            stimulus parameters (require 'nReps' and 'threshold' for spikes)
        dendVoltage : array (default: None)
            Voltage recorded at some point in the dendrite for this cell

        Returns
        -------
        Nothing

        """
        if not self.Params["plotFlag"]:
            return
        # if self.Params['Parallel']:
        #     return  # no plots when doing parallel runs...
        nReps = self.Params["nReps"]
        threshold = self.Params["threshold"]
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
            "stimInfo",
            "spikeTimes",
            "inputSpikeTimes",
            "somaVoltage",
            "time",
            "stimWaveform",
            "stimTimebase",
        ]

        try:
            for rk in requiredKeys:
                assert rk in k  # make sure keys are there
            res_mode = "syn"
        except:
            res_mode = "reps"

        results = {}
        # results with be a dict with params, stiminfo, and trials as keys
        print("\n*** analysis_filewriter\n")

        results['basename'] = self.Params["simulationFilename"]
        results["Params"] = self.Params  # include all the parameters of the run too
        # clean up Parame to remove PosixPath from filenames
        for v in [
            "initIVStateFile",
            "initANStateFile",
            "simulationFilename",
            "hocfile",
        ]:
            results["Params"][v] = str(results["Params"][v])
        results['runInfo'] = None  # lazy, all info is in "params" 
        results['modelPars'] = copy.deepcopy(self.post_cell.status)
        results["Results"] = result
        results["mode"] = res_mode
        fout = self.Params[
            "simulationFilename"
        ]  # base name created by make-filename - do not overwrite
        if len(tag) > 0:
            fout = Path(str(fout).replace("_all_", "_" + tag + "_"))
        if self.Params["tagstring"] is not None:
            fnb = fout.stem
            fp = fout.parent
            fnb = str(fnb) + "_" + self.Params["tagstring"] + ".p"
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
            print(len(self.hg.colorMap))
            print(len(list(hf.sec_groups.keys())))
            for i, s in enumerate(hf.sec_groups.keys()):
                print(s)
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
    curdir = Path.cwd()
    model = ModelRun()  # create instance of the model

    parser = argparse.ArgumentParser(
        description="Simulate activity in a reconstructed model cell",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )
    parser.add_argument(
        dest="cell", action="store", default=None, help="Select the cell (no default)"
    )
    parser.add_argument(
        "--type",
        "-T",
        dest="cellType",
        action="store",
        default="Bushy",
        choices=model.cellChoices,
        help="Define the cell type (default: Bushy)",
    )
    parser.add_argument(
        "--model",
        "-M",
        dest="modelName",
        action="store",
        default="XM13",
        choices=model.modelNameChoices,
        help="Define the model type (default: XM13)",
    )
    parser.add_argument(
        "--modeltype",
        dest="modelType",
        action="store",
        default="II",
        choices=model.modelTypeChoices,
        help="Define the model type (default: XM13)",
    )
    parser.add_argument(
        "--sgcmodel",
        type=str,
        dest="SGCmodelType",
        action="store",
        default="cochlea",
        choices=model.SGCmodelChoices,
        help="Define the SGC model type (default: Zilany)",
    )
    parser.add_argument(
        "--protocol",
        "-P",
        dest="runProtocol",
        action="store",
        default="runIV",
        choices=model.protocolChoices,
        help="Protocol to use for simulation (default: IV)",
    )
    parser.add_argument(
        "-H",
        "--defaulthoc",
        action="store_true",
        dest="usedefaulthoc",
        default=True,
        help="Use default hoc file for this cell",
    )
    parser.add_argument(
        "--hocfile",
        dest="hocfile",
        action="store",
        default=None,
        help='hoc file to use for simulation (default is the selected "cell".hoc)',
    )
    parser.add_argument(
        "-F",
        "--full",
        dest="fullhocfile",
        action="store_true",
        default=False,
        help='Use "full" hoc file as in "VCN_c02_Full.hoc instead of VCN_c02.hoc")',
    )

    parser.add_argument(
        "--inputpattern",
        "-i",
        type=str,
        dest="inputPattern",
        action="store",
        default=None,
        help="cell input pattern to use (substitute) from cell_config.py",
    )
    parser.add_argument(
        "--stimulus",
        "-s",
        type=str,
        dest="soundtype",
        action="store",
        default="tonepip",
        choices=model.soundChoices,
        help="Define the stimulus type (default: tonepip)",
    )
    parser.add_argument(
        "--check",
        "-/",
        action="store_true",
        default=False,
        dest="checkcommand",
        help="Only check command line for valid input; do not run model",
    )
    parser.add_argument(
        "-C",
        "--configfile",
        type=str,
        default=None,
        dest="configfile",
        help="Read a formatted configuration file (JSON, TOML) for commands",
    )

    # lowercase options are generally parameter settings:
    parser.add_argument(
        "-d",
        "--dB",
        type=float,
        default=30.0,
        dest="dB",
        help="Set sound intensity dB SPL (default 30)",
    )
    parser.add_argument(
        "-f",
        "--frequency",
        type=float,
        default=4000.0,
        dest="F0",
        help="Set tone frequency, Hz (default 4000)",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.1,
        dest="pip_duration",
        help="Set sound stimulus duration (sec; default 0.1)",
    )
    parser.add_argument(
        "-r", "--reps", type=int, default=1, dest="nReps", help="# repetitions"
    )
    parser.add_argument(
        "--seed", type=int, default=1, dest="seed", help="AN starting seed"
    )
    parser.add_argument(
        "-S",
        "--SRType",
        type=str,
        default="HS",
        dest="SRType",
        choices=model.SRChoices,
        help=("Specify SR type (from: %s)" % model.SRChoices),
    )
    parser.add_argument(
        "--synapsetype",
        type=str,
        default="multisite",
        dest="ANSynapseType",
        choices=model.ANSynapseChoices,
        help=("Specify AN synapse type (from: %s)" % model.ANSynapseChoices),
    )
    parser.add_argument(
        "--depression",
        type=int,
        default=0,
        dest="ANSynapticDepression",
        choices=[0, 1],
        help=(
            "Specify AN depression flag for multisite synapses (from: %s)" % str([0, 1])
        ),
    )

    parser.add_argument(
        "--fmod",
        type=float,
        default=20.0,
        dest="fmod",
        help="Set SAM modulation frequency",
    )
    parser.add_argument(
        "--dmod",
        type=float,
        default=100.0,
        dest="dmod",
        help="Set SAM modulation depth (in percent)",
    )
    parser.add_argument(
        "--S2M",
        type=float,
        default=0,
        dest="signalToMasker",
        help="Signal to Masker ratio (dB)",
    )
    parser.add_argument(
        "--cmmrmode",
        type=str,
        default="CMR",
        dest="CMMRmode",
        choices=model.cmmrModeChoices,
        help=("Specify mode (from: %s)" % model.cmmrModeChoices),
    )

    parser.add_argument(
        "--spirou",
        type=str,
        dest="spirou",
        action="store",
        default="all",
        choices=model.spirouChoices,
        help="Specify spirou experiment type.... ",
    )
    # Morphology
    parser.add_argument(
        "--soma_inflate",
        type=float,
        dest="soma_inflation",
        action="store",
        default=1.0,
        help="Specify factor by which to inflate soma AREA",
    )
    parser.add_argument(
        "--soma_autoinflate",
        action="store_true",
        dest="soma_autoinflate",
        default=False,
        help="Automatically inflate soma based on table",
    )
    parser.add_argument(
        "--dendrite_inflate",
        type=float,
        dest="dendrite_inflation",
        action="store",
        default=1.0,
        help="Specify factor by which to inflate total dendritic AREA",
    )
    parser.add_argument(
        "--dendrite_autoinflate",
        action="store_true",
        dest="dendrite_autoinflate",
        default=False,
        help="Automatically inflate dendrite area based on table",
    )
    parser.add_argument(
        "--dendrite_from_soma",
        action="store_true",
        dest="dendrite_fromsoma",
        default=False,
        help="Automatically inflate dendrite area based on soma inflation",
    )
    parser.add_argument(
        "--ASA_from_soma",
        action="store_true",
        dest="ASA_fromsoma",
        default=False,
        help="Automatically inflate dendrite area based on soma inflation",
    )

    parser.add_argument(
        "--tagstring",
        type=str,
        default=None,
        dest="tagstring",
        help="Add a tag string to the output filename to distinguish it",
    )
    parser.add_argument(
        "-a",
        "--AMPAScale",
        type=float,
        default=1.0,
        dest="AMPAScale",
        help="Set AMPAR conductance scale factor (default 1.0)",
    )
    parser.add_argument(
        "--allmodes",
        action="store_true",
        default=False,
        dest="all_modes",
        help=("Force run of all modes (CMR, CMD, REF) for stimulus configuration."),
    )
    parser.add_argument(
        "--sequence",
        type=str,
        default="",
        dest="sequence",
        help=("Specify a sequence for the primary run parameters"),
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        default=False,
        dest="plotFlag",
        help="Plot results as they are generated - requires user intervention... ",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        dest="nWorkers",
        help='Number of "workers" for parallel processing (default: 4)',
    )
    parser.add_argument(
        "--noparallel",
        action="store_false",
        default=True,
        dest="Parallel",
        help="Use parallel or not (default: True)",
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        default=False,
        dest="auto_initialize",
        help="Force auto initialization if reading the state fails in initialization",
    )
    parser.add_argument(
        "--saveall",
        action="store_true",
        default=False,
        dest="save_all_sections",
        help="Save data from all sections in model",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="Print out extra stuff for debugging",
    )

    # Parser arguments for gif noise generator:
    # parameters for generating a noise signal to generating GIF model of cell
    parser.add_argument(
        "--gifi",
        type=float,
        default=0.0,
        dest="gif_i0",
        help="Set Noise for GIF current level (default 0 nA)",
    )
    parser.add_argument(
        "--gifsigma",
        type=float,
        default=0.2,
        dest="gif_sigma",
        help="Set Noise for GIF variance (default 0.2 nA)",
    )
    parser.add_argument(
        "--giffmod",
        type=float,
        default=0.2,
        dest="gif_fmod",
        help="Set Noise for GIF fmod (default 0.2 Hz)",
    )
    parser.add_argument(
        "--giftau",
        type=float,
        default=3.0,
        dest="gif_tau",
        help="Set Noise for GIF tau (default 0.3 ms)",
    )
    parser.add_argument(
        "--gifdur",
        type=float,
        default=10.0,
        dest="gif_dur",
        help="Set Noise for GIF duration (default 10 s)",
    )
    parser.add_argument(
        "--gifskew",
        type=float,
        default=0.0,
        dest="gif_skew",
        help="Set Noise for GIF to have skewed distribution (0 = normal)",
    )

    # parser.add_argument('-p', '--print', action="store_true", default=False, dest = 'print_info',
    #     help='Print extra information during analysis')
    #  parser.add_argument('-l', '--list', action="store_true", default=False, dest = 'list_results',
    #     help='List results to screen')

    args = parser.parse_args()
    if args.hocfile is not None:
        args.usedefaulthoc = False

    if args.configfile is not None:
        config = None
        if args.configfile is not None:
            if ".json" in args.configfile:
                # The escaping of "\t" in the config file is necesarry as
                # otherwise Python will try to treat is as the string escape
                # sequence for ASCII Horizontal Tab when it encounters it
                # during json.load
                config = json.load(open(args.configfile))
            elif ".toml" in args.configfile:
                config = toml.load(open(args.configfile))

        vargs = vars(args)  # reach into the dict to change values in namespace
        for c in config:
            if c in args:
                print("c: ", c)
                vargs[c] = config[c]

    # now copy into the model.Params structure
    for key, value in list(vars(args).items()):
        if key in list(model.Params.keys()):
            model.Params[key] = value
        else:
            raise ValueError(f"Parameter {k:s} is not in Params key list")

    # now add the configs if there are any to add.
    # if config is not None:
    #     for k in config.keys():
    #         if k in list(model.Params.keys()):
    #             model.Params[k] = config[k]
    #         else:
    #             raise ValueError(f'Parameter {k:s} is not in Params key list')
    print(
        json.dumps(model.Params, indent=4)
    )  # pprint doesn't work well with ordered dicts

    # print(model.Params['cell'])
    # if model.Params['hocfile'] == None: # just use the matching hoc file
    #     model.Params['hocfile'] = model.Params['cell'] + '.hoc'
    # allargs = ''
    # parsedargs, unparsed = parser.parse_known_args()
    # for parg in vars(parsedargs):
    #     try:
    #         allargs = allargs + ('{0:s}: {1:s}\n'.format(parg, getattr(parsedargs, parg)))
    #     except:
    #         pass
    if model.Params["checkcommand"]:
        model.Params["commandline"] = " ".join(sys.argv)
        print(
            json.dumps(model.Params, indent=4)
        )  # pprint doesn't work well with ordered dicts
        print(model.Params["commandline"])

    if not model.Params["checkcommand"]:
        model.Params["Icmds"] = "[-2.0, 2.0, 0.5]"
        model.run_model()  # then run the model

    # if sys.flags.interactive == 0:
    #     pg.QtGui.QApplication.exec_()


if __name__ == "__main__":
    main()
