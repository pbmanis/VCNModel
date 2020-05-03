import sys
import time
from pathlib import Path
import argparse
from dataclasses import dataclass, field
from typing import Union, Dict, List
import dataclasses
from collections import OrderedDict
import numpy as np
import json
import toml
from pprint import PrettyPrinter

"""
Define data structures used for:
    Command line parsing
    Specifying model parameters (general)
    Specifying runInfo (run instance) parameters

"""


@dataclass(frozen=True)
class CmdChoices:
    parser = None
    cellChoices = ["Bushy", "TStellate", "DStellate"]
    modelNameChoices = [
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
    modelTypeChoices = ["II", "II-I", "I-II", "I-c", "I-t", "II-o"]
    SGCmodelChoices = [
        "Zilany",
        "cochlea",
    ]  # cochlea is python model of Zilany data, no matlab, JIT computation; Zilany model creates matlab instance for every run.
    cmmrModeChoices = ["CM", "CD", "REF"]  # comodulated, codeviant, reference
    SRChoices = [
        "LS",
        "MS",
        "HS",
        "fromcell",
    ]  # AN SR groups (assigned across all inputs)
    protocolChoices = [
        "initIV",
        "testIV",
        "runIV",
        "initandrunIV",
        "initAN",
        "runANPSTH",
        "runANIO",
        "runANSingles",
        "gifnoise",
    ]
    soundChoices = ["tonepip", "noise", "stationaryNoise", "SAM", "CMMR"]
    speciesChoices = ["mouse", "guineapig"]
    spirouChoices = ["all", "max=mean", "all=mean", "removelargest", "largestonly"]
    ANSynapseChoices = ["simple", "multisite"]

    srname = ["LS", "MS", "HS"]  # runs 0-2, not starting at 0


# set up the Params data structure. This should hold everything that might be modified
# or need to be known when running the CmdChoices.
#


@dataclass
class Params:

    setup: bool = False  # true once we have setup the cell and filenames
    cellID: str = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
    cell: object = None  # model instance (neuron/hoc)
    AMPAScale: float = 1.0  # Use the default scale for AMPAR conductances
    ANSynapseType: str = "simple"  # or multisite
    ANSynapticDepression: int = 0  # depression calculation is off by default
    initIVStateFile: Union[str, Path, None] = None  # 'IVneuronState_%s.dat'
    initANStateFile: Union[str, Path, None] = None  # 'ANneuronState_%s.dat'
    simulationFilename: Union[str, Path, None] = None
    shortSimulationFilename: Union[str, Path, None] = None
    simPath: Union[str, Path, None] = None
    hocfile: Union[str, Path, None] = None
    usedefaulthoc: bool = False
    cellType: str = CmdChoices.cellChoices[0]
    modelName: str = CmdChoices.modelNameChoices[0]
    modelType: str = CmdChoices.modelTypeChoices[0]
    SGCmodelType: str = CmdChoices.SGCmodelChoices[0]
    species: str = CmdChoices.speciesChoices[0]

    # cell specific parameters related to geometry
    fullhocfile: bool = False  # use the "full" hoc file (cellname_Full.hoc)
    dt: float = 0.025
    celsius: float = 37  # set the temperature.
    Ra: float = 150.0  # ohm.cm
    soma_inflation: float = 1.0  # factor to multiply soma section areas by
    soma_autoinflate: bool = False  #
    dendrite_inflation: float = 1.0
    dendrite_autoinflate: bool = False
    dendrite_fromsoma: bool = False
    ASA_inflation: float = 1.0
    ASA_fromsoma: bool = False
    lambdaFreq: float = 2000.0  # Hz for segment number
    # spontaneous rate (group, in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
    srnames = ["LS", "MS", "HS"]  # runs 0-2, not starting at 0    # same as CmcChoices
    SRType: str = CmdChoices.SRChoices[2]
    inputPattern:  Union[
        str,
        None,
    ] = None  # ID of cellinput pattern (same as cellID): for substitute input patterns.
    synno:  Union[int, None] = None  # for selection of single synaptic inputs

    lastfile: Union[Path, str, None] = None

    seed: int = 100  # always the same start - avoids lots of recomutation
    # in production, vary this or set to None for random values
    all_modes: bool = False

    # general control parameters
    plotFlag: bool = False
    auto_initialize: bool = False
    nWorkers: int = 8
    Parallel: bool = True
    verbose: bool = False
    save_all_sections: bool = False
    commandline: str = ""  # store command line on run
    checkcommand: bool = False
    configfile: Union[str, None] = None
    tagstring:  Union[str, None] = None
    initialization_time: float = 50.


# claas ModelPar:
#
#     axon:bool = False
#     cellClass:[str, None]: None
#     dendrites:bool= False,
#     hillock:bool= False,
#     initialsegment:bool= False,
#     modelName:[str, None]= None,
#     modelType:str= 'II',
#     morphology': None,
#     myelinatedaxon': False,
#     na': None,
#     name': None,
#     pumps': False,
#     soma': True,
#     species': 'mouse',
#     temperature': 34.0,
#     ttx': False,
#     unmyelinatedaxon': False,

# runinfo parameters are filled by generate_run at the initialization of a single trace run
def definj():
    return{
        "pulse": np.linspace(-1.0, 2.00, 16, endpoint=True)
    }

def defstarts():
    return [0.1]
def defemptydict():
    return {}
def defemptylist():
    return []

@dataclass
class RunInfo:
    # general and filenames
    folder: str = "Simulations"
    fileName: str = "Normal"
    runProtocol: str = CmdChoices.protocolChoices[
        2
    ]  # testIV is default because it is fast and should be run often
    runName: str = "Run"
    manipulation: str = "Canonical"
    preMode: str = "cc"
    postMode: str = "cc"
    TargetCellType: str = ""  # celltype, # valid are "Bushy", "Stellate", "MNTB"
    electrodeSection:Union[object, str, None] = None  # electrodeSection
    electrodeSectionName:str = None
    dendriticElectrodeSection:Union[object, str, None] = None  # dendriticElectrodeSection,
    dendriticSectionDistance:float = 100.0  # microns.

    nReps: int = 1
    seeds: list = field(default_factory = defemptylist)
    # current injection parameters
    sequence: str = ""  # sequence for run - may be string [start, stop, step]
    nStim: int = 1
    stimFreq: float = 200.0  # hz
    stimInj: dict = field(default_factory = definj)
      # iRange,  # nA, a list of test levels for current clamp
    stimDur: float = 100.0  # msec
    stimDelay: float = 5.0  # msec
    stimPost: float = 3.0  # msec
    vnStim: float = 1
    vstimFreq: float = 200.0  # hz
    vstimInj: float = 50  # mV amplitude of step (from holding)
    vstimDur: float = 50.0  # msec
    vstimDelay: float = 2.0  # msec
    vstimPost: float = 3.0  # msec
    vstimHolding: float = -60  # holding, mV

    # sound parameters
    initialization_time: float = 50.0  # nominal time to let system settle, in msec
    run_duration: float = 0.35  # in sec
    soundtype: str = "SAM"  # or 'tonepip'
    pip_duration: float = 0.1  # duration in seconds for tone pip
    pip_start: list = field(default_factory = defstarts)  # start (delay to start of pip)
    pip_offduration: float = 0.05  # time after pip ends to keep running
    
    Fs: float = 100e3  # cochlea/zilany model rate
    F0: float = 16000.0  # stimulus frequency
    dB: float = 30.0  # in SPL
    RF: float = 2.5e-3  # rise-fall time
    fmod: float = 20  # hz, modulation if SAM
    dmod: float = 0.0  # percent if SAM
    threshold: float = -35.0  # spike threshold, mV
    signalToMasker: float = 0.0
    CMMRmode: str = "CM"
    Spirou: str = "all"

    # gif model parameters
    # parameters for generating a noise signal to generating GIF model of cell
    gif_i0: float = 0.0  # base current level
    gif_sigma: float = 0.5  # std of noise
    gif_fmod: float = 0.2  # mod in Hz
    gif_tau: float = 3.0  # tau, msec
    gif_dur: float = 10.0  # seconds
    gif_skew: float = 0.0  # a in scipy.skew
    runTime: object = time.asctime()  # store date and time of run
    inFile: Union[str, Path, None] = None
    # 'ANFiles/AN10000Hz.txt', # if this is not None, then we will use these spike times...
    inFileRep: int = 1  # which rep to use (or array of reps)
    spikeTimeList: dict = field(default_factory = defemptydict)  # Dictionary of spike times
    v_init: float = -61.0  # from Rothman type II model - not appropriate in all cases
    useSaveState: bool = True  # useSavedState,  # use the saved state.


def build_parser():
    """
    Set up the command line parser for running the model
    This include both verbs for running different types of protocols, and
    parameters that define the cell scaffold (structure, channels) that is being run
    The parameters in this file should also appear in the Params data structure (dict) for the class
    """

    parser = argparse.ArgumentParser(
        description="Simulate activity in a reconstructed model cell",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )
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
        choices=CmdChoices.cellChoices,
        help="Define the cell type (default: Bushy)",
    )
    parser.add_argument(
        "--model",
        "-M",
        dest="modelName",
        action="store",
        default="XM13",
        choices=CmdChoices.modelNameChoices,
        help="Define the model type (default: XM13)",
    )
    parser.add_argument(
        "--modeltype",
        dest="modelType",
        action="store",
        default="II",
        choices=CmdChoices.modelTypeChoices,
        help="Define the model type (default: XM13)",
    )
    parser.add_argument(
        "--sgcmodel",
        type=str,
        dest="SGCmodelType",
        action="store",
        default="cochlea",
        choices=CmdChoices.SGCmodelChoices,
        help="Define the SGC model type (default: Zilany)",
    )
    parser.add_argument(
        "--protocol",
        "-P",
        dest="runProtocol",
        action="store",
        default="runIV",
        choices=CmdChoices.protocolChoices,
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
        choices=CmdChoices.soundChoices,
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
        choices=CmdChoices.SRChoices,
        help=("Specify SR type (from: %s)" % CmdChoices.SRChoices),
    )
    parser.add_argument(
        "--synapsetype",
        type=str,
        default="multisite",
        dest="ANSynapseType",
        choices=CmdChoices.ANSynapseChoices,
        help=("Specify AN synapse type (from: %s)" % CmdChoices.ANSynapseChoices),
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
        choices=CmdChoices.cmmrModeChoices,
        help=("Specify mode (from: %s)" % CmdChoices.cmmrModeChoices),
    )

    parser.add_argument(
        "--spirou",
        type=str,
        dest="spirou",
        action="store",
        default="all",
        choices=CmdChoices.spirouChoices,
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
    return parser


def getCommands():
    parser = build_parser()
    args = parser.parse_args()

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
                print('getting toml')
                config = toml.load(open(args.configfile))

        vargs = vars(args)  # reach into the dict to change values in namespace

        for c in config:
            if c in args:
                # print("Getting parser variable: ", c)
                vargs[c] = config[c]
            else:
                raise ValueError(f"config variable {c:s} does not match with comand parser variables")

        print("All configuration file variables read OK")
    # now copy into the Param dataclass
    params = Params()
    runinfo = RunInfo()
    parnames = dir(params)
    runnames = dir(runinfo)
    for key, value in vargs.items():
        if key in parnames:
            # print('key: ', key)
            # print(str(value))
            exec(f"params.{key:s} = {value!r}")
        elif key in runnames:
            exec(f"runinfo.{key:s} = {value!r}")
    
    return(args, params, runinfo)


if __name__ == "__main__":
    getCommands()
