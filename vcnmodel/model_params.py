"""
Define data structures used for:
    
    *Command line parsing
    *Specifying model parameters (general)
    *Specifying runInfo (run instance) parameters

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

"""

import time
from pathlib import Path
import argparse
from dataclasses import dataclass, field
from typing import Union
import numpy as np
import json
import toml



display_orient_cells = {
    "VCN_c02": [140.0, 0.0, -144.0],
    "VCN_c06": [140.0, -59.0, -12.0],
    "VCN_c05": [140.0, -46.0, 121.0],
    "VCN_c09": [140.0, -74.0, 18.0],
    "VCN_c11": [140.0, -2.0, -181.0],
    "VCN_c10": [140.0, 5.0, -35.0],
    "VCN_c13": [140.0, -22.0, 344.0],
    "VCN_c17": [140.0, -158.0, 39.0],
    "VCN_c30": [140.0, -134.0, -181.0],
}


@dataclass(frozen=True)
class CmdChoices:
    parser = None
    cellChoices = ["Bushy", "TStellate", "DStellate"]
    modelNameChoices = [
        "XM13",
        "XM13_nacncoop",
        "XM13A_nacncoop",
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
    ]
    # cochlea is python model of Zilany data, no matlab, JIT computation;
    # Zilany model creates matlab instance for every run.
    cmmrModeChoices = ["CM", "CD", "REF"]  # comodulated, codeviant, reference
    SRChoices = [
        "LS",
        "MS",
        "HS",
        "fromcell",
        "mixed1",
    ]  # AN SR groups (assigned across all inputs)
    dendriteChoices = [
        "normal",
        "passive",
        "active",
        "allpassive"
    ]
    dendriteExptChoices = [
        "default",
        "Full",
        "AxonOnly",
        "NoDend",
        "NoDistal",
        "NoUninnervated",
    ]
    axonExptChoices = [
        "default",
        "standardized"
    ]
    protocolChoices = [
        "initIV",
        "Zin",
        "runIV",
        "initandrunIV",
        "initVC",
        "runIVSpikeThreshold",
        "runVC",
        "initAN",
        "runANPSTH",
        "runANIO",
        "runANSingles",
        "runANThreshold",
        "gifnoise",
    ]
    soundChoices = ["tonepip", "noise", "stationaryNoise", "SAM", "CMMR"]
    speciesChoices = ["mouse", "guineapig"]
    SpirouChoices = [
        "all",
        "max=mean",
        "all=mean",
        "removelargest",
        "removetwolargest",
        "largestonly",
        "twolargest",
        "threelargest",
        "fourlargest",
    ]
    ANSynapseChoices = ["simple", "multisite"]

    srname = ["LS", "MS", "HS"]  # runs 0-2, not starting at 0

    displayModeChoices = ["None", "vm", "sec-type", "mechanism"]
    displayStyleChoices = ["cylinders", "graph", "volume", "surface"]
    channelChoices = ["klt", "kht", "ihvcn", "nacncoop", "nacn", "najsr"]

# set up the Params data structure. This should hold everything that might be modified
# or need to be known when running the CmdChoices.
#


@dataclass
class Params:

    setup: bool = False  # true once we have setup the cell and filenames
    cellID: Union[str, None] = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
    cell: object = None  # model instance (neuron/hoc)
    AMPAScale: float = 1.0  # Use the default scale for AMPAR conductances
    ANSynapseType: str = "simple"  # or multisite
    ANSynapticDepression: int = 0  # depression calculation is off by default
    SynapseConfiguration: Union[dict, None] = None
    synapseConfig: Union[dict, None] = None
    initStateFile: Union[str, Path, None] = None  # 'IVneuronState_%s.dat'
    simulationFilename: Union[str, Path, None] = None
    shortSimulationFilename: Union[str, Path, None] = None
    simPath: Union[str, Path, None] = None
    hocfile: Union[str, Path, None] = None
    meshInflate: bool = True  # use the hoc file that has been inflated to the mesh area
    usedefaulthoc: bool = False
    cellType: str = CmdChoices.cellChoices[0]
    modelName: str = CmdChoices.modelNameChoices[0]
    dataTable: str = ""
    dendriteMode: str = CmdChoices.dendriteChoices[0]
    dendriteExpt: str = CmdChoices.dendriteExptChoices[0]
    axonExpt: str = CmdChoices.axonExptChoices[0]
    modelType: str = CmdChoices.modelTypeChoices[0]
    SGCmodelType: str = CmdChoices.SGCmodelChoices[0]
    species: str = CmdChoices.speciesChoices[0]
    displayStyle: str = CmdChoices.displayStyleChoices[0]
    displayMechanism: str = CmdChoices.channelChoices[0]
    displayMode: str = CmdChoices.displayModeChoices[0]

    # cell specific parameters related to geometry
    fullhocfile: bool = False  # use the "full" hoc file (cellname_Full.hoc) (obselete)
    hocstring: str=""  # a string containing the hoc file that was read and used at the time of the simulation
    dtIC: float = 0.025  # ok.
    dtVC: float = 0.005  # voltage clamp; need shorter steop size for transient measure
    celsius: float = 37  # set the temperature.
    Ra: float = 150.0  # ohm.cm
    soma_inflation: float = 1.0  # factor to multiply soma section areas by
    soma_autoinflate: bool = False
    dendrite_inflation: float = 1.0
    dendrite_autoinflate: bool = False
    dendrite_fromsoma: bool = False
    ASA_inflation: float = 1.0
    ASA_fromsoma: bool = False
    lambdaFreq: float = 2000.0  # Hz for segment number
    area_adjustment_method: str = "pt3d"
    # spontaneous rate (group, in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
    srnames = ["LS", "MS", "HS", "mixed1"]  # runs 0-2, not starting at 0    # same as CmcChoices
    SRType: str = CmdChoices.SRChoices[2]
    inputPattern: Union[
        str, None,
    ] = None  # ID of cellinput pattern (same as cellID): for substitute input patterns.
    synno: Union[int, None] = None  # for selection of single synaptic inputs

    lastfile: Union[Path, str, None] = None

    seed: int = 100  # always the same start - avoids lots of recomputation
    # in production, vary this or set to None for random values
    all_modes: bool = False

    # general control parameters
    plotFlag: bool = False
    auto_initialize: bool = False
    nWorkers: int = 8
    Parallel: bool = True
    verbose: bool = False
    save_all_sections: bool = False
    commandline: str = ""  # store command line on run, all parser args
    commands: str = ""  # store command line on run, the actual given commands
    checkcommand: bool = False
    testsetup: bool = False
    configfile: Union[str, None] = None
    tagstring: Union[str, None] = None
    initialization_time: float = 200.0


# runinfo parameters are filled by generate_run at the initialization of a single trace run
def definj():
    return {"pulse": np.linspace(-1.0, 2.00, 16, endpoint=True)}


def defVCsteps():
    """
    Steps are relative to HOLDING, which is typicall -60 mV
    """
    start = -20.
    stop = 120.
    step = 10.
    npts = int((stop-start)/step)+1
    return {"pulse": np.linspace(start, stop, npts, endpoint=True)}


def defstarts():
    return [0.1]


def defemptydict():
    return {}


def defemptylist():
    return []

def defdatadict():
    return {'ChannelCompartments': "", 'ChannelData': ""}

@dataclass
class RunInfo:
    # general and filenames
    folder: str = "Simulations"
    fileName: str = "Normal"
    runProtocol: str = CmdChoices.protocolChoices[
        2
    ]  # tZin is default because it is fast and should be run often
    runName: str = "Run"
    manipulation: str = "Canonical"
    preMode: str = "CC"
    postMode: str = "CC"
    TargetCellType: str = ""  # celltype, # valid are "Bushy", "Stellate", "MNTB"
    electrodeSection: Union[object, str, None] = None  # electrodeSection
    electrodeSectionName: Union[str, None] = None
    dendriticElectrodeSection: Union[
        object, str, None
    ] = None  # dendriticElectrodeSection,
    dendriticSectionDistance: float = 100.0  # microns.
    tableData: dict = field(default_factory=defdatadict)  # the table data itself, which might change between runs...

    nReps: int = 1
    seeds: list = field(default_factory=defemptylist)
    # current injection parameters
    sequence: str = ""  # sequence for run - may be string [start, stop, step]
    nStim: int = 1
    stimFreq: float = 200.0  # hz
    stimInj: dict = field(default_factory=definj)
    stimVC: dict = field(default_factory=defVCsteps)
    # iRange,  # nA, a list of test levels for current clamp
    stimDur: float = 100.0  # msec
    stimDelay: float = 5.0  # msec
    stimPost: float = 3.0  # msec
    vnStim: float = 1
    vstimFreq: float = 200.0  # hz
    vstimInj: float = 50  # mV amplitude of step (from holding)
    vstimDur: float = 100.0  # msec
    vstimDelay: float = 5.0  # msec
    vstimPost: float = 25.0  # msec
    vstimHolding: float = -65.0  # holding, mV

    # sound parameters
    initialization_time: float = 50.0  # nominal time to let system settle, in msec
    run_duration: float = 0.35  # in sec
    soundtype: str = "SAM"  # or 'tonepip'
    noise_seed: int = 9  # noise generator seed value
    pip_duration: float = 0.1  # duration in seconds for tone pip
    pip_start: list = field(default_factory=defstarts)  # start (delay to start of pip)
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
    spikeTimeList: dict = field(
        default_factory=defemptydict
    )  # Dictionary of spike times
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
        "--dendritemode",
        dest="dendriteMode",
        default="normal",
        choices=CmdChoices.dendriteChoices,
        help="Choose dendrite table (normal, active, passive)",
    )

    parser.add_argument(
        "--hocfile",
        dest="hocfile",
        action="store",
        default=None,
        help='hoc file to use for simulation (default is the selected "cell".hoc)',
    )

    parser.add_argument(
        "--nomeshInflate",
        dest="meshInflate",
        action="store_false",
        default=False,
        help="use uninflated hoc file, not mesh-inflated file (default: False, uses mesh inflated file)",
    )

    parser.add_argument(
        "-D",
        "--dendriteexpt",
        dest="dendriteExpt",
        default='default',
        choices=CmdChoices.dendriteExptChoices,
        help="Choose dendrite experiment (default, Full, NoDend, NoDistal, NoUninnervated)",
    )
    parser.add_argument(
        "-A",
        "--axonexpt",
        dest="axonExpt",
        default='default',
        # choices=CmdChoices.axonExptChoices,
        help="Choose dendrite/axon experiment (default, standardized)",
    )
    parser.add_argument(
        "--datatable",
        type=str,
        dest="dataTable",
        action="store",
        default="",
        help="Specify the data table for this run",
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
        "--style",
        dest="displayStyle",
        default=CmdChoices.displayStyleChoices[0],
        choices=CmdChoices.displayStyleChoices,
        help="Render cell with neuronvis with style: {str(CmdChoices.displayStyleChoices):s}",
    )
    parser.add_argument(
        "--mechanism",
        dest="displayMechanism",
        default=CmdChoices.channelChoices[0],
        help=f"Render channel mechanisms: {str(CmdChoices.channelChoices):s}",
    )
    parser.add_argument(
        "--displaymode",
        dest="displayMode",
        default=CmdChoices.displayModeChoices[0],
        choices=CmdChoices.displayModeChoices,
        help=f"Render cell with neuronvis : set mode to one of {str(CmdChoices.displayModeChoices):s}",
    )

    parser.add_argument(
        "--displayscale",
        dest="displayscale",
        action="store_true",
        default=False,
        help="use display scale and orientation from table for generating renderings",
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
        "--testsetup",
        action="store_true",
        default=False,
        help="Test all setup, but do not run simulations"
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
    # parser.add_argument(
    #     "--duration",
    #     type=float,
    #     default=0.1,
    #     dest="pip_duration",
    #     help="Set sound stimulus duration (sec; default 0.1)",
    # )
    parser.add_argument(
        "-r", "--reps", type=int, default=1, dest="nReps", help="# repetitions"
    )
    parser.add_argument(
        "--seed", type=int, default=1, dest="seed", help="AN starting seed"
    )
    parser.add_argument(
        "--noise_seed", type=int, default=1, dest="noise_seed", help="Noise generator starting seed"
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
        "--pip_start",
        type=float,
        default=0.1,
        dest="pip_start",
        help="Set delay to onset of acoustic stimulus",
    )
    parser.add_argument(
        "--pip_duration",
        type=float,
        default=0.1,
        dest="pip_duration",
        help="Set duration of acoustic stimulus",
    )
    parser.add_argument(
        "--pip_offduration",
        type=float,
        default=0.05,
        dest="pip_offduration",
        help="Time to continue simulation AFTER sound ends",
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
        "--Spirou",
        type=str,
        dest="Spirou",
        action="store",
        default="all",
        choices=CmdChoices.SpirouChoices,
        help="Specify Spirou experiment type.... ",
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
        "--hold",
        type=float,
        action="store",
        dest="vstimHolding",
        default=-80.0,
        help="Holding voltage in VClamp (mV) (default: -80 mV)",
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


def getCommands(toml_dir='.'):
    parser = build_parser()
    args = parser.parse_args()

    if args.configfile is not None:
        config = None
        if args.configfile is not None:
            cfile = Path(toml_dir, args.configfile)
            if ".json" in str(args.configfile):
                # The escaping of "\t" in the config file is necesarry as
                # otherwise Python will try to treat is as the string escape
                # sequence for ASCII Horizontal Tab when it encounters it
                # during json.load
                config = json.load(open(cfile))
                print(f"Reading JSON configuration file: {str(cfile):s}")
            elif ".toml" in str(args.configfile):
                print(f"Reading TOML configuration file: {str(cfile):s}")
                config = toml.load(open(cfile))

        vargs = vars(args)  # reach into the dict to change values in namespace

        for c in config:
            if c in args:
                # print("Getting parser variable: ", c)
                vargs[c] = config[c]
            else:
                raise ValueError(
                    f"config variable {c:s} does not match with comand parser variables"
                )

        print("   ... All configuration file variables read OK")
    else:
        vargs = vars(args)
    # now copy into the Param dataclass
    params = Params()
    runinfo = RunInfo()
    parnames = dir(params)
    runnames = dir(runinfo)
    for key, value in vargs.items():
        if key in parnames:
            # print('key: ', key)
            # print(".  ", str(value))
            exec(f"params.{key:s} = {value!r}")
        elif key in runnames:
            exec(f"runinfo.{key:s} = {value!r}")

    return (args, params, runinfo)


if __name__ == "__main__":
    a, p, r = getCommands()
    print(p)
