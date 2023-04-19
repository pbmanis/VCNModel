"""
Define some data structures that are used in multiple analysis locations
"""
from dataclasses import dataclass, field
import numpy as np

def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


def def_empty_dict():
    return {}

@dataclass()
class PData:
    """
    data class for some parameters that control what we read
    """
    gradeA: list = field(default_factory=def_empty_list)
    default_modelName: str = "XM13_nacncoop"
    soma_inflate: bool = True
    dend_inflate: bool = True
    basepath: str = ""  # config["baseDataDirectory"]
    renderpath: str = ""  # " str(Path(self.config["codeDirectory"], "Renderings"))
    revcorrpath: str = ""
    thiscell: str = ""

@dataclass()
class RevCorrPars:
    ntrials: int = 1  # number of trials
    ninputs: int = 1  # number of inputs for the cell under study
    algorithm: str = "RevcorrSPKS"  # save the algorithm name
    # clip trace to avoid end effects
    min_time: float = 10.0  # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time: float = 250.0  # this window needs to be at least as long as maxwin
    binw: float = 0.1  # binwidth, ms
    minwin: float = (
        -5.0
    )  # start of revcorr display/calculation, ms relative to reference event
    maxwin: float = 2.5  # end of display/calculation relative to reference event
    amax: float = 0.0  # peak amplitude
    si: dict = field(default_factory=def_empty_dict)  # stimulus parameters
    ri: dict = field(default_factory=def_empty_dict)  # runInfon parameters


@dataclass()
class RevCorrData:
    C: list = field(default_factory=def_empty_list)  # using Brian 1.4 correlation
    CB: list = field(
        default_factory=def_empty_list
    )  # using elephant/ binned correlation
    CBT: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    TC: float = 0.0  # list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)  # spike times
    tx: np.array = field(default_factory=def_empty_np)  #
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(
        default_factory=def_empty_np
    )  # timebase for average revcorr
    sv_all: list = field(
        default_factory=def_empty_list
    )  # np.array = field(default_factory=def_empty_np)  #
    sv_avg: np.array = field(default_factory=def_empty_np)
    sv_trials: list = field(default_factory=def_empty_list)
    sites: np.array = field(default_factory=def_empty_np)
    nsp_avg: int = 0
    npost_spikes: int = 0
    npre_spikes: int = 0
    mean_pre_intervals: np.array = field(default_factory=def_empty_np)
    mean_post_intervals: float = 0.0
    max_coin_rate: float = 0
    participation: np.array = field(default_factory=def_empty_np)
    s_pair: float = 0.0
    ynspike: np.array = field(default_factory=def_empty_np)
    pre_w: list = field(default_factory=def_empty_list)
    pre_st: list = field(default_factory=def_empty_list)


@dataclass()
class SpikeData:
    """
    Data class to hold information about each spike
    """

    trial: int = -1  # the trial the spike came from
    time_index: int = 0  # index into the time array for this spike
    dt: float = 0.025  # sample interval, msec
    start_win: float = -5.0
    end_win: float = 5.0
    waveform: np.array = None  # the waveform of this postspike, clipped
    prespikes: np.array = None  # time indices to pre spikes


@dataclass()
class Patterns:
    """
    Hold results for different input spike patterns relative to the output
    spikes. This structure is used by spike_pattern_analysis
    """
    selected_inputs: field(default_factory=def_empty_list)
    event_spike_pattern: field(default_factory=def_empty_list)
    event_pair_wise: field(default_factory=def_empty_list)
    filter_table: field(default_factory=def_empty_list)
    filter_table2: field(default_factory=def_empty_list)
    all_spikes: field(default_factory=def_empty_list)
    n_post: int = 0
    n_pre: int = 0
    name: str = ""
