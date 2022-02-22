import dataclasses
import typing
from dataclasses import dataclass, field
from typing import List, Tuple, Union

import numpy as np

from vcnmodel.util import trace_calls as TRC

"""
Compute the reverse correlation between two spike trains
For all of the spikes in st1, compute the time difference
with each spike in st2, find where the difference is within a
time window we are interested in, then get and return the indices
for all of those spikes.
"""

def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


def def_empty_dict():
    return {}


@dataclass
class RevCorrPars:
    ntrials: int = 1  # number of trials
    ninputs: int = 1  # number of inputs for the cell under study
    algorithm: str = "RevcorrSPKS"  # save the algorithm name
    # clip trace to avoid end effects
    min_time: float = 10.0  # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time: float = 250.0  # this window needs to be at least as long as maxwin
    binw: float = 0.1  # binwidth, ms
    minwin: float = -5.0  # start of revcorr display/calculation, ms relative to reference event
    maxwin: float = 2.5  # end of display/calculation relative to reference event
    amax: float = 0.0  # peak amplitude
    si: dict = field(default_factory=def_empty_dict)  # stimulus parameters
    ri: dict = field(default_factory=def_empty_dict)  # runInfon parameters


@dataclass
class RevCorrData:
    C: list = field(default_factory=def_empty_list)  # using Brian 1.4 correlation
    CB: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    CBT: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    TC: float = 0.0  # list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)  # spike times
    tx: np.array = field(default_factory=def_empty_np)  # 
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(default_factory=def_empty_np)  # timebase for average revcorr
    sv_all: np.array = field(default_factory=def_empty_np)  # 
    sv_avg: np.array = field(default_factory=def_empty_np)
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

from numba import jit


@jit(
    parallel=True, cache=True,
)
def nb_revcorr(st1, st2, binwidth, corrwindow):
    xds = np.zeros(int((corrwindow[1] - corrwindow[0]) / binwidth))

    # [None]*len(st1)
    for i, sp in enumerate(st1):
        diff = st2 - sp
        v = np.where((corrwindow[0] <= diff) & (diff <= corrwindow[1]))
        # print('v: ', v)
    #  if len(v) > 0:
    #      iv = [int(vx/binwidth) for vx in v]
    #      xds[iv] = 1
    #      print('xds: ', xds)# for d in diff:
    #     if corrwindow[0] <= d <= corrwindow[1]:
    #         xds[i] = d
    # print(np.array(xds))
    # print(np.array(xds).ravel().shape)
    return xds, len(st1)  # return the n postsynaptic spikes


def reverse_correlation(
    st1: Union[np.ndarray, List] = None,
    st2: Union[np.ndarray, List] = None,
    binwidth: float = 0.1,
    corrwindow: Union[List, Tuple] = [
        -5.0,
        1.0,
    ],  # time window to examine correlation relative to st1
) -> (np.ndarray, int):

    
    if st1 is None or st2 is None:
        raise ValueError(
            "coincident_spikes_correlation: reference and comparator must be defined"
        )

    xds = np.zeros(int((corrwindow[1] - corrwindow[0]) / binwidth))
    for i, sp in enumerate(st1):
        diff = st2 - sp
        v = diff[np.where((corrwindow[0] < diff) & (diff < corrwindow[1]))]
        # print('v: ', v)
        if len(v) > 0:
            indxs = [int(vx / binwidth) for vx in v]
            xds[indxs] = xds[indxs] + 1

    return xds, len(st1)  # return the n postsynaptic spikes