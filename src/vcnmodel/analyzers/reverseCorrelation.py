import typing
from typing import Union, List, Tuple
import numpy as np

"""
Comopute the reverse correlation between two spike trains
For all of the spikes in st1, compute the time difference
with each spike in st2, find where the difference is within a
time window we are interested in, then get and return the indices
for all of those spikes.
"""


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