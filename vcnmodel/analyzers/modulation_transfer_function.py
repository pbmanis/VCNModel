""" modulation_transfer_function.py - calculate rate modulation transfer function

Calculates rate modulation transfer function from a spike train, for the specified frequency

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import numpy as np
import vcnmodel.analyzers.flatten_spike_array as VAFlatten

def modulation_transfer_function(spikes, freq:float=None, time_window:tuple=(0,1.0), nreps: int = 0):
    """
    Calculate for the specified frequency

    Parameters
    ----------
    spikes : Spike train, in sec.
        If the data comes from repeated trials, the spike train needs to be flattened into a 1d array before
        calling.
    freq : Stimulus frequency in Hz
    time_window: window of time in seconds for the measurement window (used to compute rate)
    nreps: number of repetitions in this data set
    extras: bool (default True):
        whether to include rMTF and entrainment in the calculation.
    Returns
    -------
        float: MTF value for this spike train
    """
    assert freq is not None
    assert nreps > 0
    assert spikes is not None
    
    spikes = VAFlatten.flatten_spike_array(spikes, time_window = time_window, isi_flag=False)
    period = 1.0 / freq
    n_spikes = len(spikes)

    # find the number of full cycles that fit in the window duration,
    # recalculate the window and get the spikes from that

    twopi_per = 2.0 * np.pi / period
    phasev = twopi_per * np.fmod(
        spikes, period
    )  # convert to time within a cycle period in radians
    rMTF = 0.0

    if n_spikes > 0:
        earliest_spike = np.min(spikes)
        zsp = spikes - earliest_spike
        n_periods = int(np.floor((time_window[1]-time_window[0]) / period))
        max_time = n_periods * period
        mtf_spikes = zsp[(zsp >= 0.0) & (zsp < max_time)]
        # print(f"^^^^^^^ period: {period:.5f} window: {window_duration:.2f}  max_time: {max_time:.2f}  nspikes: {VSR.n_spikes:d} mtf_spikes: {mtf_spikes.shape[0]:d}")
        rMTF = mtf_spikes.shape[0] / max_time / nreps

    return rMTF