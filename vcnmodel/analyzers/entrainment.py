""" entrainment.py - calculate spike entrainment for a periodic stimulus

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Distributed under MIT/X11 license. See license.txt for more infomation. 

Entrainment:

From Rudnick and Hemmert, 2017:

"The entrainment index p is calculated using the definition given by Joris et al. (1994).
First an inter-spike interval histogram is constructed. Then the number of intervals k 
belonging to the first maximum (between 0.5 and 1.5 of the stimulus period) 
is divided by the total number of intervals n as given by:
p = k/n  ( Equation 6 )
The values of EI are between 0 and 1. 
A value of 1 indicates that the neuron fired a minimum of one spike during each 
stimulus period. Entrainment falls to 0 when the neuron can no longer fire action
potentials for adjacent stimulus periods."

The above is the same as the definition given in Ashida et al. PLoS 2019.
The original definition was from Joris et al. 1994; in that paper "1/CF" was used
as the interval and window. Entrainment was also described in Rothman and Young, 1996.
Note in both Joris et al., 1994 and Rothman and Young, 1995, the window used
to measure the spikes at the "fundamental" (not "subharmonic") intervals is not
explicitely or unambiguously stated. 

In this implementation, we exclude any spikes with an interval less
than half of the modulation or frequency period. The means that cells that fire
more than once per cycle with a short ISI are not penalized (only the first spike
of the interval counts). Another way of stating this is that 
"""

import numpy as np
import vcnmodel.analyzers.flatten_spike_array as VAFlatten

def entrainment(spikes, freq:float=None, time_window:tuple=(0., 10), minisi: float = 0.):
    """Compute spike entrainment (see notes above)

    Args:
        spikes (_type_): a 1d- array or list of spike times (if list, each array is a repetition)
        freq (float, optional): _description_. Defaults to None.
        window_duration (float, optional): _description_. Defaults to 1.0.
        nreps (int, optional): _description_. Defaults to 0.
        minisi (_type_, optional): _description_. Defaults to 0..

    Returns:
        _type_: _description_
    """
    assert spikes is not None
    assert freq is not None
    assert len(time_window) == 2
    
    n_spikes = len(spikes)
    if n_spikes == 0:
        return 0.0
    
    # find the number of full cycles that fit in the window duration,
    # recalculate the window and get the spikes from that
    period = 1.0 / freq
    twopi_per = 2.0 * np.pi / period

    entrainment = 0.0

    spikes = VAFlatten.flatten_spike_array(mtf_spikes, time_window = time_window,
                    isi_flag=True)
    earliest_spike = np.min(spikes)
    zsp = spikes - earliest_spike
    n_periods = int(np.floor(np.diff(time_window) / period))
    max_time = n_periods * period
    mtf_spikes = zsp[(zsp >= 0.0) & (zsp < max_time)]
    # print(f"^^^^^^^ period: {period:.5f} window: {window_duration:.2f}  max_time: {max_time:.2f}  nspikes: {VSR.n_spikes:d} mtf_spikes: {mtf_spikes.shape[0]:d}")

    # now calculate entrainment
    spike_isis = np.diff(mtf_spikes)  # ISI from spikes in full cycles
    spike_isis = spike_isis[spike_isis > minisi/freq]  # only include intervals within each trial (across trials will be negative), and remove closely aligned spikes
    fhalf = 0.5/freq
    fonehalf = 1.5/freq
    # print(fhalf, fonehalf, freq)
    n_entrain_total = len(spike_isis)
    entrained = spike_isis[(spike_isis >= fhalf) & (spike_isis < fonehalf)]
    if n_entrain_total == 0:
        entrainment = 0.0
    else:
        entrainment = len(entrained)/float(n_entrain_total)
    
    return entrainment
