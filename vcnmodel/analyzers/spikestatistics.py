"""
Spike statistics
----------------
In all functions below, spikes is a sorted list of spike times

These functions are from Brian 1.4, and are not included in Brian 2.
From 2017.
See the public archive on github at:
    https://github.com/brian-team/brian
        at brian/brian/tools/statistics.py (renamed here to be more 
        descriptive of the functions)

"""
from __future__ import print_function
import numpy as np
import quantities as pq

from operator import itemgetter

__all__ = [
    "firing_rate",
    "CV",
    "correlogram",
    "autocorrelogram",
    "CCF",
    "ACF",
    "CCVF",
    "ACVF",
    "group_correlations",
    "sort_spikes",
    "total_correlation",
    "vector_strength",
    "gamma_factor",
    "get_gamma_factor_matrix",
    "get_gamma_factor",
    "spike_triggered_average",
]

# First-order statistics
def firing_rate(spikes):
    """
    Rate of the spike train.
    """
    if len(spikes) < 2:
        return np.nan
    return (len(spikes) - 1) / (spikes[-1] - spikes[0])


def CV(spikes):
    """
    Coefficient of variation.
    """
    if spikes == []:
        return np.nan
    ISI = diff(spikes)  # interspike intervals
    return std(ISI) / mean(ISI)


# Second-order statistics
def correlogram(T1, T2, width=20.0*pq.ms, bin_width=1.0*pq.ms, T=None):
    """
    Returns a cross-correlogram with lag in [-width,width] and given bin_width size.
    T is the total duration (optional) and should be greater than or
    equal ti the duration of T1 and T2.
    The result is in Hz (rate of coincidences in each bin).

    N.B.: units are discarded, but all must match in msec
    TODO: optimise?
    """
    if (T1 == []) or (T2 == []):  # empty spike train
        return np.nan
    # Remove units
    width = float(width)
    T1 = np.array(T1)
    T2 = np.array(T2)
    i = 0
    j = 0
    n = int(np.ceil(width / bin_width)) # Histogram length
    l = []
    for t in T1:
        while i < len(T2) and T2[i] < t - width:  # other possibility use searchsorted
            i += 1
        while j < len(T2) and T2[j] < t + width:
            j += 1
        l.extend(T2[i:j] - t)
    H, _ = np.histogram(l, bins=np.arange(2 * n + 1) * bin_width - n * bin_width) #, new = True)

    # Divide by time to get rate
    if T is None and len(T1) > 0:
        T = np.max((T1[-1], T2[-1]))*pq.ms - np.min((T1[0], T2[0]))*pq.ms
    if len(T1) == 0:  # handle empty spike train
        T = 0*pq.ms
    # Windowing function (triangle)
    W = np.zeros(2 * n)
    W[:n] = T - bin_width * np.arange(n - 1, -1, -1)
    W[n:] = T - bin_width * np.arange(n)

    return H / W


def autocorrelogram(T0, width=20.0, bin_width=1.0, T=None):
    """
    Returns an autocorrelogram with lag in [-width,width] and given bin_width size.
    T is the total duration (optional) and should be greater than the duration of T1 and T2.
    The result is in Hz (rate of coincidences in each bin).

    N.B.: units are discarded.
    """
    return correlogram(T0, T0, width, bin_width, T)


def CCF(T1, T2, width=20.0, bin_width=1.0, T=None):
    """
    Returns the cross-correlation function with lag in [-width,width] and given bin_width size.
    T is the total duration (optional).
    The result is in Hz**2:
    CCF(T1,T2)=<T1(t)T2(t+s)>

    N.B.: units are discarded.
    """
    return correlogram(T1, T2, width, bin_width, T) / bin_width


def ACF(T0, width=20, bin_width=1, T=None):
    """
    Returns the autocorrelation function with lag in [-width,width] and given bin_width size.
    T is the total duration (optional).
    The result is in Hz**2:
    ACF(T0)=<T0(t)T0(t+s)>

    N.B.: units are discarded.
    """
    return CCF(T0, T0, width, bin_width, T)


def CCVF(T1, T2, width=20, bin_width=1, T=None):
    """
    Returns the cross-covariance function with lag in [-width,width] and given bin_width size.
    T is the total duration (optional).
    The result is in Hz**2:
    CCVF(T1,T2)=<T1(t)T2(t+s)>-<T1><T2>

    N.B.: units are discarded.
    """
    return CCF(T1, T2, width, bin_width, T) - firing_rate(T1) * firing_rate(T2)


def ACVF(T0, width=20, bin_width=1, T=None):
    """
    Returns the autocovariance function with lag in [-width,width] and given bin_width size.
    T is the total duration (optional).
    The result is in Hz**2:
    ACVF(T0)=<T0(t)T0(t+s)>-<T0>**2

    N.B.: units are discarded.
    """
    return CCVF(T0, T0, width, bin_width, T)


def spike_triggered_average(
    spikes, stimulus, max_interval, dt, onset=None, display=False
):
    """
    Spike triggered average reverse correlation. 
    spikes is an array containing spike times
    stimulus is an array containing the stimulus
    max_interval (second) is the duration of the averaging window
    dt (second) is the sampling period
    onset (second) before which the spikes are discarded. Note: it will be at least as long as max_interval
    display (default=False) display the number of spikes processed out of the total number
    output the spike triggered average and the corresponding time axis
    """
    stimulus = stimulus.flatten()
    if onset < max_interval or onset == None:
        onset = max_interval
    nspikes = len(spikes)
    sta_length = int(max_interval/dt)
    spike_triggered_ensemble=np.zeros((nspikes,sta_length))
    time_axis = linspace(0, max_interval, sta_length)
    onset = float(onset)
    for ispike, spike in enumerate(spikes):
        if display == True:
            print("sta: spike #", ispike, " out of :", nspikes)
        if spike > onset:
            spike = int(spike / dt)
            # print stimulus[spike-sta_length:spike].shape
            spike_triggered_ensemble[ispike, :] = stimulus[spike - sta_length : spike]
            ispike += 1

    return sum(spike_triggered_ensemble, axis=0)[::-1] / (nspikes - 1), time_axis


def total_correlation(T1, T2, width=20, T=None):
    """
    Returns the total correlation coefficient with lag in [-width,width].
    T is the total duration (optional).
    The result is a real (typically in [0,1]):
    total_correlation(T1,T2)=int(CCVF(T1,T2))/rate(T1)
    
    Modified: width has neg and pos parts. 
    """
    if (T1 == []) or (T2 == []):  # empty spike train
        return np.nan
    # Remove units
    if not np.isscalar(width):
        twidth = sum(width)
        nwidth = width[0]
        pwidth = width[1]
    else:
        nwidth = pwidth = twidth = width
    T1 = np.array(T1)
    T2 = np.array(T2)
    # Divide by time to get rate
    if T is None and len(T1) > 0:
        T = max(T1[-1], T2[-1]) - min(T1[0], T2[0])
    i = 0
    j = 0
    x = 0
    for t in T1:
        while i < len(T2) and T2[i] < t - nwidth:  # other possibility use searchsorted
            i += 1
        while j < len(T2) and T2[j] < t + pwidth:
            j += 1
        x += sum(
            1.0 / (T - abs(T2[i:j] - t))
        )  # counts coincidences with windowing (probabilities)
    #    return float(x / firing_rate(T1)) - float(firing_rate(T2) *2  * twidth)
    return float(x / firing_rate(T1)) - float(firing_rate(T2) * twidth)


def sort_spikes(spikes):
    """
    Sorts spikes stored in a (i,t) list by time.
    """
    spikes = sorted(spikes, key=itemgetter(1))
    return spikes


def group_correlations(spikes, delta=None):
    """
    Computes the pairwise correlation strength and timescale of the given pool of spike trains.
    spikes is a (i,t) list and must be sorted.
    delta is the length of the time window, 10*pq.ms by default.
    """
    aspikes = np.array(spikes)
    N = aspikes[:, 0].max() + 1  # neuron count
    T = aspikes[:, 1].max()  # total duration
    spikecount = np.zeros(N)
    tauc = np.zeros((N, N))
    S = np.zeros((N, N))
    if delta is None:
        delta = 10  # size of the window
    windows = (
        -2 * delta * ones(N)
    )  # windows[i] is the end of the window for neuron i = lastspike[i}+delta
    for i, t in spikes:
        sources = t <= windows  # neurons such that (i,t) is a target spike for them
        if sum(sources) > 0:
            indices = nonzero(sources)[0]
            S[indices, i] += 1
            delays = t - windows[indices] + delta
            #            print i, t, indices, delays
            tauc[indices, i] += delays
        spikecount[i] += 1
        windows[i] = t + delta  # window update

    tauc /= S

    S = S / tile(spikecount.reshape((-1, 1)), (1, N))  # normalize S
    rates = spikecount / T
    S = S - tile(rates.reshape((1, -1)), (N, 1)) * delta

    S[isnan(S)] = 0.0
    tauc[isnan(tauc)] = 0.0

    return S, tauc


# Phase-locking properties
def vector_strength(spikes, period):
    """
    Returns the vector strength of the given train
    """
    return abs(np.mean(np.exp(np.array(spikes) * 1j * 2 * np.pi / period)))


# Normalize the coincidence count of two spike trains (return the gamma factor)
def get_gamma_factor(
    coincidence_count, model_length, target_length, target_rates, delta
):
    NCoincAvg = 2 * delta * target_length * target_rates
    norm = 0.5 * (1 - 2 * delta * target_rates)
    gamma = (coincidence_count - NCoincAvg) / (norm * (target_length + model_length))
    return gamma


# Normalize the coincidence matrix between a set of  trains (return the gamma factor matrix)
def get_gamma_factor_matrix(
    coincidence_matrix, model_length, target_length, target_rates, delta
):

    target_lengthMAT = tile(target_length, (len(model_length), 1))
    target_rateMAT = tile(target_rates, (len(model_length), 1))
    model_lengthMAT = tile(model_length.reshape(-1, 1), (1, len(target_length)))
    NCoincAvg = 2 * delta * target_lengthMAT * target_rateMAT
    norm = 0.5 * (1 - 2 * delta * target_rateMAT)

    # print  target_rateMAT
    print(coincidence_matrix)
    # print NCoincAvg
    # print (norm * (target_lengthMAT + model_lengthMAT))
    gamma = (coincidence_matrix - NCoincAvg) / (
        norm * (target_lengthMAT + model_lengthMAT)
    )
    gamma = triu(gamma, 0) + triu(gamma, 1).T
    return gamma


# Gamma factor


def gamma_factor(source, target, delta, normalize=True, dt=None):
    """
    Returns the gamma precision factor between source and target trains,
    with precision delta.
    source and target are lists of spike times.
    If normalize is True, the function returns the normalized gamma factor 
    (less than 1.0), otherwise it returns the number of coincidences.
    dt is the precision of the trains, by default it is defaultclock.dt
    
    Reference
    * R. Jolivet et al., 'A benchmark test for a quantitative assessment of simple neuron models',
    Journal of Neuroscience Methods 169, no. 2 (2008): 417-424.
    """

    source = np.array(source)
    target = np.array(target)
    target_rate = firing_rate(target) * Hz

    if dt is None:
        delta_diff = delta
    else:
        source = np.array(rint(source / dt), dtype=int)
        target = np.array(rint(target / dt), dtype=int)
        delta_diff = int(rint(delta / dt))

    source_length = len(source)
    target_length = len(target)

    if target_length == 0 or source_length == 0:
        return 0

    if source_length > 1:
        bins = 0.5 * (source[1:] + source[:-1])
        indices = digitize(target, bins)
        diff = abs(target - source[indices])
        matched_spikes = diff <= delta_diff
        coincidences = sum(matched_spikes)
    else:
        indices = [
            amin(abs(source - target[i])) <= delta_diff for i in xrange(target_length)
        ]
        coincidences = sum(indices)

    # Normalization of the coincidences count
    #    NCoincAvg = 2 * delta * target_length * target_rate
    #    norm = .5*(1 - 2 * target_rate * delta)
    #    gamma = (coincidences - NCoincAvg)/(norm*(source_length + target_length))

    # TODO: test this
    gamma = get_gamma_factor(
        coincidences, source_length, target_length, target_rate, delta
    )

    if normalize:
        return gamma
    else:
        return coincidences


if __name__ == "__main__":

    import matplotlib.pyplot as mpl

    print(vector_strength([1.1, 1, 0.9], 2))

    N = 100000
    T1 = np.cumsum(np.random.random(N) * 10)
    T2 = np.cumsum(np.random.random(N) * 10)
    duration = T1[int(N / 2)]  # Cut so that both spike trains have the same duration
    T1 = T1[T1 < duration]
    T2 = T2[T2 < duration]
    print(firing_rate(T1))
    C = CCVF(T1, T2, bin_width=1)
    print(total_correlation(T1, T2))
    mpl.plot(C)
    mpl.show()
