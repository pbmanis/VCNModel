""" Compute the cv and regularity according to Young et al., J. Neurophys,
    60: 1, 1988.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2020 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import warnings  # suppress nan of empty array on 2d-np arrays.
from dataclasses import dataclass

import numpy as np


def isi_cv(
    splist: list,
    binwidth: float = 0.001,
    reftime: float = 0.0,
    t0: float = 0.0,
    t1: float = 0.30,
    tgrace: float = 0.0,
    ARP: float = 0.0007,
):
    """Compute the cv and regularity according to Young et al., J. Neurophys,
        60: 1, 1988.
        Analysis is limited to isi's starting at or after t0 but before t1, and
        ending completely before t1 + tgrace(to avoid end effects). t1 should
        correspond to the the end of the stimulus Version using a list of numpy
        arrays for cvisi The ARP is the "Absolute Refractory Period" (see
        Rothman et al. 1993; Rothman and Manis, 2003c) The default value used is
        typically 0.7 msec. this is used as a correction factor for comparision
        with those papers and others that use this value.

    Parameters
    ----------
    splist : list
        Spike times, as a list on a per-trial basis. Units: seconds_
    binwidth : float, optional
        bin width for computing CV, by default 0.001
    reftime : float, optional
        reference time for stimulus onset, by default 0.0
    t0 : float, optional
        start of analysis window for CV, by default 0.0
    t1 : float, optional
        end of analysis window for CV, by default 0.30
    tgrace : float, optional
        grace period for CV, by default 0.0
    ARP : float, optional
        Absolute Refractory period, by default 0.0007

    Returns
    -------
    tuple
        cvisi time bins, the full 2d array,
        cvt, mean and standard deviation.
    """
    cvisit = np.arange(0, t1 + tgrace, binwidth)  # build time bins
    bincount = np.zeros_like(cvisit).astype(int)
    cvisi = np.nan * np.zeros((len(splist), len(cvisit)))  # and matching CV array
    for trial, spkt in enumerate(splist):  # for all the traces
        spkt = np.array(spkt) - reftime
        in_win = np.where((spkt >= t0) & (spkt < t1))[0]
        spkt = spkt[in_win]
        if len(in_win) < 2:
            continue
        isib = np.floor(np.array(spkt) / binwidth).astype(
            int
        )  # indices for spike times
        isii = np.diff(spkt)  # associated intervals
        for j, spike in enumerate(spkt):  # loop over spikes
            if j >= len(spkt) - 1:  # second interval would be out of bounds
                continue
            if (
                (spike < t0) or (spike > t1 + tgrace) or (spkt[j + 1] > t1 + tgrace)
            ):  # first spike of the pair is outside the window or second spike is outside
                continue
            spkindex = int(np.floor(spike / binwidth))
            cvisi[trial, spkindex] = isii[j]
            bincount[spkindex] += 1
    # Suppress some warnings that are can occur here because sometimes the array is empty...
    # but that is ok
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="Mean of empty slice")
        warnings.filterwarnings(
            action="ignore", message="Degrees of freedom <= 0 for slice"
        )
        cvmn = np.nanmean(cvisi, axis=0) - ARP
        cvsd = np.nanstd(cvisi, axis=0)
    cvt = cvisit

    return cvisit, cvisi, cvt, cvmn, cvsd


def firing_rate(spikes):
    """
    Rate of the spike train, single trial
    """
    if len(spikes) < 2:
        return np.nan
    return (len(spikes) - 1) / (spikes[-1] - spikes[0])

def mean_firing_rate(spikes):
    if len(spikes) < 2:
        return np.nan
    ISI = np.diff(spikes)  # interspike intervals
    return 1.0 / np.mean(ISI)

def CV(spikes):
    """
    Coefficient of variation.
    """
    if spikes == []:
        return np.nan
    ISI = np.diff(spikes)  # interspike intervals
    return np.std(ISI) / np.mean(ISI)


if __name__ == "__main__":

    """
    Test the calculation
    """
    import matplotlib.pyplot as mpl

    N = 100
    reps = 50
    rng = np.random.default_rng(1)
    spikes = rng.exponential(scale=5.0, size=(reps, N))
    T1 = [np.cumsum(s) for s in spikes]
    print(np.mean(T1[0]))
    print(f"Firing rate (first trial): {firing_rate(T1[0]):.3f}")
    for i in range(len(T1)):
        print(
            f"max time for trial: {i:d} = {np.max(T1[i]):.2f} ms   nspikes = {len(T1[i]):d}"
        )
    cvisit, cvisi, cvt, cvm, cvs = isi_cv(T1, binwidth=2, t0=10, t1=210, tgrace=25)
    print("len of isis: ", len(cvisi), len(cvisit))
    print("mean: ", np.mean(cvs / cvm))
    mpl.plot(cvt, cvs, "bo-")
    mpl.plot(cvt, cvm, "rs-")
    mpl.show()
