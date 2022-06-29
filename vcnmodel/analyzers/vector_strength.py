""" vector_strength.py - calculate standard vector strength

Calculates vector strength and related parameters from a spike train, for the specified frequency

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import dataclasses
from dataclasses import dataclass, field

import matplotlib.pyplot as mpl
import numpy as np
import scipy.stats


@dataclass
class VSResult:
    vs: float = np.nan  # vector strength
    n_spikes: int = 0  # number of spikes
    Rayleigh: float = np.nan  # coefficient
    pRayleigh: float = 1.0  # p value for NOT flat
    circ_phase: float = np.nan  # mean phase (circularized phase. radians)
    circ_phaseMean: float = np.nan
    circ_phaseSD: float = np.nan  # SD in phase
    circ_timeSD: float = np.nan  # SD in time units
    frequency: float = np.nan  # modulation or VS measure frequency

    # carrier_frequency: float = np.nan  # carrier frequency
    # dmod: float = 0.  # modulation depth
    # dB : float = 0. # stimulus intensity
    # n_inputs : int = 0  # number of inputs to cell
    # an_vs : float = np.nan  # vs for an input (placeholder)
    # an_circ_phaseMean: float=np.nan
    # an_circ_phaseSD: float=np.nan
    # an_circ_timeSD: float = np.nan


class VectorStrength:
    def __init__(self):
        pass

    def vector_strength(self, spikes, freq=None):
        """
        Calculate vector strength and related parameters from a spike train, for the specified frequency

        Parameters
        ----------
        spikes : Spike train, in sec.
            If the data comes from repeated trials, the spike train needs to be flattened into a 1d array before
            calling.
        freq : Stimulus frequency in Hz

        Returns
        -------
            A dataclass of type VSResult
        """
        assert freq is not None
        VSR = VSResult()
        amax = 0.0
        if spikes is None:
            print("NO spikes")

        VSR.frequency = freq
        period = 1.0 / freq
        VSR.n_spikes = len(spikes)
        twopi_per = 2.0 * np.pi / period
        print(np)
        print(np.fmod)
        phasev = twopi_per * np.fmod(
            spikes, period
        )  # convert to time within a cycle period in radians
        VSR.circ_phase = phasev
        sumcos = np.sum(np.cos(phasev))
        sumsin = np.sum(np.sin(phasev))
        mean_phase = np.arctan2(sumsin, sumcos)

        # print('circ phase, mean: ', VSR.circ_phase, mean_phase)
        sumc2 = sumcos * sumcos
        sums2 = sumsin * sumsin

        if VSR.n_spikes == 0:
            return VSR
        VSR.vs = (1.0 / VSR.n_spikes) * np.sqrt(
            sumc2 + sums2
        )  # standard vector strength computation
        # print('VSR nsp @ modf: ', VSR.frequency,  VSR.n_spikes, 'VS: ', VSR.vs, 'sumc2, s2: ', sumc2, sums2)
        VSR.Rayleigh = (
            VSR.n_spikes * VSR.vs * VSR.vs
        )  # Raleigh coefficient (Ashida et al, 2010 and many others: note some papers report 2nvs^2)
        VSR.pRayleigh = np.exp(
            -VSR.n_spikes * VSR.vs * VSR.vs
        )  # p value for n > 50 (see Ashida et al. 2010).
        VSR.circ_phaseMean = scipy.stats.circmean(phasev)
        VSR.circ_phaseSD = scipy.stats.circstd(phasev)  # radians
        VSR.circ_timeSD = VSR.circ_phaseSD / (2 * np.pi * freq)  # convert to time
        self.spikes = spikes
        # index_row = self.parent.selected_index_rows[0]
        # selected = self.parent.table_manager.get_table_data(
        #     index_row
        # )
        return VSR

    def print_VS(self, d):

        self.VS_colnames = f"Cell,Configuration,carrierfreq,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,AN_phase,AN_phasesd,maxArea,ninputs"

        line += f"{d.carrier_frequency:.1f},{d.carrier_frequency:.1f},{d.dmod:.1f},{d.dB:.1f},"
        line += f"{d.vs:.4f},"
        line += f"{d.n_spikes:d},"
        line += f"{d.circ_phaseMean:.4f},"
        line += f"{d.circ_phaseSD:.4f},"
        line += f"{d.Rayleigh:.4f},"
        line += f"{d.pRayleigh:.4e},"
        line += f"{d.an_vs:.4f},"
        line += f"{d.an_circ_phaseMean:.4f},"
        line += f"{d.an_circ_phaseSD:.4f},"
        line += f"{d.amax:.4f},"
        line += f"{d.n_inputs:d}"
        print(line)

    def vs_test__print(self, d):
        print(f"{'Vector Strength':>24s}: {d['r']:12.4f}")
        print(f"{'Spike Count':>24s}: {d['n']:12d}")
        print(
            f"{'mean phase (radians)':>24s}: {scipy.stats.circmean(d['ph'])/(2*np.pi):12.4f}"
        )
        print(f"{'SD phase':>24s}: {scipy.stats.circstd(d['ph'])/(2.*np.pi):12.4f}")
        print(f"{'Rayleigh':>24s}: {d['R']:12.4f}")
        print(f"{'p value':>24s}: {d['p']:12.4e}")
        print(f"{'d':>24s}: {d['d']:12.4f}")

    def compute_vs(self, freq=100.0, nsp=1000.0, sd=0.0, rg=None):
        if sd <= 1:
            sdn = (
                (2.0 / freq) + sd * rg.standard_normal(size=nsp) + np.pi
            )  # just express in temrs of time
        else:
            sdn = (1.0 / freq) * rg.uniform(size=nsp)  # uniform across the interval
        spikes = (
            np.cumsum(np.ones(nsp) * (1.0 / freq)) + sdn
        )  # locked spike train with jitter
        phsp = 2 * np.pi * freq * np.fmod(spikes, 1.0 / freq)
        return (spikes, phsp)


def plot_ticks(ax, data, color):
    ax.plot(np.tile(data, (2, 1)), [np.zeros_like(data), np.ones_like(data)], color)


def test_vs():
    from numpy.random import default_rng

    rg = default_rng(12345)
    freq = 100.0
    nsp = 1000
    sd = 0.002  # in seconds, unles > 1 then is uniform
    x = np.array(
        [
            0.0,
            0.00005,
            0.0001,
            0.0002,
            0.0003,
            0.0005,
            0.00075,
            0.001,
            0.002,
            0.003,
            0.004,
            0.0050,
            0.0075,
            2,
        ]
    )
    y = np.zeros_like(x)
    ph = np.zeros_like(x)
    vs = np.zeros_like(x)
    fig, ax = mpl.subplots(3, 1)
    ax = ax.ravel()
    VSC = VectorStrength()
    for i, sd in enumerate(x):
        spikes, phsp = VSC.compute_vs(freq, nsp, sd, rg)
        vsd = VSC.vector_strength(spikes, freq=freq)
        y[i] = vsd.circ_timeSD
        vs[i] = vsd.vs
        ph[i] = vsd.circ_phaseSD
    x[x == 2] = 0.020
    ax[0].plot(x, y, "o")
    ax[0].plot([0, np.max(x)], [0, np.max(x)], "--k", alpha=0.5)
    ax[1].plot(x, vs, "x-")
    ax[2].plot(x, ph, "s-")
    mpl.show()

    fig, ax = mpl.subplots(2, 1)
    ax = ax.ravel()
    plot_ticks(ax[0], phsp, "r-")
    plot_ticks(ax[1], spikes, "b-")
    mpl.show()
    print("mean isi: ", np.mean(np.diff(spikes)), "stim period: ", 1.0 / freq)
    vsd = VSC.vector_strength(spikes, freq=freq)
    print("VS: ", vsd.vs, "SD sec: ", vsd.circ_timeSD)
    print("circ_phase min max: ", np.min(vsd.circ_phase), np.max(vsd.circ_phase))
    bins = np.linspace(0, 2 * np.pi, 36)
    mpl.hist(vsd.circ_phase, bins)
    mpl.show()
