#!/usr/bin/env python
"""
Spike time tiling calculation, from:
Cutts, C.S. and Eglen, S.J., "Detecting Pairwise Correlations in Spike Trains:
An Objective Comparison of Methods and Application to the Study of
Retinal Waves", The Journal of Neuroscience, October 22, 2014,
34(43):14288-14303

Implementation by P.B. Manis, Ph.D., UNC Chapel Hill
November, 2017

The cython version does not work corectly. The python version, surpisingly, is
the fastest (compared to cython and numba)

"""
import functools
import numpy as np
import matplotlib.pyplot as mpl
from numba import jit
import pyximport
from typing import Union
import time


# The cython version is not correct. Since the python version
# is actually faster, we ditched it.
from vcnmodel.analyzers import sttc_cython
# pyximport.install()


def time_func(func):
    """
    A decorator to show execution time of functions.
    Place inside (after) winprint if using
    Output is to terminal.
    """

    @functools.wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        #print(f"Starting : {func.__name__!r}")
        start_time = time.perf_counter()  # 1
        value = func(self, *args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"F{func.__name__!r} itook {run_time:.4f} secs")
        return value

    return wrapper_timer


# @time_func
@jit(
    nopython=True, cache=True,
)
def nb_tiletimes(times, st, rate, itile):
    npts = times.shape[0]
    tiles = np.zeros(npts)
    #    itile = int(self.tilewindow/rate)
    for i in range(len(st)):
        ix0 = int(st[i] / rate) - itile
        if ix0 < 0:
            ix0 = 0
        ix1 = int(st[i] / rate) + itile
        if ix1 > len(tiles):
            ix1 = len(tiles)
        tiles[ix0:ix1] = 1
    ta = np.sum(tiles) / npts
    return ta, tiles


class STTC:
    def __init__(self, seed=0, engine="python"):
        np.random.seed(0)
        self.engine = engine

        # self.set_spikes(times, rate, st1, st2, dt)

    def set_spikes(self, rate, st1, st2, tilewindow):
        self.tilewindow = tilewindow
        self.sample_rate = rate
        self.st1 = np.array(st1, dtype=float)
        self.st2 = np.array(st2, dtype=float)
        npts = np.max([np.max(st1), np.max(st2)]) + self.tilewindow * 2
        npts = int(npts / self.sample_rate)
        self.times = np.arange(npts) * self.sample_rate

    def calc_ccf_sttc(self, corrwindow=[-5.0e-3, 1.0e-3], binwidth=0.2e-3):
        deltaT = np.arange(corrwindow[0], corrwindow[1], binwidth)
        self.ccf = np.zeros_like(deltaT)
        self.original_st2 = self.st2.copy()
        for i, t in enumerate(deltaT):
            self.st2 = self.original_st2 + t
            self.ccf[i] = self.calc_sttc()
        return self.ccf, deltaT

    def calc_sttc(self, tw=None, engine="python"):

        if tw is None:
            tw = self.tilewindow
        else:
            self.tilewindow = tw
        self.st1 = np.sort(self.st1)
        self.st2 = np.sort(self.st2)
        itile = int(self.tilewindow / self.sample_rate)
        if self.engine == "python":
            ta, tatiles = self.tiletimes_python(
                self.times, self.st1, self.sample_rate, itile
            )
            tb, tbtiles = self.tiletimes_python(
                self.times, self.st2, self.sample_rate, itile
            )
        elif self.engine == "numba":
            ta, tatiles = nb_tiletimes(self.times, self.st1, self.sample_rate, itile)
            tb, tbtiles = nb_tiletimes(self.times, self.st2, self.sample_rate, itile)
        # elif self.engine == "cython":
        #     ta, tatiles = self.tile_times_cython(
        #         self.times, self.st1, self.sample_rate, itile
        #     )
        #     tb, tbtiles = self.tile_times_cython(
        #         self.times, self.st2, self.sample_rate, itile
        #     )
        else:
            raise ValueError("Engine is not recognized - check call to calc_sttc")
        pa = self.proportion(self.st1, tbtiles)
        pb = self.proportion(self.st2, tatiles)
        #        print('pa, pb: ', pa, pb)
        sttc = 0.5 * (((pa - tb) / (1 - pa * tb)) + ((pb - ta) / (1 - pb * ta)))
        self.tatiles = tatiles
        self.tbtiles = tbtiles
        return sttc

    @time_func
    def tile_times_cython(
        self, times: np.array, st: Union[list, np.array], rate: float, itile: int
    ):
        ta = 0.0
        tiles = np.zeros_like(times)

        sttc_cython.sttc_cython(
            np.array(times),  # data array (input)
            np.array(st),
            self.sample_rate,
            itile,
            ta,
            tiles,
        )
        return ta, tiles

    @time_func
    def tiletimes_python(self, times, spike_times, rate, itile):
        """
        Compute the total time tiled in spike train st
        """
        ta = 0.0
        tiles = np.zeros(times.shape[0])

        for st in spike_times:
            ix0 = int(st / rate) - itile
            if ix0 < 0:
                ix0 = 0
            ix1 = int(st / rate) + itile
            if ix1 > len(tiles):
                ix1 = len(tiles)
            tiles[ix0:ix1] = 1
        ta = np.sum(tiles) / len(tiles)
        return ta, tiles

    def proportion(self, st, tile):
        ist = [int(np.floor(sti / self.sample_rate)) for sti in st]
        # ist = [x for x in ist if x < len(tile)]
        #  print(len(tile), self.sample_rate, np.max(ist))
        p = np.sum(tile[ist[:-1]]) / len(st)
        return p

    def tests(
        self,
        i, j, axes,
        distribution="exp",
        pdelete=0.0,
        independent=True,
        dither=0.0,
        tilewindow=1e-2,
        spikerate=20.0,
        duration=1.0,
    ):

        assert distribution in ["exp", "exponential", "poisson", "regular"]
        print("setting up: ", distribution, independent)
        samplerate = 0.0001  # seconds
        nspikes = int(spikerate*duration)
        print(f"nspikes: {nspikes:d}")
        if distribution in ["exp", "exponential"]:
            print("exp")
            st1 = np.random.exponential(1.0 / spikerate, nspikes)
            st1 = np.cumsum(st1)
        elif distribution == "regular":
            print("regular")
            st1 = np.linspace(
                1. / samplerate, duration, int(spikerate*duration)
            )
        elif distribution == "poisson":
            print("poisson")
            st1 = np.random.poisson(1.0 / spikerate, nspikes)
            st1 = np.cumsum(st1)

        if independent:
            print("indepdent", independent)
            if distribution == "regular":
                st2 = np.linspace(
                    int(1.0 / samplerate), duration, int(spikerate*duration)
                ) +0.005
            elif distribution in ["exp", "exponential"]:
                st2 = np.random.exponential(1.0 / spikerate, nspikes,)
                st2 = np.cumsum(st2)
            elif distribution in ["poisson"]:
                st2 = np.random.poisson(1.0 / spikerate, nspikes)
                st2 = np.cumsum(st2)
        else:
            print("identical", independent)
            st2 = st1.copy()
        if pdelete != 0:
            st2 = np.random.choice(
                st2, int((1.0 - pdelete) * st1.shape[0]), replace=False
            )
        if dither > 0:
            print("dither")
            st2 = st2 + np.random.randn(len(st2)) * dither
        else:
            print("no dither")
        if not independent and dither == 0.0:
            assert np.array_equal(st1, st2)
        else:
            assert not np.array_equal(st1, st2)
        print("passed independent test of either equality or inequality")
        print('len st1, st2: ', len(st1), len(st2), np.max(st1), np.max(st2))

        self.set_spikes(samplerate, st1, st2, tilewindow=tilewindow)
        sttc_ccf, wins_ccf = self.calc_ccf_sttc(corrwindow=[0, tilewindow], binwidth=0.001)
        sttc = self.calc_sttc()
        print(
            "# of spikes in spike train 1: {0:d}, in spike train 2: {1:d} ".format(
                st1.shape[0], st2.shape[0]
            )
        )
        print("STTC value: {0:.3f} ".format(sttc))
        # self.plot_sttc(st1, st2, sttc_ccf, wins_ccf, i, j, axes)

    def plot_sttc(self, st1, st2, sttc_ccf, wins_ccf, i, j, axes):
        axes[i, j].eventplot(
            [st1,  st2],
            linelengths=[0.75, 0.75],
            lineoffsets=[0, 1],
            colors=["b", "m"],
        )

        axes[i+1, j].plot(wins_ccf, sttc_ccf, 'k-')


if __name__ == "__main__":

    distributions = ["regular", "poisson", "exponential"]
    fig, axes = mpl.subplots(6, 3)
    for i, dist in enumerate(distributions):
        for j, eng in enumerate(['python', 'numba']):
            S = STTC(seed=0, engine=eng)
            S.tests(
                i*2, j, axes,
                distribution=dist,
                pdelete=0.0,
                independent=True,
                dither=0.001,
                tilewindow=0.020,
                spikerate=20,
                duration=1.0,
            )
            if i == 0:
                axes[i, j].set_title(eng)
            if j == 0:
                axes[i*2, j].set_ylabel(dist)

    mpl.show()
