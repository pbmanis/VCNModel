#!/usr/bin/env python
"""
Spike time tiling calculation, from:
Cutts, C.S. and Eglen, S.J., "Detecting Pairwise Correlations in Spike Trains:
An Objective Comparison of Methods and Application to the Study of
Retinal Waves", The Journal of Neuroscience, October 22, 2014, 34(43):14288-14303

Implementation by P.B. Manis, Ph.D., UNC Chapel Hill
November, 2017

"""
import functools
import numpy as np
import matplotlib.pyplot as mpl
from numba import jit
import pyximport
from typing import Union
import time
pyximport.install()
from vcnmodel.analyzers import sttc_cython

def time_func(func):
    """
    A decorator to show execution time of functions.
    Place inside (after) winprint if using
    Output is to terminal.
    """

    @functools.wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        print(f"Starting : {func.__name__!r}")
        start_time = time.perf_counter()  # 1
        value = func(self, *args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer


@time_func
@jit(nopython=True, cache=True,)
def nb_tiletimes(times, st, rate, itile):
    npts = times.shape[0]
    tiles = np.zeros(npts)
#    itile = int(self.tilewindow/rate)
    for i in range(len(st)):
        ix0 = int(st[i]/rate)-itile
        if ix0 < 0:
            ix0 = 0
        ix1 = int(st[i]/rate)+itile
        if ix1 > len(tiles):
            ix1 = len(tiles)
        tiles[ix0:ix1] = 1
    ta = np.sum(tiles)/npts
    return ta, tiles



class STTC():
    def __init__(self, seed=0, engine="python"):
        np.random.seed(0)
        self.engine = engine
        
        #self.set_spikes(times, rate, st1, st2, dt)
    
    def set_spikes(self, rate, st1, st2, tilewindow):
        self.tilewindow = tilewindow
        self.sample_rate = rate
        self.st1 = np.array(st1, dtype=float)
        self.st2 = np.array(st2, dtype=float)
        npts = np.max([np.max(st1), np.max(st2)]) + self.tilewindow*2
        npts = int(npts/self.sample_rate)
        print('npts: ', npts)
        self.times = np.arange(npts)*self.sample_rate
    
    def calc_ccf_sttc(self, corrwindow=[-5., 1.], binwidth=0.1):
        deltaT = np.arange(corrwindow[0], corrwindow[1], binwidth)
        self.ccf = np.zeros_like(deltaT)
        self.original_st2 = self.st2.copy()
        for i, t in enumerate(deltaT):
            self.st2 = self.original_st2 + t
            self.ccf[i] = self.calc_sttc()
        return self.ccf
        
    def calc_sttc(self, tw=None, engine="python"):
        
        if tw is None:
            tw = self.tilewindow
        else:
            self.tilewindow = tw
        st1 = np.sort(self.st1)
        st2 = np.sort(self.st2)
        itile =  int(self.tilewindow/self.sample_rate)
        if self.engine == "python":
            ta, tatiles = self.tiletimes_python(self.times, st1, self.sample_rate, itile)
            tb, tbtiles = self.tiletimes_python(self.times, st2, self.sample_rate, itile)
        elif self.engine == "numba":
            ta, tatiles = nb_tiletimes(self.times, st1, self.sample_rate, itile)
            tb, tbtiles = nb_tiletimes(self.times, st2, self.sample_rate, itile)
        elif self.engine == "cython":
            print('st1')
            ta, tatiles = self.tile_times_cython(self.times, st1, self.sample_rate, itile)
            print('st2')
            tb, tbtiles = self.tile_times_cython(self.times, st2, self.sample_rate, itile)
        else:
            raise ValueError("Engine is not recognized - check call to calc_sttc")
        pa = self.proportion(st1, tbtiles)
        pb = self.proportion(st2, tatiles)
#        print('pa, pb: ', pa, pb)
        sttc = 0.5*(((pa-tb)/(1-pa*tb)) + ((pb-ta)/(1-pb*ta)))
        self.tatiles = tatiles
        self.tbtiles = tbtiles
        return sttc
    
    @time_func
    def tile_times_cython(self, times:np.array, st:Union[list, np.array], rate:float, itile:int):
        ta = 0.
        tiles = np.zeros_like(times)
        print('tiles: ', tiles[:10])
        print('st: ', st[:10])
        sttc_cython.sttc_cython(
                np.array(times), # data array (input)
                np.array(st),
                self.sample_rate,
                itile,
                ta,
                tiles
                )
        return ta, tiles

    @time_func
    def tiletimes_python(self, times, st, rate, itile):
        """
        Compute the total time tiled in spike train st
        """
        ta = 0.
        tiles = np.zeros(times.shape[0])

        for i in range(len(st)):
            ix0 = int(st[i]/rate)-itile
            if ix0 < 0:
                ix0 = 0
            ix1 = int(st[i]/rate)+itile
            if ix1 > len(tiles):
                ix1 = len(tiles)
            tiles[ix0:ix1] = 1
        ta = np.sum(tiles)/len(tiles)
        return (ta, tiles)
    
    def proportion(self, st, tile):
        ist = [int(np.floor(sti/self.sample_rate)) for sti in st]
        # ist = [x for x in ist if x < len(tile)]
        print(len(tile), self.sample_rate, np.max(ist))
        p = np.sum(tile[ist[:-1]])/len(st)
        return p

    def tests(self, distribution='exp', pdelete=0., independent=True, dither=0., tilewindow=1.0, nspikes=100):
        
        assert distribution in ['exp', 'exponential', 'poisson', 'regular']
        print('setting up: ', distribution, independent)
        samplerate = 0.1 # ms
        spikerate = 0.001 # firing rate
        if distribution in ['exp', 'exponential']:
            print("exp")
            st1 = np.random.exponential(1./spikerate, nspikes)
            st1 = np.cumsum(st1)
        elif distribution == 'regular':
            print("regular")
            st1 = np.linspace(int(10./samplerate),
                int(9000./samplerate), int(10./samplerate))
        elif distribution == 'poisson':
            print("poisson")
            st1 = np.random.poisson(1./spikerate, nspikes)
            st1 = np.cumsum(st1)
        
        if independent:
            print("indepdent", independent)
            st2 = np.random.exponential(1./spikerate, nspikes, dtype=float)
            st2 = np.cumsum(st2)
        else:
            print("identical", independent)
            st2 = st1.copy()
        if pdelete != 0:
            st2 = np.random.choice(st2,
                    int((1.0-pdelete)*st1.shape[0]), replace=False)
        if dither > 0:
            print("dither")
            st2 = st2 + np.random.randn(len(st2))*dither
        else:
            print("no dither")
        if not independent:
            assert(np.array_equal(st1, st2))
        else:
            assert(not np.array_equal(st1,st2))
        print("passed indepdent test of either equality or inequality")
#        print('len st1, st2: ', len(st1), len(st2), np.max(st1), np.max(st2))
        self.set_spikes(samplerate, st1, st2, tilewindow=tilewindow)
        sttc = self.calc_sttc()
        print('# of spikes in spike train 1: {0:d}, in spike train 2: {1:d} '.format(st1.shape[0], st2.shape[0]))
        print('STTC value: {0:.3f} '.format(sttc))
        self.plot_sttc(st1, st2)

    
    def plot_sttc(self, st1, st2):
        mpl.figure()
        st1x = np.repeat(st1, 3)
        st1x[2::3] = np.NaN
        st1y = 1 + 0*st1x
        st1y[1::3] = st1y[0::3] + 0.25
        
        st2x = np.repeat(st2, 3)
        st2x[2::3] = np.NaN
        st2y = 1.2 + 0*st2x
        st2y[1::3] = st2y[0::3] + 0.25
        mpl.plot(st1x, st1y, 'b-', linewidth=1.5)
        mpl.plot(st2x, st2y, 'r-', linewidth=1.5)
        txa = np.argwhere(self.tatiles == 0.)
        txb = np.argwhere(self.tbtiles == 0.0)
        ta = self.tatiles
        ta[txa] = np.nan
        tb = self.tbtiles
        tb[txb] = np.nan
        mpl.plot(self.times, ta*0.9, 'b-')
        mpl.plot(self.times, tb*0.95, 'r-')
        
        # mpl.ylim([0., 2.])
        # #mpl.show()
    

if __name__ == '__main__':
    
    S = STTC(seed=0, engine="python")
    S.tests(distribution='regular', pdelete=0.0, independent=False, dither=.0,
        tilewindow=2.0, nspikes=1000)
    S.tests(distribution='poisson', pdelete=0.0, independent=False, dither=.0,
        tilewindow=2.0, nspikes=1000)
    S.tests(distribution='exponential', pdelete=0.0, independent=False, dither=.0,
        tilewindow=2.0, nspikes=1000)
    mpl.show()
