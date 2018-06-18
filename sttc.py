#!/usr/bin/env python
from __future__ import print_function

"""
Spike time tiling calculation, from:
Cutts, C.S. and Eglen, S.J., "Detecting Pairwise Correlations in Spike Trains: An Objective Comparison of Methods and Application to the Study of Retinal Waves", The Journal of Neuroscience, October 22, 2014, 34(43):14288-14303

Implementation by P.B. Manis, Ph.D., UNC Chapel Hill
November, 2017

"""

import numpy as np
import matplotlib.pyplot as mpl
from numba import jit

@jit(nopython=True, cache=True,)
def nb_tiletimes(time, st, rate, itile):
    npts = time.shape[0]
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
    def __init__(self, seed=0):
        np.random.seed(0)
        
        #self.set_spikes(time, rate, st1, st2, dt)
    
    def set_spikes(self, rate, st1, st2, tilewindow):
        self.tilewindow = tilewindow
        self.sample_rate = rate
        self.st1 = st1
        self.st2 = st2
        npts = np.max([np.max(st1), np.max(st2)]) + self.tilewindow*2
        npts = int(npts/self.sample_rate)
        print('npts: ', npts)
        self.time = np.arange(npts)*self.sample_rate
    
    def calc_sttc(self, tw=None):
        
        if tw is None:
            tw = self.tilewindow
        else:
            self.tilewindow = tw
        st1 = np.sort(self.st1)
        st2 = np.sort(self.st2)
        ta, tatiles = self.tiletimes(self.time, self.sample_rate, st1)
        tb, tbtiles = self.tiletimes(self.time, self.sample_rate, st2)
        pa = self.proportion(st1, tbtiles)
        pb = self.proportion(st2, tatiles)
#        print('pa, pb: ', pa, pb)
        sttc = 0.5*(((pa-tb)/(1-pa*tb)) + ((pb-ta)/(1-pb*ta)))
        self.tatiles = tatiles
        self.tbtiles = tbtiles
        return sttc
    
    def tiletimes(self, time, rate, st):
        """
        Compute the total time tiled in spike train st
        """
        itile = int(self.tilewindow/rate)
        #ta, tiles = nb_tiletimes(time, st, rate, itile)  # not faster with numba... 
        ta = 0.
        tiles = np.zeros(time.shape[0])

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
        ist = [int(sti/self.sample_rate) for sti in st]
        p = np.sum(tile[ist])/len(st)
        return p

    def tests(self, distribution='exp', pdelete=0., independent=True, dither=0., tilewindow=1.0):
        
        assert distribution in ['exp', 'exponential', 'poisson', 'regular']
        samplerate = 0.1 # ms
        spikerate = 0.001 # firing rate
        nspikes = 100 # number of spikes to test
        if distribution in ['exp', 'exponential']:
            st1 = np.random.exponential(1./spikerate, nspikes)
            st1 = np.cumsum(st1)
        elif distribution == 'regular':
            st1 = np.linspace(int(10./samplerate),
                int(9000./samplerate), int(10./samplerate))
        elif distribution == 'poisson':
            st1 = np.random.poisson(1./spikerate, nspikes)
            st1 = np.cumsum(st1)
        
        if independent:
            st2 = np.random.exponential(1./spikerate, nspikes)
            st2 = np.cumsum(st1)
        else:
            st2 = st1
        st2 = np.random.choice(st2,
                    int((1.0-pdelete)*st1.shape[0]), replace=False)
        if dither > 0:
            st2 = st2 + np.random.randn(len(st2))*dither
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
        mpl.plot(self.time, ta*0.9, 'b-')
        mpl.plot(self.time, tb*0.95, 'r-')
        
        mpl.ylim([0., 2.])
        #mpl.show()
    

if __name__ == '__main__':
    
    S = STTC(seed=0)
    S.tests(distribution='regular', pdelete=0.3, independent=False, dither=2.0,
        tilewindow=2.0)
    
    
        
        
    
        
        