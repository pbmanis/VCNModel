#!python
#cython: language_level=3
cimport cython
import numpy as np
from libc.stdio cimport printf


def sttc_cython(
        double[:] time, 
        double[:] spike_times, # data array (input)
        double rate, # sample rate single float input
        long int itile,
        double ta,
        double[:] tiles, # result tiles (output, must preallocate)
        ):
        cdef long int npts, i, ix0, ix1

        npts = time.shape[0]  # number of trials

        for i in range(len(spike_times)):
            ix0 = int(spike_times[i]/rate)-itile
            if ix0 < 0:
                ix0 = 0
            ix1 = int(spike_times[i]/rate)+itile
            if ix1 > len(tiles):
                ix1 = len(tiles)
            tiles[ix0:ix1] = 1
        ta = np.sum(tiles)/npts
        return ta, tiles
