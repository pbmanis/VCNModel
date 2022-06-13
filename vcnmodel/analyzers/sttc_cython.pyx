#!python
#cython: language_level=3
import numpy as np
cimport numpy as np
cimport cython
np.import_array()

#from libc.stdio cimport printf
ctypedef np.int_t INT_t
ctypedef np.double_t DOUBLE_t

#@cython.boundscheck(False)
#@cython.wraparound(False)
def sttc_cython(
        np.ndarray[DOUBLE_t, ndim=1] time,
        np.ndarray[DOUBLE_t, ndim=1] spike_times, # data array (input)
        double rate, # sample rate single float input
        long int itile,
        double ta,
        np.ndarray[INT_t, ndim=1] tiles, # result tiles (output, must preallocate)
        ):
        cdef long int npts
        cdef long int ix0, ix1
        cdef double st


        npts = time.shape[0]  # number of

        for st in spike_times:
            ix0 = int(st/rate)-itile
            if ix0 < 0:
                ix0 = 0
            ix1 = int(st/rate)+itile
            if ix1 > int(len(tiles)):
                ix1 = int(len(tiles))
            tiles[ix0:ix1] += 1
#            for i in range(ix0, ix1):
#                tiles[i] = tiles[i] + 1
        ta = np.sum(tiles)/npts
        return ta, tiles
