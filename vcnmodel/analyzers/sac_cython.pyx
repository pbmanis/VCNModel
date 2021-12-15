#!python
#cython: language_level=3
cimport cython
import numpy as np
from libc.stdio cimport printf

def sac_cython(
        double[:,:] X, # data array (input) - spike event times, seconds
        long[:] event_lengths, # lengths of events in the array X
        double twin, # time window single float input, seconds
        double binw, # bin width, seconds
        double[:] yout, # result SAC (output)
        long ns, # total spike count
        ):
        cdef long int N, n, m, i, j
        cdef double xd, t0, t1, x0
        
        N = X.shape[0]  # number of trials
        t0 = -twin - (binw / 2.0)
        t1 = twin + (binw / 2.0)
        for m in range(N):  # for each trial
            for n in range(m+1, N):  # against each other trial, except:
                if m == n:
                    continue  # skip identical trains
                for i in range(int(event_lengths[n])):  # cross correlate all spikes in Ym
                    x0 = X[n,i]  # major speedup pulling this up to here!
                    for j in range(int(event_lengths[m])):
                        xd = X[m,j] - x0
                        if (xd > t0) & (xd < t1):
                            yout[ns] = xd  # store
                            ns = ns + 1  # keep track of storage location*/
        return yout, ns

#
# compute cross-correlation like SAC, between two different spike trains
#
def xac_cython(
        double[:,:] X, # data array (input)
        double[:,:] Y, # second array (input)
        long[:] event_lengths_x, # lengths of events in the array X
        long[:] event_lengths_y, # lengths of events in the array Y
        double twin, # time win single float input
        double binw,
        double[:] yout, # result XAC (output)
        long ns, # spike count #2 (output)
        ):
        cdef long int N, M, n, m, i, j
        cdef double xd, x0, t0, g1

        N = X.shape[0]  # number of trials
        M = Y.shape[0]  # number of trials
        t0 = -twin - (binw / 2.0)
        t1 = twin + (binw / 2.0)
        for n in range(N):  # for each trial
            for m in range(M):  # against each other trial, except:
                for i in range(int(event_lengths_x[n])):  # cross correlate all spikes in Yi
                    # look at all differences from i'th spike in m'th trial to spikes in nth trial
                    x0 = X[n][i]
                    for j in range(int(event_lengths_y[m])):
                        xd = Y[m][j] - x0
                        if (xd > t0) & (xd < t1):
                            yout[ns] = xd  # store result
                            ns = ns + 1  # keep track of storage location*/
        return yout, ns