#!python
#cython: language_level=3
cimport cython
import numpy as np
from libc.stdio cimport printf

def sac_cython(
        double[:,:] X, # data array (input)
        long[:] event_lengths, # lengts of events in the array X
        double twin, # time win single float input
        double binw, 
        double delay,
        double dur,
        double[:] yout, # result SAC (output)
        long ns, # total spike count
        ):
        cdef long int N, maxn, n, m, i, j, k
        cdef double xd

        N = X.shape[0]  # number of trials

        for m in range(N):  # for each trial
            for n in range(N):  # against each other trial, except:
                if m == n:
                    continue  # skip identical trains
                for i in range(int(event_lengths[m])):  # cross correlate all spikes in Ym
                    for j in range(int(event_lengths[n])):
                        if (X[m, i] <= delay) or (X[m, i] > (delay + dur)):
                            continue
                        if (X[n, j] <= delay) or (X[n, j] > (delay + dur)):
                            continue
                        xd = X[n,j] - X[m,i]
                        if (xd >= (-twin - (binw / 2.0))) & (xd <= ( twin + (binw / 2.0))):
                            yout[ns] = xd  # store
                            ns = ns + 1  # keep track of storage location*/
        return yout, ns

#
# compute cross-correlation like SAC, between two different spike trains
#
def xac_cython(
        double[:,:] X, # data array (input)
        double[:,:] Y, # second array (input)
        double twin, # time win single float input
        double binw, 
        double delay,
        double dur,
        double[:] yout, # result XAC (output)
        long[:] spcount1, # spike count, int, output
        long[:] spcount2, # spike count, int, output
        long ns, # spike count #2 (output)
        ):
        cdef long int N, M, maxn, n, m, i, j, k
        cdef double xd, x
        cdef double[:,:] Xa, Ya

        N = X.shape[0]  # number of trials
        M = Y.shape[0]  # number of trials

        Xa = np.nan*np.zeros_like(X)
        Ya = np.nan*np.zeros_like(Y)
        # first we trim out spikes outside the delay:delay+dur time window
        for i in range(N):
            k = 0  # count spikes in each trial
            for j in range(len(X[i])):
                if (X[i,j] >= delay) and (X[i,j] < (delay + dur)):
                    Xa[i,j] = X[i,j]
                    k += 1
            spcount1[i] = k

        for i in range(M):
            k = 0  # count spikes in each trial
            for j in range(len(Y[i])):
                if (Y[i,j] >= delay) and (Y[i,j] < (delay + dur)):
                    Ya[i,j] = Y[i,j]
                    k += 1
            spcount2[i] = k

        for n in range(N):  # for each trial
            for m in range(M):  # against each other trial, except:
                for i in range(len(Ya[m])):  # cross correlate all spikes in Yi
                    # look at all differences from i'th spike in m'th trial to spikes in nth trial
                    for j in range(len(Xa[n])):
                        xd = Xa[n][j] - Ya[m][i]
                        if (xd >= (-twin - (binw / 2.0))) & (xd <= ( twin + (binw / 2.0))):
                            yout[ns] = xd  # store result
                            ns = ns + 1  # keep track of storage location*/
        return yout, spcount1, spcount2, ns