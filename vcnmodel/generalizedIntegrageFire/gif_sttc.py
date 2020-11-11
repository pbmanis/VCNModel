#!/usr/bin/env/python
"""'
compute STTC for different windows from the data in the specified file
Compares gif vs input model sttcs.

"""
from __future__ import print_function

import sys
import numpy as np
import matplotlib.pyplot as mpl
import pickle
import sttc
from numba import jit

SC = sttc.STTC()  # create instance
if len(sys.argv) == 1:
    cellno = 19
else:
    cellno = int(sys.argv[1])
      
@jit(nopython=True, cache=True)
def nb_spikethr(V, threshold, ref_ind):    

    spks = []

    t=0
    while (t<len(V)-1) :
    
        if (V[t] >= threshold and V[t-1] <= threshold) :
            spks.append(t)
            t += ref_ind
        t += 1
                
    spks = np.array(spks)
    return spks
      

gn = 0.0
fig, ax = mpl.subplots(nrows=3, ncols=1)
fig.set_size_inches(w=5., h=8., forward=True)


def compute(gn, ax, cell=19, hist=False, trace=False):

    fn = 'GFIT_VCN_c%02d_mGBC_gn=%6.3f_gifnoise.p' % (cell, gn)
    #fn = 'GFIT_Original_gifnoise.p'

    h = open(fn, 'rb')
    d = pickle.load(h)
    h.close()

    # data now in d
    st_m = d['V']
    st_e = d['ExpData']
    rate = d['gifpars']['dt']
    print('Rate: ', rate)
    tilewindow = 1.0
    refract = 1.5
    threshold = -20.
    ref_ind = int(refract/rate)      
    stexpt = nb_spikethr(st_e, threshold, ref_ind)*rate
    stmodel = nb_spikethr(st_m, threshold, ref_ind)*rate
    print('gn: %f  erate: %f  mrate: %f' % (gn, np.mean(np.diff(stexpt)), np.mean(np.diff(stmodel))))
    if trace:
        nx = len(d['ExpData'])
        ax[0].plot(d['time'][:nx], d['ExpData'], 'k-', linewidth=0.75)
        ax[0].plot(d['time'], d['V'], 'r-', linewidth=0.33)
        ax[0].set_xlim([0, 7500.])
    nbins = int(np.max(rate*len(st_m))/100.)
    # print nbins
    #x = np.array([stexpt, stmodel])
    x = [stexpt, stmodel]# print x
    # print x.shape
    if hist:
        ax[1].hist(x, nbins, histtype='step', stacked=True, fill=False, label=['Expt', 'Model'])
        ax[1].set_title('rates')
        ax[1].legend()


    SC.set_spikes(rate, stexpt, stmodel, tilewindow)
    tilerange = np.logspace(np.log10(0.1), np.log10(10), 20) # [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0]
    st = np.zeros(len(tilerange))
    for i, tw in enumerate(tilerange):
        st[i] = SC.calc_sttc(tw=tw)
     
    print('STTC : ' % st)
    stx = stexpt + np.random.randn(len(stexpt))*0.15  # jitter to show effect of timing on STTC.
    SC.set_spikes(rate, stexpt, stx, tilewindow)
    stxw = np.zeros(len(tilerange))
    for i, tw in enumerate(tilerange):
        stxw[i] = SC.calc_sttc(tw=tw)
    
    ax[2].plot(tilerange, st, 'o-', label='gn %.2f' % gn)
    #ax[1].plot(tilerange, stxw, 's-')

gns = [0., 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
# gns = [0., 0.5, 10.0] # for testing code... 
for i in range(len(gns)):
    if i == len(gns)-1:
        hist = True
        trace = True
    else:
        hist = False
        trace = False
    compute(gns[i], ax, cell=cellno, hist=hist, trace=trace)
ax[2].legend(fontsize=10)
mpl.savefig('gif_sttc.pdf')
#mpl.show()
#print st1
# SC.plot_sttc(st1, st2) # just plot first 100
