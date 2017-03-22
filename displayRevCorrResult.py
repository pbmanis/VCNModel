"""
displayresult shows the results of a model_run. Just enter the filename in the fn field
"""

from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
from cnmodel.util import sound
import sac_campagnola as SAC
import os.path
import pycircstat as PCS
import pylibrary.PlotHelpers as PH
import pylibrary.Utility as pu
import spikestatistics as SPKS
import cell_config as SC
from collections import OrderedDict
from cycler import cycler
import random
from matplotlib import rc
rc('text', usetex=False)
import seaborn
#plt.style.use('seaborn-muted')


def norm(p, n):
    pmin = np.min(p)
    pmax = np.max(p)
    return (p[n] - pmin) / float(pmax - pmin)
    
# pattern_list = ['c18', 'c08', 'c09', 'no']
# patterns =  pattern_list #[0] # ['no']
# fn = {}
# fn[pattern_list[0]] = 'AN_Result_VCN_c18_delays_N020_040dB_4000.0_MS.p'
# fn[pattern_list[1]] = 'AN_Result_VCN_c18_VCN_c08_delays_N020_040dB_4000.0_MS.p'
# fn[pattern_list[2]] = 'AN_Result_VCN_c18_VCN_c09_delays_N020_040dB_4000.0_MS.p'
# fn[pattern_list[3]] = 'AN_Result_VCN_c18_delays_N020_040dB_4000.0_FM10.0_DM000_MS.p'

def clean_spiketimes(spikeTimes, mindT=0.7):
    """
    Clean up spike time array, removing all less than mindT
    spikeTimes is a 1-D list or array
    mindT is difference in time, same units as spikeTimes
    If 1 or 0 spikes in array, just return the array
    """
    if len(spikeTimes) > 1:
        dst = np.diff(spikeTimes)
        st = np.array(spikeTimes[0])  # get first spike
        sok = np.where(dst > mindT)
        st = np.append(st, [spikeTimes[s+1] for s in sok])
        # print st
        spikeTimes = st
    return spikeTimes

def plot_revcorr(p, ax, respike=False):
    d={}
    basepath = 'VCN_Cells/VCN_{0:3s}/Simulations/AN'.format(p)
    h = open(os.path.join(basepath, fn[p]))
    d[p] = pickle.load(h)
    h.close()
    seaborn.set_style('ticks')
    syninfo = SC.VCN_Inputs['VCN_{0:3s}'.format(p)]

#    print ('# cells: ', d[p].keys())
    if not respike:
        st = d[p]['spikeTimes']
    else:
        st = {}
        dt = (d[p]['time'][1]-d[p]['time'][0])
        for k in d[p]['somaVoltage'].keys():
            st[k] = pu.findspikes(d[p]['time'], d[p]['somaVoltage'][k], -20, dt=dt, mode='peak')
            st[k] = clean_spiketimes(st[k])
#    print (st[0])
    ndet1 = 0
    for n in st:
        ndet1 = ndet1 + len(st[n])   
    ndet0 = 0
    for n in d[p]['spikeTimes']:
        ndet0 = ndet0 + len(d[p]['spikeTimes'][n])
        
    print('Detected %d and %d spikes', ndet0, ndet1)
    an = d[p]['inputSpikeTimes']
    print('cell: ', p)
    print ('syninfo: ', syninfo)
    print('ans in syninfo: ', len(syninfo))
    print ('# an inputs: ', len(an[0]))
    print ("# trials: ", len(st))
    binw = 0.2
    width = 5.
    tx = np.arange(-width, 0, binw)
    #colors = build_colors(len(an[0]))
    #ax.set_prop_cycle( cycler('color', colr_list))
    ninputs = len(syninfo[1])
    sites = np.zeros(ninputs)
    amax = 0.
    for isite in range(ninputs): # precompute areas
        area = syninfo[1][isite][0]
        if area > amax:
            amax = area
        sites[isite] = int(np.around(area*SC.synperum2))
        
    
    for isite in range(ninputs):  # for each ANF input (get from those on first trial)
        nt = 0
        for trial in range(len(st)):  # sum across trials
#            print ('trials; ', trial, len(st), len(an), isite, len(an[trial]))
            ann = an[trial][isite]
            if trial == 0:
                C = SPKS.correlogram(st[trial], ann, width=width, bin=binw, T=None)
            else:
                C = C + SPKS.correlogram(st[trial], ann, width=width, bin=binw, T=None)
            nt = nt + 1
        nc = len(C)
        color = plt.cm.viridis(norm(sites, isite))
        ax.plot(tx, C[0:nc/2], color=color, label=('Input {0:2d} N={1:3d}'.format(isite, int(sites[isite]))))
    #plt.legend()    
#    PH.nice_plot(ax, spines=['left', 'bottom'], position=0)
    seaborn.despine(ax=ax)
    PH.adjust_spines(ax, distance=0)
#    ax.set_ylim(0, 0.20)
    ax.set_ylabel('R', fontsize=10)
    ax.set_xlabel('T (ms)', fontsize=10)
    
    ax.set_title('VCN_{0:3s} [{1:d}-{2:d}] Amax={3:.1f}'.format(p, int(np.min(sites)), int(np.max(sites)), amax), y=0.9, x=0.02,
        horizontalalignment='left', fontsize=8)


def plot_SAC():
    #fig, ax = plt.subplots(len(pattern_list)+1, 3)
    fig = plt.figure()

    gs_outer = gridspec.GridSpec(4, 1)  # 4 rows, one column
    gs_inner = []
    ax = {}
    for i, outer in enumerate(gs_outer):
        gs_inner.append(gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=outer))
        ax1 = plt.Subplot(fig, gs_inner[-1][:-1, 0])
        fig.add_subplot(ax1)
        ax2 = plt.Subplot(fig, gs_inner[-1][3, 0])
        fig.add_subplot(ax2)
    
        ax3 = plt.Subplot(fig, gs_inner[-1][:, 1])
        fig.add_subplot(ax3)
        ax4 = plt.Subplot(fig, gs_inner[-1][:, 2])
        fig.add_subplot(ax4)
        ax5 = plt.Subplot(fig, gs_inner[-1][:, 3])
        fig.add_subplot(ax5)
        ax[i] = [ax1, ax2, ax3, ax4, ax5]
        for a in ax[i]:
            PH.nice_plot(a)
        #PH.noaxes(ax1, 'x')

    #fig.tight_layout()
    # plt.show()
    # exit()


    # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
    #print 'data keys: ', d.keys()
    #print len(d['time'])
    fmod = d[patterns[0]]['stimInfo']['fmod']  # modulation frequency
    sacht = 2.5
    sacmax = 100
    phaseht = 15

    #print d['spikeTimes']
    for j, pattern in enumerate(patterns):
        w = d[patterns[j]]['stimWaveform'][0][0].generate()
        t = d[patterns[j]]['stimWaveform'][0][0].time
        ax[j][1].plot(t, w)  # stimulus underneath
        for i, st in enumerate(d[pattern]['spikeTimes'].keys()):
            ax[j][0].plot(d[pattern]['spikeTimes'][st]/1000., i*np.ones(len( d[pattern]['spikeTimes'][st])), 'o', markersize=2.5, color='b')
            inputs = len(d[pattern]['inputSpikeTimes'][st])
            # for j in range(inputs):
            #     t = (d['ref']['inputSpikeTimes'][st][j])/1000.
            #     y = (i+0.1+j*0.05)*np.ones(len(t))
            #     ax[0,].plot(t, y, 'o', markersize=2.5, color='r')
    #for i, st in enumerate(d['ref']['somaVoltage'].keys()):
    #    ax[2,].plot(d['ref']['time'], d['ref']['somaVoltage'][st], color='k')
        ax[j][0].set_xlim((0, 0.8))
        ax[j][1].set_xlim((0, 0.8))


    sac = SAC.SAC()
    xf = []
    dt = 0.1
    u = 2.0*np.pi/(1000.0/fmod)
    for j, pattern in enumerate(patterns):
        X = []
        R = {}
        pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                'dur': 1000., 'ddur': 200.}
        pars_x =     pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                'dur': 1000., 'ddur': 200.}
        tsynch = np.arange(0., 1000., dt)
        x_spl = []
        y_spl = []
        Nspikes = 0
        for i, st in enumerate(d[pattern]['spikeTimes'].keys()):

            X.append(d[pattern]['spikeTimes'][st])
            xw = np.zeros(tsynch.shape[0])
    #        print [int(xx/dt) for xx in X[-1]]
            Nspikes += np.shape(X[-1])[0]
            # calculate vector strength as well.
            x_spl.append(np.cos(np.array(X[-1])*u))
            y_spl.append(np.sin(np.array(X[-1])*u))
        x_spl = np.hstack(x_spl)
        y_spl = np.hstack(y_spl)
        vs = np.sqrt(np.sum(x_spl)**2 + np.sum(y_spl)**2)/Nspikes
        th = np.arctan2(y_spl, x_spl)
        p, z = PCS.rayleigh(th, w=None, d=None, axis=None)
        if j == 0:
            X0 = X  # save the X
        yh, bins = sac.SAC_asm(X, pars)
        print ('mean vs: for %s at fMod = %.2f:  %f angle: %f' % (pattern, fmod, vs, th.mean()))
        print ('rayleigh: p=%f  z=%f (p > 0.05 means data is uniformly distributed)' % (p, z))
        ax[j][2].bar(bins[:-1], yh)
        ax[j][2].set_xlim((-sacmax, sacmax))
        ax[j][2].set_ylim((0, sacht))
        phasebinnum = 101
        phasebins = np.linspace(-0.5, 0.5, num=phasebinnum)
        thhist, thbins = np.histogram(th/np.pi, bins=phasebins, density=False)
        ax[j][3].bar(thbins[:-1]+0.5, thhist, width=1.0/phasebinnum)
        ax[j][3].set_ylim((0, phaseht))

        if j == 0:
            continue
        ycc, ccbins = sac.XAC(X0, X, pars_x, trialbased=True)
        ax[j][4].bar(ccbins[:-1], ycc)
        ax[j][4].set_xlim((-sacmax, sacmax))
    
        # now cross-correlate data...
    

    #PH.nice_plot(ax.flatten().tolist())
    plt.show()

if __name__ == '__main__':
    gbcs = [  8,     9,   17,   18,    19, 20, 21, 22]
    thrs = [-20.,  -30., -35., -25.,                 ]
    SR = 'MS'
    fn = {}
    for g in gbcs:
        fn['c{0:02d}'.format(g)] = 'AN_Result_VCN_c{0:02d}_delays_N050_040dB_4000.0_{1:2s}.p'.format(g, SR)


    patterns = fn.keys()
    print ('patterns: ', patterns)
    p = 'c09'
    # fig, ax = plt.subplots(2, 4, figsize=(10,6))
    lmar = 0.1
    rmar = -0.01
    tmar = 0.04
    hspace = 0.04
    vspace = 0.1
    ncols = 4
    nrows = len(gbcs)/ncols
    gwid = (1.0-lmar-rmar)/ncols - hspace
    wid = gwid-hspace
    ght = (1.0-2*tmar)/nrows
    ht = ght - vspace
    yp = np.linspace(tmar, 1.0-tmar-ght, nrows)[::-1]
    lp = np.linspace(lmar, 1.0-lmar-gwid, ncols)
    plots  = ['rcorr']
    sizer = OrderedDict([])
    k = 0
    
    for ir, r in enumerate(range(nrows)):
        for ic, c in enumerate(range(ncols)):
            pname = ('VCN_c{0:02d}'.format(gbcs[k]))
            sizer[pname] = [lp[ic], wid, yp[ir], ht]
            k = k + 1
    # sizer = OrderedDict([('VCN_c08', [l1, wid, yp[0], ht]), ('VCN_c09', [l1, wid, yp[1], ht]),
    #         ('VCN_c17', [l1, wid, yp[2], ht]), ('VCN_c18', [l1, wid, yp[3], ht]),
    #         ('VCN_c19', [l2, wid, yp[0], ht]), ('VCN_c20', [l2, wid, yp[1], ht]),
    #         ('VCN_c21', [l2, wid, yp[2], ht]), ('VCN_c22', [l2, wid, yp[3], ht]),
    # ])  # dict elements are [left, width, bottom, height] for the axes in the plot.
    gr = [(a, a+1, 0, 1) for a in range(0, len(sizer.keys()))]   # just generate subplots - shape does not matter
    axmap = OrderedDict(zip(sizer.keys(), gr))
    P = PH.Plotter(rcshape=sizer, label=False, figsize=(10 , 6))
    P.figure_handle.suptitle('Reverse Correlations, AN types={0:2s}'.format(SR))
    for l in P.axlabels:
        l.set_text(' ')
    for i, g in enumerate(sizer.keys()):
        p = g[4:]
        plot_revcorr(p, P.axdict[g], respike=True)
    plt.show()
    exit()
    