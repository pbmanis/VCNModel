


import sys
import os.path
import pickle
#import neuronvis.sim_result as sr

import pylibrary.Utility as pu  # access to spike finder routine

import time

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import pylibrary.pyqtgraphPlotHelpers as pgh
import analyze_run as ar
import calyxPlots as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pylibrary.PlotHelpers as PH
from collections import OrderedDict

baseName = 'VCN_Cells'
cell = 'VCN_c08'
#filename = 'AN_Result_VCN_c08_delays_N050_040dB_4000.0_HS.p'
#filename = 'AN_Result_VCN_c08_delays_N005_040dB_4000.0_MS.p'

#synfile_template = 'AN_Result_VCN_c18_reparented755V2Syn%03d_N005_040dB_4000.0_ 3.p'
#synfile_template = 'AN_Result_VCN_c18_reparented755V2_Syn%03d_N002_040dB_4000.0_HS.p'
synfile_template = 'AN_Result_%s_Syn%03d_N005_040dB_4000.0_%2s.p'

def readFile(filename):
    f = open(filename, 'r')
    d = pickle.load(f)
    f.close()
    if isinstance(d, list):
        stimInfo = d[0]
        spikeTimes = d[1]
        inputSpikeTimes = d[2]
    else:
        stimInfo = d['stimInfo']
        spikeTimes = d['spikeTimes']
        inputSpikeTimes = d['inputSpikeTimes']

    return(spikeTimes, inputSpikeTimes, stimInfo, d)


def ANfspike(spikes, stime, nReps):
    nANF = len(spikes[0])
    sl1 = np.full(nReps*nANF, np.nan)
    sl2 = np.full(nReps*nANF, np.nan)
    for r in range(nReps):
        for i in range(nANF):
            # print spikes[r][i]
            # print stime
            rs = np.where(spikes[r][i] > stime)[0]
            # print spikes[r][i][rs[0]]
            # return
            if len(spikes[r][i]) > 0:
                sl1[r*nANF+i] = spikes[r][i][rs[0]]
            if len(spikes[r][i]) > 1:
                sl2[r*nANF+i] = spikes[r][i][rs[1]]
#    print fsl
    print 'Auditory nerve: '
    print '   mean First spike latency:  %8.3f ms stdev: %8.3f  (N=%3d)' % (np.nanmean(sl1), np.nanstd(sl1),
        np.count_nonzero(~np.isnan(sl1)))
    print '   mean Second spike latency: %8.3f ms stdev: %8.3f  (N=%3d)' % (np.nanmean(sl2), np.nanstd(sl2),
        np.count_nonzero(~np.isnan(sl2)))

def getFirstSpikes(spikes, stime, nReps):
    sl1 = np.full(nReps, np.nan)
    sl2 = np.full(nReps, np.nan)
    for r in range(nReps):
        rs = np.where(spikes[r] > stime)[0]  # get spike times post stimulus onset
        if len(spikes[r]) > 0 and len(rs) > 0:
            sl1[r] = spikes[r][rs[0]]
        if len(spikes[r]) > 1 and len(rs) > 0:
            sl2[r] = spikes[r][rs[1]]  # second spikes
    return(sl1, sl2)    

def CNfspike(spikes, stime, nReps):
    # print spikes
    print 'stime: ', stime
    sl1, sl2 = getFirstSpikes(spikes, stime, nReps)

    print 'Cochlear Nucleus Bushy Cell: '
    print '   mean First spike latency:  %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl1), np.nanstd(sl1),
        np.count_nonzero(~np.isnan(sl1)))
    print '   mean Second spike latency: %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl2), np.nanstd(sl2), 
        np.count_nonzero(~np.isnan(sl2)))
    return(sl1, sl2)


def plot1(spikeTimes, inputSpikeTimes, stimInfo):
    win = pgh.figure(title='AN Inputs')
    layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
    # flatten spike times
    spt = []
    nReps = stimInfo['nReps']
    stime = stimInfo[ 'pip_start'][0]*1000.
    for r in range(nReps):
        spt.extend(spikeTimes[r])

    ANFs = [0]
    anf = []
    for r in range(nReps):
        for n in range(len(ANFs)):
            anf.extend(inputSpikeTimes[r][ANFs[n]])  # pick the ANF
    #print stimInfo
    srGroups = ['', 'LS', 'MS', 'HS']
    (anpsth, anbins) = np.histogram(anf, bins=np.linspace(0., 250., 250), density=False)
    anrate = (anpsth*1e3)/nReps  # convert spikes/msec/50 trials to spikes/sec per trial in each bin
    curve = pg.PlotCurveItem(anbins, anrate, stepMode=True, fillLevel=0, brush=(0, 0, 0, 255))
    layout.getPlot(1).addItem(curve)
    layout.getPlot(1).setLabel('left', 'spikes/sec (1 msec bins)')
    if 'SR' in stimInfo.keys():
        sr =stimInfo['SR']
    else:
        sr = 'HS'
    layout.getPlot(1).setTitle('AN (1 fiber: %.1f kHz, %ddb SPL, %sR, %d Reps)' % 
            (stimInfo['F0']/1000., stimInfo['dB'], sr, stimInfo['nReps']))
    (CNpsth, CNbins) = np.histogram(spt, bins=np.linspace(0., 250., 250), density=False)
    CNrate = (CNpsth*1e3)/nReps    
    curve = pg.PlotCurveItem(CNbins, CNrate, stepMode=True, fillLevel=0, brush=(0, 0, 255, 255))
    layout.getPlot(0).addItem(curve)

    (sl1, sl2) = getFirstSpikes(spikeTimes, stime, nReps)  # get first and secons in respond to stimulus only
    (CNpsthFS, CNbinsFS) = np.histogram(sl1, bins=np.linspace(0., 250., 250), density=False)
    CNrateFS = (CNpsthFS*1e3)/nReps
    curveFS = pg.PlotCurveItem(CNbinsFS, CNrateFS, stepMode=True, fillLevel=0,
        brush=(255, 0, 0, 255), pen=None)
    layout.getPlot(0).addItem(curveFS)

    (CNpsthSecS, CNbinsSecS) = np.histogram(sl2, bins=np.linspace(0., 250., 250), density=False)
    CNrateSecS = (CNpsthSecS*1e3)/nReps
    curveSecS = pg.PlotCurveItem(CNbinsSecS, CNrateSecS, stepMode=True, fillLevel=0,
        brush=(255, 0, 255, 255), pen=None)
    layout.getPlot(0).addItem(curveSecS)
    
    layout.getPlot(0).setLabel('left', 'spikes/sec (1 msec bins)')
    layout.getPlot(0).setTitle('VCN Cell 18')

    pgh.show()


def plotPSTH(infile):
    spikeTimes, inputSpikeTimes, stimInfo, d = readFile(infile)
    nReps = stimInfo['nReps']
    starttime = 1000.*stimInfo['pip_start'][0]
    stimdur = stimInfo['pip_duration']
    print 'AN: '
    ANfspike(inputSpikeTimes, starttime, nReps)
    print 'CN: '
    CNfspike(spikeTimes, starttime, nReps)
    plot1(spikeTimes, inputSpikeTimes, stimInfo)


def get_dimensions(n, pref='height'):
    nopt = np.sqrt(n)
    inoptw = int(nopt)
    inopth = int(nopt)
    while inoptw*inopth < n:
        if pref == 'width':
            inoptw += 1
            if inoptw * inopth > (n-inopth):
                inoptw -= 1
                inopth += 1
        else:
            inopth += 1
            if inoptw * inopth > (n-inoptw):
                inopth -= 1
                inoptw += 1
            
    return(inopth, inoptw)


def plotSingles(inpath, SR='HS'):
    nInputs = 6
    nrow, ncol = get_dimensions(nInputs, pref='width')
    fig, ax = plt.subplots(nInputs, 1, figsize=(4,6))
    # gs = gridspec.GridSpec(nrow, ncol)
    # ax = OrderedDict()
    # for i, g in enumerate(gs):
    #     ax[i] = plt.Subplot(fig, g)
    # win = pgh.figure(title='AN Inputs')
    # win.resize(300*ncol, 125*nrow)
    # layout = pgh.LayoutMaker(cols=ncol,rows=nrow, win=win, labelEdges=True, ticks='talbot')
    for i in range(nInputs):
        fname = synfile_template % (cell, i, SR)
        fname = os.path.join(inpath, fname)
        print fname
        spikeTimes, inputSpikeTimes, stimInfo, d = readFile(fname)
        nReps = stimInfo['nReps']
        # print 'nreps: ', nReps
        pl = ax[i]
#        PH.nice_plot(pl)
        sv = d['somaVoltage']
        tm = d['time']
        # print sv
        # print tm
        for j in range(nReps):
            # print 'j: ', j
            # print tm.shape
            # print sv[j].shape
            pl.plot(tm-100., sv[j])
           # print dir(pl)
            pl.set_ylabel('mV')
            pl.set_ylim(-65., 25.)
            pl.set_xlim(-50., 100.)
#            pl.plot(tm, sv[j], pen=pg.mkPen(pg.intColor(j, nReps), width=1.0))
#            pl.setLabel('left', 'mV')
#            pl.setYRange(-65, -5)
            # PH.do_talbotTicks(pl, ndec=0)
            #PH.do_talbotTicks(pl, ndec=0)
#    pgh.show()
    plt.draw()
    plt.show()

def readIVFile(filename):
    print 'Reading IV results file: %s', filename
    f = open(filename, 'r')
    d = pickle.load(f)
    f.close()
    tr = {}
    print d.keys()
    for i in range(len(d['Results'])):
       dr = d['Results'][i]
       k = dr.keys()[0]
       tr[k] = dr[k]
    arun = ar.AnalyzeRun(tr)  # set up the analysis for this data
    arun.IV()  # now analyze the data
    print('IV Summary: ')
    print('  Rinss = {:8.1f} Mohm   Rinpk = {:8.1f} Mohm'.format(arun.IVResult['Rinss'], arun.IVResult['Rinpk']))
#    print arun.IVResult.keys()
    order = np.argsort(arun.IVResult['I'])
    utau = u'\u03C4'
    print(u'   {:^7s} {:^8s} {:^9s}{:^9}'.format('I (nA)', u'\u03C4 (ms)', 'Vss (mV)', 'Vmin (mV)')) 
    print arun.IVResult
    for k in order[::-1]:
        tf = arun.IVResult['taus'].keys()[k]
        print('  {:6.1f}  {:7.2f} {:9.1f} {:9.1f}'.format(arun.IVResult['I'][k], 
            arun.IVResult['taus'][tf]['tau'].value,  arun.IVResult['Vss'][k], arun.IVResult['Vmin'][k]))
    return(arun.IVResult, d, tr)


    
def plotIV(infile):
    # win = pgh.figure(title='AN Inputs')
    # win.resize(600, 400)
    # layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
    # layout.gridLayout.setRowStretchFactor(0, 12)
    # layout.gridLayout.setRowStretchFactor(1, 1)

    fig, ax = plt.subplots(2, 1, figsize=(6, 4))
    ivr, d, tr = readIVFile(infile)
    mons = d['runInfo']['electrodeSection']
    # blueline = pg.mkPen('b')
    # blkline = pg.mkPen('k')
    # redline = pg.mkPen('r')
    taufit = ivr['taufit']
    for k, i in enumerate(d['runInfo']['stimInj']):
        ax[0].plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticV'], 'k-')
        ax[1].plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticI'], 'b-')
        # layout.getPlot((0,0)).plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticV'], pen=blkline)
        # layout.getPlot((1,0)).plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticI'], pen=blueline)
          #plot_run(infile, tr[i], init=False)
    for k in taufit[0].keys():
        ax[0].plot(taufit[0][k], taufit[1][k], 'r-')
 #        layout.getPlot((0,0)).plot(taufit[0][k], taufit[1][k], pen=redline)
    ax[0].set_xlim(0, 200.)
    ax[0].set_xlabel('T (ms)')
    ax[1].set_xlim(0., 200.)
    ax[1].set_xlabel('T (ms)')
    PH.nice_plot(ax[0])
    PH.nice_plot(ax[1])
    # layout.getPlot((0,0)).setXRange(0, 200.)
    # layout.getPlot((0,0)).setLabels(left='V (mV)', bottom='T (ms)')
    # layout.getPlot((1,0)).setXRange(0, 200.)
    # layout.getPlot((0,0)).setLabels(left='I (nA)', bottom='T (ms)')
    # layout.getPlot((0,0)).setTitle(infile)
    #    cplts.plotFit(1, ivr['taufit'][0], ivr['taufit'][1],  c='r')
#    cplts.plotFit(1, ivr['ihfit'][0], ivr['ihfit'][1],  c='b')
    # pgh.show()


if len(sys.argv) > 1:
    analysis = sys.argv[1]
    SR = 'HS'
    if len(sys.argv) > 2:
        SR = sys.argv[2]
else:
    exit()

if analysis == 'psth':
    infile = os.path.join(baseName, cell, 'Simulations/AN', filename)
    plotPSTH(infile)
elif analysis == 'iv':
    infile = os.path.join(baseName, cell, 'Simulations/IV', '%s' % cell + '.p')
    plotIV(infile)
elif analysis == 'singles':
    inpath = os.path.join(baseName, cell, 'Simulations/AN')
    plotSingles(inpath, SR)
plt.show()




