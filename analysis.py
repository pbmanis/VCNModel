


import sys
import os.path
import pickle
#import neuronvis.sim_result as sr
import argparse
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
import seaborn
from matplotlib import rc
rc('text', usetex=False)
import cell_config as CFG

baseName = 'VCN_Cells'

filename_template = 'AN_Result_{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'

synfile_template = 'AN_Result_{0:s}_Syn{1:03d}_N{2:03d}_040dB_4000.0_{3:2s}.p'
excsynfile_template = 'AN_Result_{0:s}_ExcludeSyn{1:03d}_N{2:03d}_040dB_4000.0_{3:2s}.p'
synIOfile_template = 'AN_Result_{0:s}_SynIO{1:03d}_N{2:03d}_040dB_4000.0_{3:2s}.p'

all_cells = ['09', '09h', '09nd', '17',   '18',    '19', '20', '21', '22']

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

def readFile(filename, cmd):
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
#    print 'cmd: ', cmd
    if cmd['respike']:
        spiketimes = {}
        for k in d['somaVoltage'].keys():
            spikeTimes[k] = pu.findspikes(d['time'], d['somaVoltage'][k], cmd['threshold'])
            spikeTimes[k] = clean_spiketimes(spikeTimes[k], 0.7)
    return(spikeTimes, inputSpikeTimes, stimInfo, d)

#def findspikes(x, v, thresh, t0=None, t1= None, dt=1.0, mode='schmitt', interpolate=False, debug=False):
    """ findspikes identifies the times of action potential in the trace v, with the
    times in t. An action potential is simply timed at the first point that exceeds
    the threshold... or is the peak. 
    4/1/11 - added peak mode
    if mode is none or schmitt, we work as in the past.
    if mode is peak, we return the time of the peak of the AP instead
    7/15/11 - added interpolation flag
    if True, the returned time is interpolated, based on a spline fit
    if False, the returned time is just taken as the data time. 
    """

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
        if len(spikes[r]) > 1 and len(rs) > 1:
            sl2[r] = spikes[r][rs[1]]  # second spikes
    sl1 = sl1[~np.isnan(sl1)]  # in case of no spikes at all
    sl2 = sl2[~np.isnan(sl2)]  # in case of no second spikes
    return(sl1, sl2)

def CNfspike(spikes, stime, nReps):
    # print spikes
    #print 'stime: ', stime
    sl1, sl2 = getFirstSpikes(spikes, stime, nReps)

    print 'Cochlear Nucleus Bushy Cell: '
    print '   mean First spike latency:  %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl1), np.nanstd(sl1),
        np.count_nonzero(~np.isnan(sl1)))
    print '   mean Second spike latency: %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl2), np.nanstd(sl2), 
        np.count_nonzero(~np.isnan(sl2)))
    return(sl1, sl2)


def plot1(spikeTimes, inputSpikeTimes, stimInfo, cmd):
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
    (CNpsth, CNbins) = np.histogram(spt, bins=np.linspace(0., 250., 251), density=False)
    CNrate = (CNpsth*1e3)/nReps    
    curve = pg.PlotCurveItem(CNbins, CNrate, stepMode=True, fillLevel=0, brush=(0, 0, 255, 255))
    layout.getPlot(0).addItem(curve)

    (sl1, sl2) = getFirstSpikes(spikeTimes, stime, nReps)  # get first and secons in respond to stimulus only
    (CNpsthFS, CNbinsFS) = np.histogram(sl1, bins=np.linspace(0., 250., 251), density=False)
    CNrateFS = (CNpsthFS*1e3)/nReps
    curveFS = pg.PlotCurveItem(CNbinsFS, CNrateFS, stepMode=True, fillLevel=0,
        brush=(255, 0, 0, 255), pen=None)
    layout.getPlot(0).addItem(curveFS)

    (CNpsthSecS, CNbinsSecS) = np.histogram(sl2, bins=np.linspace(0., 250., 251), density=False)
    CNrateSecS = (CNpsthSecS*1e3)/nReps
    curveSecS = pg.PlotCurveItem(CNbinsSecS, CNrateSecS, stepMode=True, fillLevel=0,
        brush=(255, 0, 255, 255), pen=None)
    layout.getPlot(0).addItem(curveSecS)
    
    layout.getPlot(0).setLabel('left', 'spikes/sec (1 msec bins)')
    layout.getPlot(0).setTitle('%s' % cmd['cell'])

    pgh.show()


def plotPSTH(infile, cmd):
    spikeTimes, inputSpikeTimes, stimInfo, d = readFile(infile, cmd)
    nReps = stimInfo['nReps']
    starttime = 1000.*stimInfo['pip_start'][0]
    stimdur = stimInfo['pip_duration']
    #print 'AN: '
    ANfspike(inputSpikeTimes, starttime, nReps)
    #print 'CN: '
    CNfspike(spikeTimes, starttime, nReps)
    plot1(spikeTimes, inputSpikeTimes, stimInfo, cmd)

def plotVm(infile, cmd):
    spikeTimes, inputSpikeTimes, stimInfo, d = readFile(infile, cmd)
    nReps = stimInfo['nReps']
    starttime = 1000.*stimInfo['pip_start'][0]
    stimdur = stimInfo['pip_duration']
    #print 'AN: '
    ANfspike(inputSpikeTimes, starttime, nReps)
    fig, ax = plt.subplots(2, sharex=True)
    spiketimes = {}
    for k in d['somaVoltage'].keys():
        spikeTimes[k] = pu.findspikes(d['time'], d['somaVoltage'][k], cmd['threshold'])
        spikeTimes[k] = clean_spiketimes(spikeTimes[k], 0.7)
        ax[0].plot(d['time'], d['somaVoltage'][k])
    ANFs = [0]
    anf = []
    for r in range(nReps):
        for n in range(len(ANFs)):
            anf.extend(inputSpikeTimes[r][ANFs[n]])  # pick the ANF
    #print stimInfo
    srGroups = ['', 'LS', 'MS', 'HS']
    (anpsth, anbins) = np.histogram(anf, bins=np.linspace(0., 250., 250), density=False)
    anrate = (anpsth*1e3)/nReps  # convert spikes/msec/50 trials to spikes/sec per trial in each bin
    ax[1].step(anbins[:-1], anpsth, where='pre')
    
    
    
    
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

def plotIO(cmd):  # plots ALL IO's in one place
    l1 = 0.11
    l2 = 0.41
    l3 = 0.71
    wid = 0.25
    ht = 0.13
    yp = [0.81, 0.62, 0.43, 0.24, 0.05]
    sizer = OrderedDict([('VCN_c09', [l1, wid, yp[0], ht]), ('VCN_c17', [l2, wid, yp[0], ht]),
                         ('VCN_c09h', [l1, wid, yp[1], ht]), ('VCN_c18', [l2, wid, yp[1], ht]),
                         ('VCN_c09nd', [l1, wid, yp[2], ht]), ('VCN_c19', [l2, wid, yp[2], ht]), 
                         #('VCN_c08', [l1, wid, yp[3], ht]), 
                         ('VCN_c20', [l3, wid, yp[0], ht]), 
                         ('VCN_c21', [l3, wid, yp[1], ht]), ('VCN_c22', [l3, wid, yp[2], ht]),
    ])  # dict elements are [left, width, bottom, height] for the axes in the plot.
    gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
    axmap = OrderedDict(zip(sizer.keys(), gr))
    P = PH.Plotter(rcshape=sizer, label=False, figsize=(5, 8), labeloffset=[0.6, 0.])
    seaborn.set_style('ticks')
    for ic, cn in enumerate(all_cells):
        cell = 'VCN_c{0:s}'.format(cn)
        inpath = os.path.join(baseName, cell, 'Simulations/AN')
        sites, celltype = CFG.makeDict(cell)
    #    print sites
        nInputs = 0
        for s in sites:
            nInputs += len(s['postlocations'].keys())
        print ('cell: ', cell, ' nInputs: ', nInputs)
        vmax = -50.

        template = synIOfile_template
        tmin = 0.
        trange = [0., 50.]
        for i in range(nInputs):
            fname = template.format(cell, i, cmds['nReps'], cmds['SR'])
            fname = os.path.join(inpath, fname)
            print fname
            try:
                spikeTimes, inputSpikeTimes, stimInfo, d = readFile(fname, cmd)
            except:
                print('Missing: ', fname)
                continue
            sv = d['somaVoltage']
            tm = d['time']
            plt.setp(P.axdict[cell].get_yticklabels(), fontsize=6)
            iof = np.zeros(len(sv.keys()))
            iofdv = np.zeros(len(sv.keys()))
            for k in sv.keys():
                iof[k] = np.max(sv[k])
                #iofdv[k] = np.max(np.diff(sv[k])/np.mean(np.diff(tm)))
            P.axdict[cell].plot(np.array(stimInfo['gSyns'])/1000., iof, 'o-', markersize=3, linewidth=0.6)
        # P.axdict[cell].set_title(cell, y=0.9, x=0.02,
        #         horizontalalignment='left', fontsize=6)
        P.axdict[cell].tick_params(direction='in', length=5., width=1., labelsize=6)
        P.axdict[cell].set_xlabel('g', fontsize=9)
        P.axdict[cell].set_ylabel('V (mV)', fontsize=9)
    seaborn.despine(P.figure_handle)
    plt.savefig('ANIO_all.pdf')
            
def plotSingles(inpath, cmd):
    seaborn.set_style('ticks')
    sites, celltype = CFG.makeDict(cmd['cell'])
#    print sites
    nInputs = 0
    for s in sites:
        nInputs += len(s['postlocations'].keys())
    #nInputs = len(sites)
    print ('cell: ', cmd['cell'], ' nInputs: ', nInputs)
    nrow, ncol = get_dimensions(nInputs, pref='width')
    fig, ax = plt.subplots(nInputs, 1, figsize=(4.75,6))
    fig.suptitle('{0:s}  SR: {1:s}'.format(cmd['cell'], cmd['SR']))
    vmax = -50.
    template = synfile_template
    tmin = -100.
    trange = [-50., 100.]

    if cmd['analysis'] == 'omit':
        template = excsynfile_template
    # if cmd['analysis'] == 'io':
    #     template = synIOfile_template
    #     tmin = 0.
    #     trange = [0., 50.]
    #     fig2, ax2 = plt.subplots(2, 2, figsize=(4.75,6))
    #     fig2.suptitle('{0:s}  SR: {1:s}'.format(cmd['cell'], cmd['SR']))
    #     axr = ax2.ravel()
        
    for i in range(nInputs):
        fname = template.format(cmds['cell'], i, cmds['nReps'], cmds['SR'])
        fname = os.path.join(inpath, fname)
        print fname
        try:
            spikeTimes, inputSpikeTimes, stimInfo, d = readFile(fname, cmd)
        except:
            print('Missing: ', fname)
            continue
        #print ('stiminfo: ', stimInfo)
        if cmd['respike']:
            spiketimes = {}
            for k in d['somaVoltage'].keys():
                dt = d['time'][1]-d['time'][0]
                spikeTimes[k] = pu.findspikes(d['time'], d['somaVoltage'][k], thresh=cmd['threshold'])
                spikeTimes[k] = clean_spiketimes(spikeTimes[k])
        nReps = stimInfo['nReps']
        pl = ax[i]
        pl.set_title('Input {0:d}: N sites: {1:d}'.format(i+1, sites[i]['nSyn']), fontsize=8, x=0.05, y=0.92,
            horizontalalignment = 'left', verticalalignment='top')
        sv = d['somaVoltage']
        tm = d['time']
        for j in range(nReps):
            m = np.max(sv[j])
            if m > vmax:
                vmax = m
        for j in range(nReps):
            pv = pl.plot(tm-tmin, sv[j], linewidth=0.8)
            pl.vlines(spikeTimes[j]-tmin, -j*2+vmax, -j*2-2+vmax,
                color=pv[0].get_c(), linewidth=1.5)  # spike marker same color as trace
            pl.set_ylabel('mV')
            pl.set_ylim(-65., vmax)
            pl.set_xlim(trange[0], trange[1])
        if pl != ax[-1]:
            pl.set_xticklabels([])
        else:
            plt.setp(pl.get_xticklabels(), fontsize=9)
            pl.set_xlabel('ms')
        # if cmd['analysis'] == 'io':
        #     plt.setp(pl.get_yticklabels(), fontsize=9)
        #     iof = np.zeros(len(sv.keys()))
        #     iofdv = np.zeros(len(sv.keys()))
        #     for k in sv.keys():
        #         iof[k] = np.max(sv[k])
        #         iofdv[k] = np.max(np.diff(sv[k])/np.mean(np.diff(tm)))
        #     axr[0].plot(stimInfo['gSyns'], iof, 'o-', markersize=4)
        #     axr[1].plot(stimInfo['gSyns'], iofdv, 's-', markersize=4)
           # axr.set_title('Input {0:d}: N sites: {1:d}'.format(i+1, sites[i]['nSyn']), fontsize=8, x=0.05, y=0.92,
        #        horizontalalignment = 'left', verticalalignment='top')
        
    for i in range(nInputs): # rescale all plots the same
        ax[i].set_ylim(-65., vmax+cmd['nReps']*3)
    ax[0].set_ylim(-65., vmax)
    
    seaborn.despine(fig)
#    seaborn.despine(fig2)


def readIVFile(filename):
    print 'Reading IV results file: %s', filename
    f = open(filename, 'r')
    d = pickle.load(f)
    f.close()
    tr = {}
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
#    print arun.IVResult
    for k in order[::-1]:
        tf = arun.IVResult['taus'].keys()[k]
        print('  {:6.1f}  {:7.2f} {:9.1f} {:9.1f}'.format(arun.IVResult['I'][k], 
            arun.IVResult['taus'][tf]['tau'].value,  arun.IVResult['Vss'][k], arun.IVResult['Vmin'][k]))
    return(arun.IVResult, d, tr)


    
def plotIV(infile, plottau=False):
    # you don't have to use an ordered dict for this, I just prefer it when debugging
    sizer = OrderedDict([('A', [0.2, 0.6, 0.3, 0.5]), ('B', [0.2, 0.6, 0.1, 0.15])
            # ('C1', [0.72, 0.25, 0.65, 0.3]), ('C2', [0.72, 0.25, 0.5, 0.1]),
            # ('D', [0.08, 0.25, 0.1, 0.3]), ('E', [0.40, 0.25, 0.1, 0.3]), ('F', [0.72, 0.25, 0.1, 0.3]),
    ])  # dict elements are [left, width, bottom, height] for the axes in the plot.
    gr = [(a, a+1, 0, 1) for a in range(0, len(sizer.keys()))]   # just generate subplots - shape does not matter
    axmap = OrderedDict(zip(sizer.keys(), gr))
    P = PH.Plotter((len(sizer.keys()), 1), axmap=axmap, label=True, figsize=(6., 4.))
#    PH.show_figure_grid(P.figure_handle)
    P.resize(sizer)  # perform positioning magic
    P.figure_handle.suptitle(cell)
    ivr, d, tr = readIVFile(infile)
    mons = d['runInfo']['electrodeSection']
    taufit = ivr['taufit']
    for k, i in enumerate(d['runInfo']['stimInj']):
        P.axdict['A'].plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticV'], 'k-', linewidth=0.75)
        P.axdict['B'].plot(tr[i]['monitor']['time'], tr[i]['monitor']['postsynapticI'], 'b-', linewidth=0.75)
    if plottau:
        for k in taufit[0].keys():
            P.axdict['A'].plot(taufit[0][k], taufit[1][k], 'r-')
    P.axdict['A'].set_xlim(0, 150.)
    P.axdict['A'].set_ylabel('V (mV)')
    P.axdict['A'].set_xticklabels([])
    
    P.axdict['B'].set_xlim(0., 150.)
    P.axdict['B'].set_xlabel('T (ms)')
    P.axdict['B'].set_ylabel('I (nA)')


def parse_cmdline():
    parser = argparse.ArgumentParser(description='Analyze protocols from a reconstructed model cell')
    parser.add_argument(dest='cell', action='store',
               default=None,
               help='Select the cell (no default)')
    parser.add_argument('-a', '--analysis', dest='analysis', action='store',
               default='iv', choices=['psth', 'iv', 'singles', 'omit', 'io', 'voltage'],
               help = 'Specify an analysis')
    parser.add_argument('-r', '--reps', type=int, default=1, dest='nReps',
               help='# repetitions in file (filename)')
    parser.add_argument('-t', '--threshold', type=float, default=-20.0, action='store',
               dest='threshold',
               help='# Spike threshold for recalculating spikes (mV). Negative numbers must be quoted with a leading space: " -30."')
    parser.add_argument('-s', '--SR', dest='SR', default=None,
               choices=['LS', 'MS', 'HS'], help='Select SR group (default is None)')
    parser.add_argument('--respike', action='store_true', default=False, dest='respike',
               help='recompute spikes from Vm waveforms (default: False)')
    args = vars(parser.parse_args())
    return args


cmds = parse_cmdline()

if cmds['analysis'] == 'psth':
    infile = os.path.join(baseName, cmds['cell'], 'Simulations/AN', filename_template.format(cmds['cell'], cmds['nReps'], cmds['SR']))
    plotPSTH(infile, cmds)
if cmds['analysis'] == 'voltage':
    infile = os.path.join(baseName, cmds['cell'], 'Simulations/AN', filename_template.format(cmds['cell'], cmds['nReps'], cmds['SR']))
    plotVm(infile, cmds)

elif cmds['analysis'] == 'iv':
    infile = os.path.join(baseName, cmds['cell'], 'Simulations/IV', '%s' % cmds['cell'] + '.p')
    plotIV(infile)
elif cmds['analysis'] in ['io']:
    plotIO(cmds)
elif cmds['analysis'] in ['singles', 'omit']:
    inpath = os.path.join(baseName, cmds['cell'], 'Simulations/AN')
    plotSingles(inpath, cmds)
plt.show()




