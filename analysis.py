


import sys
import os.path
import pickle
#import neuronvis.sim_result as sr

import pylibrary.Utility as pu  # access to spike finder routine

import time

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import pylibrary.pyqtgraphPlotHelpers as pgh

import numpy as np

filename = 'AN_Result_VCN_c18_reparented755N050_040dB_4000.0_MS.p'


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


starttime = 1000.*stimInfo['pip_start'][0]
sstartdur = stimInfo['pip_dur']
nReps = stimInfo['nReps']



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


def CNfspike(spikes, stime, nReps):
    # print spikes
    # print
    
    sl1 = np.full(nReps, np.nan)
    sl2 = np.full(nReps, np.nan)
    for r in range(nReps):
        rs = np.where(spikes[r] > stime)[0]  # get spike times post stimulus onset
        if len(spikes[r]) > 0:
            sl1[r] = spikes[r][rs[0]]
        if len(spikes[r]) > 1:
            sl2[r] = spikes[r][rs[1]]  # second spikes
    print 'Cochlear Nucleus Bushy Cell: '
    print '   mean First spike latency:  %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl1), np.nanstd(sl1),
        np.count_nonzero(~np.isnan(sl1)))
    print '   mean Second spike latency: %8.3f ms stdev: %8.3f (N=%3d)' % (np.nanmean(sl2), np.nanstd(sl2), 
        np.count_nonzero(~np.isnan(sl2)))

print 'AN: '
ANfspike(inputSpikeTimes, starttime, nReps)
print 'CN: '
CNfspike(spikeTimes, starttime, nReps)

# plot PSTH's


win = pgh.figure(title='AN Inputs')
layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
# flatten spike times
spt = []
for r in range(nReps):
    spt.extend(spikeTimes[r])

oneANF = inputSpikeTimes[0]
anf = []
for r in range(len(oneANF)):
    anf.extend(oneANF[r])
#print stimInfo

(anpsth, anbins) = np.histogram(anf, bins=np.linspace(0., 250., 250), density=False)
anrate = (anpsth*1e3)/nReps  # convert spikes/msec/50 trials to spikes/sec per trial in each bin
curve = pg.PlotCurveItem(anbins, anrate, stepMode=True, fillLevel=0, brush=(0, 0, 0, 255))
layout.getPlot(1).addItem(curve)
layout.getPlot(1).setLabel('left', 'spikes/sec (1 msec bins)')
layout.getPlot(1).setTitle('AN (1 fiber: %.1f kHz, %ddb SPL, %sR, %d Reps)' % (stimInfo['F0']/1000., stimInfo['dB'], 'MS', stimInfo['nReps']))
(CNpsth, CNbins) = np.histogram(spt, bins=np.linspace(0., 250., 250), density=False)
CNrate = (CNpsth*1e3)/nReps
curve = pg.PlotCurveItem(CNbins, CNrate, stepMode=True, fillLevel=0, brush=(0, 0, 255, 255))
layout.getPlot(0).addItem(curve)
layout.getPlot(0).setLabel('left', 'spikes/sec (1 msec bins)')
layout.getPlot(0).setTitle('VCN Cell 18')

pgh.show()





