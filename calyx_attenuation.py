#!/usr/bin/env python

"""
Compute attenuation along the processes and make a plot

"""
import neuron as h
from neuron import *
import pylibrary.PlotHelpers as PH
import os
import numpy as np
import matplotlib.pyplot as mpl
import cellInitialization as CI
# import here so we can parse commands more quickly 
# (and without neuron garbage)
from neuronvis.hoc_reader import HocReader
from neuronvis.hoc_viewer import HocViewer
import neuronvis.hoc_graphics

h.load_file('stdrun.hoc')  # don't forget!

hoc_file = os.path.join('MorphologyFiles', 'P30_calyx.hoc')
hoc = HocReader(hoc_file)
electrodesite = hoc.sections['sections[26]']
relectrodesite = hoc.sections['sections[292]']
istim = h.iStim(0.5, sec=electrodesite)
CI.initModel(hoc, electrodeSite=electrodesite)
clampV = -65
tdelay = 1
tstep = [0.5, 10]
#view = HocViewer(hoc)
print 'hoc file: ', hoc_file
#print 'hoc: ', hoc.sections
#print dir(hoc.sections['sections[0]'])
for s in hoc.sections.keys():
    hoc.sections[s].Ra = 200.

#for s in hoc.sections.keys():
#    print s, hoc.sections[s].allseg(0.5).v
vec={}
for var in ['time', 'V', 'V2', 'IChan', 'Vcmd']:
    vec[var] = h.Vector()

secvec = {}
for s in hoc.sections.keys():
    secvec[s] = h.Vector() # one for each section!
    
vcPost = h.SEClamp(0.5, sec=electrodesite)
vcPost.dur1 = tdelay
vcPost.amp1 = clampV
vcPost.dur2 = tstep[0]
vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
vcPost.dur3 = tstep[1]
vcPost.amp3 = clampV
vcPost.rs = 1e-6

vcPost.amp2 = -60

istim.delay = 0.5
istim.dur = 200  # these actually do not matter...
istim.iMax = 1.

hoc.h.tstop = 10
hoc.h.dt = 0.01  # force small time step. cvode is probably off.
hoc.h.celsius = 37.0
hoc.h.finitialize(-50)
vec['IChan'].record(vcPost._ref_i, sec=electrodesite)
vec['V'].record(electrodesite()._ref_v, sec=electrodesite)
vec['V2'].record(relectrodesite()._ref_v, sec=relectrodesite)
for s in hoc.sections.keys():
    secvec[s].record(hoc.sections[s]()._ref_v, sec=hoc.sections[s])
vec['time'].record(h._ref_t)
hoc.h.run(25)
# for s in hoc.sections.keys():
#     print s, hoc.sections[s].allseg(0.5).v
tx = np.array(vec['time'])
vx = np.array(vec['V'])
v2 = np.array(vec['V2'])
alls = np.zeros((len(hoc.sections.keys()), len(v2)))

secdist = np.zeros((len(hoc.sections.keys()), 2))

electrodesite.push()  # make this the currently accessed section
hoc.h.distance()
for i, s in enumerate(hoc.sections.keys()):
    print i, s 
    d = hoc.h.distance(hoc.sections[s](0.5).x, sec=hoc.sections[s])
    secdist[i,0] = d
    alls[i,:] = np.array(secvec[s])
    secdist[i,1] = np.max(alls[i,:])
# for i in range(len(vec['time'])):
#     print 't: %f v: %f' % (tx[i], vx[i])
f, ax = mpl.subplots(2)
#ax[0].plot(tx, vx)
#ax[0].plot(tx, v2)

for i in range(alls.shape[0]):
    ax[0].plot(tx, alls[i])
ax[0].set_ylim([-66, -58.])
mpl.xlabel('Time (ms)')
mpl.ylabel('V (mV)')
PH.nice_plot(ax[0])   

ax[1].scatter(secdist[:,0], secdist[:,1], s=15, c='r', alpha=0.5)    
ax[1].set_ylim([-66, -58.])
ax[1].set_xlim([-5., 105.])
PH.refline(ax[1], -65)
mpl.xlabel('Distance (um)')
mpl.ylabel('max V (mV)')
PH.nice_plot(ax[1])

mpl.show()

