"""
A toy model of a soma plus a tapering dendrite, with synaptic inputs
scattered along the dendrite. This is to investigate how dendritic
electrotonus and tapering affect the shape of the EPSC at the soma
as recorded under voltage clamp
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as mpl
import neuron
import nrn
from neuron import h
from cnmodel.util import custom_init

print( "<< Creating a model cell with tapering dendrites >>")
soma = h.Section(name="Soma")  # one compartment of about 29000 um2
soma.nseg = 1
soma.L = 20
soma.diam = 20
soma.insert('leak')
dends = []
ndendsec = 25
#den_dia = np.ones(ndendsec) * 5
den_dia = np.linspace(2, 0.1, ndendsec)  # taper
print ('dend dia: ', den_dia)
den_len = 10.
den = h.Section(name='dend_0')
dends.append(den)
dends[-1].connect(soma)
for i in range(ndendsec):
    dends.append(h.Section(cell=soma))
    dends[-1].connect(dends[-2])
    dends[-1].insert('leak')
    dends[-1].diam = den_dia[i]
    dends[-1].L = den_len
    dends[-1].Ra = 350.
#h.topology()

nsyn = 10
syns = []
vs = []
nc = []
vx = []
tipi = 20.
stimtimes = np.arange(0, nsyn*tipi, tipi) + 5.
#stimtimes = np.zeros(nsyn) + 5.
print( stimtimes)
print( 'build syns')
for i in range(nsyn):
    idend = i*int(ndendsec/nsyn)
    print('idend: ', idend)
    syns.append(h.Exp2Syn(dends[idend](0.5)))
    syns[-1].e = 0
    syns[-1].tau1 = 0.2
    syns[-1].tau2 = 2
    vs.append(h.VecStim())
    print (stimtimes[i])
    vx.append(h.Vector([stimtimes[i]]))
    vs[-1].play(vx[-1])
    nc.append(h.NetCon(vs[-1], syns[-1]))
    nc[-1].weight[0] = 0.01 *(1/np.exp(-i/2.))

print ('build clamp')
vstim = h.SEClamp(0.5, soma) # set up a single-electrode clamp
vstim.dur1 = 10.0
vstim.amp1 = -60
vstim.dur2 = 10.0
vstim.amp2 = -60.0
vstim.dur3 = np.max(stimtimes)+tipi
vstim.amp3 = -60.0
vstim.rs = 0.01

print( 'build vectors')
res = {}
res['v_soma'] = h.Vector()
res['v_soma'].record(soma(0.5)._ref_v)
res['i_inj'] = h.Vector()
res['i_inj'].record(vstim._ref_i)
res['time'] = h.Vector()
res['time'].record(h._ref_t)
res['delec'] = h.Vector()
res['delec'].record(dends[-1](0.5)._ref_v)

print( 'init')
h.dt = 0.025

custom_init(v_init=-65)
h.finitialize()
h.tstop = np.max(stimtimes) + tipi
h.t = 0.
while h.t < h.tstop:
        h.fadvance()
# h.batch_save() # save nothing
# h.batch_run(h.tstop, h.dt, "dend.dat")
fig, ax = mpl.subplots(4, 1, figsize=(5, 10))

ax = ax.ravel()
ax[0].plot(res['time'], res['i_inj'])
for i in range(nsyn):
    it0 = np.argmin(np.fabs(np.array(res['time']) - (5+i*tipi)))
    it1 = np.argmin(np.fabs(np.array(res['time']) - (5+i*tipi + tipi/2.)))
    print( it0, it1)
    y = np.array(res['i_inj'])[it0:it1]
    y = y - y[0]
    ax[1].plot(np.array(res['time'])[it0:it1] - (5+i*tipi), y)
    ax[2].plot(np.array(res['time'])[it0:it1] - (5+i*tipi), y/np.min(y))
    rx = np.array(res['time'])[it0:it1] - np.array(res['time'])[it0]
    if i == 0:
        dvsum = y
    else:
        dvsum += y
ax[1].plot(rx, (dvsum - dvsum[0])/(-np.min(dvsum - dvsum[0])), 'k-', linewidth=2.0)
ax[3].plot(np.array(res['time']), np.array(res['delec']), 'g-')
mpl.show()

    
    