"""
Check NMDA_kampa mechanism (voltage clamp)

Make a single site version of the multisite synapse, set AMPA to 0. 
Go voltage clamp, vary V and watch currents

"""
from __future__ import print_function


from collections import OrderedDict
import numpy as np
from neuron import h
import cnmodel.util as util
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.synapses import GluPSD
import matplotlib.pyplot as mpl
from lmfit import Model
from lmfit.models import ExponentialModel

# define some model equations and lmfit models
def modwoodhull(x, amp, ko, delta):
    """
    Modified woodhull eq from Pliss, Yang and Xu-Friedman, 2009
    """
    mg = 1.0
    F = h.FARADAY
    R = h.R
    z = 2.0
    T = 273.16 + 34.
    vRev = 0.
    return amp * (x-vRev)/(1.0 + (mg/ko)*np.exp(-delta*z*F*x/(R*T)))

wmodel = Model(modwoodhull)
wmodel.set_param_hint('delta', min=0.0, max=1.0)

def boltz(x, amp, vh, vr):
    """
    Simple Boltzmann function
    """
    vRev = 0.
    return amp * (x-vRev)/(1.0 + np.exp((vh-x)/vr))

bmodel = Model(boltz)


tstop = 100.
pre_cell = cells.SGC.create()
post_cell = cells.Bushy.create(debug=True, ttx=True)
n_synapses = 1
dt = 0.025

synapses = []
for i in range(n_synapses):
    synapses.append(pre_cell.connect(post_cell))


clampV = -60.  # mV
vccontrol = h.VClamp(0.5, sec=post_cell.soma)
vccontrol.dur[0] = 10.0
vccontrol.amp[0] = clampV
vccontrol.dur[1] = 100.0
vccontrol.amp[1] = clampV
vccontrol.dur[2] = 20.0
vccontrol.amp[2] = clampV

#
# set up stimulation of the presynaptic axon/terminal
#

istim = h.iStim(0.5, sec=pre_cell.soma)
stim = {
    'NP': 1,
    'Sfreq': 100.0,
    'delay': 10.0,
    'dur': 0.5,
    'amp': 10.0,
    'PT': 0.0,
    'dt': dt,
}
(secmd, maxt, tstims) = util.make_pulse(stim)

if tstop is None:
    tstop = len(secmd) * dt

istim.delay = 0
istim.dur = 1e9 # these actually do not matter...
istim.iMax = 0.0

# istim current pulse train
i_stim_vec = h.Vector(secmd)
i_stim_vec.play(istim._ref_i, dt, 0)

# create hoc vectors for each parameter we wish to monitor and display
synapse = synapses[0]

all_psd = []
for syn in synapses:
    # collect all PSDs across all synapses
    all_psd.extend(syn.psd.all_psd)
    for p in syn.psd.ampa_psd:
        p.gmax = 0.
    post_cell.make_psd(syn.terminal, 'multisite')

# print 'terminal: ', dir(synapses[0].terminal)
# print 'synapses: ', dir(synapses[0])
print('synapse psd: ', synapses[0].psd.ampa_psd[0].gmax)
print('v shift: ', synapses[0].psd.nmda_psd[0].vshift)

#print 'all psd: ', all_psd    
#
#
# Run simulation
#

temp = 34.0
h.tstop = tstop # duration of a run
h.celsius = temp
h.dt = dt
iterations = 1
results = {}
vcmds = range(-100, 60, 10)
f, ax = mpl.subplots(3, 1)
ax = np.ravel(ax)
print(dir(post_cell.soma(0.5)))
ivdat = np.zeros((2, len(vcmds)))
for i, vc in enumerate(vcmds):
    vccontrol.amp[1] = vc
    results['v_pre'] = h.Vector()
    results['v_pre'].record(pre_cell.soma(0.5)._ref_v)
    results['t'] = h.Vector()
    results['t'].record(h._ref_t)
    results['i_soma'] = h.Vector()
    results['i_soma'].record(all_psd[0]._ref_i)
    results['cmd'] = h.Vector()
    results['cmd'].record(post_cell.soma(0.5)._ref_v)
    util.custom_init()
    while h.t < tstop:
        h.fadvance()
        #h.run()
    if (vc - 40.) < 1.0:
        i40 = i
        itrace = np.array(results['i_soma'])
    ax[0].plot(results['t'], results['i_soma'])
    ax[1].plot(results['t'], results['cmd'])
    tpos = int(14./dt)
    ivdat[0,i] = results['cmd'][tpos]
    ivdat[1,i] = results['i_soma'][tpos]
    
ivdat[1,:] = ivdat[1,:]/ivdat[1,i40]  # normalize to 40 mV data

ax[2].plot(ivdat[0,:], ivdat[1,:], marker='s', linestyle='None')
    

# do some fitting

# exp decay
ipk = np.argmax(itrace)
t = np.arange(0, (len(itrace)-ipk)*dt, dt)
gmodel = ExponentialModel(independent_vars=['x']) #, prefix='', missing=None, name=None, **kwargs)
result = gmodel.fit(itrace[ipk:], x=t, amplitude=0.02, decay=20.)
print(result.fit_report())
ax[0].plot(t+ipk*dt, result.best_fit, 'r-')

# IV curve...
wresult = wmodel.fit(ivdat[1,:], x=ivdat[0,:]/1000., amp=1, ko=1.9, delta=0.5)
print(wresult.fit_report())
ax[2].plot(ivdat[0,:], wresult.best_fit, 'r--')

bresult = bmodel.fit(ivdat[1,:], x=ivdat[0,:]/1000., amp=1, vh=0., vr=0.02)
print(bresult.fit_report())
ax[2].plot(ivdat[0,:], bresult.best_fit, 'b--')

mpl.show()
