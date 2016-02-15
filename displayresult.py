"""
displayresult shows the results of a model_run. Just enter the filename in the fn field
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from cnmodel.util import sound
import sac_campagnola as SAC
import os.path
import pycircstat as PCS

patterns = ['c18', 'c08', 'c09']
#patterns = ['c18']
fn = {}
fn[patterns[0]] = 'AN_Result_VCN_c18_delays_N020_040dB_4000.0_MS.p'
fn[patterns[1]] = 'AN_Result_VCN_c18_VCN_c08_delays_N020_040dB_4000.0_MS.p'
fn[patterns[2]] = 'AN_Result_VCN_c18_VCN_c09_delays_N020_040dB_4000.0_MS.p'

basepath = 'VCN_Cells/VCN_c18/Simulations/AN'

d = {}

for p in patterns:
    h = open(os.path.join(basepath, fn[p]))
    d[p] = pickle.load(h)
    h.close()

# d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
#print 'data keys: ', d.keys()
#print len(d['time'])
fmod = d[patterns[0]]['stimInfo']['fmod']  # modulation frequency

w = d[patterns[0]]['stimWaveform'][0][0].generate()
t = d[patterns[0]]['stimWaveform'][0][0].time
fig, ax = plt.subplots(4, 3)
ax[0,0].plot(t, w)  # stimulus at top
#print d['spikeTimes']
for j, pattern in enumerate(patterns):
    for i, st in enumerate(d[pattern]['spikeTimes'].keys()):
        ax[1+j,0].plot(d[pattern]['spikeTimes'][st]/1000., i*np.ones(len( d[pattern]['spikeTimes'][st])), 'o', markersize=2.5, color='b')
        inputs = len(d[pattern]['inputSpikeTimes'][st])
        # for j in range(inputs):
        #     t = (d['ref']['inputSpikeTimes'][st][j])/1000.
        #     y = (i+0.1+j*0.05)*np.ones(len(t))
        #     ax[0,].plot(t, y, 'o', markersize=2.5, color='r')
#for i, st in enumerate(d['ref']['somaVoltage'].keys()):
#    ax[2,].plot(d['ref']['time'], d['ref']['somaVoltage'][st], color='k')


sac = SAC.SAC()
xf = []
dt = 0.1
u = 2.0*np.pi/(1000.0/fmod)
for j, pattern in enumerate(patterns):
    X = []
    R = {}
    pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
            'baseper': 1.0, 'stimdur': 1000., 'delay': 0.,
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
    yh, bins = sac.SAC_asm(X, pars)
    print 'mean vs: for %s at fMod = %.2f:  %f angle: %f' % (pattern, fmod, vs, th.mean())
    print 'rayleigh: p=%f  z=%f (p > 0.05 means data is uniformly distributed)' % (p, z)

    ax[1+j, 1].bar(bins[:-1], yh)
plt.show()