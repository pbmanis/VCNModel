"""
displayresult shows the results of a model_run. Just enter the filename in the fn field

runtype defines IV or AN
modeltype is needed to pull the correct model data
for example:
dispalyresult.py 14  should show the IV for 14 for the model in modeltype

"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
from cnmodel.util import sound
import sac_campagnola as SAC
import pycircstat as PCS
import pylibrary.PlotHelpers as PH

findfiles = False
modeltype = 'XM13'
runtype = 'IV'
if sys.argv[1] not in ['find']:
    patterns = sys.argv[1:]
    patterns = ['c%02d' % int(p) for p in patterns]
else:
    print 'searching'
    patterns = ['c%02d' % i for i in [8, 9, 17, 18, 19, 20, 21, 22]]
    findfiles = True
fn = {}
bp = {}
for p in patterns:
    if runtype is not 'IV':
        fn[p] = 'AN_Result_VCN_{0:s}_delays_N020_040dB_4000.0_FM50.0_DM100_MS.p'.format(p)
        bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
    else:
        fn[p] = 'VCN_{0:s}_{1:s}.p'.format(p, modeltype)
        bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/IV'.format(p)

d = {}

# simply look for files
if findfiles:
    for p in patterns:
        try:
            h = open(os.path.join(bp[p], fn[p]))
            print '       found %s' % p
            h.close()
        except:
            print 'did not find file for %s' % p
    exit(1)

if runtype == 'IV':
    fig, ax = plt.subplots(2,1)
    a = ax.ravel()
    d = {}
    for p in patterns:
        h = open(os.path.join(bp[p], fn[p]))
        d[p] = pickle.load(h)
        h.close()
    dt = 0.1
    print ('pattern: ', patterns)
    for j, pattern in enumerate(patterns):
        for k in range(len(d[pattern]['Results'])):
            k0 = d[pattern]['Results'][k]
            rkey = k0.keys()[0]
            v = k0[rkey]['monitor']['postsynapticV']
            a[0].plot(k0[rkey]['monitor']['time'], v)
            inj = k0[rkey]['monitor']['postsynapticI']
            a[1].plot(k0[rkey]['monitor']['time'], inj)
            


    plt.show()
    exit(0)
    #fig, ax = plt.subplots(len(patterns)+1, 3)

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
    PH.noaxes(ax1, 'x')

#fig.tight_layout()
# plt.show()
# exit()

for p in patterns:
    h = open(os.path.join(bp[p], fn[p]))
    d[p] = pickle.load(h)
    h.close()

# d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
#print 'data keys: ', d.keys()
#print len(d['time'])
fmod = d[patterns[0]]['stimInfo']['fmod']  # modulation frequency
sacht = 2.5
sacmax = 100
phaseht = 15

f, ax2 = plt.subplots(1)

dt = 0.1
for j, pattern in enumerate(patterns):
    for i in d[pattern]['somaVoltage'].keys():
        v = d[pattern]['somaVoltage'][i]
        ax2.plot(d[pattern]['time'], v)
#    plt.show()
#    print np.array(d[patterns[j]]['stimWaveform'][0])
#    w = d[patterns[j]]['stimWaveform'][0]
#    t = np.linspace(0., len(d[patterns[j]]['stimWaveform'][0])*dt, len(d[patterns[j]]['stimWaveform'][0])) 
#    ax[j][1].plot(t, w)  # stimulus underneath

    for i, st in enumerate(d[pattern]['spikeTimes'].keys()):
        ax[j][0].plot(d[pattern]['spikeTimes'][st]/1000., i*np.ones(len( d[pattern]['spikeTimes'][st])), 'o', markersize=2.5, color='b')
        inputs = len(d[pattern]['inputSpikeTimes'][st])
        # for j in range(inputs):
        #     t = (d['ref']['inputSpikeTimes'][st][j])/1000.
        #     y = (i+0.1+j*0.05)*np.ones(len(t))
        #     ax[0,].plot(t, y, 'o', markersize=2.5, color='r')
#for i, st in enumerate(d['ref']['somaVoltage'].keys()):
#    ax[2,].plot(d['ref']['time'], d['ref']['somaVoltage'][st], color='k')
#    ax[j][0].set_xlim((0, 0.25))
#    ax[j][1].set_xlim((0, 0.25))


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
    print 'mean vs: for %s at fMod = %.2f:  %f angle: %f' % (pattern, fmod, vs, th.mean())
    print 'rayleigh: p=%f  z=%f (p > 0.05 means data is uniformly distributed)' % (p, z)
    ax[j][2].bar(bins[:-1], yh)
    ax[j][2].set_xlim((-sacmax, sacmax))
#    ax[j][2].set_ylim((0, sacht))
    phasebinnum = 101
    phasebins = np.linspace(-0.5, 0.5, num=phasebinnum)
    thhist, thbins = np.histogram(th/np.pi, bins=phasebins, density=False)
    ax[j][3].bar(thbins[:-1]+0.5, thhist, width=1.0/phasebinnum)
#    ax[j][3].set_ylim((0, phaseht))

    if j == 0:
        continue
    ycc, ccbins = sac.XAC(X0, X, pars_x, trialbased=True)
    ax[j][4].bar(ccbins[:-1], ycc)
    ax[j][4].set_xlim((-sacmax, sacmax))
    
    # now cross-correlate data...
    

#PH.nice_plot(ax.flatten().tolist())
plt.show()