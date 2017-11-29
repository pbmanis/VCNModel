#!/usr/bin/python
"""
displayresult shows the results of a model_run. Just enter the filename in the fn field

runtype defines IV or AN
modeltype is needed to pull the correct model data
for example:
dispalyresult.py 14  should show the IV for 14 for the model in modeltype

AN_Result_VCN_c09_Syn006_N001_030dB_16000.0_MS.p


"""
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
import argparse
from collections import OrderedDict
#from cnmodel.util import sound
import sac_campagnola as SAC
import pycircstat as PCS
import pylibrary.PlotHelpers as PH
import GIF_fit as GFit

import warnings
warnings.filterwarnings("ignore")

class DisplayResult():
    def __init__(self):
        self.findfiles = False

        self.modeltypes = ['mGBC', 'XM13', 'RM03']
        self.runtypes = ['AN', 'an', 'IO', 'IV', 'iv', 'gifnoise']
        self.modetypes = ['find', 'singles', 'IO', 'multi']
        self.Params = OrderedDict()
        self.fn = {}
        self.bp = {}
        self.Params['patterns'] = ''
        self.Params['modeltype'] = 'XM13'
        self.Params['runtype'] = 'IV'
        self.Params['modetypes'] = 'singles'
        self.Params['threshold'] = 0.
        self.findfiles = False
        
    def file_setup(self):
        if self.Params['modetype'] in ['multi']:
            # multiple cells
            self.Params['patterns'] = ['c%s' % p for p in self.Params['cell']]
            self.findfiles = True
        else:
            self.Params['patterns'] = ['c%s' % p for p in self.Params['cell']] #i for i in [8, 9, 17, 18, 19, 20, 21, 22]]
            self.findfiles = False

        print self.Params['patterns']
        for p in self.Params['patterns']:
            if self.Params['runtype'] in ['AN'] and self.Params['modetype'] not in ['IO']:
                self.fn[p] = 'AN_Result_VCN_{0:s}_{1:s}_delays_N020_040dB_4000.0_FM50.0_DM100_MS.p'.format(p, self.Params['modeltype'])
                self.fn[p] = 'AN_Result_VCN_{0:s}_{1:s}_Syn{2:03d}_N001_030dB_16000.0_MS.p'.format(p, self.Params['modeltype'], self.Params['synno'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
            
            elif self.Params['runtype']in ['AN'] and self.Params['modetype'] in ['IO']:
                self.fn[p] = 'AN_Result_VCN_{0:s}_{1:s}_SynIO{2:03d}_N001_030dB_16000.0_MS.p'.format(p, self.Params['modeltype'], self.Params['synno'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
        
            elif self.Params['runtype'] in ['IV']:
                self.fn[p] = 'VCN_{0:s}_{1:s}.p'.format(p, self.Params['modeltype'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/IV'.format(p)
            
            elif self.Params['runtype'] in ['gifnoise']:
                self.fn[p] = 'VCN_{0:s}_{1:s}_gifnoise.p'.format(p, self.Params['modeltype'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/Noise'.format(p)
        self.lookupfiles()

    def lookupfiles(self):  # simply look for files
        if self.findfiles:
            for p in self.Params['patterns']:
                try:
                    h = open(os.path.join(self.bp[p], self.fn[p]))
                    print '       found %s' % p
                    h.close()
                except:
                    print 'did not find file for %s' % p, self.fn[p]
                    exit(1)

    def show(self):
        self.file_setup()
        if self.Params['runtype'] == 'IV':
            self.show_IV()
            #fig, ax = plt.subplots(len(self.Params['patterns'])+1, 3)
        if self.Params['runtype'] == 'AN':
            self.show_AN()
        if self.Params['runtype'] == 'gifnoise':
            self.show_gifnoise()

    def show_AN(self):    
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
        d = {}
        for p in self.Params['patterns']:
            h = open(os.path.join(self.bp[p], self.fn[p]))
            d[p] = pickle.load(h)
            h.close()

        # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
        #print 'data keys: ', d.keys()
        #print len(d['time'])

        for j, pattern in enumerate(self.Params['patterns']):
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
                i#nputs = len(d[pattern]['inputSpikeTimes'][st])
                # for j in range(inputs):
                #     t = (d['ref']['inputSpikeTimes'][st][j])/1000.
                #     y = (i+0.1+j*0.05)*np.ones(len(t))
                #     ax[0,].plot(t, y, 'o', markersize=2.5, color='r')
        #for i, st in enumerate(d['ref']['somaVoltage'].keys()):
        #    ax[2,].plot(d['ref']['time'], d['ref']['somaVoltage'][st], color='k')
        #    ax[j][0].set_xlim((0, 0.25))
        #    ax[j][1].set_xlim((0, 0.25))

    def show_SAC(self):
        sac = SAC.SAC()
        f, ax = plt.figure()
        d = {}
        for p in self.Params['patterns']:
            h = open(os.path.join(self.bp[p], self.fn[p]))
            d[p] = pickle.load(h)
            h.close()

        fmod = d['patterns'[0]]['stimInfo']['fmod']  # modulation frequency
        #sacht = 2.5
        sacmax = 100
        #phaseht = 15

        f, ax2 = plt.subplots(1)
        u = 2.0*np.pi/(1000.0/fmod)
        for j, pattern in enumerate(self.Params['patterns']):
            X = []
            pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                    'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                    'dur': 1000., 'ddur': 200.}
            pars_x =     pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                    'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                    'dur': 1000., 'ddur': 200.}
            x_spl = []
            y_spl = []
            Nspikes = 0
            for i, st in enumerate(d[pattern]['spikeTimes'].keys()):

                X.append(d[pattern]['spikeTimes'][st])
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

    def show_IV(self):
            fig, ax = plt.subplots(2,1)
            a = ax.ravel()
            d = {}
            for p in self.Params['patterns']:
                h = open(os.path.join(self.bp[p], self.fn[p]))
                d[p] = pickle.load(h)
                h.close()
            print ('pattern: ',self.Params['patterns'])
            for j, pattern in enumerate(self.Params['patterns']):
                for k in range(len(d[pattern]['Results'])):
                    k0 = d[pattern]['Results'][k]
                    rkey = k0.keys()[0]
                    v = k0[rkey]['monitor']['postsynapticV']
                    a[0].plot(k0[rkey]['monitor']['time'], v, linewidth=0.5)
                    inj = k0[rkey]['monitor']['postsynapticI']
                    a[1].plot(k0[rkey]['monitor']['time'], inj)

            plt.show()

    def show_gifnoise(self):
            fig, ax = plt.subplots(2,1)
            plt.suptitle('gifnoise compare')
            dt_beforespike = 3.
            ax = ax.ravel()
            d = {}
            for p in self.Params['patterns']:
                h = open(os.path.join(self.bp[p], self.fn[p]))
                d[p] = pickle.load(h)
                h.close()
            for j, pattern in enumerate(self.Params['patterns']):
                for k in range(len(d[pattern]['Results'])):
                    k0 = d[pattern]['Results'][k]
                    rkey = k0.keys()[0]
                    data = k0[rkey]['monitor']
                    v = data['postsynapticV']
                    inj = data['postsynapticI']
                    dt = data['time'][1]-data['time'][0]
                    GF = GFit.GIFFitter(path=self.bp[pattern], dt=dt)
                    GF.GIF.gl      = 100.        # nS, leak conductance
                    GF.GIF.C       = 20.     # nF, capacitance
                    GF.GIF.El      = -65.0            # mV, reversal potential

                    GF.GIF.Vr      = -55.0            # mV, voltage reset
                    GF.GIF.Tref    = 0.2             # ms, absolute refractory period

                    GF.GIF.Vt_star = -45.0            # mV, steady state voltage threshold VT*
                    GF.GIF.DV      = 10.             # mV, threshold sharpness
                  #GF.set_templates(aec=None, train=train_template, test=None)
                    GF.Exp.addTrainingSetTrace(v, 1e-3, inj, 1e-9, np.max(data['time']), 'Array')
                    #GF.set_files_test(AECtrace=None, trainingset=1008, testsets = range(1009, 1018), filetype='Igor')
                    #GF.set_timescales(ts=[0.1, 0.5, 2., 5.0, 10.0, 20.0])
                    print('Starting Parameters: ')
                    GF.GIF.printParameters()
                    GF.fit(threshold=self.Params['threshold'],
                            ax=ax[0], current=inj, beforeSpike=dt_beforespike)
                    print('Final Parameters: ')
                    GF.GIF.printParameters()
                    print ('v: ', v[0])
                    print('len inj, dt: ', len(inj), GF.dt, len(inj)*GF.dt)
                    (time, Vs, I_a, V_t, S) = GF.GIF.simulate(inj, v[0])  # simulate response to current trace I with starting voltage V0
                    # a[0].plot(time, Vs, 'r-', linewidth=0.75)
    #                 a[0].plot(data['time'], v, 'b-', linewidth=0.5)
                    ax[1].plot(data['time'], inj, 'k-', linewidth=0.5)
                    #GF.GIF.plotParameters()
            plt.show()
            exit(0)


if __name__ == '__main__':
    
    DR = DisplayResult()
    parser = argparse.ArgumentParser(description='Display model results')
    parser.add_argument(dest='cell', action='store', nargs='+', type=str,
                   default=None,
                   help='Select the cell(s) (no default)')
    parser.add_argument('-p', '--protocol', dest='runtype', action='store',
                    default='IV', help=('Select the run type (default: IV) from: %s' % DR.runtypes))
    parser.add_argument('-m', '--modetype', dest='modetype', action='store',
                    default='multi', help=('Select the mode type  (default: multi) from: %s' % DR.modetypes))
    parser.add_argument('-M', '--modeltype', dest='modeltype', action='store',
                    default='XM13', help=('Select the model type (default XM13) from: %s '% DR.modeltypes))
    parser.add_argument('-n', '--synno', dest='synno', action='store',
                    default=None, help='Select the synno')
    parser.add_argument('-t', '--threshold', dest='threshold', type=float, action='store',
                    default=0., help=('Spike detection threshold, mV'))
                   
    args = vars(parser.parse_args())
#    model.print_modelsetup()

    for k in args.keys():
        DR.Params[k] = args[k]
    
    DR.show()
        