#!/usr/bin/python
"""
displayresult shows the results of a model_run. Just enter the filename in the fn field

runtype defines IV or AN
modeltype is needed to pull the correct model data
for example:
dispalyresult.py 14  should show the IV for 14 for the model in modeltype

AN_Result_VCN_c09_Syn006_N001_030dB_16000.0_MS.p

# this version reads the .prm file matching that in model_run


"""
from __future__ import print_function
import os
from pathlib import Path
import pickle
import argparse
import numpy as np
import pycircstat as PCS
from collections import OrderedDict
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as mpl

from cnmodel.util import sound
import vcnmodel.sac_campagnola as SAC
import pylibrary.plotting.plothelpers as PH 
from cnmodel.util import vector_strength
# import GIF_fit as GFit
import vcnmodel.vspfile

# import warnings
# warnings.filterwarnings("ignore")

class DisplayResult():
    def __init__(self):
        self.findfiles = False

        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelNameChoices = ['XM13', 'XM13_nacn', 'RM03', 'mGBC', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr', 'XM13_nacncoop']
        self.modelTypeChoices = ['II', 'II-I', 'I-II', 'I-c', 'I-t', 'II-o']
        self.SGCmodelChoices = ['Zilany', 'cochlea']  # cochlea is python model of Zilany data, no matlab, JIT computation; Zilany model creates matlab instance for every run.
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initAN', 'runANPSTH', 'runANIO', 'runANSingles', 'runANOmitOne', 'gifnoise']
        self.soundChoices = ['tonepip', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        self.speciesChoices = ['mouse', 'guineapig']
        self.spirouChoices = ['all', 'max=mean', 'all=mean']
        self.ANSynapseChoices = ['simple', 'multisite']

        self.modetypes = ['find', 'singles', 'IO', 'multi']
        self.Params = OrderedDict()
        self.fn = {}  # filename
        self.bp = {}  # base path
        self.ID = ''
        self.tag = 'delays'
        self.Params['modeltype'] = 'singles'
        self.Params['modelName'] = 'XM13'
        self.Params['modelType'] = 'II'
        self.Params['runProtocol'] = 'AN'
        self.Params['threshold'] = 0.
        self.Params['depression'] = 0
        self.Params['SR'] = 'MS'
        self.Params['pdf'] = None
        self.Params['soma_inflation'] = 1.0
        self.Params['dendrite_inflation'] = 1.0
        self.Params['noplot'] = False
        self.findfiles = False
        self.vspfile = 'VS.p'

    def make_filename(self, tag=''):
        print(self.Params)

        addarg=DR.Params['spirou']
        if addarg == 'all':
            addarg = 'all'
        elif addarg == 'all=mean':
            addarg = '_allmean'
        elif addarg == 'max=mean':
            addarg = '_mean'
        p = self.Params['cell']
        if 'AN' in self.Params['runProtocol']:
            simmode = 'AN'
        else:
            simmode = 'IV'
        outPath = Path('/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells', f"{p:s}", 'Simulations/', simmode)
        self.bp[p] = 'VCN_Cells/{0:s}/Simulations/AN'.format(p)
        # print(('bp: ', self.Params['cell'], self.bp))
        # if self.Params['soundtype'] in ['SAM', 'sam']:
        #     ofile = os.path.join(self.bp[p], 'AN_Result_' + self.Params['cell'] + '_%s_%s_%s_N%03d_%03ddB_%06.1f_FM%03.1f_DM%03d_%2s%s' %
        #         (tag, self.Params['modelName'], self.Params['ANSynapseType'], self.Params['nReps'],
        #         int(self.Params['dB']), self.Params['F0'],
        #         self.Params['fmod'], int(self.Params['dmod']), self.Params['SRType'], addarg) + '.p')
        #
        # else:
        #     ofile = os.path.join(self.bp[p], 'AN_Result_' + self.Params['cell'] + '_%s_%s_N%03d_%03ddB_%06.1f_%2s_%s' % (
        #         self.Params['modelName'], tag,
        #          self.Params['nReps'],
        #         int(self.Params['dB']), self.Params['F0'], self.Params['SRType'], addarg) + '.p')
        fn = f"AN_Result_{p:s}_{self.Params['modelName']:s}_{self.Params['modelType']:s}"
        if self.Params['inflationfactor'] != 1.0:
            fn += f"_soma={self.Params['inflationfactor']:.3f}"
        fn += f"_{addarg:s}_{self.Params['ANSynapseType']:s}"
        fn += f"_{self.Params['nReps']:03d}_{self.Params['soundtype']:s}"
        fn += f"_{int(self.Params['dB']):03d}dB_{self.Params['F0']:06.1f}"
        if self.Params['soundtype'] in ['SAM', 'sam']:
            fn += f"_{self.Params['fmod']:03.1f}_{int(self.Params['dmod']):03d}"
            fn += f"_{self.Params['SR']:2s}.p"
            ofile = Path(outPath, fn)
        else:
            fn += f"_{self.Params['SR']:2s}.p"
            ofile = Path(outPath, fn)
        print('outpuath: ', outPath)
        self.fn[p] = ofile
        self.filename = ofile
                
    def file_setup(self):
        if self.Params['modetype'] in ['multi']:
            # multiple cells
            self.Params['patterns'] = ['c%s' % p for p in self.Params['cell']]
            self.findfiles = True
        else:
            self.Params['patterns'] = ['c%s' % p for p in self.Params['cell']] #i for i in [8, 9, 17, 18, 19, 20, 21, 22]]
            self.findfiles = False

        print(self.Params['patterns'])
        for p in self.Params['patterns']:
            if self.Params['runProtocol'] in ['AN'] and self.Params['modetype'] not in ['IO']:
                self.fn[p] = self.make_filename()
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
            
            elif self.Params['runProtocol']in ['AN'] and self.Params['modetype'] in ['IO']:
                self.fn[p] = 'AN_Result_VCN_{0:s}_{1:s}_SynIO{2:03d}_N001_030dB_16000.0_MS.p'.format(p, self.Params['modeltype'], self.Params['synno'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
        
            elif self.Params['runProtocol'] in ['IV']:
                self.fn[p] = 'VCN_{0:s}_{1:s}.p'.format(p, self.Params['modeltype'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/IV'.format(p)
            
            elif self.Params['runProtocol'] in ['gifnoise']:
                self.fn[p] = 'VCN_{0:s}_{1:s}_gifnoise.p'.format(p, self.Params['modeltype'])
                self.bp[p] = 'VCN_Cells/VCN_{0:s}/Simulations/Noise'.format(p)
        self.lookupfiles()
        print(('Files: ', self.fn))
        
    def lookupfiles(self):  # simply look for files
        if self.findfiles:
            for p in self.Params['patterns']:
                try:
                    h = open(os.path.join(self.bp[p], self.fn[p]))
                    print('       found %s' % p)
                    h.close()
                except:
                    print('did not find file for %s' % p, self.fn[p])
                    exit(1)

    def show(self, pattern):
        self.Params['patterns'] = pattern
        print(self.Params['patterns'])
        #self.file_setup()
        print(('runtype: ', self.Params['runProtocol']))
        if self.Params['runProtocol'] == 'IV':
            self.show_IV()
            #fig, ax = mpl.subplots(len(self.Params['patterns'])+1, 3)
        if self.Params['runProtocol'] in ['AN', 'runANPSTH']:
            self.show_AN()
        if self.Params['runProtocol'] == 'gifnoise':
            self.show_gifnoise()

    def show_AN(self, filename=None):
        print('Show an')
        # bulid plot layout
        d = {}
        # for p in self.Params['patterns']:  # for each cell/run get the data
        #     h = open(self.fn[p], 'rb')
        #     d[p] = pickle.load(h)
        #     h.close()
        if filename is not None:
            with open(filename, 'rb') as fh:
                print('opening specific file: ', filename)
                d['c09'] = pickle.load(fh)
            self.Params['patterns'] = ['c09']
        sizer = OrderedDict([('A', {'pos': [0.08, 0.4, 0.71, 0.22]}), ('B', {'pos': [0.55, 0.4, 0.71, 0.22]}),
                             ('C', {'pos': [0.08, 0.4, 0.39, 0.22]}), ('D', {'pos': [0.55, 0.4, 0.39, 0.22]}),
                             ('E', {'pos': [0.08, 0.4, 0.07, 0.22]}), ('F', {'pos': [0.55, 0.4, 0.07, 0.22]}),
        ])  # dict elements are [left, width, bottom, height] for the axes in the plot.
        n_panels = len(sizer.keys())
        gr = [(a, a+1, 0, 1) for a in range(0, n_panels)]   # just generate subplots - shape does not matter
        axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.Plotter((n_panels, 1), axmap=axmap, label=True, figsize=(8., 6.))
        P.resize(sizer)  # perform positioning magic

        # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
        sk = sizer.keys()
        si = d[self.Params['patterns'][0]]['Params']
        totaldur = si['pip_start'] + np.max(si['pip_start']) + si['pip_duration'] + si['pip_offduration']
        for j, pattern in enumerate(self.Params['patterns']): # for all cells in the "pattern"

            print('pattern: ', pattern)
            # print(d[pattern])
            print('len d[pattern][trials]: ', len(d[pattern]['trials']))
            print(d[pattern]['trials'].keys())
            allt = []
            for i in range(len(d[pattern]['trials'])):  # for all trails in the measure.
                trial = d[pattern]['trials'][i]
                w = trial['stimWaveform']
                stb = trial['stimTimebase']
                # print(trial['somaVoltage'])
                P.axdict['A'].plot(trial['time']/1000., trial['somaVoltage'], linewidth=0.5)
                print(stb)
                print(w)
                # P.axdict['B'].plot(stb, w, linewidth=0.5)  # stimulus underneath
                P.axdict['C'].plot(trial['spikeTimes']/1000., i*np.ones(len( trial['spikeTimes'])),
                        'o', markersize=2.5, color='b')
                allt.append(trial['spikeTimes'][0]/1000.)
                inputs = len(trial['inputSpikeTimes'])
                for k in range(inputs):
                    tk = trial['inputSpikeTimes'][k]/1000.
                    y = (i+0.1+k*0.05)*np.ones(len(tk))
                    P.axdict['E'].plot(tk, y, '|', markersize=2.5, color='r', linewidth=0.5)
            # allt = np.array(allt).ravel()
            # # print(allt)
            # P.axdict['F'].hist(allt, range=(0., totaldur))
            # print('plotted hist)')
        
        
        for a in ['A', 'B', 'C', 'E', 'F']:  # set some common layout scaling
            P.axdict[a].set_xlim((0., totaldur))
            P.axdict[a].set_xlabel('T (s)')
        P.axdict['A'].set_ylabel('mV', fontsize=8)
        P.axdict['D'].set_xlabel('Phase', fontsize=8)
        P.axdict['C'].set_ylabel('Trial', fontsize=8)
        P.axdict['E'].set_ylabel('Trial, ANF', fontsize=8)
        P.axdict['B'].set_title('Stimulus', fontsize=9)
        P.axdict['E'].set_title('ANF Spike Raster', fontsize=9)
        P.axdict['C'].set_title('Bushy Spike Raster', fontsize=9)
        P.axdict['F'].set_title('PSTH', fontsize=9)
        # combine all spikes into one array, plot PSTH
        data = d[pattern]
        allst = []
        for trial in data['trials']:
            # print (trial)
            allst.extend(data['trials'][trial]['spikeTimes']/1000.)
        allst = np.array(allst)
        allst = np.sort(allst)
            # combine all spikes into one array, plot PSTH
            # the histogram of the data
        P.axdict['F'].hist(allst, 100, normed=1, facecolor='blue', alpha=0.75)

        if si['soundtype'] == 'SAM':  # calculate vs and plot histogram
            # print('allst: ', allst)
            phasewin = [data['Params']['pip_start'][0] + 0.25*data['Params']['pip_duration'], 
                data['Params']['pip_start'][0] + data['Params']['pip_duration']]
            # print (phasewin)
            spkin = allst[np.where(allst > phasewin[0])]
            spikesinwin = spkin[np.where(spkin <= phasewin[1])]

            # set freq for VS calculation
            if si['soundtype'] == 'tone':
                f0 = data['Params']['F0']
                print ("Tone: f0=%.3f at %3.1f dbSPL, cell CF=%.3f" % 
                        (data['Params']['f0'], data['Params']['dB'], data['Params']['F0'] ))
            if si['soundtype'] == 'SAM':
                f0 = data['Params']['fmod']
                tstring = ("SAM Tone: f0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f" %
                     (data['Params']['F0'], data['Params']['dB'], data['Params']['fmod'], data['Params']['dmod'], data['Params']['F0']))
                print(tstring)
                # P.figure_handle.suptitle(tstring + '\n' + self.fn[self.Params['patterns'][0]].replace('_', '\_'), fontsize=9)
           # print('spikes: ', spikesinwin)
            print((' Sound type: ', si['soundtype']))
            if len(spikesinwin) < 10:
                vs = 0.
                print('AN Vector Strength: Insufficient spike count, n = %d' % (len(spikesinwin)))
            else:
                vs = vector_strength(spikesinwin*1000., f0)  # vs expects spikes in msec
                print('AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (f0, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n']))
                vspfile.add_data(self.Params['patterns'][0], self.Params['spirou'], vs['r'])
                # fh = open(self.vspfile, 'rb')  # set file handle to write  - open, and append
                # vsp = pickle.load(fh)  # get vs data for all cells (outside this file)
                # fh.close
                # vsp[self.Params['patterns'][0]][self.Params['spirou']] = vs['r']  # update entry for this cell
                # fh = open(self.vspfile, 'wb')  # set file handle to write  - open, which will make an empty file
                # pickle.dump(vsp, fh)  # save it.
                # fh.close()
            # print(vs['ph'])
                P.axdict['D'].hist(vs['ph'], bins=2*np.pi*np.arange(30)/30.)
                P.axdict['D'].set_xlim((0., 2*np.pi))
                P.axdict['D'].set_title('Phase (VS = {0:.3f})'.format(vs['r']), fontsize=9, horizontalalignment='center')
                    
        # make figure output filename
        if self.Params['inflationfactor'] != 1.0:
            t = 'Inflated'
        else:
            t = 'Original'
        figname = str(self.fn) # f"{t:s}_{self.Params['patterns'][0]:s}_{self.Params['spirou']:s}.pdf"
        # P.figure_handle.suptitle(figname.replace('_', '\_'))
        mpl.savefig('figname_1.pdf')  # rasterized to 300 dpi is ok for documentation.
        mpl.show()
        # mpl.close()



    def show_SAC(self):
        sac = SAC.SAC()
        f, ax = mpl.figure()
        d = {}
        for p in self.Params['patterns']:
            h = open(os.path.join(self.bp[p], self.fn[p]))
            d[p] = pickle.load(h)
            h.close()

        fmod = d['patterns'[0]]['stimInfo']['fmod']  # modulation frequency
        #sacht = 2.5
        sacmax = 100
        #phaseht = 15

        f, ax2 = mpl.subplots(1)
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
            print('mean vs: for %s at fMod = %.2f:  %f angle: %f' % (pattern, fmod, vs, th.mean()))
            print('rayleigh: p=%f  z=%f (p > 0.05 means data is uniformly distributed)' % (p, z))
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
        mpl.show()

    def show_IV(self):
            fig, ax = mpl.subplots(2,1)
            a = ax.ravel()
            d = {}
            for p in self.Params['patterns']:
                h = open(os.path.join(self.bp[p], self.fn[p]))
                d[p] = pickle.load(h)
                h.close()
            print(('pattern: ',self.Params['patterns']))
            for j, pattern in enumerate(self.Params['patterns']):
                for k in range(len(d[pattern]['Results'])):
                    k0 = d[pattern]['Results'][k]
                    rkey = k0.keys()[0]
                    v = k0[rkey]['monitor']['postsynapticV']
                    a[0].plot(k0[rkey]['monitor']['time'], v, linewidth=0.5)
                    inj = k0[rkey]['monitor']['postsynapticI']
                    a[1].plot(k0[rkey]['monitor']['time'], inj)

            mpl.show()

    # def show_gifnoise(self):
    #         fig, ax = mpl.subplots(2,1)
    #         mpl.suptitle('gifnoise compare')
    #         dt_beforespike = 3.
    #         ax = ax.ravel()
    #         d = {}
    #         for p in self.Params['patterns']:
    #             h = open(os.path.join(self.bp[p], self.fn[p]))
    #             d[p] = pickle.load(h)
    #             h.close()
    #         for j, pattern in enumerate(self.Params['patterns']):
    #             for k in range(len(d[pattern]['Results'])):
    #                 k0 = d[pattern]['Results'][k]
    #                 rkey = k0.keys()[0]
    #                 data = k0[rkey]['monitor']
    #                 v = data['postsynapticV']
    #                 inj = data['postsynapticI']
    #                 dt = data['time'][1]-data['time'][0]
    #                 GF = GFit.GIFFitter(path=self.bp[pattern], dt=dt)
    #                 GF.GIF.gl      = 100.        # nS, leak conductance
    #                 GF.GIF.C       = 20.     # nF, capacitance
    #                 GF.GIF.El      = -65.0            # mV, reversal potential
    #
    #                 GF.GIF.Vr      = -55.0            # mV, voltage reset
    #                 GF.GIF.Tref    = 0.2             # ms, absolute refractory period
    #
    #                 GF.GIF.Vt_star = -45.0            # mV, steady state voltage threshold VT*
    #                 GF.GIF.DV      = 10.             # mV, threshold sharpness
    #               #GF.set_templates(aec=None, train=train_template, test=None)
    #                 GF.Exp.addTrainingSetTrace(v, 1e-3, inj, 1e-9, np.max(data['time']), 'Array')
    #                 #GF.set_files_test(AECtrace=None, trainingset=1008, testsets = range(1009, 1018), filetype='Igor')
    #                 #GF.set_timescales(ts=[0.1, 0.5, 2., 5.0, 10.0, 20.0])
    #                 print('Starting Parameters: ')
    #                 GF.GIF.printParameters()
    #                 GF.fit(threshold=self.Params['threshold'],
    #                         ax=ax[0], current=inj, beforeSpike=dt_beforespike)
    #                 print('Final Parameters: ')
    #                 GF.GIF.printParameters()
    #                 print(('v: ', v[0]))
    #                 print(('len inj, dt: ', len(inj), GF.dt, len(inj)*GF.dt))
    #                 (time, Vs, I_a, V_t, S) = GF.GIF.simulate(inj, v[0])  # simulate response to current trace I with starting voltage V0
    #                 # a[0].plot(time, Vs, 'r-', linewidth=0.75)
    # #                 a[0].plot(data['time'], v, 'b-', linewidth=0.5)
    #                 ax[1].plot(data['time'], inj, 'k-', linewidth=0.5)
    #                 #GF.GIF.plotParameters()
    #         mpl.show()
    #         exit(0)


if __name__ == '__main__':

    DR = DisplayResult()
    
    parser = argparse.ArgumentParser(description='Display model results')
    parser = argparse.ArgumentParser(description='Simulate activity in a reconstructed model cell',
                    argument_default=argparse.SUPPRESS,
                    fromfile_prefix_chars='@')
    parser.add_argument('--cell', dest='cell', action='store',
                   default=None,
                   help='Select the cell (no default)')
    parser.add_argument('--type', '-T', dest='cellType', action='store',
                   default='Bushy', choices=DR.cellChoices,
                   help='Define the cell type (default: Bushy)')
    parser.add_argument('--model', '-M', dest='modelName', action='store',
                   default='XM13', choices=DR.modelNameChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--modeltype', dest='modelType', action='store',
                   default='II', choices=DR.modelTypeChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--sgcmodel', type=str, dest='SGCmodelType', action='store',
                   default='Zilany', choices=DR.SGCmodelChoices,
                   help='Define the SGC model type (default: Zilany)')
    parser.add_argument('--protocol', '-P', dest='runProtocol', action='store',
                   default='runIV', choices=DR.protocolChoices,
                   help='Protocol to use for simulation (default: IV)')
    parser.add_argument('-H', action='store_true', dest='usedefaulthoc',
                  help='Use default hoc file for this cell')
    parser.add_argument('--hocfile', dest='hocfile', action='store',
                  default=None,
                  help='hoc file to use for simulation (default is the selected "cell".hoc)')
    parser.add_argument('--inputpattern', '-i', type=str, dest='inputPattern', action='store',
                  default=None,
                  help='cell input pattern to use (substitute) from cell_config.py')
    parser.add_argument('--stimulus', '-s', type=str, dest='soundtype', action='store',
                   default='tonepip', choices=DR.soundChoices,
                   help='Define the stimulus type (default: tonepip)')
    parser.add_argument('--check', '-/', action='store_true', default=False, dest='checkcommand',
                   help='Only check command line for valid input; do not run model')
    parser.add_argument('--noplot', action='store_false', default=False, dest='noplot',
                   help='Suppress plots')
    parser.add_argument('--inflation', type=float, dest='inflationfactor', default=1.0,
                help='select inflation factor')
                   
    # lowercase options are generally parameter settings:
    parser.add_argument('-d', '--dB', type=float, default=30., dest='dB',
        help='Set sound intensity dB SPL (default 30)')
    parser.add_argument('-f', '--frequency', type=float, default=4000., dest='F0',
        help='Set tone frequency, Hz (default 4000)')
    parser.add_argument('--duration', type=float, default=0.1, dest='pip_duration',
        help='Set sound stimulus duration (sec; default 0.1)')
    parser.add_argument('-r', '--reps', type=int, default=1, dest = 'nReps',
        help='# repetitions')
    parser.add_argument('-S', '--SRType', type=str, default='HS', dest = 'SRType',
        choices=DR.SRChoices,
        help=('Specify SR type (from: %s)' % DR.SRChoices))
    parser.add_argument('--synapsetype', type=str, default='multisite', dest = 'ANSynapseType',
        choices=DR.ANSynapseChoices,
        help=('Specify AN synpase type (from: %s)' % DR.ANSynapseChoices))
    parser.add_argument('--depression', type=int, default=0, dest = 'ANSynapticDepression',
        choices=[0, 1],
        help=('Specify AN depression flag for multisite synapses (from: %s)' % str([0, 1])))
        
    parser.add_argument('--fmod', type=float, default=20, dest = 'fmod',
        help='Set SAM modulation frequency')
    parser.add_argument('--dmod', type=float, default=100., dest = 'dmod',
        help='Set SAM modulation depth (in percent)')
    parser.add_argument('--S2M', type=float, default=0, dest = 'signalToMasker',
        help='Signal to Masker ratio (dB)')
    parser.add_argument('--cmmrmode', type=str, default='CMR', dest = 'CMMRmode',
        choices=DR.cmmrModeChoices,
        help=('Specify mode (from: %s)' % DR.cmmrModeChoices))

    parser.add_argument('--soma-inflate', type=float, dest='soma_inflation', action='store', default=1.0,
            help='Specify factor by which to inflate soma AREA')
    parser.add_argument('--soma-autoinflate', action='store_true', dest='soma_autoinflate', default=False,
            help='Automatically inflate soma based on table')
    parser.add_argument('--dendrite-inflate', type=float, dest='dendrite_inflation', action='store', default=1.0,
            help='Specify factor by which to inflate total dendritic AREA')
    parser.add_argument('--dendrite-autoinflate', action='store_true', dest='dendrite_autoinflate', default=False,
            help='Automatically inflate dendrite area based on table')
                
    parser.add_argument('--spirou', type=str, dest='spirou', action='store', default='all',
            choices = DR.spirouChoices,
            help='Specify spirou experiment type.... ')
    parser.add_argument('-n', '--synno', dest='synno', action='store',
                    default=None, help='Select the synno')
    parser.add_argument('-t', '--threshold', dest='threshold', type=float, action='store',
                    default=0., help=('Spike detection threshold, mV'))
    parser.add_argument('--last', action='store_true',  dest='lastfile', default=False,
                help='Just the file in lastfile.txt')
    parser.add_argument('-F', '--filename', type=str, default=None, dest='filename',
                    help="just use a filename")

    args = vars(parser.parse_args())
#    model.print_modelsetup()
    print(args)


    for k in args.keys():
        DR.Params[k] = args[k]
        print('k: ', k)

    print(args['lastfile'])
    if args['filename']:
        DR.filename = args['filename']
        DR.show_AN(DR.filename)

        exit()

    elif args['lastfile']:
        lf = Path('lastfile.txt').read_text()
        DR.filename = lf
        print('reading: ', DR.filename)
        
    DR.make_filename(tag='delays')
    if os.path.isfile(DR.filename):
        print(('FOUND ofile: ', DR.filename))
        DR.show(pattern=[DR.Params['cell']])
    else:
        print(('NOT FOUND: ', DR.filename))
    #DR.show()
        