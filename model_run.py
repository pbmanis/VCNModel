__author__ = 'pbmanis'


"""
celltype options:
Bushy_RM03
Bushy_XM13
Calyx
Stellate
MNTB
L23pyr (sort of... not a very good rendition)

"""

import sys
import os.path
import pickle
import time
import argparse
from collections import OrderedDict
import pprint
import json

#import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer
from neuronvis.hoc_reader import HocReader
import neuronvis.hoc_graphics as hoc_graphics
from channel_decorate import ChannelDecorate
from generate_run import GenerateRun
import cellInitialization as cellInit

from nrnlibrary.protocols import protocol
from nrnlibrary import cells
from nrnlibrary import synapses
from nrnlibrary.util import get_anspikes
from nrnlibrary.util import sound
import nrnlibrary.util as nu

import pylibrary.Utility as pu  # access to spike finder routine
import pyqtgraph.multiprocess as mproc

from neuron import h
import neuron

try:
    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui
    import pylibrary.pyqtgraphPlotHelpers as pgh
    HAVE_PG = True
except:
	HAVE_PG = False

if HAVE_PG:
    import render

import numpy as np

verbose = False
showCell = False


        #infile = 'L23pyr.hoc'
        #infile = 'LC_nmscaled_cleaned.hoc'
        #infile = 'Calyx-68cvt2.hoc'
        #infile = 'Calyx-S53Acvt3.hoc'
        #infile = 'wholeThing_cleaned.hoc'
        #infile = 'MNTB_Cell2_cleaned.hoc'
        #infile = 'VCN_Dend.hoc'
        #infile = 'somaOnly.hoc'

class ModelRun():
    def __init__(self, args=None):

        # use v2 files for model with rescaled soma
        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelChoices = ['XM13', 'RM03', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr']
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initAN', 'runANPSTH', 'runANSingles']
        self.soundChoices = ['tone', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        
        # IV_neuronStateFile = 'an_neuronstateV2.dat'
        # make_init = False
        # test_init = False
        # IV_mode = False
        #
        # # AN
        # AN_neuronStateFile = 'aniv_neuronstateV2.dat'
        # make_ANIntialConditions = False
        # ANPSTH_mode = True
        # ANSingles = False
        self.srname = ['**', 'LS', 'MS', 'HS']  # runs 1-3, not starting at 0
        
        self.Params = OrderedDict()
        
        self.Params['initIVStateFile'] = 'an_neuronstateV2.dat'
        self.Params['initANStateFile'] = 'aniv_neuronstateV2.dat'
        self.Params['infile ']= ''

        self.Params['cellType'] = self.cellChoices[0]
        self.Params['modelType'] = self.modelChoices[0]
        self.Params['SRType'] = self.SRChoices[2]
        self.Params['SR'] = self.Params['SRType']  # actually used SR this might be cell-defined, rather than command defined
        self.Params['runProtocol'] = self.protocolChoices[2]  # testIV is default because it is fast and should be run often

        self.Params['run_duration'] = 0.25 # in sec
        self.Params['pip_duration'] = 0.1
        self.Params['pip_start'] = [0.1]
        self.Params['Fs'] = 100e3
        self.Params['F0'] = 4000.
        self.Params['dB'] = 40.
        self.Params['RF'] = 2.5e-3
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params['threshold'] = -20
        self.Params['plotFlag'] = False
        
    def printModelSetup(self):
        for p in self.Params.keys():
            print '%18s = ' % p, self.Params[p]
        print '-----------'

    def set_celltype(self, cellType):
        if cellType not in self.cellChoices:
            print 'Celltype must be one of: %s. Got: %s', (', '.join(self.cellChoices), cellType)
            exit()
        self.Params['cellType'] = cellType


    def set_modeltype(self, modelType):
        if modelType not in self.modelChoices:
            print 'Model type must be one of: %s. Got: %s ' % (', '.join(self.modelChoices), modelType)
            exit()
        self.Params['modelType'] = modelType


    def set_SR(self, SRType):
        if self.SRType not in self.SRChoices:
            print 'SR type must be one of: %s. Got: %s ' % (', '.join(self.SRChoices), SRType)
            exit()
        self.Params['SRType'] = SRType


    def set_starttime(self, starttime):
        self.Params['StartTime'] = starttime



    def runModel(self, parMap={}):
        if verbose:
            print 'runModel entry'
        if 'id' in parMap.keys():
            self.idnum = parMap['id']
        else:
            self.idnum = 9999

        if parMap == {}:
            self.plotFlag = True
        filename = os.path.join('MorphologyFiles/', self.Params['infile'])
        print 'reading input file: %s' % filename
        self.hf = HocReader(filename)
        self.hf.h.celsius = 38.
        self.hf.h.Ra = 150.
        print 'Ra is: ', self.hf.h.Ra
        print 'Temp is: ', self.hf.h.celsius
        for group in self.hf.sec_groups.keys():
            g = self.hf.sec_groups[group]
            for section in list(g):
                self.hf.get_section(section).Ra = self.hf.h.Ra
                if verbose:
                    print 'section: ', section
                    print 'Ra: ', self.hf.get_section(section).Ra

        electrodeSection = list(self.hf.sec_groups['soma'])[0]
        self.electrodeSite = self.hf.get_section(electrodeSection)
        self.electrodeSection = 'soma'
        self.hg = hoc_graphics
        self.get_hoc_file(filename)
        sg = self.hf.sec_groups['soma']
        self.distances(self.hf.get_section(list(sg)[0]).name()) # make distance map from soma
        if verbose:
            print 'Parmap in runModel: ', parMap
        self.cd = ChannelDecorate(self.hf, celltype=self.Params['cellType'], modeltype=self.Params['modelType'],
                             parMap=parMap)

        # self.hf.h.topology()
        # self.cd.channelValidate(self.hf, verify=False)

        #return
        # for group in self.hf.sec_groups.keys():
        #     g = self.hf.sec_groups[group]
        #     for section in list(g):
        #         secinfo = self.hf.get_section(section)
        # handle the following protocols:
        # ['initIV', 'initAN', 'runIV', 'runANPSTH', 'runANSingles']
        
        if self.Params['runProtocol'] == 'initIV':
            if verbose:
                print 'initIV'
            self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag)
            cellInit.getInitialConditionsState(self.hf, tdur=3000., 
                filename=self.Params['initIVStateFile'], electrodeSite=self.electrodeSite)
            print 'Ran to get initial state for %f msec' % self.hf.h.t
            return

        if self.Params['runProtocol'] == 'testIV':
            if verbose:
                print 'test_init'
            self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag, )
            cellInit.testInitialConditions(self.hf, filename=self.Params['initIVStateFile'],
                electrodeSite=self.electrodeSite)
            #self.R.testRun()
            return  # that is ALL, never make init and then keep running.

        if self.Params['runProtocol'] == 'runANPSTH':
            if verbose:
                print 'ANPSTH'
            self.ANRun(self.hf)

        if self.Params['runProtocol'] == 'initAN':
            if verbose:
                print 'Init AN'
            self.ANRun(self.hf, make_ANIntialConditions=True)

        if self.Params['runProtocol'] == 'runANSingles':
            if verbose:
                print 'ANSingles'
            self.ANRun_singles(self.hf)
            
        if self.Params['runProtocol'] == 'runIV':
            if verbose:
                print 'iv_mode'
            self.IVRun(parMap)

        if showCell:
            self.render = HocViewer(self.hf)
            cylinder=self.render.draw_cylinders()
            cylinder.set_group_colors(self.section_colors, alpha=0.8, mechanism=['nav11', 'gbar'])


    def IVRun(self, parMap={}):
        if verbose:
            print 'generateRun'
        print 'IVRun begins'
        self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag)

        if verbose:
            print 'doRun'
        self.R.doRun(self.Params['infile'], parMap, save='monitor', restoreFromFile=True)
        if verbose:
            print '  doRun completed'
            print self.R.IVResult
        #basename = self.R.saveRuns(self.R.results)
        #self.R.arun.saveIVResult(basename)
        isteps = self.R.IVResult['I']
        if verbose:
            for k, i in enumerate(self.R.IVResult['tauih'].keys()):
                print 'ih: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['tauih'][i]['tau'].value)
                print '        dV : %f' % self.R.IVResult['tauih'][i]['a'].value
            for k, i in enumerate(self.R.IVResult['taus'].keys()):
                print 'i: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['taus'][i]['tau'].value)
                print '       dV : %f' % (self.R.IVResult['taus'][i]['a'].value)

            print 'Nspike, Ispike: ', self.R.IVResult['Nspike'], self.R.IVResult['Ispike']
            print 'Rinss: ', self.R.IVResult['Rinss']
            print 'Vm: ', np.mean(self.R.IVResult['Vm'])
        print 'ivresult keys: ', self.R.IVResult['taus'].keys()
        if len(self.R.IVResult['taus'].keys()) == 0:
            taum_mean = 0.
            tauih_mean = 0.
        else:
            taum_mean = np.mean([self.R.IVResult['taus'][i]['tau'].value for k, i in enumerate(self.R.IVResult['taus'].keys())])
            tauih_mean = np.mean([self.R.IVResult['tauih'][i]['tau'].value for k, i in  enumerate(self.R.IVResult['tauih'].keys())])
        #print 'taum_mean: ', taum_mean
        #print 'tauih_mean: ', tauih_mean
        # construct dictionary for return results:
        self.IVSummary = {'basefile': self.R.basename,
                          'parmap': parMap, 'ID': self.idnum,
                          'Vm': np.mean(self.R.IVResult['Vm']),
                          'Rin': self.R.IVResult['Rinss'],
                          'taum': taum_mean, 'tauih': tauih_mean,
                          'spikes': {'i': self.R.IVResult['Ispike'], 'n': self.R.IVResult['Nspike']},
                          }

        #print 'ivsummary: ', self.IVSummary
        print 'model_run::runModel: write summary for file set = ', self.R.basename
        return self.IVSummary

    #synconfig consists of a list of tuples.
    # Each element in the list corresponds to one terminal and all of it's active zones
    # each tuple consists of N sites (calculated from area * average synapses/um2)
    # delay, and SR
    # the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
    # and the diameter. The brances are not included. distances are in microns
    
    VCN_c18_synconfig = [[int(216.66*0.65), 0., 2, 49.2, 1.222],
                         [int(122.16*0.65), 0., 2, 82.7, 1.417], 
                         [int(46.865*0.65), 0., 2, 67.3, 1.309],
                         [int(84.045*0.65), 0., 2, 22.4, 1.416],
                         [int(80.27*0.65),  0., 2, 120.3, 0.687],
                        ]
    # VCN_c18_synconfig_original = [(int(216.66*0.65), 0., 2), (int(122.16*0.65), 0., 2),
    #     (int(46.865*0.65), 0., 2), (int(84.045*0.65), 0., 2), (int(2.135*0.65), 0, 2), (int(3.675*0.65), 0, 2), (int(80.27*0.65), 0, 2)]
    test_synconfig = [(80, 0, 2)]  # if switching configs, need to re-run with AN_InitialConditions to set up the
    # save/restore state. Any change in input configuration requires a new "state" 


    def ANRun(self, hf, verify=False, seed=0, synapseConfig=VCN_c18_synconfig, make_ANIntialConditions=False):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        """

        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        for i, s in enumerate(synapseConfig):
            s[1] = s[3]*0.001/3.0
            print 'delay for input %d is %8.4f msec' % (i, s[1])

        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV

        stimInfo = self.Params
        # {'Morphology': self.Params['infile'], 'synapseConfig': synapseConfig,
        #             'runDur': self.run_duration, 'pip_dur': self.Params['pip_duration'], 'pip_start': self.Params['pip_start'],
        #             'run_duration': self.run_duration,
        #             'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
        #             'cellType': self.Params['cellType'], 'modelType': self.Params['modelType'], 'nReps': nReps, 'threshold': threshold}

        preCell, postCell, synapse, self.electrodeSite = self.configureCell(hf, synapseConfig, stimInfo)

        # see if we need to save the cell state now.
        if make_ANIntialConditions:
            print 'getting initial conditions for AN'
            cellInit.getInitialConditionsState(hf, tdur=3000., 
                filename=self.Params['initANStateFile'], electrodeSite=self.electrodeSite)
            cellInit.testInitialConditions(hf, filename=self.Params['initANStateFile'],
                electrodeSite=self.electrodeSite)
            return

        seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        print 'AN Seeds: ', seeds
        stimInfo['seeds'] = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        celltime = []
        
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.

        nWorkers = self.Params['nWorkers']
        TASKS = [s for s in range(nReps)]
        tresults = [None]*len(TASKS)

        # run using pyqtgraph's parallel support
        with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
            for j, x in tasker:
                tresults = self.singleANRun(hf, j, synapseConfig,
                    stimInfo, seeds, preCell, postCell, self.an_setup_time)
                tasker.results[j] = tresults
        # retreive the data
        for j, N in enumerate(range(nReps)):
           
            celltime.append(tresults[j]['time']) # (self.time)
            spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'],
                    threshold, t0=0., t1=stimInfo['run_duration']*1000., dt=1.0, mode='peak')
            inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
            somaVoltage[N] = np.array(tresults[j]['Vsoma'])

        # for j, N in enumerate(range(nReps)):
        #     print 'Rep: %d' % N
        #
        #     preCell, postCell = self.singleANRun(hf, j, synapseConfig, stimInfo, seeds, preCell, postCell, self.an_setup_time)
        #
        #     celltime.append(self.time)
        #     spikeTimes[N] = pu.findspikes(self.time, self.Vsoma, threshold, t0=0., t1=hf.h.tstop, dt=1.0, mode='peak')
        #     inputSpikeTimes[N] = [preCell[i]._spiketrain for i in range(len(preCell))]
        #     somaVoltage[N] = np.array(self.Vsoma)

        total_elapsed_time = time.time() - self.start_time
#        total_run_time = time.time() - run_time
        print "Total Elapsed Time = %8.2f min (%8.0fs)" % (total_elapsed_time/60., total_elapsed_time)
        print "Total Setup Time = %8.2f min (%8.0fs)" % (self.setup_time/60., self.setup_time)
        print "Total AN Calculation Time = %8.2f min (%8.0fs)" % (self.an_setup_time/60., self.an_setup_time)
        print "Total Neuron Run Time = %8.2f min (%8.0fs)" % (self.nrn_run_time/60., self.nrn_run_time)
        
        result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes, 
            'somaVoltage': somaVoltage, 'time': np.array(celltime[0])}
        
        self.analysis_filewriter(self.Params['infile'], result, tag='delays')
        self.plotAN(np.array(celltime[0]), result['somaVoltage'], result['stimInfo'])


    def ANRun_singles(self, hf, verify=False, seed=0, synapseConfig=VCN_c18_synconfig):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        This routine is special - it runs nReps for each synapse, turning off all of the other synapses
        by setting the synaptic conductance to 0 (to avoid upsetting initialization)
        
        """

        self.start_time = time.time()

        nReps = self.Params['nReps']
        threshold = -20. # spike threshold, mV

        stimInfo = self.Params
        # {'Morphology': self.Params['infile'], 'synapseConfig': synapseConfig,
        #             'runDur': self.run_duration, 'pip_dur': self.Params['pip_duration'], 'pip_start': self.Params['pip_start'],
        #             'run_duration': self.run_duration,
        #             'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
        #             'cellType': self.Params['cellType'], 'modelType': self.Params['modelType'], 'nReps': nReps, 'threshold': threshold}

        preCell, postCell, synapse, self.electrodeSite = self.configureCell(hf, synapseConfig, stimInfo)

        # see if we need to save the cell state now.
        if make_ANIntialConditions:
            print 'getting initial conditions for AN'
            cellInit.getInitialConditionsState(hf, tdur=3000., 
                filename=self.Params['initANStateFile'], electrodeSite=self.electrodeSite)
            cellInit.testInitialConditions(hf, filename=self.Params['initANStateFile'],
                electrodeSite=self.electrodeSite)
            return
        nSyns = len(synapseConfig)
        seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        print 'AN Seeds: ', seeds
        stimInfo['seeds'] = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        celltime = []
        
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.
        # get the gMax's
        gMax = np.zeros(nSyns)
        for i, s in enumerate(synapse):
            for p in s.psd.ampa_psd:
                gMax[i] = p.gmax
        print 'synapse gMax: ', gMax
        
        for k in range(nSyns):
            # only enable gsyn on the selected input
            for i, s in enumerate(synapse):
                for p in s.psd.ampa_psd:
                    if i != k:
                        p.gmax = 0.  # disable all
                    else:
                        p.gmax = gMax[i]  # except the shosen one

            nWorkers = self.Params['nWorkers']
            TASKS = [s for s in range(nReps)]
            tresults = [None]*len(TASKS)

            # run using pyqtgraph's parallel support
            with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
                for j, x in tasker:
                    tresults = self.singleANRun(hf, j, synapseConfig,
                        stimInfo, seeds, preCell, postCell, self.an_setup_time)
                    tasker.results[j] = tresults
            # retreive the data
            for j, N in enumerate(range(nReps)):
#                print 'Rep: %d' % N

#                preCell, postCell = self.singleANRun(hf, j, synapseConfig, stimInfo, seeds, preCell, postCell)
                
                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'], threshold, t0=0.,
                        t1=stimInfo['run_duration']*1000, dt=1.0, mode='peak')
                inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]['Vsoma'])

            total_elapsed_time = time.time() - self.start_time
    #        total_run_time = time.time() - run_time
            print "Total Elapsed Time = %8.2f min (%8.0fs)" % (total_elapsed_time/60., total_elapsed_time)
            print "Total Setup Time = %8.2f min (%8.0fs)" % (self.setup_time/60., self.setup_time)
            print "Total AN Calculation Time = %8.2f min (%8.0fs)" % (self.an_setup_time/60., self.an_setup_time)
            print "Total Neuron Run Time = %8.2f min (%8.0fs)" % (self.nrn_run_time/60., self.nrn_run_time)
        
            result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes, 
                'somaVoltage': somaVoltage, 'time': np.array(tresults[j]['time'])}
        
            self.analysis_filewriter(self.Params['infile'], result, tag='Syn%03d' % k)
        #self.plotAN(np.array(self.time), result['somaVoltage'], result['stimInfo'])


    def configureCell(self, hf, synapseConfig, stimInfo):
        sg = hf.sec_groups['soma']
        postCell = cells.Generic.create(soma=hf.get_section(list(sg)[0]))

        self.cd.channelValidate(hf, verify=False)

        preCell = []
        synapse = []
        # reconfigure syanpses to set the spont rate group
        stimInfo['SR'] = stimInfo['SRType']
        for i, syn in enumerate(synapseConfig):
            if stimInfo['SRType'] is 'fromcell':  # use the one in the table
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=syn[2]))
                stimInfo['SR'] = self.srname[syn[2]] # use and report value from table
            else:
                try:
                    srindx = self.srname.index(stimInfo['SRType'])
                    print 'retrieved index %d with SR type %s' % (srindx, stimInfo['SRType'])
                except:
                    raise ValueError('SR type "%s" not found in Sr type list' % stimInfo['SRType'])
                    
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=syn[srindx]))  # override
            synapse.append(preCell[-1].connect(postCell, pre_opts={'nzones':syn[0], 'delay':syn[1]}))
        for i, s in enumerate(synapse):
            s.terminal.relsite.Dep_Flag = 0  # turn off depression computation
            #print dir(s.psd.ampa_psd[0])

        #  ****** uncomment here to adjust gmax.

            # for p in s.psd.ampa_psd:
            #     p.gmax = p.gmax*2.5
            #
            #print 'Depression flag for synapse %d = %d ' % (i, s.terminal.relsite.Dep_Flag)

        # checking the gmax for all synapses
        # for s in synapse:
        #     psds = s.psd.ampa_psd
        #     for p in psds:
        #         #print p.gmax
        #         p.gmax = p.gmax
        # for s in synapse:
        #     psds = s.psd.ampa_psd
        #     for p in psds:
        #         #print p.gmax
        #         pass
        electrodeSection = list(hf.sec_groups['soma'])[0]
        electrodeSite = hf.get_section(electrodeSection)
        return (preCell, postCell, synapse, electrodeSite)


    def singleANRun(self, hf, j, synapseConfig, stimInfo, seeds, preCell, postCell, an_setup_time):
        cellInit.restoreInitialConditionsState(hf, electrodeSite=None, filename=self.Params['initANStateFile'])
        stim=[]
        # make independent inputs for each synapse
        ANSpikeTimes=[]
        an0_time = time.time()
        nrn_run_time = 0.
        for i, syn in enumerate(synapseConfig):
            stim.append(sound.TonePip(rate=stimInfo['Fs'], duration=stimInfo['run_duration'], 
                                  f0=stimInfo['F0'], dbspl=stimInfo['dB'],
                                  ramp_duration=stimInfo['RF'], pip_duration=stimInfo['pip_duration'],
                                  pip_start=stimInfo['pip_start']))
            print stim[-1]
            nseed = seeds[j, i]
            preCell[i].set_sound_stim(stim[-1], seed=nseed)  # generate spike train, connect to terminal
            ANSpikeTimes.append(preCell[i]._spiketrain)
        an_setup_time += (time.time() - an0_time)
        nrn_start = time.time()
        Vsoma = hf.h.Vector()
        rtime = hf.h.Vector()
        Vsoma.record(postCell.soma(0.5)._ref_v, sec=postCell.soma)
        rtime.record(hf.h._ref_t)
        print '...running '
        hf.h.finitialize()
        hf.h.tstop = stimInfo['run_duration']*1000.
        hf.h.t = 0.
        hf.h.batch_save() # save nothing
        hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
        # old - don't d this! #hf.h.finitialize()
        # hf.h.run()
        nrn_run_time += (time.time() - nrn_start)
        print '...done'
        anresult = {'Vsoma': np.array(Vsoma), 'time': np.array(rtime), 'ANSpikeTimes': ANSpikeTimes}
        return anresult
        

    def plotAN(self, celltime, somaVoltage, stimInfo):
        nReps = stimInfo['nReps']
        threshold = stimInfo['threshold']
        win = pgh.figure(title='AN Inputs')
        layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
        for j, N in enumerate(range(len(somaVoltage))):
            layout.plot(0, celltime, somaVoltage[N], pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.plot(0, [np.min(celltime), np.max(celltime)], [threshold, threshold], pen=pg.mkPen((0.5, 0.5, 0.5), width=0.5))
            dvdt = np.diff(somaVoltage[N])
            layout.plot(1, celltime[:-1], dvdt, pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.getPlot(0).setXLink(layout.getPlot(1))

        pgh.show()


    def analysis_filewriter(self, filebase, result, tag=''):
        k = result.keys()
        requiredKeys = ['stimInfo', 'spikeTimes', 'inputSpikeTimes', 'somaVoltage', 'time']
        for rk in requiredKeys:
            assert rk in k

        stimInfo = result['stimInfo']
        fname = os.path.splitext(self.Params['infile'])[0]
        f = open('AN_Result_' + fname + '_%s_N%03d_%03ddB_%06.1f_%2s' % (tag, stimInfo['nReps'],
                int(stimInfo['dB']), stimInfo['F0'], stimInfo['SR']) + '.p', 'w')
        pickle.dump(result, f)
        f.close()        

    def distances(self, section):
        self.hf.distanceMap = {}
        self.hf.h('access %s' % section) # reference point
        d = self.hf.h.distance()
        for si in self.hf.sections.keys():
            self.hf.h('access %s' % si)
            self.hf.distanceMap[si] = self.hf.h.distance(0.5) # should be distance from first point


    def get_hoc_file(self, infile):
        if self.hf.file_loaded is False:
            exit()
        self.section_list = self.hf.get_section_prefixes()
#        print 'section groups: ', self.hf.sec_groups.keys()
        self.hf.sec_groups.keys()
        if len(self.hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(self.hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.colorMap[i]
#        else: # single section name, assign colors to SectionList types:
#        self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
#            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

        (v, e) = self.hf.get_geometry()
        self.clist = []

        for si in self.hf.sections: # self.section_list[s]:
            self.hf.h('access %s' % si)
            sr = self.hf.h.SectionRef()
            n1 = self.hf.h.cas().name()
            if sr.has_parent() == 1:
                x=sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])


if __name__ == "__main__":
    model = ModelRun(sys.argv[1:])  # create instance of the model
    pp = pprint.PrettyPrinter(indent=4, width=60)
    
    parser = argparse.ArgumentParser(description='Simulate activity in a reconstructed model cell')
    parser.add_argument(dest='infile', action='store',
                   default='None', 
                   help='Select the input file (no default)')
    parser.add_argument('--type', '-T', dest='cellType', action='store',
                   default='Bushy', choices=model.cellChoices,
                   help='Define the cell type (default: Bushy)')
    parser.add_argument('--model', '-M', dest='modelType', action='store',
                   default='XM13', choices=model.modelChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--protocol', '-P', dest='runProtocol', action='store',
                   default='IV', choices=model.protocolChoices,
                   help='Protocol to use for simulation (default: IV)')
                  
    # lowercase options are generally parameter settings:

    parser.add_argument('-r', '--reps', type=int, default=1, dest = 'nReps',
        help='# repetitions')
    parser.add_argument('--modf', type=float, default=10, dest = 'modFreq', 
        help='Set SAM modulation frequency')
    parser.add_argument('--moddepth', type=float, default=100., dest = 'modDepth', 
        help='Set SAM modulation depth (in percent)')
    parser.add_argument('--S2M', type=float, default=0, dest = 'signalToMasker',
        help='Signal to Masker ratio (dB)')
    parser.add_argument('--cmmrmode', type=str, default='CMR', dest = 'CMMRmode',
        choices=model.cmmrModeChoices, 
        help=('Specify mode (from: %s)' % model.cmmrModeChoices))
    parser.add_argument('--allmodes', action="store_true", default = False, dest = 'all_modes',
        help=('Force run of all modes (CMR, CMD, REF) for stimulus configuration.'))
    parser.add_argument('-S', '--SRType', type=str, default='HS', dest = 'SRType',
        choices=model.SRChoices,
        help=('Specify SR type (from: %s)' % model.SRChoices))
    parser.add_argument('--sequence', type=str, default='[1,2,5]', dest = 'sequence',
            help=('Specify a sequence for the primary run parameters'))
    parser.add_argument('--plot',  action="store_true", default=False, dest = 'showPlot', 
            help='Plot results as they are generated - requires user intervention... ')
type=int,  default=4, dest = 'nWorkers', 
            help='Number of "workers" for parallel processing (default: 4)')

    # parser.add_argument('-p', '--print', action="store_true", default=False, dest = 'print_info',
    #     help='Print extra information during analysis')
    #  parser.add_argument('-l', '--list', action="store_true", default=False, dest = 'list_results',
    #     help='List results to screen')


    args=vars(parser.parse_args())   
#    model.printModelSetup()
    for k in args.keys():
        model.Params[k] = args[k]
    print(json.dumps(model.Params, indent=4))  # pprint doesn't work well with ordered dicts

    
    model.runModel() # then run the model
    if showCell:
        QtGui.QApplication.instance().exec_()
