__author__ = 'pbmanis'

"""
model_run.py

Run a model based on a hoc structure, decorating the structure with ion channels and synapses.

Directory structure:
VCN_Cells/   # top level for data
    cell_ID/    # VCN_c18, for example (first argument in call should be this directory name)
        MorphologyFiles/  # location for swc and hoc files
            VCN_c18_755V2.hoc (cell body scaled version)
            VCN_c18_755.hoc  (unscaled version)
            misc.swc  (various swc files that were translated to hoc files)
        InitializationFiles/  # location for Init files
            IV/  # initialization for just the IV with no synaptic input
                VCN_c18_755v2.ninit  # different base structures
                VCN_c18755.ninit
            AN/  # initialization for synaptic inputs
                (ditto)  # different base structures of input arrangements
        Simulations/
            IV/  # results from IV simulations
            AN/  # results from AN simulations
         
Usage:
Simulate activity in a reconstructed model cell

positional arguments:
  cell                  Select the cell (no default)

optional arguments:
  -h, --help            show this help message and exit
  --type {Bushy,TStellate,DStellate}, -T {Bushy,TStellate,DStellate}
                        Define the cell type (default: Bushy)
  --model {XM13,RM03,XM13PasDend,Calyx,MNTB,L23Pyr}, -M {XM13,RM03,XM13PasDend,Calyx,MNTB,L23Pyr}
                        Define the model type (default: XM13)
  --sgcmodel {Zilany, cochlea}, -M {Zilany, cochlea}
  --protocol {initIV,testIV,runIV,initAN,runANPSTH,runANSingles}, -P {initIV,testIV,runIV,initAN,runANPSTH,runANSingles}
                        Protocol to use for simulation (default: IV)
  --inputpattern        Cell ID to use for input pattern if substituting inputs
  --hoc INFILE, -H INFILE
                        hoc file to use for simulation (default is the
                        selected "cell".hoc)
  -r NREPS, --reps NREPS
                        # repetitions
  --modf MODFREQ        Set SAM modulation frequency
  --moddepth MODDEPTH   Set SAM modulation depth (in percent)
  --S2M SIGNALTOMASKER  Signal to Masker ratio (dB)
  --cmmrmode {CM,CD,REF}
                        Specify mode (from: ['CM', 'CD', 'REF'])
  --allmodes            Force run of all modes (CMR, CMD, REF) for stimulus
                        configuration.
  -S {LS,MS,HS,fromcell}, --SRType {LS,MS,HS,fromcell}
                        Specify SR type (from: ['LS', 'MS', 'HS', 'fromcell'])
  --sequence SEQUENCE   Specify a sequence for the primary run parameters
  --plot                Plot results as they are generated - requires user
                        intervention...
  --workers NWORKERS    Number of "workers" for parallel processing (default:
                        4)

ex:

"""

import sys
import os.path
import os
import errno
import pickle
import time
import argparse
from collections import OrderedDict
import pprint
import json

#import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer

import neuronvis.hoc_graphics as hoc_graphics
from channel_decorate import ChannelDecorate
from generate_run import GenerateRun
import cellInitialization as cellInit

from cnmodel import cells
from cnmodel.util import sound

import cochlea
import thorns

import pylibrary.Utility as pu  # access to spike finder routine
import pyqtgraph.multiprocess as mproc


try:
    import pyqtgraph as pg
    #from pyqtgraph.Qt import QtGui
    from PyQt4 import QtCore, QtGui
    import pylibrary.pyqtgraphPlotHelpers as pgh
    HAVE_PG = True
except:
	HAVE_PG = False

# if HAVE_PG:
#     import render

import numpy as np

import cell_config

verbose = False
showCell = False




class ModelRun():
    def __init__(self, args=None):

        # use v2 files for model with rescaled soma
        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelChoices = ['XM13', 'RM03', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr']
        self.SGCmodelChoices = ['Zilany', 'cochlea']  # cochlea is python model of Zilany data, no matlab, JIT computation
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initAN', 'runANPSTH', 'runANSingles']
        self.soundChoices = ['tonepip', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        self.speciesChoices = ['mouse', 'guineapig']
        

        self.srname = ['**', 'LS', 'MS', 'HS']  # runs 1-3, not starting at 0
        self.cellID = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
        self.Params = OrderedDict()

        self.Params['initIVStateFile'] = 'IVneuronState.dat'
        self.Params['initANStateFile'] = 'ANneuronState.dat'
        self.Params['infile ']= None

        self.Params['cellType'] = self.cellChoices[0]
        self.Params['modelType'] = self.modelChoices[0]
        self.Params['SGCmodelType'] = self.SGCmodelChoices[0]
        self.Params['species'] = self.speciesChoices[0]
        self.Params['SRType'] = self.SRChoices[2]
        self.Params['SR'] = self.Params['SRType']  # actually used SR this might be cell-defined, rather than command defined
        self.Params['inputPattern'] = None # ID of cellinput pattern (same as cellID): for substitute input patterns.

        self.Params['runProtocol'] = self.protocolChoices[2]  # testIV is default because it is fast and should be run often
        self.Params['nReps'] = 1
        self.Params['run_duration'] = 0.8 # in sec
        self.Params['soundtype'] = 'SAM'  # or 'tonepip'
        self.Params['pip_duration'] = 0.1
        self.Params['pip_start'] = [0.1]
        self.Params['Fs'] = 100e3
        self.Params['F0'] = 4000.
        self.Params['dB'] = 40.
        self.Params['RF'] = 2.5e-3
        self.Params['fmod'] = 20 # hz, modulation if SAM
        self.Params['dmod'] = 0 # percent if SAM
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params['threshold'] = -20
        self.Params['plotFlag'] = False
        self.Params['auto_initialize'] = False
        self.Params['nWorkers'] = 4
        self.baseDirectory = 'VCN_Cells'
        self.morphDirectory = 'Morphology'
        self.initDirectory = 'Initialization'
        self.simDirectory = 'Simulations'


    def printModelSetup(self):
        """
        Print out all of the parameters in the model
        """
        for p in self.Params.keys():
            print ('{:18s} = {:12}'.format(p, self.Params[p]))
        print '-----------'


    def set_celltype(self, cellType):
        """
        Set the cell Type, as requested. The cell type must be in the cell choices
        
        Parameters
        ----------
        cellType : string
            The type of the cell that will be the basis for the model
        
        Returns
        -------
            Nothing
        
        """
        if cellType not in self.cellChoices:
            print 'Celltype must be one of: %s. Got: %s', (', '.join(self.cellChoices), cellType)
            exit()
        self.Params['cellType'] = cellType


    def set_modeltype(self, modelType):
        """
        Set the model Type, as requested. The model type must be in the model choices
        
        Parameters
        ----------
        modelType : string
            The type of the model that will be used (condutance settings)
        
        Returns
        -------
            Nothing
        
        """
        if modelType not in self.modelChoices:
            print 'Model type must be one of: %s. Got: %s ' % (', '.join(self.modelChoices), modelType)
            exit()
        self.Params['modelType'] = modelType


    def set_SR(self, SRType):
        """
        Set the SR, overriding SR in the cell_config file. The SR type must be in the SR choices
        
        Parameters
        ----------
        SRType : string
            The SR type that will be used for AN fibers.
        
        Returns
        -------
            Nothing
        
        """
        if self.SRType not in self.SRChoices:
            print 'SR type must be one of: %s. Got: %s ' % (', '.join(self.SRChoices), SRType)
            exit()
        self.Params['SRType'] = SRType

    def set_starttime(self, starttime):
        self.Params['StartTime'] = starttime

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def runModel(self, parMap={}):
        """
        Main entry routine for running all models
        
        Parameters
        ----------
        parMap : dict (default: empty)
            A dictionary of parameters, passed to models that are run (not used).
        
        Returns
        -------
            Nothing
        
        """
        if verbose:
            print 'runModel entry'
        if 'id' in parMap.keys():
            self.idnum = parMap['id']
        else:
            self.idnum = 9999
        self.cellID = os.path.splitext(self.Params['cell'])[0]

        if parMap == {}:
            self.plotFlag = True
        ivinitfile = os.path.join(self.baseDirectory, self.cellID, 
                            self.initDirectory, self.Params['initIVStateFile'])
        ivinitdir = os.path.join(self.baseDirectory, self.cellID, 
                            self.initDirectory)
        self.mkdir_p(ivinitdir) # confirm existence of that file
        filename = os.path.join(self.baseDirectory, self.cellID, self.morphDirectory, self.Params['infile'])
        
        if self.Params['cellType'] in ['Bushy', 'bushy']:
            print 'runModel creating a bushy cell: '
            self.postCell = cells.Bushy.create(morphology=filename, decorator=ChannelDecorate,
                    species=self.Params['species'], 
                    modelType=self.Params['modelType'], )
        elif self.Params['cellType'] in ['tstellate', 'TStellate']:
            print 'runModel creating a t-stellate cell: '
            self.postCell = cells.TStellate.create(morphology=filename, decorator=ChannelDecorate,
                    species=self.Params['species'], 
                    modelType=self.Params['modelType'], ) 
        elif self.Params['cellType'] in ['dstellate', 'DStellate']:
            print 'runModel creating a D-stellate cell: '
            self.postCell = cells.DStellate.create(morphology=filename, decorator=ChannelDecorate,
                    species=self.Params['species'], 
                    modelType=self.Params['modelType'], ) 
        else:
            raise ValueError("cell type %s not implemented" % (self.Params['cellType']))

                      
        # print 'PostCell: ', dir(self.postCell)
       #  print 'PostCell hr: ', dir(self.postCell.hr)

        self.postCell.hr.h.celsius = 38.
        self.postCell.hr.h.Ra = 150.
        print 'Ra is: ', self.postCell.hr.h.Ra
        print 'Temp is: ', self.postCell.hr.h.celsius
        for group in self.postCell.hr.sec_groups.keys():
            g = self.postCell.hr.sec_groups[group]
            for section in list(g):
                self.postCell.hr.get_section(section).Ra = self.postCell.hr.h.Ra
                if verbose:
                    print 'section: ', section
                    print 'Ra: ', self.postCell.hr.get_section(section).Ra

        electrodeSection = list(self.postCell.hr.sec_groups['soma'])[0]
        self.electrodeSite = self.postCell.hr.get_section(electrodeSection)
        self.electrodeSection = 'soma'
        self.hg = hoc_graphics
        self.get_hoc_file(self.postCell.hr)
        #sg = self.postCell.hr.sec_groups['soma']
        #self.distances(self.postCell.get_section(list(sg)[0]).name()) # make distance map from soma
#        print self.postCell.hr.distanceMap
        dendriticElectrode = list(self.postCell.hr.sec_groups['dendrite'])[0]
        self.dendriticElectrodeSite = self.postCell.hr.get_section(dendriticElectrode)
        self.dendriticElectrodeSection = 'dendrite'
        if verbose:
            print 'Parmap in runModel: ', parMap
            self.postCell.hr.h.topology()

        # handle the following protocols:
        # ['initIV', 'initAN', 'runIV', 'runANPSTH', 'runANSingles']
        
        if self.Params['runProtocol'] == 'initIV':
            if verbose:
                print 'initIV'
            ivinitfile = os.path.join(self.baseDirectory, self.cellID, 
                                self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(self.postCell.hr, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection, 
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=self.postCell.irange,
                             plotting = HAVE_PG and self.plotFlag)
            cellInit.getInitialConditionsState(self.postCell, tdur=3000., 
                filename=ivinitfile, electrodeSite=self.electrodeSite)
            print 'Ran to get initial state for %f msec' % self.postCell.hr.h.t
            return

        if self.Params['runProtocol'] == 'testIV':
            if verbose:
                print 'test_init'
            ivinitfile = os.path.join(self.baseDirectory, self.cellID, 
                                self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(self.postCell, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection, 
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=self.postCell.irange,
                             plotting = HAVE_PG and self.plotFlag, )
            cellInit.testInitialConditions(self.postCell, filename=ivinitfile,
                electrodeSite=self.electrodeSite)
            #self.R.testRun()
            return  # that is ALL, never make init and then keep running.

        if self.Params['runProtocol'] == 'runANPSTH':
            if verbose:
                print 'ANPSTH'
            self.ANRun(self.postCell)

        if self.Params['runProtocol'] == 'initAN':
            if verbose:
                print 'Init AN'
            self.ANRun(self.postCell, make_ANIntialConditions=True)

        if self.Params['runProtocol'] == 'runANSingles':
            if verbose:
                print 'ANSingles'
            self.ANRun_singles(self.postCell)
            
        if self.Params['runProtocol'] == 'runIV':
            if verbose:
                print 'iv_mode'
            self.IVRun(parMap)

        if showCell:
            self.render = HocViewer(self.postCell)
            cylinder=self.render.draw_cylinders()
            cylinder.set_group_colors(self.section_colors, alpha=0.8, mechanism=['nav11', 'gbar'])


    def IVRun(self, parMap={}):
        """
        Main entry routine for running all IV (current-voltage relationships with somatic electrode)
        
        Parameters
        ----------
        parMap : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).
        
        Returns
        -------
            summary : dict
                A summary of the results, including the file, parMap, resting input resistance,
                time constant, and spike times
        
        """
        print 'IvRun: IVRun begins'
        if verbose:
            print 'IVRun: calling generateRun'
        self.R = GenerateRun(self.postCell.hr, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection,
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=self.postCell.irange,
                             plotting = HAVE_PG and self.plotFlag)
        ivinitfile = os.path.join(self.baseDirectory, self.cellID, 
                                self.initDirectory, self.Params['initIVStateFile'])
        self.R.runInfo.folder = os.path.join('VCN_Cells', self.cellID, self.simDirectory, 'IV')
        if verbose:
            print 'IVRun: calling doRun'
        self.R.doRun(self.Params['infile'], parMap=parMap, save='monitor', restoreFromFile=True, initfile=ivinitfile,
            workers=self.Params['nWorkers'])
        if verbose:
            print '   doRun completed'
        isteps = self.R.IVResult['I']
        if verbose:
            print 'IVRun: Results summary: '
            for k, i in enumerate(self.R.IVResult['tauih'].keys()):
                print '   ih: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['tauih'][i]['tau'].value)
                print '           dV : %f' % self.R.IVResult['tauih'][i]['a'].value
            for k, i in enumerate(self.R.IVResult['taus'].keys()):
                print '   i: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['taus'][i]['tau'].value)
                print '          dV : %f' % (self.R.IVResult['taus'][i]['a'].value)

            print '   Nspike, Ispike: ', self.R.IVResult['Nspike'], self.R.IVResult['Ispike']
            print '   Rinss: ', self.R.IVResult['Rinss']
            print '   Vm: ', np.mean(self.R.IVResult['Vm'])
        if len(self.R.IVResult['taus'].keys()) == 0:
            taum_mean = 0.
            tauih_mean = 0.
        else:
            taum_mean = np.mean([self.R.IVResult['taus'][i]['tau'].value for k, i in
                enumerate(self.R.IVResult['taus'].keys())])
            tauih_mean = np.mean([self.R.IVResult['tauih'][i]['tau'].value for k, i in
                enumerate(self.R.IVResult['tauih'].keys())])
        # construct dictionary for return results:
        self.IVSummary = {'basefile': self.R.basename,
                          'parmap': parMap, 'ID': self.idnum,
                          'Vm': np.mean(self.R.IVResult['Vm']),
                          'Rin': self.R.IVResult['Rinss'],
                          'taum': taum_mean, 'tauih': tauih_mean,
                          'spikes': {'i': self.R.IVResult['Ispike'], 'n': self.R.IVResult['Nspike']},
                          }
#        print 'model_run::runModel::IVRun write summary for file set = ', self.R.basename
        return self.IVSummary

    def checkForANStateFile(self):
        print 'state file: ', self.Params['initANStateFile']
        print 'cellid: ', self.Params['cell']
        print 'base: ', self.baseDirectory
        print 'init: ', self.initDirectory
        statefile = os.path.join(self.baseDirectory, self.Params['cell'], 
                            self.initDirectory, self.Params['initANStateFile'])
        return(os.path.isfile(statefile))

    def ANRun(self, postCell, verify=False, seed=0, make_ANIntialConditions=False):
        """
        Establish AN inputs to soma, and run the model.
        Requires a synapseConfig list of dicts from cell_config.makeDict()
        each list element represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])

        Parameters
        ----------
        postCell : Cell object
            Provides access to cell class, as well as neuron and hoc file information
        
        verify : boolean (default: False)
            Flag to control printing of various intermediate results
        
        seed : int (default: 0)
            A random number seed to use when generating AN inputs (not used)
        
        make_ANIntitialConditions : bool (default : False)
            Flag to control whether the initial conditions need to be recomputed and stored.
            
        Returns
        -------
            Nothing
        

        """
        if self.Params['inputPattern'] is not None:
            fromdict = self.Params['inputPattern']
            print 'Cell id: %s  using input pattern: %s' % (self.cellID, fromdict)
        else:
            fromdict = self.cellID
            
        synapseConfig, celltype = cell_config.makeDict(fromdict)
        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        # for i, s in enumerate(synapseConfig):
        #     s[1] = s[3]*0.001/3.0
        #     print 'delay for input %d is %8.4f msec' % (i, s[1])

        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV

        stimInfo = self.Params
        # {'Morphology': self.Params['infile'], 'synapseConfig': synapseConfig,
        #             'runDur': self.run_duration, 'pip_dur': self.Params['pip_duration'], 'pip_start': self.Params['pip_start'],
        #             'run_duration': self.run_duration,
        #             'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
        #             'cellType': self.Params['cellType'], 'modelType': self.Params['modelType'], 'nReps': nReps, 'threshold': threshold}

        preCell, synapse, self.electrodeSite = self.configureCell(postCell, synapseConfig, celltype, stimInfo)

        # see if we need to save the cell state now.
        if make_ANIntialConditions:
            print 'getting initial conditions for AN'
            aninitfile = os.path.join(self.baseDirectory, self.cellID, 
                                self.initDirectory, self.Params['initANStateFile'])
            cellInit.getInitialConditionsState(postCell.hr, tdur=3000., 
                filename=aninitfile, electrodeSite=self.electrodeSite, reinit=self.Params['auto_initialize'])
            cellInit.testInitialConditions(postCell.hr, filename=aninitfile,
                electrodeSite=self.electrodeSite)
            return

        seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        print 'AN Seeds: ', seeds
        stimInfo['seeds'] = seeds  # keep the seed values too.
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        stimWaveform = {}
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.

        nWorkers = self.Params['nWorkers']
        TASKS = [s for s in range(nReps)]
        tresults = [None]*len(TASKS)

        # run using pyqtgraph's parallel support
        with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
            for j, x in tasker:
                tresults = self.singleANRun(postCell.hr, j, synapseConfig,
                    stimInfo, seeds, preCell, postCell, self.an_setup_time)
                tasker.results[j] = tresults
        # retreive the data
        for j, N in enumerate(range(nReps)):
           
            celltime.append(tresults[j]['time']) # (self.time)
            spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'],
                    threshold, t0=0., t1=stimInfo['run_duration']*1000., dt=1.0, mode='peak')
            inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
            somaVoltage[N] = np.array(tresults[j]['Vsoma'])
            dendriteVoltage[N] = np.array(tresults[j]['Vdend'])
            stimWaveform[N] = np.array(tresults[j]['stim']) # save the stimulus
            

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
            'somaVoltage': somaVoltage, 'dendriteVoltage': dendriteVoltage, 'stimWaveform': stimWaveform, 'time': np.array(celltime[0])}
        
        self.analysis_filewriter(self.Params['cell'], result, tag='delays')
        self.plotAN(np.array(celltime[0]), result['somaVoltage'], result['stimInfo'],
            dendVoltage=result['dendriteVoltage'])


    def ANRun_singles(self, postCell, verify=False, seed=None):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        This routine is special - it runs nReps for each synapse, turning off all of the other synapses
        by setting the synaptic conductance to 0 (to avoid upsetting initialization)

        Parameters
        ----------
        hf : hoc_reader object
            Access to neuron and file information
        
        verify : boolean (default: False)
            Flag to control printing of various intermediate results
        
        seed : int (default: None)
            A random number seed to use when generating AN inputs.
            Allows different runs with different spike trains, but also
            to recreate a particular run.

        Returns
        -------
            Nothing
        
        """

        self.start_time = time.time()
        synapseConfig, celltype = cell_config.makeDict(self.cellID)
        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV

        stimInfo = self.Params
        # {'Morphology': self.Params['infile'], 'synapseConfig': synapseConfig,
        #             'runDur': self.run_duration, 'pip_dur': self.Params['pip_duration'], 'pip_start': self.Params['pip_start'],
        #             'run_duration': self.run_duration,
        #             'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
        #             'cellType': self.Params['cellType'], 'modelType': self.Params['modelType'], 'nReps': nReps, 'threshold': threshold}

        preCell, synapse, self.electrodeSite = self.configureCell(postCell, synapseConfig, celltype, stimInfo)


        nSyns = len(synapseConfig)
        if seed is None:
            seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        else:
            seeds = np.random.randint()
        print 'AN Seeds: ', seeds
        stimInfo['seeds'] = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = True
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

            tresults = [None]*nReps

            if parallel:
                nWorkers = self.Params['nWorkers']
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
                    for j, x in tasker:
                        tresults = self.singleANRun(postCell.hr, j, synapseConfig,
                            stimInfo, seeds, preCell, postCell, self.an_setup_time)
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    tresults[j] = self.singleANRun(postCell.hr, j, synapseConfig,
                            stimInfo, seeds, preCell, postCell, self.an_setup_time)
            for j, N in enumerate(range(nReps)):
#                print 'Rep: %d' % N

#                preCell, postCell = self.singleANRun(hf, j, synapseConfig, stimInfo, seeds, preCell, postCell)
                
                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'], threshold, t0=0.,
                        t1=stimInfo['run_duration']*1000, dt=1.0, mode='peak')
                inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]['Vsoma'])
                dendriteVoltage[N] = np.array(tresults[j]['Vdend'])

            total_elapsed_time = time.time() - self.start_time
    #        total_run_time = time.time() - run_time
            print "Total Elapsed Time = %8.2f min (%8.0fs)" % (total_elapsed_time/60., total_elapsed_time)
            print "Total Setup Time = %8.2f min (%8.0fs)" % (self.setup_time/60., self.setup_time)
            print "Total AN Calculation Time = %8.2f min (%8.0fs)" % (self.an_setup_time/60., self.an_setup_time)
            print "Total Neuron Run Time = %8.2f min (%8.0fs)" % (self.nrn_run_time/60., self.nrn_run_time)
        
            result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes, 
                'somaVoltage': somaVoltage, 'dendriteVoltage': dendriteVoltage, 'time': np.array(tresults[j]['time'])}
        
            self.analysis_filewriter(self.Params['cell'], result, tag='Syn%03d' % k)
        self.plotAN(np.array(result['time']), result['somaVoltage'], result['stimInfo'])

    def configureCell(self, thisCell, synapseConfig, celltype, stimInfo):
        """
        Configure the cell. This routine builds the cell in Neuron, adds presynaptic inputs
        as described in the synapseConfig, and configures those according to parameters in 
        stiminfo.
        
        Parameters
        ----------
        thisCell : hoc_reader object
            Access to neuron and file information
        
        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.
        
        celltype : string
            A string describing the cell type. Determines how the channels are populated
            and how the cell is built, even if it comes from a .hoc structure.
            
        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)

        Returns
        -------
            preCell : list
                A list of the preCell hoc objects attached to the synapses
            postCell : cells object
                The target cell
            synapse : list
                A list of the synapses that were created and connected to this cell
            electrodeSite : Neuron hoc object
                The section and location where the electrode is attached to the cell.
        
        
        """
        
        debug = False
        if debug:
            print 'hf.sec_groups : ', thisCell.sec_groups.keys()
        #self.createCell(hf, celltype, self.Params['modelType'], self.Params['species'])
        if debug:
            print thisCell.print_all_mechs()


        preCell = []
        synapse = []
        # reconfigure syanpses to set the spont rate group
        stimInfo['SR'] = stimInfo['SRType']
        for i, syn in enumerate(synapseConfig):
            # print 'SYN: ', syn
            if stimInfo['SRType'] is 'fromcell':  # use the one in the table
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=syn['SR']))
                stimInfo['SR'] = self.srname[syn[2]] # use and report value from table
            else:
                try:
                    srindex = self.srname.index(stimInfo['SRType'])
                    print 'retrieved index %d with SR type %s' % (srindex, stimInfo['SRType'])
                except:
                    raise ValueError('SR type "%s" not found in Sr type list' % stimInfo['SRType'])
                    
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=srindex))  # override
            # print 'precell: ', preCell
            synapse.append(preCell[-1].connect(thisCell, pre_opts={'nzones':syn['nSyn'], 'delay':syn['delay2']}))
        for i, s in enumerate(synapse):
            s.terminal.relsite.Dep_Flag = 0  # turn off depression computation
            #print dir(s.psd.ampa_psd[0])

        #  ****** uncomment here to adjust gmax.

            for p in s.psd.ampa_psd:
                p.gmax = p.gmax*2.5

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
        electrodeSection = list(thisCell.hr.sec_groups['soma'])[0]
        electrodeSite = thisCell.hr.get_section(electrodeSection)
        return (preCell, synapse, electrodeSite)

    def set_dbspl(self, signal, dbspl):
        """Scale the level of `signal` to the given dB_SPL."""
        p0 = 20e-6
        rms = np.sqrt(np.sum(signal**2) / signal.size)

        scaled = signal * 10**(dbspl / 20.0) * p0 / rms

        return scaled


    def singleANRun(self, hf, j, synapseConfig, stimInfo, seeds, preCell, postCell, an_setup_time):
        """
        Perform a single run with all AN input on the target cell turned off except for input j.
        
        Parameters
        ----------
        hf : hoc_reader object
            Access to neuron and file information
        
        j : int
            The input that will be active in this run
             
        synapseConfig : dict
            A dictionary with information about the synapse configuration to use. 
            
        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)
        
        preCell : list
            A list of the preCell hoc objects attached to the synapses
        
        postCell : cells object
            The target cell
        
        an_setup_time : time object
            
        Returns
        -------
        anresult : dict
            A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'
        
        """
        filename = os.path.join(self.baseDirectory, self.cellID,
                    self.initDirectory, self.Params['initANStateFile'])
        cellInit.restoreInitialConditionsState(hf, electrodeSite=None, filename=filename)
        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        if isinstance(stimInfo['pip_start'], list):
            pips = stimInfo['pip_start']
        else:
            pips = [stimInfo['pip_start']]
        if stimInfo['soundtype'] == 'tonepip':
            stim = sound.TonePip(rate=stimInfo['Fs'], duration=stimInfo['run_duration'], 
                              f0=stimInfo['F0'], dbspl=stimInfo['dB'],
                              ramp_duration=stimInfo['RF'], pip_duration=stimInfo['pip_duration'],
                              pip_start=pips)
        elif stimInfo['soundtype'] == 'SAM':
            stim = sound.SAMTone(rate=stimInfo['Fs'], duration=stimInfo['run_duration'], 
                              f0=stimInfo['F0'], dbspl=stimInfo['dB'],
                              ramp_duration=stimInfo['RF'], fmod=stimInfo['fmod'],
                                   dmod=stimInfo['dmod'],
                              pip_duration=stimInfo['pip_duration'],
                              pip_start=pips)
        else:
            raise ValueError('StimInfo sound type %s not implemented' % stimInfo['soundtype'])

        for i, syn in enumerate(synapseConfig):
            nseed = seeds[j, i]

            if self.Params['SGCmodelType'] in ['Zilany']:
#                print 'running with Zilany'
                preCell[i].set_sound_stim(stim, seed=nseed, simulator='matlab')  # generate spike train, connect to terminal
            
            elif self.Params['SGCmodelType'] in ['cochlea']:
#                print 'running with cochlea, j=%d' % j
                wf = self.set_dbspl(stim.generate(), stimInfo['dB'])
                stim._sound = wf
                preCell[i].set_sound_stim(stim, seed=nseed, simulator='cochlea')  # generate spike train, connect to terminal

            else:
                raise ValueError('SGC model type type %s not implemented' % self.Params['SGCmodelType'])
                
            ANSpikeTimes.append(preCell[i]._spiketrain)
        an_setup_time += (time.time() - an0_time)
        nrn_start = time.time()
        Vsoma = hf.h.Vector()
        Vdend = hf.h.Vector()
        rtime = hf.h.Vector()
        dendsite = postCell.all_sections['dendrite'][-1]
        
        Vsoma.record(postCell.soma(0.5)._ref_v, sec=postCell.soma)
        Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        rtime.record(hf.h._ref_t)
#        print '...running '
        hf.h.finitialize()
        hf.h.tstop = stimInfo['run_duration']*1000.
        hf.h.t = 0.
        hf.h.batch_save() # save nothing
        hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
        # old - don't di this! #hf.h.finitialize()
        # hf.h.run()
        nrn_run_time += (time.time() - nrn_start)
#        print '...done'
        anresult = {'Vsoma': np.array(Vsoma), 'Vdend': np.array(Vdend), 'time': np.array(rtime), 'ANSpikeTimes': ANSpikeTimes, 'stim': stim}
        return anresult


    def plotAN(self, celltime, somaVoltage, stimInfo, dendVoltage=None):
        """
        """
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
        if dendVoltage is not None:
            for j, N in enumerate(range(len(dendVoltage))):
                layout.plot(0, celltime, dendVoltage[N], pen=pg.mkPen(pg.intColor(N, nReps)),
                    style=QtCore.Qt.DashLine)            
        pgh.show()

    def analysis_filewriter(self, filebase, result, tag=''):
        k = result.keys()
        requiredKeys = ['stimInfo', 'spikeTimes', 'inputSpikeTimes', 'somaVoltage', 'time']
        for rk in requiredKeys:
            assert rk in k

        stimInfo = result['stimInfo']
        outPath = os.path.join('VCN_Cells', self.cellID, self.simDirectory, 'AN')
        self.mkdir_p(outPath) # confirm that output path exists
        ID = self.cellID
        if self.Params['inputPattern'] is not None:
            ID += '_%s' % self.Params['inputPattern']
        if stimInfo['soundtype'] in ['SAM', 'sam']:
            f = open(os.path.join(outPath, 'AN_Result_' + ID + '_%s_N%03d_%03ddB_%06.1f_FM%03.1f_DM%03d_%2s' %
                (tag, stimInfo['nReps'],
                int(stimInfo['dB']), stimInfo['F0'],
                stimInfo['fmod'], int(stimInfo['dmod']), stimInfo['SR']) + '.p'), 'w')
        else:
            f = open(os.path.join(outPath, 'AN_Result_' + ID + '_%s_N%03d_%03ddB_%06.1f_%2s' % (tag, stimInfo['nReps'],
                int(stimInfo['dB']), stimInfo['F0'], stimInfo['SR']) + '.p'), 'w')
        pickle.dump(result, f)
        f.close()
        print '**** Analysis wrote output file ****\n    %s' % f        

    def get_hoc_file(self, hf):
        if hf.file_loaded is False:
            exit()
        self.section_list = hf.get_section_prefixes()
#        print 'section groups: ', self.hf.sec_groups.keys()
        hf.sec_groups.keys()
        if len(hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.colorMap[i]
#        else: # single section name, assign colors to SectionList types:
#        self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
#            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

        (v, e) = hf.get_geometry()
        self.clist = []

        for si in hf.sections: # self.section_list[s]:
            hf.h('access %s' % si)
            sr = hf.h.SectionRef()
            n1 = hf.h.cas().name()
            if sr.has_parent() == 1:
                x=sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])


if __name__ == "__main__":
    curdir = os.getcwd()
    model = ModelRun(sys.argv[1:])  # create instance of the model
    pp = pprint.PrettyPrinter(indent=4, width=60)
    
    parser = argparse.ArgumentParser(description='Simulate activity in a reconstructed model cell')
    parser.add_argument(dest='cell', action='store',
                   default=None, 
                   help='Select the cell (no default)')
    parser.add_argument('--type', '-T', dest='cellType', action='store',
                   default='Bushy', choices=model.cellChoices,
                   help='Define the cell type (default: Bushy)')
    parser.add_argument('--model', '-M', dest='modelType', action='store',
                   default='XM13', choices=model.modelChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--sgcmodel', dest='SGCmodelType', action='store',
                   default='Zilany', choices=model.SGCmodelChoices,
                   help='Define the SGC model type (default: Zilany)')
    parser.add_argument('--protocol', '-P', dest='runProtocol', action='store',
                   default='IV', choices=model.protocolChoices,
                   help='Protocol to use for simulation (default: IV)')
    parser.add_argument('--hoc', '-H', dest='infile', action='store',
                  default=None,
                  help='hoc file to use for simulation (default is the selected "cell".hoc)')
    parser.add_argument('--inputpattern', dest='inputPattern', action='store',
                  default=None,
                  help='cell input pattern to use (substitute) from cell_config.py')
    parser.add_argument('--stimulus', dest='soundtype', action='store',
                   default='tonepip', choices=model.soundChoices,
                   help='Define the stimulus type (default: tonepip)')

                  
    # lowercase options are generally parameter settings:

    parser.add_argument('-r', '--reps', type=int, default=1, dest = 'nReps',
        help='# repetitions')
    parser.add_argument('--fmod', type=float, default=20, dest = 'fmod', 
        help='Set SAM modulation frequency')
    parser.add_argument('--depth', type=float, default=100., dest = 'dmod', 
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
    parser.add_argument('--workers', type=int,  default=4, dest = 'nWorkers', 
            help='Number of "workers" for parallel processing (default: 4)')
    parser.add_argument('--auto-intialize', action="store_true", default=False, dest='auto_initialize',
            help='Force auto initialization if reading the state fails in initialization')

    # parser.add_argument('-p', '--print', action="store_true", default=False, dest = 'print_info',
    #     help='Print extra information during analysis')
    #  parser.add_argument('-l', '--list', action="store_true", default=False, dest = 'list_results',
    #     help='List results to screen')


    args = vars(parser.parse_args())   
#    model.printModelSetup()

    for k in args.keys():
        model.Params[k] = args[k]
    print model.Params['cell']
    if model.Params['infile'] == None: # just use the hoc file
        model.Params['infile'] = model.Params['cell'] + '.hoc'
    print(json.dumps(model.Params, indent=4))  # pprint doesn't work well with ordered dicts
    
    model.runModel() # then run the model
    if showCell:
        QtGui.QApplication.instance().exec_()
    try:
        exit()
    except:
        pass
