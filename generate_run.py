__author__ = 'pbmanis'
"""
GenerateRun is a class that sets up for a run after the cell has been decorated with channels.
It requires the celltype, and a section where the electrode will be inserted.
The stimulus in current clamp can consist of single pulses or pulse trains
Code for reading externally generated spike trains from a file is also included.
Methods:
    doRun(filename) will execute the run. The resulting plot will have the filename text at the top.
    doRun calls 3 private methods, _prepareRun, _initializeRun, and _executeRun, which respectively,
    set up the stimuli, initialize the state of the model, and then run the model, generating a plot
February 2014, Paul B. Manis UNC Chapel Hill
"""

import os
import numpy as np
import nrnlibrary.makestim as makestim
from pylibrary.Params import Params
import CalyxPlots as cp
import time
import csv
import pickle
import pyqtgraph as pg

verbose = False # use this for testing.

class GenerateRun():
    def __init__(self, hf, celltype=None, electrodeSection = 'soma[0]'):
        self.run_initialized = False
        self.hf = hf # get the reader structure and the hoc pointer object locally
        # use the Params class to hold program run information and states
        # runInfo holds information that is specific the the run - stimuli, conditions, etc
        # DO NOT add .p to filename...
        self.runInfo = Params(fileName='Normal', runName='Run', manipulation="Canonical", folder="Simulations",
                              preMode="cc",
                              postMode='cc',
                              TargetCellType=celltype, # valid are "Bushy", "Stellate", "MNTB"
                              electrodeSection=electrodeSection,
                              celsius=37,  # set the temperature.
                              nStim=1,
                              stimFreq=200.,  # hz
                              stimInj=2.0,  # nA
                              stimDur=10.0,  # msec
                              stimDelay=2.0,  # msec
                              stimPost=3.0,  # msec
                              vnStim=1,
                              vstimFreq=200.,  # hz
                              vstimInj=50,  # mV amplitude of step (from holding)
                              vstimDur=50.0,  # msec
                              vstimDelay=2.0,  # msec
                              vstimPost=3.0,  # msec
                              vstimHolding=-60, # holding, mV
                              runTime=time.asctime(),  # store date and time of run
                              inFile=None,
                              # 'ANFiles/AN10000Hz.txt', # if this is not None, then we will use these spike times...
                              inFileRep=1,  # which rep to use (or array of reps)
                              spikeTimeList={},  # Dictionary of spike times
                              v_init = -63.9,
        )

        if self.runInfo.inFile is None:
            self.runInfo.tstop = 1000.0 * (
                self.runInfo.nStim - 1) / self.runInfo.stimFreq + self.runInfo.stimPost + self.runInfo.stimDelay

        else:
            print 'reading spike train'
            maxt = 0.
            with open(self.runInfo.inFile, 'r') as csvfile:
                spks = csv.reader(csvfile, delimiter=',')
                for i, row in enumerate(spks):
                    if i == 0:  # ture first line
                        print row
                        maxt = float(row[1]) * 1000
                        reps = int(row[0])
                        continue
                    if int(row[0]) in self.runInfo.spikeTimeList.keys():
                        self.runInfo.spikeTimeList[int(row[0])].append(float(row[1]) * 1000.)
                    else:
                        self.runInfo.spikeTimeList[int(row[0])] = [float(row[1]) * 1000.]
                self.runInfo.tstop = maxt

        self.monitor = {} # standard monitoring
        self.allsecVec = {} # all section monitoring
        self.mons = {}
        self.filename=None

        if verbose:
            print 'runinfo initialization done'


    def doRun(self, filename=None):
        self.filename = filename
        self.hf.update() # make sure channels are all up to date
        self._prepareRun() # build the recording arrays
        print 'run preparae'
        self._initRun()  # this also sets nseg in the axon - so do it before setting up shapes
        print 'run initialized'
        self._executeRun() # now you can do the run

    def _prepareRun(self):
        """
        (private method)
        Control a single run of the model with updated display of the voltages, etc.
        Inputs: None
        Outputs: None
        Actions: optionally displays the results
        Side Effects: A number of class variables are created and modified, mostly related to the
        generation of stimuli and monitoring of voltages and currents
        """
        for var in self.hf.sections: # get morphological components
            self.allsecVec[var] = self.hf.h.Vector()
            self.allsecVec[var].record(self.hf.sections[var](0.5)._ref_v, sec=self.hf.sections[var])

        for var in ['time', 'postsynapticV', 'postsynapticI', 'i_stim0', 'v_stim0']: # get standard stuff
            self.monitor[var] = self.hf.h.Vector()

        self.hf.h.celsius = self.runInfo.celsius

        self.clist = {}  #  color list

        # make fake cell right here for testing...
        # soma = self.hf.h.Section()
        # soma.L = 20.
        # soma.diam = 20.
        # soma.insert('hh')
        #
        #electrodeSite = soma
        print self.runInfo.electrodeSection
        self.electrodeSite = self.hf.sections[self.runInfo.electrodeSection]
        if self.runInfo.postMode in ['vc', 'vclamp']:
            print 'vclamp'
            # Note to self (so to speak): the hoc object returned by this call must have a life after
            # # the routine exits. Thus, it must be "self." Same for the IC stimulus...
            self.vcPost = self.hf.h.SEClamp(0.5, sec=self.electrodeSite) #self.hf.sections[electrodeSite])
            self.vcPost.dur1 = 2
            self.vcPost.amp1 = self.runInfo.vstimHolding
            self.vcPost.dur2 = 1e9
            self.vcPost.amp2 = self.runInfo.vstimHolding  # just a tiny step to keep the system honest
            self.vcPost.dur3 = 10.0
            self.vcPost.amp3 = self.runInfo.vstimHolding
            self.vcPost.rs = 1e-6
            stim = {}
            stim['NP'] = self.runInfo.vnStim
            stim['Sfreq'] = self.runInfo.vstimFreq  # stimulus frequency
            stim['delay'] = self.runInfo.vstimDelay
            stim['dur'] = self.runInfo.vstimDur
            stim['amp'] = self.runInfo.vstimInj
            stim['PT'] = 0.0
           # print self.hf.h.soma[0]
            (secmd, maxt, tstims) = makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            secmd = secmd + self.runInfo.vstimHolding # add holding
            self.monitor['v_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['v_stim0'].play(self.vcPost._ref_amp2, self.hf.h.dt, 0, sec=self.electrodeSite)
            self.monitor['postsynapticV'].record(self.electrodeSite(0.5)._ref_v, sec=self.electrodeSite)
            self.monitor['postsynapticI'].record(self.vcPost._ref_i, sec=self.electrodeSite)
            self.mons = ['postsynapticI', 'v_stim0']
        elif self.runInfo.postMode in ['cc', 'iclamp']:
            print 'iclamp'
            stim = {}
            stim['NP'] = self.runInfo.nStim
            stim['Sfreq'] = self.runInfo.stimFreq  # stimulus frequency
            stim['delay'] = self.runInfo.stimDelay
            stim['dur'] = self.runInfo.stimDur
            stim['amp'] = self.runInfo.stimInj
            stim['PT'] = 0.0
            (secmd, maxt, tstims) = makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrodeSite)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor['i_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['i_stim0'].play(self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrodeSite)
            self.monitor['postsynapticI'].record(self.icPost._ref_i, sec=self.electrodeSite)
            self.monitor['postsynapticV'].record(self.electrodeSite(0.5)._ref_v, sec=self.electrodeSite)
            self.mons = ['postsynapticV', 'postsynapticI' ]
        else:
            print 'generate_run.py, mode %s  unknown' % self.runInfo.postMode
            return
        self.hf.h.tstop = maxt
        self.monitor['time'].record(self.hf.h._ref_t)


    def testRun(self, title='testing...'):
        self.hf.h.tstop = 10
        self.hf.h.finitialize()
        self.hf.h.run()
        pg.mkQApp()
        pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
        pl.setTitle(title)


    def _initRun(self):
        """
        (private method)
        Model initialization procedure:
        Set RMP to the resting RMP of the model cell.
        Make sure nseg is large enough in the proximal axon.
        Initialize the leak, and stabilize the pumps, etc.
        Outputs: None.
        Action: Initializes the NEURON state to begin a run.
        Does not instantiate recording or stimulating.
        """
        if self.runInfo is None:
            raise Exception('GenerateRun: initRun has no runInfo')

        if self.runInfo.postMode in ['vc', 'vclamp']:
            self.hf.h.finitialize(self.runInfo.vstimHolding)
            self.run_initialized = True
            return

        # First we set e_leak so that the rmp in each segment is the same
        self.hf.h.finitialize(self.runInfo.v_init)
        # starting way back in time
        self.hf.h.t = -1e10
        dtsav = self.hf.h.dt
        self.hf.h.dt = 1e9  # big time steps for slow process
        temp = self.hf.h.cvode.active()
        if (temp != 0):
            self.hf.h.cvode.active(0)  # turn cvode off (note, in this model it will be off because one of the mechanisms is not compatible with cvode at this time
        while (self.hf.h.t < -1e9):
            self.hf.h.fadvance()

        if (temp != 0):
            self.hf.h.cvode.active(1)
        self.hf.h.dt = dtsav
        self.hf.h.t = 0
        if (self.hf.h.cvode.active()):
            self.hf.h.cvode.re_init()
        else:
           self.hf.h.fcurrent()
        self.hf.h.frecord_init()
        self.hf.h.finitialize(self.hf.h.v_init)
        self.run_initialized = True


    def _executeRun(self, testPlot=False):
        """
        (private mmethod)
        After prepare run and initialization, this routine actually calls the run method in hoc
        assembles the data, saves it to disk and plots the results.
        Inputs: flag to put up a test plot....
        """
        assert self.run_initialized == True
        print 'executeRun: Running for: ', self.hf.h.tstop
        print 'V = : ', self.electrodeSite.v
        self.hf.h.run()
        if testPlot:
            pg.mkQApp()
            pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
            if self.filename is not None:
                pl.setTitle('%s' % self.filename)
            else:
                pl.setTitle('executeRun, no filename')
        print 'run done'
        np_monitor = {}
        for k in self.monitor.keys():
            np_monitor[k] = np.array(self.monitor[k])

        np_allsecVec = {}
        for k in self.allsecVec.keys():
            np_allsecVec[k] = np.array(self.allsecVec[k])
        self.runInfo.clist = self.clist
        results = Params(Sections=self.hf.sections.keys(), vec=np_allsecVec,
                         monitor=np_monitor,
                         distanceMap = self.hf.distanceMap,
        )
        print("Run done\n")
        self.plotRun(results)


    def saveRun(self, results):
        """
        Save the results to disk. Results must be a Param structure, which we turn into
         a dictionary...
        """
        dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        if os.path.exists(self.runInfo.folder) is False:
            os.mkdir(self.runInfo.folder)
        fn = os.path.join(self.runInfo.folder, self.runInfo.fileName + dtime + '.p')
        pfout = open(fn, 'wb')
        pickle.dump({'runInfo': self.runInfo.todict(),
                     'modelPars': [],
                     'Results': results.todict()}, pfout)
        pfout.close()

        # if recordSection:
        #     fns = os.path.join(self.runInfo.folder, self.runInfo.fileName + '_swellings_' + dtime + '.p')
        #     srsave = sr.SimulationResult()
        #     srsave.save(fns, data=np.array(self.vswel),
        #             time=np.array(self.vec['time']),
        #             hoc_file = 'MorphologyFiles/' + self.modelPars.topofile,
        #             section_map=secarray,
        #     )
            # pfout = open(fns, 'wb')
            # pickle.dump({'swellings': secarray,
            #              'time': np.array(self.vec['time']),
            #              'data': np.array(self.vswel)}, pfout)
            # pfout.close()
            # plot = pg.plot()
            # vs = np.array(self.vswel)
            # ts = np.array(self.vec['time'])
            # for i in range(len(self.swellAxonMap)):
            #     plot.plot(ts, vs[i])


    def plotRun(self, results):
        cplts = cp.CalyxPlots(title=self.filename)
        cplts.plotResults(results.todict(), self.runInfo.todict(), somasite=self.mons)
        cplts.show()
