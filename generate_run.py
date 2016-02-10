__author__ = 'pbmanis'
"""
GenerateRun is a class that sets up for a run after the cell has been decorated with channels.
It requires the celltype, and a section where the electrode will be inserted.
The stimulus in current clamp can consist of single pulses or pulse trains. cd is the
channelDecorator, which is also used to set the current range level for IV's.
Code for reading externally generated spike trains from a file is also included.
Methods:
    doRun(filename) will execute the run. The resulting plot will have the filename text at the top.
    doRun calls 3 private methods, _prepareRun, _initializeRun, and _executeRun, which respectively,
    set up the stimuli, initialize the state of the model, and then run the model, generating a plot
February 2014, Paul B. Manis UNC Chapel Hill
"""

import os
import numpy as np
import cnmodel.makestim
from pylibrary.Params import Params
import calyxPlots as cp
import analyze_run as ar
import cellInitialization as cellInit
import time
import csv
import pickle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import pyqtgraph.multiprocess as mproc
import errno



verbose = False # use this for testing.


class GenerateRun():
    def __init__(self, hf, idnum=0, celltype=None,
                 electrodeSection=None,
                 dendriticElectrodeSection=None,
                 starttime=None,
                 iRange=[-1, 0, 1],
                 plotting=False,
                 saveAllSections=False,
                 useSavedState=True):

        self.run_initialized = False
        self.plotting = plotting

        self.hf = hf # get the reader structure and the hoc pointer object locally
        self.basename = 'basenamenotset'
        self.idnum = idnum
        self.startTime = starttime
        self.saveAllSections = saveAllSections
        # use the Params class to hold program run information and states
        # runInfo holds information that is specific the the run - stimuli, conditions, etc
        # DO NOT add .p to filename...
        self.runInfo = Params(folder="Simulations", fileName='Normal', runName='Run',
                              manipulation="Canonical",
                              preMode="cc",
                              postMode='cc',
                              TargetCellType=celltype, # valid are "Bushy", "Stellate", "MNTB"
                              electrodeSection=electrodeSection,
                              dendriticElectrodeSection=dendriticElectrodeSection,
                              dendriticSectionDistance=100., # microns.
                              celsius=37,  # set the temperature.
                              nStim=1,
                              stimFreq=200.,  # hz
                              stimInj=iRange,  # nA, a list of test levels for current clamp
                              stimDur=100.0,  # msec
                              stimDelay=5.0,  # msec
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
                              v_init = -61.0,  # from Rothman type II model - not appropriate in all cases
                              useSaveState = useSavedState,  # use the saved state.
        )

        electrodeSection = list(self.hf.sec_groups[self.runInfo.electrodeSection])[0]
        self.electrodeSite = self.hf.get_section(electrodeSection)
        dend_sections = list(self.hf.sec_groups[self.runInfo.dendriticElectrodeSection])
        # invert the mappint of sections
        # first, just get the dendrite components, making a dendrite submap
        dendDistMap = {}
        for k in dend_sections:
            dendDistMap[k] = self.hf.distanceMap[k]
        revmap = dict((v, k) for k, v in dendDistMap.iteritems())
            
        # now find the distlist key that corresponds to the closest value to our desired value
        (num, dist) = min(enumerate(revmap), key=lambda x: abs(x[1]-self.runInfo.dendriticSectionDistance))
        print('Monitoring dendrite section number {:d}, at {:6.2f} microns from soma: '.format(num, dist))
        
        
        dendriticElectrodeSection = list(self.hf.sec_groups[self.runInfo.dendriticElectrodeSection])[num]
        self.dendriticElectrodeSite = self.hf.get_section(dendriticElectrodeSection)

        if self.runInfo.inFile is None:
            self.runInfo.tstop = 1000.0 * (
                self.runInfo.nStim - 1) / self.runInfo.stimFreq + self.runInfo.stimPost + self.runInfo.stimDelay

        else:
            print 'reading spike train'
            maxt = 0.
            with open(self.runInfo.inFile, 'r') as csvfile:
                spks = csv.reader(csvfile, delimiter=',')
                for i, row in enumerate(spks):
                    if i == 0:  # first line
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

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def makeFileName(self, filename=None, subdir=None):
        self.runInfo.filename = filename
        if self.startTime is None:
            self.dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        else:  # convert time object
            self.dtime = dtime = time.strftime("%y.%m.%d-%H.%M.%S", time.localtime(self.startTime))
        self.mkdir_p(self.runInfo.folder)
        #f os.path.exists(self.runInfo.folder) is False:
        #    os.mkdir(self.runInfo.folder)
        if subdir is not None:
            folder = os.path.join(self.runInfo.folder, subdir)
        else:
            folder = self.runInfo.folder
        if os.path.exists(folder) is False:
            os.mkdir(folder)
        #self.basename = os.path.join(folder, filename + self.dtime)
        self.basename = os.path.join(folder, filename)

    def _prepareRun(self, inj=None):
        """
        (private method)
        Control a single run of the model with updated display of the voltages, etc.
        Inputs: inj: override for current injection
        Outputs: None
        Actions: optionally displays the results
        Side Effects: A number of class variables are created and modified, mostly related to the
        generation of stimuli and monitoring of voltages and currents
        """
        if verbose:
            print '_prepareRun'
        for group in self.hf.sec_groups.keys(): # get morphological components
            if not self.saveAllSections:  # just save soma sections
                if group.rsplit('[')[0] == 'soma':
                    self.allsecVec['soma'] = self.hf.h.Vector()
                    section = list(self.hf.sec_groups[group])[0]
                    sec = self.hf.get_section(section)
                    self.allsecVec['soma'].record(sec(0.5)._ref_v, sec=sec)
                    break  # only save the FIRST occurance.
            else:
                g = self.hf.sec_groups[group]
                for section in list(g):
                    sec = self.hf.get_section(section)
                    self.allsecVec[sec.name()] = self.hf.h.Vector()
                    self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)

        for var in ['time', 'postsynapticV', 'dendriteV', 'postsynapticI', 'i_stim0', 'v_stim0']: # get standard stuff
            self.monitor[var] = self.hf.h.Vector()

        self.hf.h.celsius = self.runInfo.celsius

        self.clist = {}  #  color list

        electrodeSection = list(self.hf.sec_groups[self.runInfo.electrodeSection])[0]
        self.electrodeSite = self.hf.get_section(electrodeSection)
        if self.runInfo.postMode in ['vc', 'vclamp']:
            # Note to self (so to speak): the hoc object returned by this call must have a life after
            # the routine exits. Thus, it must be "self." Same for the IC stimulus...
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
            (secmd, maxt, tstims) = cnmodel.makestim.makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.stim = stim
            secmd = secmd + self.runInfo.vstimHolding # add holding
            self.monitor['v_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['v_stim0'].play(self.vcPost._ref_amp2, self.hf.h.dt, 0, sec=self.electrodeSite)
            self.monitor['postsynapticV'].record(self.electrodeSite(0.5)._ref_v, sec=self.electrodeSite)
            self.monitor['dendriteV'].record(self.dendriticElectrodeSite(0.5)._ref_v, sec=self.dendriticElectrodeSite)
            self.monitor['postsynapticI'].record(self.vcPost._ref_i, sec=self.electrodeSite)
            self.mons = ['postsynapticI', 'v_stim0']
        elif self.runInfo.postMode in ['cc', 'iclamp']:
            stim = {}
            stim['NP'] = self.runInfo.nStim
            stim['Sfreq'] = self.runInfo.stimFreq  # stimulus frequency
            stim['delay'] = self.runInfo.stimDelay
            stim['dur'] = self.runInfo.stimDur
            if inj is not None:
                stim['amp'] = inj
            else:
                stim['amp'] = self.runInfo.stimInj[0]
            stim['PT'] = 0.0
            (secmd, maxt, tstims) = cnmodel.makestim.makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.stim = stim
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrodeSite)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor['i_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['i_stim0'].play(self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrodeSite)
            self.monitor['postsynapticI'].record(self.icPost._ref_i, sec=self.electrodeSite)
            self.monitor['postsynapticV'].record(self.electrodeSite(0.5)._ref_v, sec=self.electrodeSite)
            self.monitor['dendriteV'].record(self.dendriticElectrodeSite(0.5)._ref_v, sec=self.dendriticElectrodeSite)
            self.mons = ['postsynapticV', 'postsynapticI', 'dendriteV' ]
        else:
            print 'generate_run.py, mode %s  unknown' % self.runInfo.postMode
            return
        self.hf.h.tstop = maxt+self.runInfo.stimDelay
        print ('PrepareRun: \n   maxt = {:8.2f}'.format(maxt))
        print ('   delay: {:9.3f} ms\n   duration: {:9.3f} ms'.format(self.runInfo.stimDelay, self.runInfo.stimDur))
        print ('   tstop: {:9.3f} ms'.format(self.hf.h.tstop))
        print ("   t:  {:8.2f}\n----------------\n".format(self.hf.h.t))
        self.monitor['time'].record(self.hf.h._ref_t)
        #self.hf.h.topology()
        #pg.show()
        #self.hf.h('access %s' % self.hf.get_section(self.electrodeSite).name())


    def doRun(self, filename=None, parMap=None, save=False, restoreFromFile=False, initfile=None, workers=4):
        if verbose:
            print('generat_run::doRun')
        (p, e) = os.path.splitext(filename)  # make sure filename is clean
        self.runInfo.filename = p  # change filename in structure, just name, no extension
        if parMap is None or len(parMap) == 0:
            self.makeFileName(filename = self.runInfo.filename) # base name pluse underscore
        else:
            mstr = '_'
            for k in parMap.keys():
                if k == 'id':
                    continue
                mstr += k + '_'
            #if 'id' in parMap.keys():
            #    mstr += 'ID%04d_' % parMap['id']
            self.makeFileName(self.runInfo.filename + mstr )

        if verbose:
            print ('genrate_run::doRun: basename is = {:s}'.format(self.basename))
        #self.hf.update() # make sure channels are all up to date
        self.results={}
        
        nLevels = len(self.runInfo.stimInj)
        nWorkers = workers
        print ('doRun: initfile = {:s}'.format(initfile))
        TASKS = [s for s in range(nLevels)]
        tresults = [None]*len(TASKS)
        runner = [None]*nLevels
        
        # run using pyqtgraph's parallel support
        with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
            for i, x in tasker:
                inj = self.runInfo.stimInj[i]
                self._prepareRun(inj=inj) # build the recording arrays
                self.run_initialized = cellInit.initModel(self.hf, mode='cc', restoreFromFile=restoreFromFile, 
                    filename=initfile)
                tresults = self._executeRun() # now you can do the run
                tasker.results[i] = {'r': tresults, 'i': inj}
        
#        self.results = results
        self.results={}
        for i in range(len(tresults)):
            self.results[tresults[i]['i']] = tresults[i]['r']
        for k, i in enumerate(self.runInfo.stimInj):
            if self.plotting:
                if k == 0:
                    self.mons = self.results[i].monitor.keys()
                    self.plotRun(self.results[i], init=True)
                else:
                    self.plotRun(self.results[i], init=False)
        if verbose:
            print ('doRun, calling IV')
        self.arun = ar.AnalyzeRun(self.results) # create an instance of the class with the data
        self.arun.IV()  # compute the IV on the data
        self.IVResult =self.arun.IVResult
        if verbose:
            print ('doRun, back from IV')
        if save == 'monitor':
            self.saveRuns('monitor')
        if verbose:
            print ('doRun: ivresult is: {:32}'.format(self.IVResult))
        if self.plotting:
            self.plotFits(1, self.IVResult['taufit'], c='r')
            self.plotFits(1, self.IVResult['ihfit'], c='b')
            self.cplts.show()


    def testRun(self, title='testing...'):
        self._prepareRun(inj=0.0)
        self.run_initialized = cellInit.initModel(self.hf, restoreFromFile=True)
        self.hf.h.t = 0.
        self.hf.h.tstop = 10
        #self.hf.h.run()
        self._executeRun()
        pg.mkQApp()
        pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
        pl.setTitle(title)
        QtGui.QApplication.instance().exec_()


    def _executeRun(self, testPlot=False):
        """
        (private mmethod)
        After prepare run and initialization, this routine actually calls the run method in hoc
        assembles the data, saves it to disk and plots the results.
        Inputs: flag to put up a test plot....
        """
        if verbose:
            print '_executeRun'
        assert self.run_initialized == True
        print 'Starting Vm at electrode site: ', self.electrodeSite.v
        
        # one way
        self.hf.h.t = 0
        """
        #while (self.hf.h.t < self.hf.h.tstop):
        #    for i=0, tstep/dt {
        #       self.hf.h.fadvance()
        # self.hf.h.run()  # calls finitialize, causes offset
        """
        self.hf.h.batch_save() # save nothing
        self.hf.h.batch_run(self.hf.h.tstop, self.hf.h.dt, "v.dat")
        print 'Finishing Vm: ', self.electrodeSite.v
        if testPlot:
            pg.mkQApp()
            pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
            if self.filename is not None:
                pl.setTitle('%s' % self.filename)
            else:
                pl.setTitle('executeRun, no filename')
#        print 'run done'
        np_monitor = {}
        for k in self.monitor.keys():
            np_monitor[k] = np.array(self.monitor[k])

        np_allsecVec = {}
        for k in self.allsecVec.keys():
            np_allsecVec[k] = np.array(self.allsecVec[k])
        self.runInfo.clist = self.clist
        results = Params(Sections=self.hf.sections.keys(),  vec=np_allsecVec,
                         monitor=np_monitor, stim=self.stim, runInfo=self.runInfo,
                         distanceMap = self.hf.distanceMap,
        )
        if verbose:
            print '    _executeRun completed'
        return results


    def saveRun(self, results):
        """
        Save the result of a single run to disk. Results must be a Param structure, which we turn into
         a dictionary...
        """
        fn = self.basename +  '.p'
        pfout = open(fn, 'wb')
        pickle.dump({'basename': self.basename,
                     'runInfo': self.runInfo.todict(),
                     'modelPars': [],
                     'Results': results.todict(),
                     'IVResults': self.IVResult},
                     pfout)
        pfout.close()


    def saveRuns(self, save=None):
        """
        Save the result of multiple runs to disk. Results is in a dictionary,
        each element of which is a Param structure, which we then turn into
        a dictionary...
        """
        fn = self.basename + '.p'
        pfout = open(fn, 'wb')
        pickle.dump({'basename': self.basename,
                     'runInfo': self.runInfo.todict(),
                     'modelPars': [],
                     'Results': [{k:x.todict()} for k,x in self.results.iteritems()]}, pfout)
                     
        pfout.close()
        return (self.runInfo.folder, self.basename) # return tuple to assemble name elsewhere


        # if recordSection:
        #     fns = os.path.join(self.runInfo.folder, self.runInfo.fileName + '_swellings_' + dtime + '.p')
        #     srsave = sr.SimulationResult()
        #     srsave.save(fns, data=np.array(self.vswel),
        #             time=np.array(self.vec['time']),
        #             hoc_file = 'MorphologyFiles/' + self.modelPars.file,
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


    def plotRun(self, results, init=True, show=False):
        if init:
            self.cplts = cp.CalyxPlots(title=self.filename)
        self.cplts.plotResults(results.todict(), self.runInfo.todict(), somasite=self.mons)
        if show:
            self.cplts.show()

    def plotFits(self, panel, x, c='g'):
        self.cplts.plotFit(panel, x[0], x[1], c)
