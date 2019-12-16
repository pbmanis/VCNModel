from __future__ import print_function
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
from pathlib import Path
import copy
from datetime import datetime
from collections import OrderedDict
from cnmodel.util.stim import make_pulse # makestim
import vcnmodel.NoiseTrainingGen as NG
#from NoiseTrainingGen.NoiseGen import generator
from pylibrary.params import Params
import vcnmodel.IVPlots as IVP
import vcnmodel.analyze_run as ar
import vcnmodel.cellInitialization as cellInit
import time
import csv
import pickle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import pyqtgraph.multiprocess as mproc
import errno
import signal


verbose = False # use this for testing.


class GenerateRun():
    def __init__(self, cell, idnum=0, celltype=None,
                 electrodeSection=None,
                 dendriticElectrodeSection=None,
                 starttime=None,
                 stimtype='IV',
                 iRange=[-1, 0, 1],
                 plotting=False,
                 saveAllSections=False,
                 useSavedState=True,
                 params=None):
        
        self.run_initialized = False
        self.plotting = plotting
        #print(dir(cell))
        self.cell = cell  # cnmodel cell instance
        self.params = params  # copied from the amin
        self.filename = params['cell']
        self.hf = cell.hr # get the reader structure and the hoc pointer object locally
        self.basename = 'basenamenotset'
        self.saveAllSections = params['save_all_sections']
        print('Save all sections: ', self.saveAllSections)
        self.idnum = idnum
        self.startTime = starttime
        # use the Params class to hold program run information and states
        # runInfo holds information that is specific the the run - stimuli, conditions, etc
        # DO NOT add .p to filename...
        ri = self.params['runInfo']
        self.params['runInfo'].TargetCellType = celltype
        self.params['runInfo'].electrodeSection = electrodeSection
        self.params['runInfo'].dendriticElectrodeSection=dendriticElectrodeSection
        self.params['runInfo'].stimInj = iRange
        self.params['runInfo'].gif_i0 = params['gif_i0']  # just copy down here.
        self.params['runInfo'].gif_sigma = params['gif_sigma']
        self.params['runInfo'].gif_fmod = params['gif_fmod']
        self.params['runInfo'].gif_tau = params['gif_tau']
        self.params['runInfo'].gif_dur = params['gif_dur']
        self.params['runInfo'].gif_skew = params['gif_skew']
        self.params['runInfo'].runTime = time.asctime()  # store date and time of run        
        self.params['runInfo'].useSaveState = useSavedState  # use the saved state.

        if stimtype == 'IV':
            self.params['runInfo'].postMode = 'cc'
        elif stimtype == 'VC':
            self.params['runInfo'].postMode = 'vc'
        elif stimtype == 'gifnoise':
            self.params['runInfo'].postMode = 'gifnoise'
            if params is None:
                raise ValueError('generate_run: gifnoise mode requires params to be passed')

        electrodeSection = list(self.hf.sec_groups[self.params['runInfo'].electrodeSection])[0]
        self.electrode_site = self.hf.get_section(electrodeSection)
        dend_sections = list(self.hf.sec_groups[self.params['runInfo'].dendriticElectrodeSection])
        # invert the mappint of sections
        # first, just get the dendrite components, making a dendrite submap
        dendDistMap = {}
        for k in dend_sections:
            dendDistMap[k] = self.hf.distanceMap[k]
        revmap = dict((v, k) for k, v in dendDistMap.items())
            
        # now find the distlist key that corresponds to the closest value to our desired value
        (num, dist) = min(enumerate(revmap), key=lambda x: abs(x[1]-self.params['runInfo'].dendriticSectionDistance))
        print('Monitoring dendrite section number {:d}, at {:6.2f} microns from soma: '.format(num, dist))
        
        dendriticElectrodeSection = list(self.hf.sec_groups[self.params['runInfo'].dendriticElectrodeSection])[num]
        self.dendriticElectrodeSite = self.hf.get_section(dendriticElectrodeSection)

        if self.params['runInfo'].inFile is None:
            self.params['runInfo'].tstop = 1000.0 * (
                self.params['runInfo'].nStim - 1) / self.params['runInfo'].stimFreq + self.params['runInfo'].stimPost + self.params['runInfo'].stimDelay

        else:
            print('Reading spike train')
            maxt = 0.
            with open(self.params['runInfo'].inFile, 'r') as csvfile:
                spks = csv.reader(csvfile, delimiter=',')
                for i, row in enumerate(spks):
                    if i == 0:  # first line
                        print(row)
                        maxt = float(row[1]) * 1000
                        reps = int(row[0])
                        continue
                    if int(row[0]) in self.params['runInfo'].spikeTimeList.keys():
                        self.params['runInfo'].spikeTimeList[int(row[0])].append(float(row[1]) * 1000.)
                    else:
                        self.params['runInfo'].spikeTimeList[int(row[0])] = [float(row[1]) * 1000.]
                self.params['runInfo'].tstop = maxt

        self.monitor = OrderedDict() # standard monitoring
        self.allsecVec = OrderedDict() # all section monitoring
        self.mons = OrderedDict()


        if verbose:
            print('Runinfo initialization done')

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def makeFileName(self, filename=None, subdir=None):
        self.params['runInfo'].filename = filename
        if self.startTime is None:
            self.dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        else:  # convert time object
            self.dtime = dtime = time.strftime("%y.%m.%d-%H.%M.%S", time.localtime(self.startTime))
        self.mkdir_p(self.params['runInfo'].folder)
        #f os.path.exists(self.params['runInfo'].folder) is False:
        #    os.mkdir(self.params['runInfo'].folder)
        if subdir is not None:
            folder = os.path.join(self.params['runInfo'].folder, subdir)
        else:
            folder = self.params['runInfo'].folder
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
            print('_prepareRun')
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
                    self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)  # recording of voltage all set up here

        for var in ['time', 'postsynapticV', 'dendriteV', 'postsynapticI', 'i_stim0', 'v_stim0']: # get standard stuff
            self.monitor[var] = self.hf.h.Vector()

        self.params['runInfo'].celsius = self.cell.status['temperature']
        self.hf.h.celsius = self.params['runInfo'].celsius
        
        self.clist = {}  #  color list

        electrodeSection = list(self.hf.sec_groups[self.params['runInfo'].electrodeSection])[0]
        self.electrode_site = self.hf.get_section(electrodeSection)
        
        if self.params['runInfo'].postMode in ['vc', 'vclamp']:
            # Note to self (so to speak): the hoc object returned by this call must have a life after
            # the routine exits. Thus, it must be "self." Same for the IC stimulus...
            self.vcPost = self.hf.h.SEClamp(0.5, sec=self.electrode_site) #self.hf.sections[electrode_site])
            self.vcPost.dur1 = 2
            self.vcPost.amp1 = self.params['runInfo'].vstimHolding
            self.vcPost.dur2 = 1e9
            self.vcPost.amp2 = self.params['runInfo'].vstimHolding  # just a tiny step to keep the system honest
            self.vcPost.dur3 = 10.0
            self.vcPost.amp3 = self.params['runInfo'].vstimHolding
            self.vcPost.rs = 1e-6
            stim = {}
            stim['NP'] = self.params['runInfo'].vnStim
            stim['Sfreq'] = self.params['runInfo'].vstimFreq  # stimulus frequency
            stim['delay'] = self.params['runInfo'].vstimDelay
            stim['dur'] = self.params['runInfo'].vstimDur
            stim['amp'] = self.params['runInfo'].vstimInj
            stim['PT'] = 0.0
            stim['dt'] = self.hf.h.dt
            (secmd, maxt, tstims) = make_pulse(stim)
            self.stim = stim
            secmd = secmd + self.params['runInfo'].vstimHolding # add holding
            self.monitor['v_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['v_stim0'].play(self.vcPost._ref_amp2, self.hf.h.dt, 0, sec=self.electrode_site)
            self.monitor['postsynapticV'].record(self.electrode_site(0.5)._ref_v, sec=self.electrode_site)
            self.monitor['dendriteV'].record(self.dendriticElectrodeSite(0.5)._ref_v, sec=self.dendriticElectrodeSite)
            self.monitor['postsynapticI'].record(self.vcPost._ref_i, sec=self.electrode_site)
            self.mons = ['postsynapticI', 'v_stim0']
            
        elif self.params['runInfo'].postMode in ['cc', 'iclamp']:
            stim = {}
            stim['NP'] = self.params['runInfo'].nStim
            stim['Sfreq'] = self.params['runInfo'].stimFreq  # stimulus frequency
            stim['delay'] = self.params['runInfo'].stimDelay
            stim['dur'] = self.params['runInfo'].stimDur
            if inj is not None:
                stim['amp'] = inj
            else:
                stim['amp'] = self.params['runInfo'].stimInj['pulse'][0]
            stim['PT'] = 0.0
            stim['dt'] = self.hf.h.dt
            (secmd, maxt, tstims) = make_pulse(stim) # cnmodel.makestim.makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.stim = stim
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrode_site)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor['i_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['i_stim0'].play(self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrode_site)
            self.monitor['postsynapticI'].record(self.icPost._ref_i, sec=self.electrode_site)
            self.monitor['postsynapticV'].record(self.electrode_site(0.5)._ref_v, sec=self.electrode_site)
            self.monitor['dendriteV'].record(self.dendriticElectrodeSite(0.5)._ref_v, sec=self.dendriticElectrodeSite)
            self.mons = ['postsynapticV', 'postsynapticI', 'dendriteV' ]
        
        elif self.params['runInfo'].postMode in ['gifnoise']:
            stim = {}
            self.stim = NG.NoiseGen()
            self.stim.generator(
                dt = self.hf.h.dt,
                i0=self.params['runInfo'].gif_i0,
                sigma0=self.params['runInfo'].gif_sigma,
                fmod=self.params['runInfo'].gif_fmod,
                tau=self.params['runInfo'].gif_tau,
                dur=self.params['runInfo'].gif_dur,
                skew=self.params['runInfo'].gif_skew,)
            maxt = 1000.*self.params['runInfo'].gif_dur
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrode_site)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor['i_stim0'] = self.hf.h.Vector(self.stim[1])
            self.monitor['i_stim0'].play(self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrode_site)
            self.monitor['postsynapticI'].record(self.icPost._ref_i, sec=self.electrode_site)
            self.monitor['postsynapticV'].record(self.electrode_site(0.5)._ref_v, sec=self.electrode_site)
            self.monitor['dendriteV'].record(self.dendriticElectrodeSite(0.5)._ref_v, sec=self.dendriticElectrodeSite)
            self.mons = ['postsynapticV', 'postsynapticI', 'dendriteV' ]
        
        else:
            print('generate_run.py, mode %s  unknown' % self.params['runInfo'].postMode)
            return
        
        self.hf.h.tstop = maxt+self.params['runInfo'].stimDelay
        print('PrepareRun: \n')
        print(f'   maxt:     {maxt:8.2f} ms')
        print(f"   delay:    {self.params['runInfo'].stimDelay:8.2f} ms")
        print(f"   duration: {self.params['runInfo'].stimDur:8.2f} ms")
        print(f'   tstop:    {self.hf.h.tstop:8.2f} ms')
        print(f"   h.t:      {self.hf.h.t:8.2f} ms")
        print(f"   h.dt      {self.hf.h.dt:8.2f} us")
        print("\n----------------\n")
        self.monitor['time'].record(self.hf.h._ref_t)
        #self.hf.h.topology()
        #pg.show()
        #self.hf.h('access %s' % self.hf.get_section(self.electrode_site).name())

    def doRun(self, filename=None, parMap=None, save=False, restore_from_file=False, initfile=None, workers=4):
        if verbose:
            print('generate_run::doRun')
        (p, e) = os.path.splitext(filename)  # make sure filename is clean
        self.params['runInfo'].filename = p  # change filename in structure, just name, no extension
        if parMap is None or len(parMap) == 0:
            self.makeFileName(filename = self.params['runInfo'].filename) # base name pluse underscore
        else:
            mstr = '_'
            for k in parMap.keys():
                if k == 'id':
                    continue
                mstr += k + '_'
            #if 'id' in parMap.keys():
            #    mstr += 'ID%04d_' % parMap['id']
            self.makeFileName(self.params['runInfo'].filename + mstr )

        if verbose:
            print('genrate_run::doRun: basename is = {:s}'.format(self.basename))
        #self.hf.update() # make sure channels are all up to date
        self.results={}
        
        if self.params['runInfo'].postMode in ['cc', 'iclamp']:
            s = self.params['runInfo'].stimInj['pulse']
            ipulses = np.arange(s[0], s[1], s[2])
        else:
            ipulses = [0]
        nLevels = len(ipulses)
        nWorkers = workers
        # print(f"doRun: initfile = {str(initfile):s}")
        TASKS = [s for s in range(nLevels)]
        tresults = [None]*len(TASKS)
        runner = [None]*nLevels
        signal.signal(signal.SIGCHLD, signal.SIG_DFL)  # might prevent OSError "no child process"
        # run using pyqtgraph's parallel support
        with mproc.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
            for i, x in tasker:
                inj = ipulses[i]
                self._prepareRun(inj=inj) # build the recording arrays
                self.run_initialized = cellInit.init_model(self.cell, mode='iclamp', restore_from_file=restore_from_file, 
                    filename=initfile)
                tr = {'r': self._executeRun(), 'i': inj} # now you can do the run
                tasker.results[i] = tr
        
        self.results = OrderedDict()

        for i in range(nLevels):
           #  print('level: ', i, '  tresults[i]["i"]: ', tresults[i]['i'])
            self.results[tresults[i]['i']] = tresults[i]['r']
        for k, i in enumerate(ipulses):
            if self.plotting:
                if k == 0:
                    self.mons = self.results[i].monitor.keys()
                    self.plotRun(self.results[i], init=True)
                else:
                    self.plotRun(self.results[i], init=False)
        if self.params['runInfo'].postMode in ['cc', 'iclamp']:
            if verbose:
                print ('doRun, calling IV')
            self.arun = ar.AnalyzeRun(self.results) # create an instance of the class with the data
            self.arun.IV()  # compute the IV on the data
            self.IVResult = self.arun.IVResult
            if verbose:
                print ('doRun, back from IV')
            if save == 'monitor':
                self.saveRuns(save='monitor')
            if verbose:
                print('doRun: ivresult is: {:32}'.format(self.IVResult))
            if self.plotting:
                self.plotFits('Soma', self.IVResult['taufit'], c='r')
                self.plotFits('Soma', self.IVResult['ihfit'], c='b')
                #print (dir(self.ivplots))
                self.ivplts.show()
        if self.params['runInfo']['postMode'] in ['gifnoise']:
            self.IVResult = None
            if save == 'monitor':
                self.saveRuns('gifnoise')

    def testRun(self, title='testing...', initfile=None):
        if initfile is None:
            raise ValueError('generate_run:testRun needs initfile name')
        self._prepareRun(inj=0.0)
        self.run_initialized = cellInit.init_model(self.cell, filename=initfile, restore_from_file=True)
        self.hf.h.t = 0.
        self.hf.h.tstop = 10
        #self.hf.h.run()
        self._executeRun(testPlot=True)
        pg.mkQApp()
        pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
        pl.setTitle(title)
        QtGui.QApplication.instance().exec_()

    def _executeRun(self, testPlot=False):
        """
        (private method)
        After prepare_run and the initialization, this routine actually calls the run method in hoc
        assembles the data, saves it to disk and plots the results.
        Inputs: flag to put up a test plot....
        """
        if verbose:
            print('_executeRun')
        assert self.run_initialized == True
        print(f'Starting Vm at electrode site: {self.electrode_site.v:6.2f}')
        
        # one way
        self.hf.h.t = 0

        self.hf.h.batch_save() # save nothing
        print (f'Temperature in run at start (deg C): {self.hf.h.celsius:6.1f}')
        self.hf.h.batch_run(self.hf.h.tstop, self.hf.h.dt, "v.dat")
        print(f'Finishing Vm: {self.electrode_site.v:6.2f}')
        self.monitor['time'] = np.array(self.monitor['time'])
        self.monitor['time'][0] = 0.
        if verbose:
            print(f"Post V: {self.monitor['postsynapticV']:.3f}")
        if testPlot:
           pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
           QtGui.QApplication.instance().exec_()
           # pg.mkQApp()
            # pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
            # if self.filename is not None:
            #     pl.setTitle('%s' % self.filename)
            # else:
            #     pl.setTitle('executeRun, no filename')
        print ('    Run finished')
        np_monitor = {}
        for k in self.monitor.keys():
            np_monitor[k] = np.array(self.monitor[k])

        np_allsecVec = OrderedDict()
        for k in self.allsecVec.keys():
            np_allsecVec[k] = np.array(self.allsecVec[k])
        self.params['runInfo'].clist = self.clist
        results = Params(Sections=list(self.hf.sections.keys()),  vec=np_allsecVec,
                         monitor=np_monitor, stim=self.stim, runInfo=self.params['runInfo'],
                         distanceMap = self.hf.distanceMap,
        )
        if verbose:
            print('    _executeRun completed')
        return results

    def saveRun(self, results):
        """
        Save the result of a single run to disk. Results must be a Param structure, which we turn into
         a dictionary...
        """
        print('Saving single run to: ')
        # fn = Path(f"{self.basename:s}_{self.cell.status['modelType']:s}")
        t = datetime.now()
        fn = f"{self.params['shortSimulationFilename']:s}_{t.strftime('%Y-%m-%d-%H%M%S-%f'):s}" + '.p'
        fn = Path(self.params['simulationFilename'].parent, fn)
        print(f"     {str(fn):s}")
        mp = copy.deepcopy(self.cell.status)
        del mp['decorator']
        # restructure data so all parameters are in the Params dict
        self.params['modelPars'] = mp  # replace the default with current values
        self.params['runInfo'] = self.params['runInfo'].todict()  # likewise for runInfo
        
        with( open(fn, 'wb')) as pfout:
            pickle.dump({'basename': self.basename,
                     'Params': self.params,
                     'Results': results.todict(),
                     'IVResults': self.IVResult},
                     pfout)
        # pfout.close()

    def saveRuns(self, save=None):
        """
        Save the result of multiple runs to disk. Results is in a dictionary,
        each element of which is a Param structure, which we then turn into
        a dictionary...
        """
        # if save is None:
        #     fn = f"{self.basename:s}{self.cell.status['modelName']:s}_{self.cell.status['modelType']:s}.p"
        # else:
        #     fn = f"{self.basename:s}{self.cell.status['modelName']:s}_{self.cell.status['modelType']:s}_{save:s}.p"
        #     # (f'{0:s}_{1:s}_{2:s}.p'.format(self.basename, self.cell.status['modelType'], save))
        t = datetime.now()
        
        fn = f"{self.params['shortSimulationFilename']:s}_{t.strftime('%Y-%m-%d-%H%M%S-%f'):s}" + '.p'
        fn = Path(self.params['simulationFilename'].parent, fn)
        print(f"WRITING DATA TO: {str(fn):s}")
        # pfout = open(fn, 'wb')
        mp = copy.deepcopy(self.cell.status)
        del mp['decorator']
        # restructure data so all parameters are in the Params dict
        self.params['modelPars'] = mp  # replace the default with current values
        folder = self.params['runInfo'].folder # get before turing into dict
        self.params['runInfo'] = self.params['runInfo'].todict()  # likewise for runInfo
        with(open(fn, 'wb')) as pfout:
            pickle.dump({'basename': self.basename,
                     'Params': self.params,
                     'Results': [{k:x.todict()} for k,x in self.results.items()]
                    }, pfout)
                     
        return (folder, self.basename) # return tuple to assemble name elsewhere

        # if recordSection:
        #     fns = os.path.join(self.params['runInfo'].folder, self.params['runInfo'].fileName + '_swellings_' + dtime + '.p')
        #     srsave = sr.SimulationResult()
        #     srsave.save(fns, data=np.array(self.vswel),
        #             time=np.array(self.vec['time']),
        #             hoc_file = 'MorphologyFiles/' + self.params['modelPars'].file,
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
            self.ivplts = IVP.IVPlots(title=self.filename, mode='mpl')
        self.ivplts.plotResults(results.todict(), self.params['runInfo'].todict(), somasite=self.mons)
        if show:
            self.ivplts.show()

    def plotFits(self, panel, x, c='g'):
        self.ivplts.plotFit(panel, x[0], x[1], c)
