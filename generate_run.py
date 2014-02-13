__author__ = 'pbmanis'

import os
import numpy as np
import nrnlibrary.makestim as makestim
from pylibrary.Params import Params
import CalyxPlots as cp
import time
import csv
import pickle

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
verbose = False

class GenerateRun():
    def __init__(self, hf):
        self.run_initialized = False
        self.hf = hf # get the reader structure and the hoc pointer object locally
        #  use the Params class to hold program run information and states
        # runInfo holds information that is specific the the run - stimuli, conditions, etc
        # DO NOT add .p to filename...
        self.runInfo = Params(fileName='Normal', runName='Run', manipulation="Canonical", folder="Simulations",
                              preMode="cc",
                              postMode='cc',
                              TargetCellName='Bushy',
                              TargetCell=None,
                              celsius=22,  # set the temperature.
                              nStim=1,
                              stimFreq=200.,  # hz
                              stimInj=2.0,  # nA
                              stimDur=2.0,  # msec
                              stimDelay=2.0,  # msec
                              stimPost=3.0,  # msec
                              vnStim=1,
                              vstimFreq=200.,  # hz
                              vstimInj=50,  # mv
                              vstimDur=4.0,  # msec
                              vstimDelay=2.0,  # msec
                              vstimPost=3.0,  # msec
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

        if verbose:
            print 'runinfo initialization done'


    def doRun(self):
        self.hf.update() # make sure channels are all up to date
        self.prepareRun()
        print 'run preparae'
        self.initRun()  # this also sets nseg in the axon - so do it before setting up shapes
        print 'run initialized'
        self.executeRun()
        print self.hf.h.tstop
        #self.hf.h.run()

    def prepareRun(self):
        """
        Control a single run of the model with updated display of the voltages, etc.
        Inputs: None
        Outputs: None
        Actions: optionally displays the results
        Side Effects: A number of class variables are created and modified, mostly related to the
        generation of stimuli and monitoring of voltages and currents
        """
        self.monitor = {} # standard monitoring
        self.allsecVec = {} # all section monitoring
        for var in self.hf.sections: # get morphological components
            self.allsecVec[var] = self.hf.h.Vector()
            self.allsecVec[var].record(self.hf.sections[var](0.5)._ref_v, sec=self.hf.sections[var])

        for var in ['time', 'postsynapticV', 'postsynapticI', 'i_stim0', 'v_stim0']: # get standard stuff
            self.monitor[var] = self.hf.h.Vector()

        self.hf.h.tstop = 50 # self.runInfo.tstop
#        self.hf.h.dt = 0.02  # force small time step. cvode is probably off.
        self.hf.h.celsius = self.runInfo.celsius
        npts = 1 + int(self.hf.h.tstop / self.hf.h.dt)  # number of points in a run

        self.clist = {}  #  color list

       # make fake cell right here...
       # electrodeSite = 'soma[0]'
        soma = self.hf.h.Section()
        soma.L = 20.
        soma.diam = 20.
        soma.insert('hh')

        electrodeSite = soma
       # print self.hf.sections
       # print 'esite: ', self.hf.sections[electrodeSite]
        if self.runInfo.postMode in ['vc', 'vclamp']:
            print 'vclamp'
            clampV = -60.0
            vcPost = self.hf.h.SEClamp(0.5, sec=electrodeSite) #self.hf.sections[electrodeSite])
            vcPost.dur1 = 2
            vcPost.amp1 = clampV
            vcPost.dur2 = 10.0
            vcPost.amp2 = clampV + 50.0  # just a tiny step to keep the system honest
            vcPost.dur3 = 10.0
            vcPost.amp3 = clampV
            vcPost.rs = 1e-6
            stim = {}
            stim['NP'] = self.runInfo.vnStim
            stim['Sfreq'] = self.runInfo.vstimFreq  # stimulus frequency
            stim['delay'] = self.runInfo.vstimDelay
            stim['dur'] = self.runInfo.vstimDur
            stim['amp'] = self.runInfo.vstimInj
            stim['PT'] = 0.0
           # print self.hf.h.soma[0]
            (secmd, maxt, tstims) = makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.monitor['v_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['v_stim0'].play(vcPost._ref_amp2, self.hf.h.dt, 0, sec=self.hf.h.axon[0]) # self.hf.sections[electrodeSite])
            self.monitor['postsynapticV'].record(electrodeSite(0.5)._ref_v, sec=electrodeSite)
            self.monitor['postsynapticI'].record(vcPost._ref_i, sec=electrodeSite)
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
            icPost = self.hf.h.iStim(0.5, sec=electrodeSite)
            icPost.delay = 2
            icPost.dur = 4  # these actually do not matter...
            icPost.iMax = 1.0
            self.monitor['i_stim0'] = self.hf.h.Vector(secmd)
            self.monitor['i_stim0'].play(icPost._ref_i, self.hf.h.dt, 0, sec=electrodeSite)
            self.monitor['postsynapticI'].record(icPost._ref_i, sec=electrodeSite)
            self.monitor['postsynapticV'].record(electrodeSite(0.5)._ref_v, sec=electrodeSite)
            self.mons = ['postsynapticV', 'postsynapticI' ]
        else:
            print 'mode unknown'
            return
        # ns = self.hf.h.Section()
        # ix = self.hf.h.iStim(0.1, sec=ns)
        # ix.delay = 2.0
        # ix.dur = 2.0
        # ix.iMax = 5.0
        # vec = self.hf.h.Vector()
        # vec.recor(ns._ref_v, sec=ns)
        self.electrodeSite = electrodeSite
        self.monitor['time'].record(self.hf.h._ref_t)
        # self.hf.h.tstop = 10
        # self.hf.h.finitialize()
        # self.hf.h.run()


    def initRun(self):
        """
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
        # First we set e_leak so that the rmp in each segment is the same
        self.hf.h.finitialize(self.runInfo.v_init)
        i = 0
        for si in self.hf.sections: # self.section_list[s]:
            self.hf.h('access %s' % si)
            sec = self.hf.h.cas()
            if i == 0:
                i += 1
            mechs = self.hf.get_mechanisms(si)
#            print mechs
            #e_newleak = sec.v + (sec.ina + sec.ik + sec.i_ihvcn) / (
            #    sec.g_leak * 1e3)  # hmmm. This is probably NOT right.
            #sec().e_leak = e_newleak

        # based on those conditions, make sure spatial grid is fine enough for our needs
        #h.Mesh_th = self.Mesh_th
        #h('mesh_init')
        #self.mesh_init()
        self.hf.h.finitialize(self.runInfo.v_init)
        #       h('axon[axonnode] ic.loc(iclocation)')

        # for name in modelPars.structureNames:
        #     if name in ['axon', 'heminode']:
        #         for sec in self.h.parentaxon:
        #             sec.nseg = 11

                    # get ready to run the system to a stable point for the calcium pump
        j = 0
        # # save the calcium pump information
        # for input in range(self.modelPars.AN_conv[runInfo.TargetCellName]):
        #     for sec in self.CalyxStruct[input]['axon']:
        #         savcore = sec().cabulk_capmp
        #         sec().cabulk_capmp = runInfo.ca_init
        #         savtau = sec().tau_capmp
        #         sec().tau_capmp = 1e-6  # make the pump go really fast

                # starting way back in time
        self.hf.h.t = -1e10
        dtsav = self.hf.h.dt
        self.hf.h.dt = 1e9  # big time steps for slow process
        temp = self.hf.h.cvode.active()
        if (temp != 0):
            self.hf.h.cvode.active(
                0)  # turn cvode off (note, in this model it will be off because one of the mechanisms is not compatible with cvode at this time
        while (self.hf.h.t < -1e9):
            self.hf.h.fadvance()

            # now restore the pump values and get the system ready for runs
        if (temp != 0):
            self.hf.h.cvode.active(1)
        self.hf.h.dt = dtsav
        self.hf.h.t = 0
        # for input in range(self.modelPars.AN_conv[runInfo.TargetCellName]):
        #     for sec in self.CalyxStruct[input]['axon']:
        #         sec().cabulk_capmp = savcore
        #         sec().tau_capmp = savtau

        if (self.hf.h.cvode.active()):
            self.hf.h.cvode.re_init()
        else:
           self.hf.h.fcurrent()
        self.hf.h.frecord_init()
        self.hf.h.finitialize(self.hf.h.v_init)
        self.run_initialized = True





    def executeRun(self):
        print 'Running for: ', self.hf.h.tstop
        print 'V = : ', self.electrodeSite.v
        self.hf.h.frecord_init()
#        self.hf.h.finitialize() # self.hf.h.v_init)
        self.hf.h.run()
        pg.mkQApp()
        pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
        pl.setTitle('testing...')

#        pg.mkQApp()
#        pg.plot(self.monitor['time'], self.monitor['postsynapticV'])# if self.run_initialized:
        #     self.hf.h.run()
        # else:
        #     print 'run not initialized!'
        #     return
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
        )
        dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        print("Run done\n")
        # if os.path.exists(self.runInfo.folder) is False:
        #     os.mkdir(self.runInfo.folder)
        # fn = os.path.join(self.runInfo.folder, self.runInfo.fileName + dtime + '.p')
        # pfout = open(fn, 'wb')
        # pickle.dump({'runInfo': self.runInfo.todict(),
        #              'modelPars': [],
        #              'Results': results.todict()}, pfout)
        # pfout.close()
        #
        # if recordSwellings:
        #     fns = os.path.join(runInfo.folder, runInfo.fileName + '_swellings_' + dtime + '.p')
        #     srsave = sr.SimulationResult()
        #     srsave.save(fns, data=np.array(self.vswel),
        #             time=np.array(self.vec['time']),
        #             hoc_file = 'MorphologyFiles/' + self.modelPars.topofile,
        #             section_map=secarray,
        #     )
        #     # pfout = open(fns, 'wb')
        #     # pickle.dump({'swellings': secarray,
        #     #              'time': np.array(self.vec['time']),
        #     #              'data': np.array(self.vswel)}, pfout)
        #     # pfout.close()
        #     # plot = pg.plot()
        #     # vs = np.array(self.vswel)
        #     # ts = np.array(self.vec['time'])
        #     # for i in range(len(self.swellAxonMap)):
        #     #     plot.plot(ts, vs[i])

        cplts = cp.CalyxPlots()
        cplts.plotResults(results.todict(), self.runInfo.todict(), somasite=self.mons)
        cplts.show()
