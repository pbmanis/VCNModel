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

import time
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

make_init = False
test_init = False
make_ANIntialConditions = False
ANPSTH_mode = False
ANSingles = True
IV_mode = False


class ModelRun():
    def __init__(self, args=None):
        if verbose:
            print args
        if isinstance(args, list) and len(args) > 0:
            self.set_celltype(args[0]) # must be string, not list...
        else:
            self.set_celltype('Bushy')
        if isinstance(args, list) and len(args) > 1:
            self.set_modeltype(args[1])
        else:
            self.set_modeltype('XM13')
        if isinstance(args, list) and len(args) > 2:
            infile = args[2]
        else:
            raise ValueError (' need a valid hoc file to read.')
        self.startTime = None
        self.plotFlag = False
        #infile = 'L23pyr.hoc'
        #infile = 'LC_nmscaled_cleaned.hoc'
        #infile = 'Calyx-68cvt2.hoc'
        #infile = 'Calyx-S53Acvt3.hoc'
        #infile = 'wholeThing_cleaned.hoc'
        #infile = 'MNTB_Cell2_cleaned.hoc'
        #infile = 'VCN_Dend.hoc'
        #infile = 'somaOnly.hoc'
        self.infile = infile
        self.srname = ['**', 'LS', 'MS', 'HS']  # runs 1-3, not starting at 0
        self.run_duration = 0.25
        self.pip_duration = 0.1
        self.pip_start = [0.1]
        self.Fs = 100e3
        self.f0 = 4000.
        self.dB = 40.
        self.RF = 2.5e-3
        self.SR = 3  # if None, DO NOT Reconfigure, otherwise, override
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        
        


    def runModel(self, parMap={}):
        if verbose:
            print 'runModel entry'
        if 'id' in parMap.keys():
            self.idnum = parMap['id']
        else:
            self.idnum = 9999

        if parMap == {}:
            self.plotFlag = True
        filename = os.path.join('MorphologyFiles/', self.infile)
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

        self.electrodeSection = 'soma'
        self.hg = hoc_graphics
        self.get_hoc_file(filename)
        sg = self.hf.sec_groups['soma']
        self.distances(self.hf.get_section(list(sg)[0]).name()) # make distance map from soma
        if verbose:
            print 'Parmap in runModel: ', parMap
        self.cd = ChannelDecorate(self.hf, celltype=self.cellType, modeltype=self.modelType,
                             parMap=parMap)

        # self.hf.h.topology()
        # self.cd.channelValidate(self.hf, verify=False)

        #return
        # for group in self.hf.sec_groups.keys():
        #     g = self.hf.sec_groups[group]
        #     for section in list(g):
        #         secinfo = self.hf.get_section(section)

        if make_init:
            if verbose:
                print 'make_init'
            self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.cellType,
                             starttime=self.startTime,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag)
            self.R.getInitialConditionsState(tdur=3000.)
            print 'Ran to get initial state for %f msec' % self.hf.h.t
            return

        if test_init:
            if verbose:
                print 'test_init'
            self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.cellType,
                             starttime=self.startTime,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag, )
            self.R.testRun()
            return  # that is ALL, never make init and then keep running.

        if ANPSTH_mode or make_ANIntialConditions:
            if verbose:
                print 'ANPSTH or make_ANInit'
            self.ANRun(self.hf)

        if ANSingles:
            if verbose:
                print 'ANSingles'
            self.ANRun_singles(self.hf)
            
        if IV_mode:
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
        self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.cellType,
                             starttime=self.startTime,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag)

        if verbose:
            print 'doRun'
        self.R.doRun(self.infile, parMap, save='monitor', restoreFromFile=True)
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


    def set_celltype(self, celltype):
        validCells = ['Bushy', 'Stellate', 'L23pyr']
        self.cellType = celltype
        if self.cellType not in validCells:
            print 'Celltype must be one of: %s. Got: %s', (', '.join(validCells), self.cellType)
            exit()


    def set_modeltype(self, modeltype):
        validModels = ['RM03', 'XM13', 'XM13Simple', 'MS']
        self.modelType = modeltype
        if self.modelType not in validModels:
            print 'Model type must be one of: %s. Got: %s ' % (', '.join(validModels), self.modelType)
            exit()


    def set_starttime(self, starttime):
        self.startTime = starttime

    VCN_c18_synconfig = [(int(216.66*0.65), 0., 2), (int(122.16*0.65), 0., 2), 
        (int(46.865*0.65), 0., 2), (int(84.045*0.65), 0., 2), (int(2.135*0.65), 0, 2), (int(3.675*0.65), 0, 2), (int(80.27*0.65), 0, 2)]

    test_synconfig = [(80, 0, 2)]  # if switching configs, need to re-run with AN_InitialConditions to set up the
    # save/restore state. Any change in input configuration requires a new "state" 


    def ANRun(self, hf, verify=False, seed=0, synapseConfig=VCN_c18_synconfig):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        """

        self.start_time = time.time()

        nReps = 1
        threshold = -20. # spike threshold, mV

        stimInfo = {'Morphology': self.infile, 'synapseConfig': synapseConfig,
                    'runDur': self.run_duration, 'pip_dur': self.pip_duration, 'pip_start': self.pip_start,
                    'run_duration': self.run_duration,
                    'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
                    'cellType': self.cellType, 'modelType': self.modelType, 'nReps': nReps, 'threshold': threshold}

        preCell, postCell, self.electrodeSite = self.configureCell(hf, synapseConfig, stimInfo)

        # see if we need to save the cell state now.
        if make_ANIntialConditions:
            print 'getting initial conditions for AN'
            cellInit.getInitialConditionsState(hf, tdur=3000., 
                filename='an_neuronstate.dat', electrodeSite=self.electrodeSite)
            cellInit.testInitialConditions(hf, filename='an_neuronstate.dat',
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
        for j, N in enumerate(range(nReps)):
            print 'Rep: %d' % N
            
            preCell, postCell = self.singleANRun(hf, j, synapseConfig, stimInfo, seeds, preCell, postCell)

            celltime.append(self.time)
            spikeTimes[N] = pu.findspikes(self.time, self.Vsoma, threshold, t0=0., t1=hf.h.tstop, dt=1.0, mode='peak')
            inputSpikeTimes[N] = [preCell[i]._spiketrain for i in range(len(preCell))]
            somaVoltage[N] = np.array(self.Vsoma)

        total_elapsed_time = time.time() - self.start_time
#        total_run_time = time.time() - run_time
        print "Total Elapsed Time = %8.2f min (%8.0fs)" % (total_elapsed_time/60., total_elapsed_time)
        print "Total Setup Time = %8.2f min (%8.0fs)" % (self.setup_time/60., self.setup_time)
        print "Total AN Calculation Time = %8.2f min (%8.0fs)" % (self.an_setup_time/60., self.an_setup_time)
        print "Total Neuron Run Time = %8.2f min (%8.0fs)" % (self.nrn_run_time/60., self.nrn_run_time)
        
        result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes, 
            'somaVoltage': somaVoltage, 'time': np.array(self.time)}
        
        self.analysis_filewriter(self.infile, result)
        self.plotAN(np.array(self.time), result['somaVoltage'], result['stimInfo'])

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

        nReps = 5
        threshold = -20. # spike threshold, mV

        stimInfo = {'Morphology': self.infile, 'synapseConfig': synapseConfig,
                    'runDur': self.run_duration, 'pip_dur': self.pip_duration, 'pip_start': self.pip_start,
                    'run_duration': self.run_duration,
                    'Fs': self.Fs, 'F0': self.f0, 'dB': self.dB, 'RF': self.RF, 'SR': self.SR,
                    'cellType': self.cellType, 'modelType': self.modelType, 'nReps': nReps, 'threshold': threshold}

        preCell, postCell, synapse, self.electrodeSite = self.configureCell(hf, synapseConfig, stimInfo)

        # see if we need to save the cell state now.
        if make_ANIntialConditions:
            print 'getting initial conditions for AN'
            cellInit.getInitialConditionsState(hf, tdur=3000., 
                filename='an_neuronstate.dat', electrodeSite=self.electrodeSite)
            cellInit.testInitialConditions(hf, filename='an_neuronstate.dat',
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
            for j, N in enumerate(range(nReps)):
                print 'Rep: %d' % N
            
                preCell, postCell = self.singleANRun(hf, j, synapseConfig, stimInfo, seeds, preCell, postCell)

                celltime.append(self.time)
                spikeTimes[N] = pu.findspikes(self.time, self.Vsoma, threshold, t0=0., t1=hf.h.tstop, dt=1.0, mode='peak')
                inputSpikeTimes[N] = [preCell[i]._spiketrain for i in range(len(preCell))]
                somaVoltage[N] = np.array(self.Vsoma)

            total_elapsed_time = time.time() - self.start_time
    #        total_run_time = time.time() - run_time
            print "Total Elapsed Time = %8.2f min (%8.0fs)" % (total_elapsed_time/60., total_elapsed_time)
            print "Total Setup Time = %8.2f min (%8.0fs)" % (self.setup_time/60., self.setup_time)
            print "Total AN Calculation Time = %8.2f min (%8.0fs)" % (self.an_setup_time/60., self.an_setup_time)
            print "Total Neuron Run Time = %8.2f min (%8.0fs)" % (self.nrn_run_time/60., self.nrn_run_time)
        
            result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes, 
                'somaVoltage': somaVoltage, 'time': np.array(self.time)}
        
            self.analysis_filewriter(self.infile, result, tag='Syn%03d' % k)
        #self.plotAN(np.array(self.time), result['somaVoltage'], result['stimInfo'])


    def configureCell(self, hf, synapseConfig, stimInfo):
        sg = hf.sec_groups['soma']
        postCell = cells.Generic.create(soma=hf.get_section(list(sg)[0]))

        self.cd.channelValidate(hf, verify=False)

        preCell = []
        synapse = []
        # reconfigure syanpses to set the spont rate group
        for i, syn in enumerate(synapseConfig):
            if stimInfo['SR'] is None:  # use the one in the table
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=syn[2]))
            else:
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=stimInfo['SR']))  # override
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


    def singleANRun(self, hf, j, synapseConfig, stimInfo, seeds, preCell, postCell):
        cellInit.restoreInitialConditionsState(hf, electrodeSite=None, filename='an_neuronstate.dat')
        self.stim=[]
        # make independent inputs for each synapse
        an0_time = time.time()
        for i, syn in enumerate(synapseConfig):
            self.stim.append(sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=self.dB,
                                  ramp_duration=self.RF, pip_duration=self.pip_duration,
                                  pip_start=self.pip_start))
            nseed = seeds[j, i]
            preCell[i].set_sound_stim(self.stim[-1], seed=nseed)  # generate spike train, connect to terminal
        self.an_setup_time += (time.time() - an0_time)
        nrn_start = time.time()
        self.Vsoma = hf.h.Vector()
        self.time = hf.h.Vector()
        self.Vsoma.record(postCell.soma(0.5)._ref_v, sec=postCell.soma)
        self.time.record(hf.h._ref_t)
        print '...running '
        hf.h.finitialize()
        hf.h.tstop = self.run_duration*1000.
        hf.h.t = 0.
        hf.h.batch_save() # save nothing
        hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
        # old - don't d this! #hf.h.finitialize()
        # hf.h.run()
        self.nrn_run_time += (time.time() - nrn_start)
        print '...done'
        return(preCell, postCell)


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

        # result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes,
        #     'somaVoltage': somaVoltage, 'time': np.array(self.time)}

    def analysis_filewriter(self, filebase, result, tag=''):
        k = result.keys()
        requiredKeys = ['stimInfo', 'spikeTimes', 'inputSpikeTimes', 'somaVoltage', 'time']
        for rk in requiredKeys:
            assert rk in k

        stimInfo = result['stimInfo']
        fname = os.path.splitext(self.infile)[0]
        f = open('AN_Result_' + fname + '%s_N%03d_%03ddB_%06.1f_%2s' % (tag, stimInfo['nReps'],
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
        print 'section groups: ', self.hf.sec_groups.keys()
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
    model = ModelRun(sys.argv[1:])
    model.runModel() # then run the model
    if showCell:
        QtGui.QApplication.instance().exec_()
