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
#import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer
from neuronvis.hoc_reader import HocReader
import neuronvis.hoc_graphics as hoc_graphics
from channel_decorate import ChannelDecorate
from generate_run import GenerateRun
from nrnlibrary.protocols import protocol
from nrnlibrary import cells
from nrnlibrary import synapses
from nrnlibrary.util import get_anspikes
from nrnlibrary.util import sound
import nrnlibrary.util as nu
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

ANPSTH_mode = False
IV_mode = True

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
                # print 'section: ', section
                # print 'Ra: ', self.hf.get_section(section).Ra

        self.electrodeSection = 'soma'
        self.hg = hoc_graphics
        self.get_hoc_file(filename)
        sg = self.hf.sec_groups['soma']
        self.distances(self.hf.get_section(list(sg)[0]).name()) # make distance map from soma
        if verbose:
            print 'Parmap in Runmodel: ', parMap
        self.cd = ChannelDecorate(self.hf, celltype=self.cellType, modeltype=self.modelType,
                             parMap=parMap)
#        self.cd.channelValidate(self.hf, verify=False)
        self.hf.h.topology()
        if ANPSTH_mode:
            self.ANinputs(self.hf)
        elif IV_mode:
            self.IVRun(parMap)

        if showCell:
            self.render = HocViewer(self.hf)
            cylinder=self.render.draw_cylinders()
            cylinder.set_group_colors(self.section_colors, alpha=0.8, mechanism=['nav11', 'gbar'])
        
    def IVRun(self, parMap={}):
        if verbose:
            print 'generateRun'
        self.R = GenerateRun(self.hf, idnum=self.idnum, celltype=self.cellType,
                             starttime=self.startTime,
                             electrodeSection=self.electrodeSection, cd=self.cd,
                             plotting = HAVE_PG and self.plotFlag)

        if verbose:
            print 'doRun'
        self.R.doRun(self.infile, parMap, save='monitor')
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
        #self.cd.channelValidate(self.hf)


    def set_celltype(self, celltype):
        self.cellType = celltype
        if self.cellType not in ['Bushy', 'Stellate', 'L23pyr']:
            print 'Celltype must be one of Bushy, Stellate or L23pyr, got: %s', self.cellType
            exit()


    def set_modeltype(self, modeltype):
        self.modelType = modeltype
        if self.modelType not in ['RM03', 'XM13', 'MS']:
            print 'Model type must be one of RM03, XM13, or MS, got: %s ' % (self.modelType)
            exit()


    def set_starttime(self, starttime):
        self.startTime = starttime


    def ANinputs(self, hf, verify=False, seed=0, synapseConfig=[(int(216.66*0.65), 0., 2), (int(122.16*0.65), 0., 2), 
    (int(46.865*0.65), 0., 2), (int(84.045*0.65), 0., 2), (int(2.135*0.65), 0, 2), (int(3.675*0.65), 0, 2), (int(80.27*0.65), 0, 2)]):
        """
        Establish AN inputs to soma
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        """
        self.run_duration = 0.25
        self.pip_duration = 0.1
        self.pip_start = [0.1]
        self.Fs = 100e3
        self.f0 = 1000.
        nReps = 1
        win = pgh.figure(title='AN Inputs')
        layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')

        # hf.h.tstop = 100.
        # hf.h.finitialize()
        # hf.h.run()
        #synapseConfig = [synapseConfig[2]]
        sg = hf.sec_groups['soma']
        dend = hf.sec_groups['dendrite']
        dendrite = hf.get_section(list(dend)[99])
        postCell = cells.Generic.create(soma=hf.get_section(list(sg)[0]))
        #postCell = cells.Bushy.create(nach='nav11')
#        print 'postcell.cell: ', postCell.all_sections
        
        preCell = []
        synapse = []
        vsoma = []
        time = []
        vdend = []
        for i, syn in enumerate(synapseConfig):
            preCell.append(cells.DummySGC(cf=self.f0, sr=syn[2]))
            synapse.append(preCell[-1].connect(postCell, pre_opts={'nzones':syn[0], 'delay':syn[1]}))
        self.pre_cell = preCell
        self.post_cell = postCell
        self.synapse = synapse
        for i, s in enumerate(synapse):
            s.terminal.relsite.Dep_Flag = 0  # turn off depression computation
            print 'Depression flag for synapse %d = %d ' % (i, s.terminal.relsite.Dep_Flag)

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
        
        seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))
        print seeds
        k = 0
        for j, N in enumerate(range(nReps)):
            print 'Rep: %d' % N
            self.stim=[]
            # make independent inputs for each synapse
            for i, syn in enumerate(synapseConfig):
                self.stim.append(sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=90,
                                      ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                                      pip_start=self.pip_start))
                nseed = seeds[j, i]
                preCell[i].set_sound_stim(self.stim[-1], seed=nseed)  # generate spike train, connect to terminal

            self.Vsoma = self.hf.h.Vector()
            self.time = self.hf.h.Vector()
            self.Vdend = self.hf.h.Vector()
            self.Vsoma.record(postCell.soma(0.5)._ref_v, sec=postCell.soma)
            self.Vdend.record(dendrite(0.5)._ref_v, sec=dendrite)
            self.time.record(self.hf.h._ref_t)
            print '...running '
            hf.h.tstop = self.run_duration*1000.
            hf.h.finitialize()
            hf.h.run()
            print '...done'
            vsoma.append(self.Vsoma)
            time.append(self.time)
            vdend.append(self.Vdend)
            dvdt = np.diff(np.array(vsoma[N]))
            sp = np.where(dvdt > 2.0)
            layout.plot(0, np.array(time[N]), np.array(vsoma[N]), pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.plot(1, np.array(time[N])[:-1], dvdt, pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.get_plot(0).setXLink(layout.get_plot(1))
        pgh.show()

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
