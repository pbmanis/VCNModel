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
import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer
from neuronvis.hoc_reader import HocReader
from neuronvis.hoc_graphics import HocGraphic
from channel_decorate import ChannelDecorate
from generate_run import GenerateRun
import time
try:
    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui
    import pylibrary.pyqtgraphPlotHelpers as pgh
    HAVE_PG = True
except:
	HAVE_PG = False

import numpy as np

verbose = False

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

        self.plotFlag = False
        #infile = 'L23pyr.hoc'
        #infile = 'LC_nmscaled_cleaned.hoc'
        #infile = 'Calyx-68cvt2.hoc'
        #infile = 'Calyx-S53Acvt3.hoc'
        infile = 'wholeThing_cleaned.hoc'
        #infile = 'somaOnly.hoc'
        self.infile = infile


    def runModel(self, parMap=None):
        if verbose:
            print 'runModel entry'
        self.hf = HocReader('MorphologyFiles/' + self.infile)
        self.electrodeSection = 'soma[0]'
        self.hg = HocGraphic(self.hf)
        self.get_hoc_file(self.infile)
        self.distances(self.electrodeSection) # make distance map from electrode site
        if verbose:
            print 'Parmap in Runmodel: ', parMap
        cd = ChannelDecorate(self.hf, celltype=self.cellType, modeltype=self.modelType, parMap=parMap)

       # self.render(['nav11', 'gbar'])
       # QtGui.QApplication.instance().exec_()
        if verbose:
            print 'generateRun'
        self.R = GenerateRun(self.hf, celltype=self.cellType,
                             electrodeSection=self.electrodeSection, cd=cd,
                             plotting = HAVE_PG and self.plotFlag)
        if verbose:
            print 'doRun'
        self.R.doRun(self.infile)
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

        taum_mean = np.mean([self.R.IVResult['taus'][i]['tau'].value for k, i in enumerate(self.R.IVResult['taus'].keys())])
        tauih_mean = np.mean([self.R.IVResult['tauih'][i]['tau'].value for k, i in  enumerate(self.R.IVResult['tauih'].keys())])
        #print 'taum_mean: ', taum_mean
        #print 'tauih_mean: ', tauih_mean
        # construct dictionary for return results:
        self.IVSummary = {'par': parMap, 'Vm': np.mean(self.R.IVResult['Vm']),
                          'Rin': self.R.IVResult['Rinss'],
                          'taum': taum_mean, 'tauih': tauih_mean,
                          'spikes': {'i': self.R.IVResult['Ispike'], 'n': self.R.IVResult['Nspike']},
                          }

        #print 'ivsummary: ', self.IVSummary
        return self.IVSummary
        #cd.channelValidate(self.hf)


    def set_celltype(self, celltype):
        self.cellType = celltype
        if self.cellType not in ['Bushy', 'Stellate', 'L23pyr']:
            print 'Celltype must be one of Bushy, Stellate or L23pyr, got: %s', self.cellType
            exit()


    def set_modeltype(self, modeltype):
        self.modelType = modeltype
        if self.modelType not in ['RM03', 'XM13', 'MS']:
            print 'Model type mist be one of RM03, XM13, or MS, got: ' % (self.modelType)
            exit()


    def get_hoc_file(self, infile):
        if self.hf.file_loaded is False:
            exit()
        self.section_list = self.hf.get_section_prefixes()
        self.hf.sec_groups.keys()
        if len(self.hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(self.hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.get_color_map(i)
        else: # single section name, assign colors to SectionList types:
            self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
                'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

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


    def distances(self, section):
        self.hf.distanceMap = {}
        self.hf.h('access %s' % section) # reference point
        d = self.hf.h.distance()
        for si in self.hf.sections.keys():
            self.hf.h('access %s' % si)
            self.hf.distanceMap[si] = self.hf.h.distance(0.5) # should be distance from first point


    def render(self, mech):
        pg.mkQApp()
        render = HocViewer(self.hf)

        type = 'line'
        if type == 'line':
            line = render.draw_graph()
            line.set_group_colors(self.section_colors, alpha=0.35)
        if type == 'surface':
            surface = render.draw_surface(resolution = 1.0)
            surface.set_group_colors(self.section_colors, alpha=0.35)
        elif type == 'cylinder':
            cylinder=render.draw_cylinders()
            cylinder.set_group_colors(self.section_colors, alpha=0.35)
        elif type == 'volume':
            volume = render.draw_volume(resolution = 1.0, max_size=1e9)
            #volume.set_group_colors(section_colors, alpha=0.35)

#        render.hr.read_hoc_section_lists(self.section_colors.keys())
        #surface = render.draw_surface()
        # cylinder = render.draw_cylinders()
        # cylinder.set_group_colors(self.section_colors, alpha=0.35, mechanism=mech)
        # line = render.draw_graph()
       # surface.set_group_colors(self.section_colors, alpha=0.35)
       # surface.set_group_colors(self.section_colors, alpha=0.35, mechanism=mech)


if __name__ == "__main__":
    model = ModelRun(sys.argv[1:])
    model.runModel() # then run the model
    #QtGui.QApplication.instance().exec_()
