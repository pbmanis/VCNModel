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
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
import numpy as np

class ModelRun():
    def __init__(self, args):
        if isinstance(args, list):
            celltype = args[0] # must be string, not list...
        else:
            celltype = args
        if celltype not in ['Bushy', 'Stellate', 'L23pyr']:
            print 'Celltype must be one of Bushy, Stellate or L23pyr'
            exit()
        modelType = 'RM03'
        if len(args) > 1:
            modelType = args[1]
        if modelType not in ['RM03', 'XM13', 'MS']:
            print 'Model type mist be one of RM03, XM13, or MS'
            exit()

        #infile = 'L23pyr.hoc'
        #infile = 'LC_nmscaled_cleaned.hoc'
        #infile = 'Calyx-68cvt2.hoc'
        #infile = 'Calyx-S53Acvt3.hoc'
        #infile = 'wholeThing_cleaned.hoc'
        infile = 'somaOnly.hoc'
        self.infile = infile
        self.hf = HocReader('MorphologyFiles/' + self.infile)
        cellType = celltype # 'Bushy_RM03' # possibly this should come from the morphology file itself...
        electrodeSection = 'soma[0]'
        self.hg = HocGraphic(self.hf)
        self.get_hoc_file(self.infile)
        cd = ChannelDecorate(self.hf, celltype=cellType, modeltype = modelType)
        self.distances(electrodeSection) # make distance map from electrode site

       # self.render(['klt', 'gbar'])
        self.R = GenerateRun(self.hf, celltype=cellType, electrodeSection=electrodeSection, cd=cd)
        self.R.doRun(self.infile)
        #cd.channelValidate(self.hf)


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
#        render.hr.read_hoc_section_lists(self.section_colors.keys())
        #print 'sec groups: ', render.hr.sec_groups

        surface = render.draw_surface()
       # line = render.draw_graph()
        surface.set_group_colors(self.section_colors, alpha=0.35)
        surface.set_group_colors(self.section_colors, alpha=0.35, mechanism=mech)


if __name__ == "__main__":
    ModelRun(sys.argv[1:])
    #QtGui.QApplication.instance().exec_()
