__author__ = 'pbmanis'

import sys
import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer
from neuronvis.hoc_reader import HocReader
from neuronvis.hoc_graphics import HocGraphic
import channel_decorate as cd
from generate_run import GenerateRun
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg

class ModelRun():
    def __init__(self, args):
        infile = 'LC_neuromantic_scaled.hoc'
        infile = 'Calyx-S53Acvt3.hoc'
        #infile = 'mainDenHOC_cleaned.hoc'
        #infile = 'somaOnly.hoc'
        self.hf = HocReader('MorphologyFiles/' + infile)
        self.hg = HocGraphic(self.hf)
        self.get_hoc_file(infile)
        cd.ChannelDecorate('Bushy', self.hf, self.section_list)

        #self.render(['ih', 'ghbar'])
        self.R = GenerateRun(self.hf)
        self.R.doRun()


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
