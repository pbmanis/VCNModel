__author__ = 'pbmanis'

import sys
import os
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
from neuronvis import HocViewer


class Render():
    def __init__(self, hf):
        self.hf = hf
        self.section_colors={'axon': 'r', 'hillock': 'r', 'initialsegment': 'orange',
             'unmyelinatedaxon': 'yellow', 'myelinatedaxon': 'white', 'dendrite': 'white',
             'soma': 'blue',
            # terminals (calyx of Held):
             'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

    def get_hoc_file(self, infile):
        if self.hf.file_loaded is False:
            exit()
        if not HAVE_PG:
            return()
        self.section_list = self.hf.get_section_prefixes()
        self.hf.sec_groups.keys()
        if len(self.hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(self.hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.get_color_map(i)
        else: # single section name, assign colors to SectionList types:
            self.section_colors={'axon': 'r', 'hillock': 'r', 'initialsegment': 'orange',
             'unmyelinatedaxon': 'yellow', 'myelinatedaxon': 'white', 'dendrite': 'white',
            'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

        (v, e) = self.hf.get_geometry()
        self.clist = []

        for si in self.hf.sections: # self.section_list[s]:
            print dir(si)
            self.hf.h('access %s' % si)
            sr = self.hf.h.SectionRef()
            n1 = self.hf.h.cas().name()
            if sr.has_parent() == 1:
                x=sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])

    def render(self, mech, rendertype='cylinder'):
        render = HocViewer(self.hf.hr.h)

        if rendertype == 'line':
            line = render.draw_graph()
            line.set_group_colors(self.section_colors, mechanism=mech)
        if rendertype == 'surface':
            surface = render.draw_surface()
            #surface.set_group_colors(self.section_colors, alpha=0.35)
            surface.set_group_colors(self.section_colors, mechanism=mech)
        elif rendertype == 'cylinder':
            cylinder = render.draw_cylinders()
            print 'setting cylinder colors for mech'
            cylinder.set_group_colors(self.section_colors,  mechanism=mech)
            print 'done coloring'
        elif rendertype == 'volume':
#            volume = render.draw_volume(resolution = 1.0, max_size=1e9)
            volume = render.draw_volume()
            #volume.set_group_colors(section_colors, alpha=0.35)



if __name__ == '__main__':
    
    pg.mkQApp()
    fn = sys.argv[1]
    filename = os.path.join('VCN_Cells', fn, 'Morphology', fn+'.hoc')
    post_cell = cells.Bushy.create(morphology=filename, decorator=Decorator,
            species='mouse',
            modelType='XM13')
    R = Render(post_cell)
    R.render(['nav11', 'gbar'], rendertype = 'surface')
    pg.show()
    pg.Qt.QtGui.QApplication.exec_()
    
    
    