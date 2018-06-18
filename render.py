#!/usr/bin/python

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
             'unmyelinatedaxon': 'yellow', 'myelinatedaxon': 'white', 'dendrite': 'green',
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
        # if len(self.hf.sec_groups) > 5: # multiple names, so assign colors to structure type
        #     self.section_colors = {}
        #     for i, s in enumerate(self.hf.sec_groups.keys()):
        #         self.section_colors[s] = self.hg.get_color_map(i)
        # else: # single section name, assign colors to SectionList types:
        #     self.section_colors={'axon': 'r', 'hillock': 'r', 'initialsegment': 'orange',
        #      'unmyelinatedaxon': 'yellow', 'myelinatedaxon': 'white', 'dendrite': 'white',
        #     'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
        #     'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

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

    def render(self, rendertype='cylinder', mech=None, colormap='viridis', backgroundcolor='w'): #'magma'):
        renderer = HocViewer(self.hf.hr.h)
        renderer.setBackcolor(backgroundcolor)
        if rendertype == 'line':
            g = renderer.draw_graph()
            g.set_group_colors(self.section_colors, mechanism=mech, colormap=colormap)
        elif rendertype == 'surface':
            g = renderer.draw_surface()
            g.set_group_colors(self.section_colors, mechanism=mech, colormap=colormap)
        elif rendertype == 'cylinder':
            g = renderer.draw_cylinders()
            g.set_group_colors(self.section_colors,  mechanism=mech, colormap=colormap)
        elif rendertype == 'volume':
            
#            volume = render.draw_volume(resolution = 1.0, max_size=1e9)
            g = renderer.draw_volume()
            g.set_group_colors(self.section_colors, mechanism=mech, alpha=1, colormap=colormap)
        else:
            print 'rendertype: ', rendertype
            raise ValueError('Render type %s not known: ' % rendertype)
        return g, renderer 


if __name__ == '__main__':
    
#    pg.mkQApp()
    rendertype = 'surface'
    mechanisms = ['klt', 'gbar']
    #mechanisms = None
    fn = sys.argv[1]
    filename = os.path.join('VCN_Cells', fn, 'Morphology', fn+'.hoc')
    post_cell = cells.Bushy.create(morphology=filename, decorator=Decorator,
            species='mouse',
            modelType='mGBC')
    #post_cell.distances()
    R = Render(post_cell)
    backgroundcolor = 'w'
    if rendertype in ['volume', 'line']:
        backgroundcolor = (125, 125, 125, 255)
    g, renderer = R.render(rendertype=rendertype, mech=mechanisms, backgroundcolor=backgroundcolor)
#    pg.show()
    if not sys.flags.interactive:
        pg.Qt.QtGui.QApplication.exec_()
    
    
    