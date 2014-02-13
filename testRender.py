__author__ = 'pbmanis'

import numpy as np
import neuronvis.sim_result as sr
from neuronvis.hoc_viewer import HocViewer
from neuronvis.hoc_reader import HocReader
import neuronvis.hoc_graphics as hoc_graphics
from pyqtgraph.Qt import QtGui
import os

#infile = 'LC_neuromantic_scaled.hoc'
#infile = 'Calyx-S53Acvt3.hoc'
# infile = 'Calyx-68cvt4.hoc'
infile = 'mainDenHOC_cleaned.hoc'
hf = HocReader('MorphologyFiles/' + infile)
if hf.file_loaded is False:
    exit()
section_list = hf.get_sections()
if len(section_list) > 1: # multiple names, so assign colors to structure type
    section_colors = {}
    for i, s in enumerate(section_list.keys()):
        section_colors[s] = hoc_graphics.colorMap[i]
else: # single section name, assign colors to SectionList types:
    # here we should find out the names of the section lists in hoc

    section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
        'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}


(v, e) = hf.get_geometry()
clist = []

for s in section_list.keys(): # all of the sections in the model
    hf.h('access %s' % s) # access the right section
    for i in range(len(section_list[s])):
        si = '%s[%d]' % (s, i)
        sn = eval('hf.h.%s' % (si))
        hf.h('access %s' % si)
        sr = hf.h.SectionRef()
#        print 'sn.name: %s' % (sn.name())
        n1 = hf.h.cas().name()
        #print 'cas1: ', n1,' is connected to: ',
        if sr.has_parent() == 1:
            x=sr.parent
            n2 = x.name()
            clist.append([n1, n2])
            #print 'cas2: ', n2
        else:
            clist.append([n1, None])
#            print 'nothing'


render = HocViewer(hf)
render.hr.read_hoc_section_lists(section_colors.keys())
#print 'sec groups: ', render.hr.sec_groups

surface = render.draw_surface()
surface.set_group_colors(section_colors, alpha=0.35)

#print clist
QtGui.QApplication.instance().exec_()
