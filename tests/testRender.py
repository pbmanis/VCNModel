"""
Render the hoc file using pyqtgraph as color cylinders by 
section type
"""

from pathlib import Path

import neuronvis.hoc_graphics as hoc_graphics
import neuronvis.sim_result as sr
import numpy as np
from vcnmodel.util.get_data_paths import get_data_paths
from neuronvis.hoc_reader import HocReader
from neuronvis.hoc_viewer import HocViewer
import pyqtgraph as pg

config = get_data_paths()

infile = 'VCN_c02_Full_MeshInflate.hoc'

hf = HocReader(Path(config["disk"], config["cellDataDirectory"], 'VCN_c02/Morphology', infile))
if hf.file_loaded is False:
    exit()

hf.h.topology()
# exit()
# print 'hf sections: ', hf.sections
named_section_colors = {'axon': 'r', 'AXON_0': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
    'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k',
    'soma': 'blue', 'dend': 'g', 'dendscaled_0': 'g', 'dendscaled_2': 'y', 'sections': 'b', 'hillock': 'magenta',
    'initseg': 'orange', 'undefined': 'k', 'dendrite': 'g', 'apical_dendrite': 'y',
    'Axon_Hillock': 'magenta', 'Axon_Initial_Segment': 'orange', 'Dendritic_Hub': 'red', 'Dendritic_Swelling': 'grey',
    'Proximal_Dendrite': 'green', 'Distal_Dendrite': 'yellow'}
section_colors = {}
if len(hf.sec_groups) > 1: # multiple names, so assign colors to structure type
    section_colors = {}
    for i, s in enumerate(hf.sec_groups.keys()):
        section_colors[s] = named_section_colors[s]
        print('section type %s assigned color %s' % (s, named_section_colors[s]))
        # else:
        #     icolor = i % 12
        #     section_colors[s] = hoc_graphics.colorMap[icolor]
else: # single section name, assign colors to SectionList types:
    # here we should find out the names of the section lists in hoc
    section_colors = named_section_colors
print('section colors:', section_colors)

(v, e) = hf.get_geometry()
clist = []

for s in hf.sections.keys(): # all of the sections in the model
    hf.h('access %s' % s) # access the right section

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

type = 'cylinder'
if type == 'surface':
    surface = render.draw_surface(resolution = 8.0)
    surface.set_group_colors(section_colors, alpha=0.7)
elif type == 'cylinder':
    cylinder=render.draw_cylinders()
    cylinder.set_group_colors(section_colors, alpha=0.8)
elif type == 'volume':
    volume = render.draw_volume(resolution = 10.0, max_size=1e9)

render.setCameraPosition(distance=4500.*0.110, elevation=60., azimuth=135.)
render.pan(dx=1350.*0.110, dy=-750.*(0.110), dz=0.0)

pos = render.cameraPosition()
center = render.opts['center']
dist = render.opts['distance']
elev = render.opts['elevation']
azim = render.opts['azimuth']

print('Pos (center, dist, elev, azim)', center, dist, elev, azim)

# for i in range(0,3):
#     print 'Soma section [%d] diameter = %6.1f microns' % (i, hf.h('soma[%d].diam' % i))
#     print 'Soma section [%d] length = %6.1f microns' % (i, hf.h('soma[%d].L' % i))

# print 'total soma: ', hf.h('soma.diam'), hf.h('soma.L')
# print 'AXON_0 diameter = %6.1f microns' % hf.h('AXON_0[0].diam')
# print 'AXON_0 length = %6.1f microns' % hf.h('AXON_0[0].L')

pg.QtWidgets.QApplication.instance().exec()
