#!/usr/bin/python

from __future__ import print_function
__author__ = 'pbmanis'

import sys
import importlib
import os
from pathlib import Path
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
from neuronvis.hoc_viewer import HocViewer
from pylibrary.tools import cprint as CP


class Render():
    def __init__(self, hf):
        self.hf = hf
        self.section_colors={'axon': 'red', 'hillock': 'r', 'initialsegment': 'orange',
             'unmyelinatedaxon': 'cyan', 'myelinatedaxon': 'white', 'dendrite': 'green',
             'soma': 'blue',
             'Axon_Hillock': 'red', 'Axon_Initial_Segment': 'lilac', 'Myelinated_Axon': 'darkgreen',
             "Unmyelinated_Axon": 'brown',
             'Proximal_Dendrite': 'purple', 'Dendritic_Hub': 'green', 'Dendritic_Swelling': 'magenta',
             'Distal_Dendrite': 'powderblue',
            # terminals (calyx of Held):
             'heminode': 'green', 'stalk':'yellow', 'branch': 'blue', 'neck': 'brown',
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
        if rendertype in ['line', 'graph']:
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
        elif rendertype == 'mpl':
            g = renderer.draw_mpl()
        elif rendertype == 'vispy':
            g = renderer.draw_vispy()

        else:
            print('rendertype: ', rendertype)
            raise ValueError('Render type %s not known: ' % rendertype)
        return g, renderer 


    # post_cell = set_table_and_celle(
    #      morphology=str(filename),
    #      decorator=Decorator,
    #      species='mouse',
    #      modelName='XM13',
    #      modelType="II",
    #      nach="nacncoop")
         
def set_table_and_cells(
                       filename:str,
                       dataTable:str,
                       morphology:str,
                       decorator:object,
                       species:str,
                       modelName:str,
                       modelType:str,
                       nach:str,
                       dendriteMode:str) -> object:
    from cnmodel import data
    dmodes = {"normal": "", "passive": "_pasdend", "active": "_actdend", "allpassive": "_allpassive"}
    changes = None
    nach = None  # uses default
    if dataTable is "":
        table_name = f"vcnmodel.model_data.data_{modelName:s}{dmodes[dendriteMode]:s}"
    else:
        table_name = f"vcnmodel.model_data.{dataTable:s}"
        CP.cprint('r', f"**** USING SPECIFIED DATA TABLE: {str(table_name):s}")
        knownmodes = ["normal", "actdend", "pasdend"]
        dendriteMode = "normal"
        for mode in knownmodes:
            if table_name.find(mode) > 0:
                dendriteMode = mode

        CP.cprint('c', f"Dendrite mode: {dendriteMode:s}")
    name_parts = modelName.split('_')
    if len(name_parts) > 1:
        nach = name_parts[1]
    else:
        nach = 'nav11'
    CHAN = importlib.import_module(table_name)
    channels = f"{name_parts[0]:s}_{nach:s}_channels"
    compartments = f"{name_parts[0]:s}_{nach:s}_channels_compartments"
    print('Channels: ', channels)
    print('Compartments: ', compartments)
    changes = data.add_table_data(
        channels,
        row_key="field",
        col_key="model_type",
        species="mouse",
        data=CHAN.ChannelData,
    )
    changes_c = data.add_table_data(
        compartments,
        row_key="parameter",
        col_key="compartment",
        species="mouse",
        model_type="II",
        data=CHAN.ChannelCompartments,
    )
    if changes is not None:
        data.report_changes(changes)
        data.report_changes(changes_c)
    
    post_cell = cells.Bushy.create(
        morphology=str(filename),
        decorator=Decorator,
        species=species,
        modelName=modelName,
        modelType=modelType,
        nach=nach,
    )
    return post_cell

def main():
    rendertype = 'cylinder' #'surface'
    mechanisms = None # ['klt', 'gbar']

    #mechanisms = None
    fn = sys.argv[1]
    mod = 'Full'
    if len(sys.argv) > 2:
         mod = sys.argv[2]
    if mod == 'Full':
        cell_dir = f"VCN_c{int(fn):02d}"
    else:
        cell_dir = f"VCN_c{int(fn):02d}_{mod:s}"
    filename = Path('../VCN-SBEM-Data', 'VCN_Cells', cell_dir,
         'Morphology', f"{cell_dir:s}.hoc")
    post_cell = set_table_and_cells(
         filename=filename,
         dataTable='data_XM13A_nacncoop_normal',
         morphology=str(filename),
         decorator=Decorator,
         species='mouse',
         modelName='XM13A_nacncoop',
         modelType="II",
         nach="nacncoop",
         dendriteMode="normal")
    #post_cell.distances()
    if len(sys.argv) > 3:
        rendertype = sys.argv[3]
    R = Render(post_cell)
    backgroundcolor = (32, 32, 32, 128)  # 'blue'
    if rendertype in ['volume', 'line', 'graph']:
        backgroundcolor = (125, 125, 125, 255)
    print(backgroundcolor, rendertype)
    g, renderer = R.render(rendertype=rendertype, mech=mechanisms, backgroundcolor=backgroundcolor)

    if rendertype == 'mpl':
        return None
    else:
        if not sys.flags.interactive:
            pg.Qt.QtGui.QApplication.exec_()

if __name__ == '__main__':
    
    g = main()
#    pg.mkQApp()

    if not sys.flags.interactive:
        pg.Qt.QtGui.QApplication.exec_()
    
    
    