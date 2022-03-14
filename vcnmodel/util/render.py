__author__ = "pbmanis"
import argparse
import importlib
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union, Tuple

from cnmodel import cells
from cnmodel.decorator import Decorator

import pyqtgraph as pg
from neuronvis.hoc_viewer import HocViewer
from pylibrary.tools import cprint as CP

"""
States sets up the display orientation for each cell 
for an initial view 
These values can be obtained by running
render celln Full vispy
when you rotate or move the cell, the new dict of the position
state will be printed. Just copy it into this dict to save it.
"""

states = {
    2: {
        "scale_factor": 148.5,
        "center": (-25.0, -9.38, -25.2),
        "fov": 45.0,
        "elevation": 13.0,
        "azimuth": -121.5,
        "roll": 0.0,
    },
    5: {
        "scale_factor": 86.58,
        "center": (41.95, 6.67, -26.02),
        "fov": 45.0,
        "elevation": -26.0,
        "azimuth": -123.5,
        "roll": 0.0,
    },
    6: {
        "scale_factor": 138.12,
        "center": (20.57, -7.60, -29.2),
        "fov": 45.0,
        "elevation": 73.0,
        "azimuth": -136.0,
        "roll": 0.0,
    },
    9: {
        "scale_factor": 131.56,
        "center": (39.42, -35.72, -30.19),
        "fov": 45.0,
        "elevation": 85.5,
        # "elevation": -6.5,
        # "elevation": -17.5,
        "azimuth": 131.5,
        # "azimuth": -119.75,
        # "azimuth": -102.25,
        "roll": 0.0,
    },
    10: {
        "scale_factor": 143.71,
        "center": (10.88, 24.80, 0.116),
        "fov": 45.0,
        "elevation": -20.6,
        "azimuth": 20.37,
        "roll": 0.0,
    },
    11: {
        "scale_factor": 128.18,
        "center": (-15.98, -0.060, -1.66),
        "fov": 45.0,
        "elevation": 22.0,
        "azimuth": -59.5,
        "roll": 0.0,
    },
    13: {
        "scale_factor": 136.16,
        "center": (59.741, 1.099, -20.94),
        "fov": 45.0,
        "elevation": -7.5,
        "azimuth": -105.5,
        "roll": 0.0,
    },
    17: {
        "scale_factor": 126.57,
        "center": (-15.88, -9.56, -2.81),
        "fov": 45.0,
        "elevation": -30.0,
        "azimuth": -45.0,
        "roll": 0.0,
    },
    18: {
        "scale_factor": 116.72,
        "center": (-27.21, 17.54, -1.56),
        "fov": 45.0,
        "elevation": 40.5,
        "azimuth": -20.0,
        "roll": 0.0,
    },
    30: {
        "scale_factor": 201.58,
        "center": (16.90, -17.64, -33.48),
        "fov": 45.0,
        "elevation": 70.0,
        "azimuth": -128.0,
        "roll": 0.0,
    },
}
sbem_sectypes = {
    # new swc mapping
    0: "Undefined",
    1: "soma",
    2: "Myelinated_Axon",
    3: "Basal_Dendrite",
    4: "Apical_Dendrite",
    5: "Custom",
    6: "Unspecified_Neurites",
    7: "Glia_Processes",
    8: "Blank",
    9: "Blank",
    10: "Axon_Hillock",
    11: "Unmyelinated_Axon",
    12: "Dendritic_Hub",
    13: "Proximal_Dendrite",
    14: "Distal_Dendrite",
    15: "Axon_Initial_Segment",
    16: "Axon_Heminode",
    17: "Axon_Node",
    18: "Dendritic_Swelling",
}


class Render:
    def __init__(self, cell_number:int, hf:object, renderer:str):
        self.cell_number = int(cell_number)
        self.renderer = renderer
        self.hf = hf
        self.section_colors = { # colors are xkcd color palette
            "axon": "spring green",
            "hillock": "red",
            "initialsegment": "lilac",
            "unmyelinatedaxon": "lilac",
            "myelinatedaxon": "spring green",
            "dendrite": "green",
            "soma": "black",
            "Axon_Hillock": "red",
            "Axon_Initial_Segment": "baby blue",
            "Myelinated_Axon": "dark red", # "spring green",
            "Unmyelinated_Axon": "lilac",
            "Proximal_Dendrite": "medium purple", # "cyan",
            "Dendritic_Hub": "royal blue",
            "Dendritic_Swelling": "gold",
            "Distal_Dendrite": "dark magenta", # "sky blue",
            # terminals (calyx of Held):
            "heminode": "green",
            "stalk": "yellow",
            "branch": "blue",
            "neck": "brown",
            "swelling": "magenta",
            "tip": "powderblue",
            "parentaxon": "orange",
            "synapse": "k",
        }

    def get_hoc_file(self)->None:
        if self.hf.file_loaded is False:
            exit()
        self.section_list = self.hf.get_section_prefixes()
        self.hf.sec_groups.keys()

        (v, e) = self.hf.get_geometry()
        self.clist = []

        for si in self.hf.sections:
            self.hf.h("access %s" % si)
            sr = self.hf.h.SectionRef()
            n1 = self.hf.h.cas().name()
            if sr.has_parent() == 1:
                x = sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])

    def render(
        self, view:str="cylinders", mechanism:Union[list, None]=None,
            colormap:Union[str, dict]="viridis",
            backgroundcolor:Union[str, list]="k",
            )-> Tuple[object, object]:
        viewer = HocViewer(self.hf.hr.h, renderer=self.renderer)
        if self.renderer not in ["mpl", "vispy"]:
            viewer.setBackcolor(backgroundcolor)
        print('render: renderer: ', self.renderer, view)
        if self.renderer == 'pyqtgraph':
            
            if view in ["line", "graph"]:
                g = viewer.draw_graph()
                g.set_group_colors(self.section_colors, mechanism=mechanism, colormap=colormap)
            elif view == "surface":
                g = viewer.draw_surface()
                g.set_group_colors(self.section_colors, mechanism=mechanism, colormap=colormap)
            elif view == "cylinders":
                g = viewer.draw_cylinders()
                g.set_group_colors(self.section_colors, mechanism=mechanism, colormap=colormap)
            elif view == "volume":
                #            volume = render.draw_volume(resolution = 1.0, max_size=1e9)
                g = viewer.draw_volume()
                g.set_group_colors(
                    self.section_colors, mechanism=mechanism, alpha=1, colormap=colormap
                )
        elif self.renderer == "mpl":
            g = viewer.draw_mpl()
        elif self.renderer == "vispy":
            g = viewer.draw_vispy(
                mechanism=mechanism,
                color=self.section_colors,
                state=states[self.cell_number],
            )

        else:
            raise ValueError("Render type %s not known: " % self.renderer)
        return g, viewer


def set_table_and_cells(
    filename: str,
    dataTable: str,
    species: str,
    modelName: str,
    modelType: str,
    nach: str,
    dendriteMode: str,
) -> object:
    from cnmodel import data

    dmodes = {
        "normal": "",
        "passive": "_pasdend",
        "active": "_actdend",
        "allpassive": "_allpassive",
    }
    changes = None
    nach = None  # uses default
    if dataTable == "":
        table_name = (
            f"vcnmodel.model_data.data_{modelName:s}{dmodes[dendriteMode]:s}"
        )
    else:
        table_name = f"vcnmodel.model_data.{dataTable:s}"
        CP.cprint("r", f"**** USING SPECIFIED DATA TABLE: {str(table_name):s}")
        knownmodes = ["normal", "actdend", "pasdend"]
        dendriteMode = "normal"
        for mode in knownmodes:
            if table_name.find(mode) > 0:
                dendriteMode = mode

        CP.cprint("c", f"Dendrite mode: {dendriteMode:s}")
    name_parts = modelName.split("_")
    if len(name_parts) > 1:
        nach = name_parts[1]
    else:
        nach = "nav11"
    CHAN = importlib.import_module(table_name)
    channels = f"{name_parts[0]:s}_{nach:s}_channels"
    compartments = f"{name_parts[0]:s}_{nach:s}_channels_compartments"
    print("Channels: ", channels)
    print("Compartments: ", compartments)
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
    else:
        CP.cprint("g", "No changes to data tables")

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
    parser = argparse.ArgumentParser(description="VCN Cell Morphology Renderer")
    parser.add_argument(
        "-n",
        type=int,
        default=0,
        dest="cellnumber",
        help="Cell by Number (2, 6, 11, etc), looks for VCN_cnn/Morphology/VCN_cnn_Full_MeshInflated.hoc",
    )
    parser.add_argument(
        "-f",
        dest="filename",
        type=str,
        default="None",
        help="Full Filename/path for rendering",
    )
    parser.add_argument(
        "-p",
        "--parts",
        dest="parts",
        type=str,
        default="Full",
        choices=["Full", "NoDend", ],
        help="select hoc files with parts removed",
    )
    parser.add_argument(
        "-r",
        "--renderer",
        type=str,
        dest="renderer",
        default="vispy",
        choices=["pyqtgraph", "mpl", "vispy"],
        help="Select render pipeline - pyqtgraph, matplotlib (mpl) or vispy",
    )
    parser.add_argument(
        "-v",
        "--view",
        type=str,
        default="cylinders",
        choices=["line", "graph", "cylinders", "volume"],
        help="view representation type - usually the cylinder default is best; others only work with pyqtgraph",
    )
    parser.add_argument(
        "-m",
        dest="mechanism",
        type=str,
        default="None",
        choices=["None", "klt", "kht", "ihvcn", "nacncoop", "ka", "kif", "kis"],
        help="Mechanism to render density",
    )
    
    args = parser.parse_args()
    if args.cellnumber > 0:
        hocfile = f"VCN_c{int(args.cellnumber):02d}_Full_MeshInflate.hoc"
        cell_dir = f"VCN_c{int(args.cellnumber):02d}"
        filename = Path("../VCN-SBEM-Data", "VCN_Cells", cell_dir, "Morphology", hocfile)
        fnum = args.cellnumber
    elif args.filename != "None":
        filename = Path(args.filename)
        fnum = filename.stem
        fnum = int(fnum[5:7])
    else:
        raise ValueError("Need either -n cell number or -f filename")

    post_cell = set_table_and_cells(
        filename=filename,
        dataTable="data_XM13A_nacncoop_normal",
        species="mouse",
        modelName="XM13A_nacncoop",
        modelType="II",
        nach="nacncoop",
        dendriteMode="normal",
    )
    # post_cell.distances()
    R = Render(fnum, post_cell, args.renderer)
    backgroundcolor = (32, 32, 32, 128)  # 'blue'
    if args.view in ["volume", "line", "graph"]:
        backgroundcolor = (125, 125, 125, 255)
    g, viewer = R.render(
        view=args.view, mechanism=args.mechanism, backgroundcolor=backgroundcolor
    )
    if R.renderer in ["mpl"]:
        return None
    elif R.renderer in ["vispy"]:  
        # vispy and pyqtgraph use qt
        viewer.close()
    elif R.renderer in ["pyqtgraph"]:
       if sys.flags.interactive == 0:
            pg.Qt.QtGui.QApplication.exec_()
        


if __name__ == "__main__":
    g = main()
