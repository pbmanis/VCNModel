import argparse
from dataclasses import dataclass, field
import datetime
from typing import Union, Type, Tuple
from pathlib import Path
import re
import toml
import pprint
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as mpl
import numpy as np
from neuron import h
from cnmodel import cells
from vcnmodel.adjust_areas import AdjustAreas
from vcnmodel import cell_config as cell_config
from neuronvis import swc_to_hoc, hocRender

pp = pprint.PrettyPrinter(indent=4)

"""
This module provides tools that take an original HOC file and replaces the axon sections
(including the AIS, hillock, and myelinated sections) with some variant of
a "standardized" axon. 
The modified version of the cell is written to a new file.
"""

def defemptylist():
    """
    For data classes, use this function to instantiate an empty list
    """
    return []

@dataclass(unsafe_hash = True)
class Point:
    """
    Defines a point in 3d space
    """

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

@dataclass
class Point3D:
    """
    Defines a pt3d for neuron, using point and diameter.
    All units are nominally in microns for NEURON
    """

    xyz: Point()
    d: float = 0.0

@dataclass
class List3DPoints:
    """
    dataclass to represent listed of pt3d
    structures in Neuron
    """

    x: list = field(default_factory=defemptylist)
    y: list = field(default_factory=defemptylist)
    z: list = field(default_factory=defemptylist)
    d: list = field(default_factory=defemptylist)

@dataclass
class AxonMorph:
    """
    Simple data class to summarize information about axon morphology
    for each of the axon parts.
    Length is the length of the part type
    Diam0 is the starting point diameter
    Diam1 is tne end point diameter
    DiamAvg is the average diameter
    """

    cellID: str = ""
    HillockLength: list = field(default_factory=defemptylist)
    HillockDiam0: list = field(default_factory=defemptylist)
    HillockDiam1: list = field(default_factory=defemptylist)
    HillockDiamAvg: list = field(default_factory=defemptylist)
    AISLength: list = field(default_factory=defemptylist)
    AISDiam0: list = field(default_factory=defemptylist)
    AISDiam1: list = field(default_factory=defemptylist)
    AISDiamAvg: list = field(default_factory=defemptylist)
    MyelLength: list = field(default_factory=defemptylist)
    MyelDiam0: list = field(default_factory=defemptylist)
    MyelDiam1: list = field(default_factory=defemptylist)
    MyelDiamAvg: list = field(default_factory=defemptylist)
    AxonLength: list = field(default_factory=defemptylist)
    AxonDiam0: list = field(default_factory=defemptylist)
    AxonDiam1: list = field(default_factory=defemptylist)
    AxonDiamAvg: list = field(default_factory=defemptylist)

@dataclass(frozen=True)
class StandardAxon:
    """
    Averaged axon
    format is a list of diam[0], diam[1], length, and # of sections
    all measurements are in microns
    The values in the table come from a run of get_axon_lengths on 7/19/2021,
    from the excel table where the averages are computed leaving out the "axon"
    for cell 21, for which only a short distance is visible.
    Measurement data is in VCN-SBEM-Data/MorphologyData/VCN/Axon_Measurements
    
    """
    #      dia1  dia2   len   # secs
    hil = [2.58, 1.89,  2.90,  3]
    ais = [1.46, 0.80, 16.93, 19]  
    axn = [1.47, 1.88, 55.47, 55]  # length excludes one with a very short length

@dataclass
class ModifiableAxon:
    """
    Averaged axon, but NOT frozen so it can be modified.
    format is a list of diam[0], diam[1], length, and # of sections
    all measurements are in microns
    The values in the table come from a run of get_axon_lengths on 7/19/2021,
    from the excel table where the averages are computed leaving out the "axon"
    for cell 21, for which only a short distance is visible.
    Measurement data is in VCN-SBEM-Data/MorphologyData/VCN/Axon_Measurements
    
    """
    #      dia1  dia2   len   # secs
    hil = [2.58, 1.89,  2.90,  3]
    ais = [1.46, 0.80, 16.93, 19]  
    axn = [1.47, 1.88, 55.47, 55]  # length excludes one with a very short length


"""
List of axons : includes all 10 gradeA Cells, plus axon
    reconstructions from other bushy cells
"""
gradeACells = [2, 5, 9, 10, 11, 13, 17, 18, 30]
additional_axons = [12, 14, 15, 16, 19, 20, 21, 22, 23, 24, 27, 29]
all_cells = gradeACells + additional_axons

"""
Regular expression definitions for parts of HOC file
We use these to find and change parts of the HOC file automatically
"""

re_find_create_sections = re.compile(
    r"(?P<section>create sections\[(?P<secno>[\d]+)\])", re.DOTALL
)
re_find_create_AH = re.compile(r"(?P<section>objref Axon_Hillock)", re.DOTALL)
re_find_create_AIS = re.compile(r"(?P<section>objref Axon_Initial_Segment)", re.DOTALL)
re_find_create_MA = re.compile(r"(?P<section>objref Myelinated_Axon)", re.DOTALL)

re_find_create_sections_position = re.compile(
    r"(?P<section>create sections)", re.DOTALL
)

re_find_accessed_sections = re.compile(
    r"(?P<section>access sections)\[(?P<secno>[\d]+)\]",  # "\\n(?P<source>Axon_Hillock\.append\(\))",
    re.DOTALL,
)
re_find_soma = re.compile(
    r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>soma\.append\(\))",
    re.DOTALL,
)

re_find_AH = re.compile(
    r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Axon_Hillock\.append\(\))",
    re.DOTALL,
)
re_find_AIS = re.compile(
    r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Axon_Initial_Segment\.append\(\))",
    re.DOTALL,
)
re_find_MA = re.compile(
    r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Myelinated_Axon\.append\(\))",
    re.DOTALL,
)

# assemble string - .VERBOSE not doing what I expect?
s_ah = r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Axon_Hillock\.append\(\))\s"
sa = r"(?P<connect0>connect sections\[[0-9]+\]\(0\)), (?P<connect1>sections\[[0-9]+\]\(1\))\s"
sb = r"(?P<sections>sections\[[0-9]+\]) {(?P<pt3ddata>[^}]+)"

re_find_AH_pt3d = re.compile(s_ah + sa + sb, re.DOTALL)

s_ais = r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Axon_Initial_Segment\.append\(\))\s"
re_find_AIS_pt3d = re.compile(s_ais + sa + sb, re.DOTALL)

s_ma = r"(?P<section>access sections)\[(?P<secno>[\d]+)\]\s(?P<source>Myelinated_Axon\.append\(\))\s"
re_find_MA_pt3d = re.compile(s_ma + sa + sb, re.DOTALL)

class MakeStandardAxon:
    """
    Create a substitute standardized axon for a cell.
    To do this we:

        1. read the original hoc file
        2. [calculate the origin and direction of the assembly of
            all parts of the axon, from hillock through the myleinated part.] # this is NOT done in the current version
            and has no effect on the electrical behavior of the model.
        3. make a standardized axon consisting of 3 sections:

            a. which has the average axon parameters (diam, length) for the hillock, AIS and myelinated part
            b. which goes in the direction of the original axon, by rotating around the hillock origin
                The standardized axon starts at the origin and is aligned along the X-axis (Y and Z are both zero).
                The standardized axon is made up of pt3d elements of 1 um in length.
            c. which matches with the cell by translating the origin to the original hillock 0-end.

        4. We then remove all of the axon sections from the oritingal hoc text, and starting with the first
            section that would have been the hillock, add the 3 new sections.
            We do this in 5 steps:

            a. Get the section number of the part of the soma that is connected to the hillock
            b. Check for all axon containing sections, noting those of
                axon_hillock, axon_initial_segment, and myelinated_axon which
                might be missing in the current reconstruction.

            c. Search and replace pt3d strings for each axon part,
                adding missing parts at the end of the list (and adjusting the
                create sections call to match)
            d. connect the prosethetic axon hillock end 0 to the soma section end 1

    The file is written back out to the appropriate Morhpology directory with the name
    VCN_cNN_Full_standardized_axon.hoc
    Modifications can be made to the "standardized_axon". One that is implemented is to change
    the length of the AIS.
    The comment at the top if the file is modified by insertion of a line:
    // Modified to use standardized axon
    or 
    // Modified to use standardized axon with AIS length = xx.x um
    
    Check the newly created sections manually to be sure that the connect statements
    link to the right parts, and that everything is formatted correctly.
    """

    def __init__(self, revised: bool = False):
        self.cconfig = cell_config.CellConfig(verbose=False, spont_mapping="HS",)
        self.revised = revised  # use the revised version with a standardized axon axon
        # find out where our files live
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {
                "cellDataDirectory": Path("../VCN-SBEM-Data", "VCN_Cells")
            }
        self.baseDirectory = self.datapaths["cellDataDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"

        self.dendrites = [
            "maindend",
            "dend",
            "dendrite",
            "primarydendrite",
            "Proximal_Dendrite",
            "secdend",
            "secondarydendrite",
            "Distal_Dendrite",
            "Dendritic_Hub",
            "Dendritic_Swelling",
            "Basal_Dendrite",
            "Apical_Dendrite",
        ]

        self.somas = [
            "soma",
            "Soma",
        ]

        self.axons = [
            "Axon_Hillock",
            "Axon_Initial_Segment",
            "Unmeylinated_Axon",
            "Myelinated_Axon",
            "Axon",
        ]
    
    def do_axon(self, cell: str, AIS_length:float = 0.):
        cname, fn = self.make_name(cell, add_standardized_axon=self.revised, AIS_length=AIS_length)
        self.get_axon_measures(cname, fn, AIS_length=AIS_length)
        

    def make_name(
        self, cell: str, add_standardized_axon: bool = False, AIS_length: float=0.,
    ) -> Tuple[str, Type[Path]]:
        """
        Make a full filename that points to the cell morphology for the "Full" file type
        and for the standardized axon if needed.
        
        The name is modified to add "standardized_axon" if this is true; if AIS_length
        is not 0, then the AIS length of the standardized axon is modified.
        """
        cname = f"VCN_c{cell:02d}"

        cell_hoc = f"{cname:s}_Full_MeshInflate"
        if add_standardized_axon and AIS_length == 0.0:
            cell_hoc += "_standardized_axon"
        if AIS_length > 0.0:
            cell_hoc += f"_AIS={AIS_length:06.2f}"
        cell_hoc += ".hoc"
        fn = Path(self.baseDirectory, cname, self.morphDirectory, cell_hoc)
        return cname, fn

    def convert_swc_hoc(
        self, cell: str,
        ) -> Tuple[str, Type[Path]]:
        """
        Wrapper to convert the axon-only swcs to hoc files
        Assumes specific filename structure, and generates the output
        as the "axononly meshinflate" hoc file. Should also render the axon.
        """
        cname = f"VCN_c{cell:02d}"
        if cell in additional_axons and cell not in gradeACells:
            cell_swc = f"{cname:s}_Axon_00000.swc"
            cell_hoc = f"{cname:s}_AxonOnly_MeshInflate.hoc"
            hocf = Path(self.baseDirectory, cname, self.morphDirectory, cell_hoc)
            swcf = Path(self.baseDirectory, cname, self.morphDirectory, cell_swc)
            if not hocf.is_file() and swcf.is_file():
                scales = {"x": 1.0, "y": 1.0, "z": 1.0, "r": 1.0, "soma": 1.0, "dend":1.0}
                if swcf.is_file():
                    s = swc_to_hoc.SWC(filename=swcf, secmap="sbem", scales=scales, verify=True)
                    s.topology = True
                    s.show_topology()
                    print(hocf)
                    s.write_hoc(hocf)
                    hocRender.Render(
                        hoc_file=hocf,
                        display_style='cylinders',
                        display_renderer='pyqtgraph',
                        display_mode='sec-type',
                        mechanism="None",
                        # alpha=args["alpha"],
 #                        verify=args["verify"],
 #                        sim_data=sim_data,
                        secmap="sbem",
                    )
            else:
                print("swc file not found: ", swcf)
                exit()

    def fill_coords(self, coords: list, sec: object, insertnan: bool = False) -> list:
        """
        Append pt3d coordinates to the list
        """
        if insertnan:
            coords.x.append(np.nan)
            coords.y.append(np.nan)
            coords.z.append(np.nan)
            coords.d.append(np.nan)
        for i in range(sec.n3d()):
            coords.x.append(sec.x3d(i))
            coords.y.append(sec.y3d(i))
            coords.z.append(sec.z3d(i))
            coords.d.append(sec.diam3d(i))
        return coords

    def get_axon_measures(self, cname: str, fn: Union[str, Path], AIS_length:float = 0.) -> AxonMorph:
        """
        Compute a number of measures of the axon sections from the
        reconstructions.

        This handles a single reconstruction.
        It also computes some factors needed to do attach the prosethetic axon.
        
        Parameters
        ----------
        cname : str
            cell name in format 'VCN_cnn'
        fn : str or path
            file name for hoc file
        """
        print("getaxonmeasures fn: ", fn)
        cconfig = cell_config.CellConfig()
        sinflateratio = cconfig.get_soma_ratio(cname)
        dinflateratio = cconfig.get_dendrite_ratio(cname)
        post_cell = cells.Bushy.create(
            morphology=str(fn),
            # decorator=Decorator,
            species="mouse",
            modelType="II",
            modelName="XM13",
        )
        self.HR = post_cell.hr
        AdjArea = AdjustAreas(method="pt3d")
        fns = str(fn)
        if not fns.find('MeshInflate'): # if mesh inflated, this has already been applied
            AdjArea.sethoc_fromCNcell(post_cell)
            # AdjArea.sethoc_fromstring(hdata=hocstruct2)
            AdjArea.cell.print_soma_info()
            AdjArea.adjust_diameters(
                sectypes=AdjArea.somas, inflationRatio=sinflateratio
            )
            AdjArea.adjust_diameters(
                sectypes=AdjArea.dendrites, inflationRatio=dinflateratio
            )
        AM = AxonMorph(cellID=cname)
        self.parent_soma = None
        sectypes = []
        ah_coords = List3DPoints()
        ais_coords = List3DPoints()
        myel_coords = List3DPoints()
        soma_coords = List3DPoints()
        dend_coords = List3DPoints()
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            sectypes.append(sectype)
            if sectype == "soma":
                soma_coords = self.fill_coords(soma_coords, sec)
            if sectype in self.dendrites:
                dend_coords = self.fill_coords(dend_coords, sec, insertnan=True)
            if sectype not in self.axons:
                continue
            if sectype == "Axon_Hillock":
                sref = h.SectionRef(sec=sec)
                self.parent_soma = sref.trueparent
                AM.HillockLength.append(sec.L)
                AM.HillockDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.HillockDiam0.append(segs[0].diam)  # end diameters for section
                AM.HillockDiam1.append(segs[-1].diam)
                ah_coords = self.fill_coords(ah_coords, sec)

            if sectype in ["Axon_Initial_Segment", "Unmyelinated_Axon"]:
                AM.AISLength.append(sec.L)
                AM.AISDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.AISDiam0.append(segs[0].diam)  # end diameters for section
                AM.AISDiam1.append(segs[-1].diam)
                ais_coords = self.fill_coords(ais_coords, sec)

            if sectype == "Myelinated_Axon":
                AM.MyelLength.append(sec.L)
                AM.MyelDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.MyelDiam0.append(segs[0].diam)  # end diameters for section
                AM.MyelDiam1.append(segs[-1].diam)
                myel_coords = self.fill_coords(myel_coords, sec)

            if sectype == "Axon":
                AM.AxonLength.append(sec.L)
                AM.AxonDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.AxonDiam0.append(segs[0].diam)  # end diameters for section
                AM.AxonDiam1.append(segs[-1].diam)
                # print('secs: ', secs)
        self.axon_origin = Point(ah_coords.x[0], ah_coords.y[0], ah_coords.z[0])
        self.all_coords = [dend_coords, soma_coords, ah_coords, ais_coords, myel_coords]
        self.coords = [dend_coords, soma_coords]
        # self.remove_old_axon()
        self.old_axon_coords = [
            ah_coords,
            ais_coords,
            myel_coords,
        ]  # save just the axon ones

        return AM

    def remove_old_axon(self) -> None:
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            if sectype in self.axons:
                h.delete_section(sec=sec)

    def remove_all_except_axon(self) -> None:
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            if sectype not in self.axons:
                h.delete_section(sec=sec)
    
    def plot_3d(
        self,
        ax: Type[mpl.axes],
        data: Type[Point],
        otherdata: Union[None, Type[Point]] = None,
    ) -> None:
        """
        plot the coordinates in data (which is expected to be a list of coordinates)
        in 3D
        To get some rough semblence of 3d sizes, there are 2 plots. One is just lines
        connecting everything, and the other is symbols (spheres) at each point with a diam
        proportional to the diameter.
        This not meant to be the detailed representation, but just enough to indicate the
        orientations of things.

        """

        markers = ["o", "o", "o", "o", "o"]
        mc = ["k", "y", "b", "r", "c"]
        for i, d in enumerate(data):
            print("plotting main data i: ", i)
            ax.scatter(
                d.x,
                d.y,
                d.z,
                color=mc[i],
                marker=markers[i],
                s=np.power(d.d, 2.0),
                alpha=0.5,
            )
            ax.plot(d.x, d.y, d.z, color=mc[i], linestyle="-", linewidth=2.0)
        if (
            otherdata is not None
        ):  # maybe a second plot to show the rotation/appended version.
            for i, d in enumerate(otherdata):
                print("plotting other data i: ", i)
                ax.scatter(
                    d.x,
                    d.y,
                    d.z,
                    color=mc[i + 2],
                    marker=markers[i],
                    s=np.power(d.d, 2.0),
                    alpha=0.5,
                )
                ax.plot(d.x, d.y, d.z, color=mc[i + 2], linestyle="-", linewidth=2.0)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.set_zlim(-100, 100)

    def concatenate(self, pt3d_lists: list, translate: bool = False) -> object:
        """
        Concatenate a list of sections, end to end along the x axis,
        starting with the first section at 0.
        Does not assign section numbers, just changes x positions
        
        Parameters
        ----------
        pt3d_lists :
            list of points to "concatenate" in the physical sense
            The lists are not adjoined, but the positions are adjusted
            so that they flow in order.
        translate: bool
            flag to force the result to be translated to an origin value
            specified in self.axon_origin (usually the Hillock 0 end)
    
        Returns
        -------
        pt3d_lists: updated point lists.
        """
        x0 = 0.0
        y0 = 0.0
        z0 = 0.0
        if translate:
            x0 = self.axon_origin.x
            y0 = self.axon_origin.y
            z0 = self.axon_origin.z
        print(x0, y0, z0, self.axon_origin)
        max_x = 0.0
        for index, pt in enumerate(pt3d_lists):
            pt.x += x0 + max_x
            pt.y += y0
            pt.z += z0
            max_x = np.max(pt.x) - np.min(pt.x)
        return pt3d_lists

    def pt3d_to_hoc(
        self, pt3d_lists: list, name: list = [], name_type: list = [], parent: list = []
    ):
        hoc_string = {}
        for i, pt in enumerate(pt3d_lists):
            nt = name_type[i]
            hoc_string[nt] = {}
            hoc_string[nt]["header"] = f"access sections[{name[i]:s}_start]" + "\n"
            hoc_string[nt]["header"] += f"{name_type[i]:s}.append()" + "\n"
            hoc_string[nt]["header"] += (
                f"connect sections[{name[i]:s}_start](0), sections[{parent[i]:s}_end](1)"
                + "\n"
            )
            hoc_string[nt]["header"] += f"sections[{name[i]:s}_start] " + "{\n"
            hoc_string[nt]["pt3d"] = ""
            for n, p in enumerate(pt.x):
                hoc_string[nt]["pt3d"] += (
                    f"\tpt3dadd({pt.x[n]:.6f}, {pt.y[n]:.6f}, {pt.z[n]:.6f}, {pt.d[n]:.6f})  // ARTIFICIAL AXON"
                    + "\n"
                )
            hoc_string[nt]["cloing"] = "}\n"
        return hoc_string

    def populate_pt3d_list(self, tab) -> Type[List3DPoints]:
        """
        Make a list of 3D points corresponding to the table for
        a section.

        This list makes sections that vary in x, but are aligned along
        y and z = 0
        
        Parameters
        ----------
        tab : the table specifying the values as:
            0 : min diameter, 1: maximum diameter, 2 : maximum x,

            3 : number of segments to populate
        """
        coords = List3DPoints()
        x_steps = np.linspace(0.0, tab[2], tab[3], endpoint=True)
        dias = np.linspace(tab[0], tab[1], tab[3], endpoint=True)
        coords.x = x_steps
        coords.y = np.zeros_like(coords.x)
        coords.z = np.zeros_like(coords.x)
        coords.d = dias
        return coords

    def populate_section(self, name, name_type, parent, tab):
        """
        Create the hoc code for a section, using the pt3dadd formatting
        
        Parameters
        ----------
        name : str
            the name of the section
        name_type: str
            the name of the section list (in Neuron)
        parent : str
            the name of the parent section type (typically in this usage, the soma)
        tab : list
            Table of values used to build a tapered section of pt3ds
        
        Returns
        -------
        hoc_string: the hoc code to build a single section 
        """
        hoc_string = f"access sections[{name:s}_start]" + "\n"
        hoc_string += f"{name_type:s}.append()" + "\n"
        hoc_string += (
            f"connect sections[{name:s}_start](0), sections[{parent:s}_end](1)" + "\n"
        )
        hoc_string += f"sections[{name:s}_start] " + "{\n"
        x_steps = np.linspace(0.0, tab[2], tab[3], endpoint=True)
        dias = np.linspace(tab[0], tab[1], tab[3], endpoint=True)
        for n in range(tab[3]):
            hoc_string += (
                f"\tpt3dadd({x_steps[n]:.6f}, 0., 0., {dias[n]:.6f})  // ARTIFICIAL AXON"
                + "\n"
            )
        hoc_string += "}\n"

        return hoc_string

    def collapse_list(self, axonslist: Type[List3DPoints]) -> Type[np.ndarray]:
        """
        Take the list of axons and generate a single array of x, y, z values
        Parameters
        ----------
        axonslist : a list of axon sections to collapse
        """
        newlist = List3DPoints()
        for a in axonslist:
            newlist.x.extend(a.x)
            newlist.y.extend(a.y)
            newlist.z.extend(a.z)
            newlist.d.extend(a.d)
        nl = np.vstack([newlist.x, newlist.y, newlist.z])
        return nl

    def compute_axon_vector(self, axonslist: Type[List3DPoints]):
        """
        Compute the mean 3d vector for the axon, going through all the points in the
        axon list.

        Vector is anchored at the 0 end of the axon hillock (axonslist[0])
        """
        zp = np.array([axonslist[0].x[0], axonslist[0].y[0], axonslist[0].z[0]])
        print("compute_axon_vector: anchor point: ", zp)
        print(axonslist)
        u = -1
        print("compute_axon_vector: ", axonslist[u].x)
        if len(axonslist[u].x) == 0:
            u = -2
        print("compute_axon_vector:", u, axonslist[u].x)
        dp = np.array(
            [axonslist[u].x[-1], axonslist[u].y[-1], axonslist[-u].z[-1]]
        )  # distal point
        vp = np.vstack([zp, dp])
        dp = vp - zp
        print("compute_axon_vector:", "distal point (relative to zp): ", dp)
        nl = self.collapse_list(axonslist)  # make original axon into 3d array
        nlr = np.mean(nl, axis=1)
        nlr = np.vstack([[0, 0, 0], nlr])
        print("compute_axon_vector:", "nlr: ", nlr)
        alv = R.align_vectors(dp, nlr)
        print("compute_axon_vector:", "alv: ", alv)

    def generate_standard_axon(self, translate: bool = True, AIS_length:float=0.) -> str:
        """
        Default standarized axon for VCN bushy cells
        values are diam[0], diam[1], length, and nsegments
        See VCN-SBEM-Data/Hillock-AIS-threshold.pzfx for the values
        entered here.
        
        If the AIS length is > 0, then we replace the AIS with one
        of the specified length in microns
        """
        # hil = [2.31, 1.87, 1.97, 3]
        # ais = [1.6, 0.82, 19.03, 19]
        # axn = [1.03, 1.49, 55.7, 55]

        p1 = self.populate_pt3d_list(StandardAxon.hil)
        print('gen std axon ais len: ', AIS_length)
        if AIS_length == 0.:
            p2 = self.populate_pt3d_list(StandardAxon.ais)
        else:
            ModifiableAxon.ais[2] = AIS_length
        p2 = self.populate_pt3d_list(ModifiableAxon.ais)
            
        p3 = self.populate_pt3d_list(StandardAxon.axn)
        print(ModifiableAxon.ais)
        axon = self.concatenate([p1, p2, p3], translate=True)
        self.new_axon = axon
        print("axon: ", axon)
        hoc_str = self.pt3d_to_hoc(
            axon,
            name=["hillock", "aix", "axon"],
            name_type=["Axon_Hillock", "Axon_Initial_Segment", "Myelinated_Axon"],
            parent=["soma", "hillock", "ais"],
        )
        return hoc_str

    def get_create_pos(self, hocf: str) -> int:
        m = re_find_create_sections_position.search(hocf)
        startpos = m.start()
        return startpos

    def verify_section_lists(self, hocf: str) -> str:
        """
        Verify that each of the axon section types is actually
        created at the top of the file
        """
        startpos = self.get_create_pos(hocf)
        m = re_find_create_AH.search(hocf)
        if m is None:
            add_str = "objref Axon_Hillock\nAxon_Hillock = new SectionList()\n"
            hocf = hocf[:startpos] + add_str + hocf[startpos:]
            print("Inserted Axon Hillock")

        m = re_find_create_AIS.search(hocf)
        if m is None:
            add_str = "objref Axon_Initial_Segment\nAxon_Initial_Segment = new SectionList()\n"
            hocf = hocf[:startpos] + add_str + hocf[startpos:]
            print("Inserted Axon Initial Segment")

        m = re_find_create_MA.search(hocf)
        if m is None:
            add_str = "objref Myelinated_Axon\nMyelinated_Axon = new SectionList()\n"
            hocf = hocf[:startpos] + add_str + hocf[startpos:]
            print("Inserted Myelinated_Axon")
        startpos = self.get_create_pos(hocf)
        m = re_find_create_MA.search(hocf)
        return hocf

    def update_create_sections(self, func: object, hocf: str, nsec: int) -> Tuple[str, int]:
        m = func.findall(hocf)
        if len(m) == 0:
            hocf, nrepl = re_find_create_sections.subn(
                f"create sections[{nsec+1:d}]", hocf
            )
            nsec += 1
            print("update section count done")
        return hocf, nsec

    def replace_points(
        self, func: object, hocf: str, hname: str, newpt3d: Type[Point3D]
    ) -> str:
        m = func.search(hocf)
        if m is None:
            # current_used_secs = []
            asec = re_find_accessed_sections.findall(hocf)
            for a in asec:
                if a == 0 or len(a) < 2:
                    break
                asec.append(int(a[1]))
            print(" no match found for : ", func.pattern)
            if hname == "Myelinated_Axon":
                d = re_find_AIS.search(hocf)
                print("Myelinated axon from AIS: ", d)
                print(d.groups("section"))

            elif hname == "Axon_Hillock":
                d = re_find_soma.search(hocf)
                print("Hillock to soma: ", d)
                print(d.groups("section"))
            elif hname == "Axon_Initial_Segment":
                d = re_find_AH.search(hocf)
                print("Ais from Hillock: ", d)
                print(d.groups("section"))
            else:
                raise ValueError("hname must match one of a limited set of names")
            previoussection = int(d.groups("section")[1])
            if self.new_max_sec not in asec:
                # generate the text:
                last_sec = int(self.new_max_sec[-1][1]) - 1
                print(last_sec)
                newsec = f"access sections[{last_sec:d}]" + "\n"
                newsec += f"{hname:s}.append()" + "\n"
                newsec += (
                    f"connect sections[{last_sec:d}](0), sections[{previoussection:d}](1)"
                    + "\n"
                )
                newsec += f"sections[{last_sec:d}] " + "{\n"
                newsec += newpt3d
                newsec += "}" + "\n"
            hocf += newsec
        else:
            # print('found m: ', m)
            # print(m.groups())
            # print(m.group('pt3ddata'))
            hocf = hocf.replace(m.group("pt3ddata"), newpt3d)  # print(dir(m))
            # print(dir(m[0]), '\n', m[0])
        return hocf

    def change_hoc(self, cell: str, hocstr: str) -> None:
        cname, fn = self.make_name(cell)
        with open(fn) as fh:
            hocf = fh.read()

        cs = re_find_create_sections.findall(hocf)
        max_sec = 0
        for s in cs:
            max_sec = max(max_sec, int(s[1]))
        self.original_max_sec = max_sec
        print("Original Max sec: ", max_sec)
        hocf = self.verify_section_lists(hocf)

        # find all of the AH, AIS and MA section blocks, modify the number of
        # created sections if needed
        funcs = [re_find_AH, re_find_AIS, re_find_MA]
        for f in funcs:
            hocf, max_sec = self.update_create_sections(f, hocf, max_sec)
        cs = re_find_create_sections.findall(hocf)
        print("Updated section count is: ", cs)
        self.new_max_sec = cs

        funcs = [re_find_AH_pt3d, re_find_AIS_pt3d, re_find_MA_pt3d]
        hname = ["Axon_Hillock", "Axon_Initial_Segment", "Myelinated_Axon"]
        for i, f in enumerate(funcs):
            hocf = self.replace_points(f, hocf, hname[i], hocstr[hname[i]]["pt3d"])
            # now, for the same set of funcs, if the
        # print('new hocf: ', hocf[:600])
        self.hocf = hocf

    def write_revised_hoc(self, cell: str, AIS_length=0.) -> None:
        cname, fn = self.make_name(cell, add_standardized_axon=True, AIS_length=AIS_length)
        new_comment = "//  Modified file: Uses a standardized bushy cell axon\n"
        new_comment += (
            "//  including Axon_Hillock, Axon_Initial_Segment and Myelinated_Axon\n"
        )
        if AIS_length > 0.:
            new_comment += (
            f"//    AIS length adjusted: {AIS_length:.2f} microns\n"
        )
        new_comment += (
            f"//     {datetime.datetime.now().isoformat(sep=' ', timespec='seconds'):s}"
        )
        new_comment += "\n"
        self.hocf = new_comment + self.hocf
        with open(fn, "w") as fh:
            fh.write(self.hocf)

def make_standard_axon(cell: str, write: bool = False) -> None:
    """
    Read and replace original axon with a "standard axon" (average)
    """
    # fig = mpl.figure(figsize=(12, 6))
    # ax1 = fig.add_subplot(121, projection="3d")
    PA = MakeStandardAxon(revised=False)
    PA.do_axon(cell)
    # PA.plot_3d(ax1, PA.all_coords)
    hocstr = PA.generate_standard_axon()
    PA.change_hoc(cell, hocstr)
    if write:
        PA.write_revised_hoc(cell)
    else:
        print(PA.hocf)

    # now display the new version
    # ax2 = fig.add_subplot(122, projection="3d")
    PB = MakeStandardAxon(revised=True)
    PB.do_axon(cell)
    # PB.plot_3d(ax2, PB.all_coords)
    # mpl.show()

def make_modified_axons(cell: str, write: bool = False, AIS_length:float = 0.) -> None:
    """
    Read and replace original axon with a "standard axon" (average)
    but with a specified length for the AIS
    """
    # fig = mpl.figure(figsize=(12, 6))
    # ax1 = fig.add_subplot(121, projection="3d")
    PA = MakeStandardAxon(revised=False)
    PA.do_axon(cell)
    # PA.plot_3d(ax1, PA.all_coords) # original cell


    # now create the new version
    hocstr = PA.generate_standard_axon(AIS_length=AIS_length)
    PA.change_hoc(cell, hocstr)
    if write:
        PA.write_revised_hoc(cell, AIS_length=AIS_length)
    else:
        pass
        # print(PA.hocf)

    # now display the new version
    # ax2 = fig.add_subplot(122, projection="3d")
    PB = MakeStandardAxon(revised=True)
    PB.do_axon(cell)
    # PB.plot_3d(ax2, PB.all_coords)
    # mpl.show()

def read_standard_axon(cell: str) -> None:
    fig = mpl.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection="3d")
    PA = MakeStandardAxon(revised=False)
    PA.do_axon(cell)
    PA.plot_3d(ax1, PA.all_coords)
    ax2 = fig.add_subplot(122, projection="3d")
    PB = MakeStandardAxon(revised=True)
    PB.do_axon(cell)
    PB.plot_3d(ax2, PB.all_coords)
    mpl.show()

if __name__ == "__main__":

    cell = 17
    for cell in [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]:
        make_standard_axon(cell, write=True)
        for aislen in [10., 12., 14., 16., 18., 20., 22., 24., 26.]:
            make_modified_axons(cell, write=True, AIS_length=aislen)
    #
    # make_standard_axon(cell, write=True)
    # for cell in  gradeACells:
        # PA = MakeStandardAxon(revised=False)
        # PA.convert_swc_hoc(cell)
        # make_standard_axon(cell, write=True)
        # read_standard_axon(cell)
