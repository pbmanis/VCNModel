"""
Adjust the areas of named processes to match
a specific value.
This class was designed to "inflate" the areas of swc/hoc
representations of reconstructed neurons to match the area
measured with hihger resolution serial-blockface EM reconstructions.
This steip is required because the swc/hoc representation is a smoothed
cylinderical version of a rather convoluted (or crenulated) cell surface,
and therefore underestimates the actual membrane area.

The procedure adjusts the areas of each section for a set of 
named processes (e.g., ['apical_dendrite', 'basal_dendrite']) by
applying a fixed scale factor to the diameter of each pt3d object
in that section. Because NEURON then readjustes the area
calculation the final inflation factor (a single, positive real number)
the inflation scale factor may differ slightly for each section. The
differendes are most prominent for very small diameter sections, but 
also appear for larger sections. Therefore, the procedure does the
adjustment iteratively until the area is inflated by the desired factor,
and changes less than a specified epsilon between two successive adjustments.

In the days of coronavirus, March 2020, Sunset Beach, NC.


"""

import numpy as np
from pprint import pprint
from pathlib import Path
from math import sin, cos
import matplotlib.pyplot as mpl

from neuron import h

import vcnmodel.cell_config as cell_config
from vcnmodel import h_reader
import pylibrary.tools.cprint as CP

cprint = CP.cprint

"""
Two Test cases: a simple soma,
and a soma with a couple of dendrites.
"""

hocstruct = """

objref dendrite
dendrite = new SectionList()

create sections[1]

access sections[0]
dendrite.append()
sections[0] {
  pt3dadd(-66.472813, 63.896850, 23.493839, 1.588886)
  pt3dadd(-67.496068, 62.617781, 23.633374, 1.191287)
  pt3dadd(-68.007696, 61.338711, 23.912444, 1.35)
  pt3dadd(-68.263510, 60.827083, 24.051979, 1.45)
  pt3dadd(-68.007696, 60.059641, 24.470583, 1.65)
  pt3dadd(-67.496068, 58.013129, 24.889188, 1.790698)
  pt3dadd(-67.240254, 56.222432, 24.889188, 2)
  pt3dadd(-67.240254, 56.222432, 26.889188, 3)
  pt3dadd(-67.240254, 56.222432, 28.889188, 2)
  pt3dadd(-67.240254, 56.222432, 30.889188, 1)
  pt3dadd(-67.240254, 56.222432, 32.889188, 0.5)
}

"""

hocstruct2 = """

objref soma
soma = new SectionList()
objref dendrite
dendrite = new SectionList()

create sections[5]

access sections[0]
soma.append()
sections[0] {
    pt3dadd(0., 0., -10., 0.5)
    pt3dadd(0., 0., -7.5, 2.5)
    pt3dadd(0., 0., -5., 8.5)
    pt3dadd(0., 0., -0., 15.0)
    pt3dadd(0., 0., 5., 9.5)
    pt3dadd(0., 0., 7.5, 7.5)
    pt3dadd(0., 0., 10., 5.0)
    pt3dadd(0., 0., 12.5, 1.0)
}

access sections[1]
dendrite.append()
connect sections[1](0), sections[0](1)
sections[1] {
  pt3dadd(0., 0., 12.5, 1.0)
  pt3dadd(3., 3., 18., 0.75)
  pt3dadd(5., 5., 24., 0.5)
  pt3dadd(7., 7., 32., 0.25)
  pt3dadd(9., 9., 45., 0.25)

}

access sections[2]
dendrite.append()
connect sections[2](0), sections[1](1)
sections[2] {
  pt3dadd(0., 9., 45, 0.25)
  pt3dadd(3., 12., 48., 0.75)
  pt3dadd(5., 18., 52., 0.75)
  pt3dadd(7., 24., 55., 0.25)
  pt3dadd(9., 30., 60., 0.25)

}

access sections[3]
dendrite.append()
connect sections[3](0), sections[0](1)
sections[3] {
  pt3dadd(0., 0., 12.5, 1.0)
  pt3dadd(3., -3., 18., 0.75)
  pt3dadd(5., -5., 24., 0.5)
  pt3dadd(7., -7., 32., 0.25)
  pt3dadd(9., -9., 45., 0.25)

}

access sections[4]
dendrite.append()
connect sections[4](0), sections[3](1)
sections[4] {
  pt3dadd(0., -9., 45, 0.25)
  pt3dadd(3., -12., 48., 0.75)
  pt3dadd(5., -18., 52., 0.75)
  pt3dadd(7., -24., 55., 0.25)
  pt3dadd(9., -30., 60., 0.25)

}
"""


class AdjustAreas:
    """
    This class encapsulates routines to adjust the areas of 
    hoc reconstructed files using pt3d data.
    """
    def __init__(self, method='pt3d'):
        self.method = method  # wheter to compute areas from pt3d data,  or segment data
        if self.method == 'pt3d':
            self.areafunc = self.segareasec
            self.adjust_dia_func = self.adjust_pt3d_dia
        elif self.method == 'segment':
            self.areafunc = self.segareasec
            self.adjust_dia_func = self.adjust_seg_dia

        else:
            raise ValueError('AdjustAReas: userpted must be one of pt3d, segment')

    def sethoc_fromCNcell(self, cell):
        self.cell = cell
        self.HR = cell.hr
        for s in self.cell.hr.h.allsec():
            s.Ra = 150.0
        S = SetNSegs(self.HR.h)
        S.set_nseg(freq=1000.0, d_lambda=0.025, minimum=3)  # perform d_lambda adjustment on structure
        
    def sethoc_fromfile(self, fn:[str, Path]):
        """
        Read a reconstruction and set the number of segments using rule
        
        Parameters
        ----------
        fn : filename (string or Path)
        """
        
        self.HR = h_reader.HocReader(str(fn))
        for s in self.HR.h.allsec():
            s.Ra = 150.0
        S = SetNSegs(self.HR.h)
        S.set_nseg(freq=1000.0, d_lambda=0.025, minimum=3)  # perform d_lambda adjustment on structure

    def sethoc_fromstring(self, hdata:str):
        """
        Read a reconstruction from a string
        and set the number of segments using rule
        
        Parameters
        ----------
        hdata : string
        """
        
        h(hdata) # neuron reads it
        self.HR = h_reader.HocReader(h)  # placeholder for the hoc_render structure - we only need to assign h
        for s in self.HR.h.allsec():
            s.Ra = 150.0
        S = SetNSegs(self.HR.h)
        S.set_nseg(freq=1000.0, d_lambda=0.025, minimum=3)


    def section_diam_test(self):
        ############### now with just diam in section
        self.HR.h(hocstruct)
        self.HR.h.topology()

        for s in self.HR.h.allsec():
            s.Ra = 150.0

        dend = self.HR.h.sections[0]

        a1 = self.HR.h.area(0.5, sec=dend)

        dend.diam = dend.diam * sf
        a2 = self.HR.h.area(0.5, sec=dend)

        print(f"sf: {sf:.2f}  a1: {a1:.3f}  a2: {a2:.3f}  a2/a1: {a2/a1:.2f}")

    def segareasec(self, sec: object, name: str = "original"):
        """
        Sum up the areas of all the _segments_ in a section

        """

        area = 0
        for i, seg in enumerate(sec.allseg()):
            # print(f"   segment: {i:d} area={seg.area():.3f}")
            area += seg.area()
        # print(f'{name:s} area: {area:.3f} ')
        return area
    
    def pt3dareasec(self, sec: object, name: str = "original"):
       """
       Sum up the areas of all the pt3d pieces in a section

       """

       area = 0
       for i in range(sec.n3d()):
           area += np.pi*sec.arc3d(i)*sec.diam3d(i)
       return area 
    


    def adjust_pt3d_dia(self, sec: object, sf: float):
        """
        Change the diamete4r of the detailed pt3d type of representations
        before they are converted to the segmented representation
        """
        h.pt3dconst(1, sec=sec)
        for i in range(sec.n3d()):
            diam = sec.diam3d(i)
            sec.pt3dchange(i, diam * sf)
    
    def adjust_seg_dia(self, sec: object, sf: float):
        h.pt3dconst(0, sec=sec)
        for seg in sec.allseg():
            seg.diam *= sf

    def adjust_segment_dia_areamatch(
        self, sec: object, sf: float, eps: float = 1e-4
    ):
        """
        This function expands the surface area of a section by a
        specific fractional amount. Because of the way that the
        areas are calculated, we focus on using the pt3d area
        sums before they are converted to neuron sections, then
        compute the areas from the  sumo of the individual segmentes
        in the section. The routine is iterative as the calculation
        as done in neuron is not always predicable.
        """

        if sf == 1.0:
            return
        section_area = self.areafunc(sec)  # starting value
        target_area = section_area * sf
        n = 0
        while np.fabs(target_area - section_area) > eps:
            ratio = target_area / section_area
            if ratio < 0.01:
                print("Inflation failed: ratio: ", ratio, section_area, target_area, n)
                exit()
            self.adjust_seg_dia(sec, ratio)
            section_area = self.areafunc(sec)
            n += 1
            # print(f"n: {n:d} section area: {section_area:.3f}  target  area: {target_area:.3f}")
        # print(f"sec: {str(sec):s}  iters: {n:d}")

    def adjust_pt3d_dia_areamatch(
        self, sec: object, sf: float, eps: float = 1e-4,

    ):
        """
        This function expands the surface area of a section by a
        specific fractional amount. Because of the way that the
        areas are calculated, we focus on using the pt3d area
        sums before they are converted to neuron sections, then
        compute the areas from the  sumo of the individual segmentes
        in the section. The routine is iterative as the calculation
        as done in neuron is not always predicable.
        """

        if sf == 1.0:
            return
        section_area = self.areafunc(sec)  # starting value
        target_area = section_area * sf
        n = 0
        while np.fabs(target_area - section_area) > eps:
            ratio = target_area / section_area
            if ratio < 0.01:
                print("Inflation failed: ratio: ", ratio, section_area, target_area, n)
                exit()
            self.adjust_dia_func(sec, ratio)
                
            # h.define_shape()
            section_area = self.areafunc(sec)

            # print(f"section area: {section_area:.3f}  target  area: {target_area:.3f}")
            n += 1
        print(f"sec: {str(sec):s}  iters: {n:d}")

    def adjust_diameters(self, sectypes:list, inflateRatio:float):

        adj_data = {  # data structure
            "Areas3D": [],
            "Areas3DInflated": [],
            "InflateRatio": inflateRatio,
            "Diams": [],
            "Lengths": [],
            "n3d": [],
            "sectypes": [],
        }
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            # print(f"sec: {str(sec):>14s} {sectype:>20s}  n={sec.nseg:>3d} L={sec.L:8.1f}  D={sec.diam:7.3f}  (L/nsec: {sec.L/sec.nseg:>6.2f})")
            if sectype in sectypes:
                adj_data["Areas3D"].append(self.areafunc(sec))
                h.pt3dconst(1, sec=sec)
                self.adjust_adj_data_dia_areamatch(sec, inflateRatio)
                adj_data["Areas3DInflated"].append(self.areafunc(sec))
                adj_data["Diams"].append(sec.diam)
                adj_data["Lengths"].append(sec.L)
                adj_data["n3d"].append(sec.n3d())
                adj_data["sectypes"].append(sectype)

        if self.method:
            print("\nChanging the areas by changing adj_data diams")
        else:
            print("\nChanging the areas by changing segment diams")
            
        print(f"inflated area: {np.sum(adj_data['Areas3DInflated']):10.2f} ")
        print(
            f"Infl factor: {np.sum(adj_data['Areas3DInflated'])/np.sum(adj_data['Areas3D']):6.3f} (should be: {adj_data['InflateRatio']:6.3f})"
        )

        return adj_data

    def plot_areas(self, pt3d: dict):

        f, ax = mpl.subplots(2, 2)
        ax = ax.ravel()
        title = str(pt3d['sectypes']).replace('_', '\_')
        # f.suptitle(f"Adjusted areas: {title:s}")
        ndend = len(pt3d["Areas3D"])
        xdend = np.arange(ndend)
        ax[0].plot(xdend, pt3d["Areas3D"], "k-", label="Original")
        ax[0].plot(
            xdend, pt3d["Areas3DInflated"], "r-", linewidth=0.5, label="Inflated"
        )
        ax[0].set_ylim((0, 100))
        ax[0].set_xlabel("Section number")
        ax[0].set_ylabel("Surface Area (um2)")
        ax[0].legend()

        ax[1].plot(
            xdend,
            np.array(pt3d["Areas3DInflated"]) / np.array(pt3d["Areas3D"]),
            linewidth=0.5,
            label="Ratio",
        )
        pt = np.argwhere(np.array(pt3d["n3d"]) <= 10)
        scp0 = ax[1].scatter(
            np.take(xdend, pt),
            np.take(pt3d["Areas3DInflated"], pt) / np.take(pt3d["Areas3D"], pt),
            c=np.take(pt3d["n3d"], pt),
            s=3,
            cmap="plasma",
            alpha=0.75,
        )
        f.colorbar(scp0, ax=ax[1])
        ax[1].plot(
            xdend,
            pt3d["InflateRatio"] * np.ones(ndend),
            linewidth=0.5,
            label="Target",
        )
        ax[1].set_ylim(1.0, 2.0)
        ax[1].set_xlabel("Section number")
        ax[1].set_ylabel("Ratio inflated/original")
        ax[1].legend()

        # ax[2].plot(xdend, np.array(secdiams), linewidth=0.5, label="Original Diams")
        # ax[2].plot(xdend, np.array(secdiamsInflated), linewidth=0.5, label="Inflated Diams")
        # ax[2].plot(xdend, np.array(secdiamsInflated)/np.array(secdiams), 'k-', linewidth=0.75, label="Ratios")
        # ax[2].set_ylim(0, 8.0)
        # ax[2].set_xlabel("dendrite section number")
        # ax[2].set_ylabel("Diameter (um)")
        # ax[2].legend()
        #
        d3 = np.array(pt3d["Areas3D"])
        ds = np.array(pt3d["Areas3DInflated"])
        alld = np.concatenate((d3, ds))
        dmin = np.min(alld)
        dmax = np.max(alld)
        pt = np.argwhere(np.array(pt3d["n3d"]) <= 10)
        # np.array([9, 16, 17, 86, 91, 92, 34, 35, 38, 43, 46, 48, 53, 54, 79])
        # scp = ax[2].scatter(np.take(d3, pt), np.take(ds, pt), s=4, c=np.take(pt3d["n3d"], pt), cmap="hsv", alpha=0.75)
        scp = ax[2].scatter(pt3d["n3d"], ds / d3, c=d3, s=3)
        f.colorbar(scp, ax=ax[2])
        unity = np.arange(dmin, dmax)
        scaled = unity * pt3d["InflateRatio"]
        unity = unity * 1.65
        # ax[2].plot(np.arange(len(unity)), unity, "c--", linewidth=0.25)
        # ax[2].plot(np.arange(len(scaled)), scaled, "k--", linewidth=0.25)
        ax[2].set_ylabel("3d Area inflated")
        ax[2].set_xlabel("3D Area before inflation")
        ax[2].set_title("Areas selected < 3 pt3d")

        scp = ax[3].scatter(d3, ds, s=4, c=pt3d["n3d"], cmap="hsv", alpha=0.75)
        f.colorbar(scp, ax=ax[3])
        unity = np.arange(dmin, dmax)
        scaled = unity * pt3d["InflateRatio"]
        unity = unity * 1.65
        ax[3].plot(np.arange(len(unity)), unity, "c--", linewidth=0.25)
        ax[3].plot(np.arange(len(scaled)), scaled, "k--", linewidth=0.25)
        ax[3].set_ylabel("3d Area inflated")
        ax[3].set_xlabel("3D Area before inflation")
        ax[2].set_title("Areas re inflation")
        mpl.savefig("../VCN-SBEM-Data/AreaCalc_Section.pdf")
        mpl.show()


class SetNSegs(object):
    def __init__(self, cell: object):
        """
        Set the number of segments needed  in each section
        for a particular cell with
        
        cell : Neuron hoc object for the cell
        """

        self.cell = cell

    def set_nseg(self, freq: float = 100.0, d_lambda: float = 0.1, minimum: int = 5):
        """
        Sets nseg in each section to an odd value so that its segments are no longer than
        d_lambda x the AC length constant at frequency freq in that section.
        The defaults are reasonable values for most models
        Be sure to specify your own Ra and cm before calling geom_nseg()

        To understand why this works,
        and the advantages of using an odd value for nseg,
        see  Hines, M.L. and Carnevale, N.T. NEURON: a tool for neuroscientists. The Neuroscientist 7:123-135, 2001.
        This is a python version of the hoc code.
        
        Parameters
        ----------
        freq : float, default=100. (Hz)
            Frequency in Hz to use in computing nseg.
        d_lambda : float, default=0.1
            fraction of AC length constant for minimum segment length
        
        """
        if self.cell is None:  # no hoc reader file, so no adjustments
            return
        for section in self.cell.allsec():
            # print(section, section.L, d_lambda, freq)
            nseg = 1 + 2 * int(
                (section.L / (d_lambda * self._lambda_f(section, frequency=freq)) + 0.9)
                / 2.0
            )
            if nseg < minimum:
                nseg = minimum  # ensure at least 3 segments per section...
            try:
                section.nseg = int(nseg)
            except:
                print("nseg: ", nseg)
                raise (ValueError)

    def _lambda_f(self, section: object, frequency: float = 100.0):
        """
        get lambda_f for the section (internal)
        
        Parameters
        ----------
        freq : float, default=100. (Hz)
            Frequency in Hz to use in computing nseg.
        section : Neuron section object
        
        Returns
        -------
        section length normalized by the length constant at freq.
        """
        self.cell("access %s" % section.name())
        if self.cell.n3d() < 2:
            return 1e-5 * np.sqrt(
                section.diam / (4.0 * np.pi * frequency * section.Ra * section.cm)
            )
        # above was too inaccurate with large variation in 3d diameter
        # so now we use all 3-d points to get a better approximate lambda
        x1 = self.cell.arc3d(0)
        d1 = self.cell.diam3d(0)
        lam = 0.001
        # print(self.cell.n3d(), x1, d1)
        for i in range(int(self.cell.n3d()) - 1):
            x2 = self.cell.arc3d(i)
            d2 = self.cell.diam3d(i)
            # print("  seg: ", i, x2, d2)
            lam = lam + ((x2 - x1) / np.sqrt(d1 + d2))
            x1 = x2
            d1 = d2
        # print("lam: ", lam)
        #  length of the section in units of lambda
        lam = (
            lam
            * np.sqrt(2.0)
            * 1e-5
            * np.sqrt(4.0 * np.pi * frequency * section.Ra * section.cm)
        )
        # print("_lambda_f: ", section.L, section.Ra, section.cm, section.nseg, lam)
        return section.L / lam




def simple_hoc():
    h(hocstruct2)
    h.topology()
    S = setnsegs(h)
    S.set_nseg(freq=1000.0, d_lambda=0.01)
    for sec in h.allsec():
        print(sec.nseg)

    # s = h.Shape()
    # s.show(0)
    # while(1):
    #     pass


def by_section_diam(fn, cname: str):
    # re-read original data
    HR2 = h_reader.HocReader(str(fn))  # h.load_file(str(fn))
    # print(dir(h))
    # exit()
    cconfig = cell_config.CellConfig()
    sinflateratio = cconfig.get_soma_ratio(cname)
    dinflateratio = cconfig.get_dendrite_ratio(cname)

    S = SetNSegs(HR2.h)
    S.set_nseg(freq=1000.0, d_lambda=0.01)

    asoma = 0.0
    aother = 0.0
    asoma_inf = 0.0
    adend_inf = 0.0
    aother_inf = 0.0
    dendAreas = []
    dendAreasInflated = []
    secdiams = []
    seclengths = []
    secdiamsInflated = []
    dendAreasSeg = []
    dendColors = []
    for secname, sec in HR2.sections.items():
        sectype = HR2.find_sec_group(secname)
        sec.nseg = 3
        # sec = HR.sections[secname]
        # print(f"sec: {str(sec):>14s} {sectype:>20s}  n={sec.nseg:>3d} L={sec.L:8.1f}  D={sec.diam:7.3f}  (L/nsec: {sec.L/sec.nseg:>6.2f})")
        if sectype in ["soma"]:
            asoma += areafunc(sec)  # h.area(0.5, sec=sec)
            sec.diam = sec.diam * sinflateratio
            asoma_inf += areafunc(sec)  # h.area(0.5, sec=sec)
            print("soma : sec, d: ", sec, sec.diam, asoma_inf)
        elif sectype in dendrites:
            dendColors.append(sec.nseg)
            dendAreas.append(areafunc(sec))  # h.area(0.5, sec=sec))
            secdiams.append(sec.diam)
            sec.diam = sec.diam * dinflateratio
            secdiamsInflated.append(sec.diam)
            seclengths.append(sec.L)
            dendAreasInflated.append(areafunc(sec))  # h.area(0.5, sec=sec))
            adend_inf += areafunc(sec)  # h.area(0.5, sec=sec)
            print(
                f"dend : {sectype:>24s} sec: {secname:s}, d: {sec.diam:6.3f} sum(A): {adend_inf:9.3f}"
            )
        else:
            aother += areafunc(sec)  # h.area(0.5, sec=sec)
            sec.diam = sec.diam * 1
            aother_inf += areafunc(sec)  # h.area(0.5, sec=sec)
            cprint(
                "r",
                f"other: {sectype:>24s} sec: {secname:s}, d: {sec.diam:6.3f}: sum(A): {aother_inf:9.3f}",
            )

    # HR2.h.define_shape()

    print("\nChanging the areas by changing section diameters")
    print(
        f"soma inflated area: {asoma_inf:10.2f} Infl factor: {asoma_inf/asoma:6.3f} (should be: {sinflateratio:6.3f})"
    )
    print(
        f"dend inflated area: {adend_inf:10.2f} Infl factor: {adend_inf/asoma:6.3f} (should be: {dinflateratio:6.3f})"
    )
    print(
        f"other inflated area: {aother_inf:10.2f} Infl factor: {aother_inf/aother:6.3f} (should be: {1.0:6.3f})"
    )


# def plot_areas(ax, xdend, a1, a2):
dendrites = [
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

somas = [
    "soma",
    "Soma",
]


def recon_hoc():
    cname = "VCN_c17"
    basepath = "/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells"
    cell = f"{cname:s}/Morphology/{cname:s}_Full.hocx"
    fn = Path(basepath, cell)
    print(fn.is_file())
    cconfig = cell_config.CellConfig()
    sinflateratio = cconfig.get_soma_ratio(cname)
    dinflateratio = cconfig.get_dendrite_ratio(cname)

    AdjArea = AdjustAreas(method='pt3d')
    AdjArea.sethoc_fromfile(fn)
    # AdjArea.sethoc_fromstring(hdata=hocstruct2)
    pt3d = AdjArea.adjust_diameters(sectypes = dendrites, inflateRatio=dinflateratio)
    AdjArea.plot_areas(pt3d)


if __name__ == "__main__":
    # simple_test()
    # neuron_example()
    recon_hoc()
