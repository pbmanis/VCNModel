from dataclasses import dataclass, field
from pathlib import Path
import toml
import pprint
import numpy as np
from cnmodel import cells
from src.vcnmodel.adjust_areas import AdjustAreas
from src.vcnmodel import cell_config as cell_config

pp = pprint.PrettyPrinter(indent=4)

allcells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
bestcells = [2, 5, 6, 9,  10, 11, 13, 17, 18, 30]


def defemptylist():
    return []


@dataclass
class AxonMorph:
    cellID: str = ''
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


class MeasureAxons(object):
    def __init__(self):
        self.cconfig = cell_config.CellConfig(
                verbose=False,
                spont_mapping='HS',
                )

        # find out where our files live
        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {"cellDataDirectory": Path("../VCN-SBEM-Data", "VCN_Cells")}
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

    def print_table(self, sl, meas):
        h_s = []
        ais_s = []
        myel_s = []
        ax_s = []
        # print(sl)
        for c in sl:
            if meas == 'Length':
                h_s.append(np.sum(sl[c].HillockLength))
                ais_s.append(np.sum(sl[c].AISLength))
                myel_s.append(np.sum(sl[c].MyelLength))
                ax_s.append(np.sum(sl[c].AxonLength))
            if meas == 'AvgDiam':
                h_s.append(np.sum(sl[c].HillockDiamAvg))
                ais_s.append(np.sum(sl[c].AISDiamAvg))
                myel_s.append(np.sum(sl[c].MyelDiamAvg))
                ax_s.append(np.sum(sl[c].AxonDiamAvg))
            if meas == 'BegDiam':
                h_s.append(np.sum(sl[c].HillockDiam0))
                ais_s.append(np.sum(sl[c].AISDiam0))
                myel_s.append(np.sum(sl[c].MyelDiam0))
                ax_s.append(np.sum(sl[c].AxonDiam0))
            if meas == 'EndDiam':
                h_s.append(np.sum(sl[c].HillockDiam1))
                ais_s.append(np.sum(sl[c].AISDiam1))
                myel_s.append(np.sum(sl[c].MyelDiam1))
                ax_s.append(np.sum(sl[c].AxonDiam1))
            print(f"{c:s}   {h_s[-1]:6.2f}", end="")
            print(f"     {ais_s[-1]:6.2f}", end="")
            print(f"     {myel_s[-1]:6.2f}", end="")
            print(f"     {ax_s[-1]:6.2f}")
        print('-'*50)
        print(f"Avg:      {np.mean(h_s):6.2f}     {np.mean(ais_s):6.2f}     {np.mean(myel_s):6.2f}", end="")
        print(f"     {np.mean(ax_s):6.2f}")
        print(f"Std:      {np.std(h_s):6.2f}     {np.std(ais_s):6.2f}     {np.std(myel_s):6.2f}", end="")
        print(f"     {np.std(ax_s):6.2f}")

    def run(self):
        sl = {}
        for cell in bestcells:
            cname = f"VCN_c{cell:02d}"
            cell_hoc = f"{cname}_Full.hoc"  # could be f"{cname}"_Full_standardized_axon.hoc""
            fn = Path(self.baseDirectory, cname, self.morphDirectory, cell_hoc)
            # print(fn.is_file())
            # h.load_file(str(fn))
            # for sec in h.sections:
            #     print(sec, sec.name())
            secs = self.get_axon_measures(cname, fn)
            sl[cname] = secs
        # pp.pprint(sl)
        print("\nLengths: \nCell       Hillock     AIS      Myel      Axon  (um)")
        self.print_table(sl, "Length")

        print("\nAvgDiam: \nCell       Hillock     AIS      Myel      Axon  (um)")
        self.print_table(sl, "AvgDiam")

        print("\nBegDiam: \nCell       Hillock     AIS      Myel      Axon  (um)")
        self.print_table(sl, "BegDiam")

        print("\nEndDiam: \nCell       Hillock     AIS      Myel      Axon  (um)")
        self.print_table(sl, "EndDiam")

    def get_axon_measures(self, cname, fn):
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
        AdjArea.sethoc_fromCNcell(post_cell)
        # AdjArea.sethoc_fromstring(hdata=hocstruct2)
        AdjArea.cell.print_soma_info()
        AdjArea.adjust_diameters(sectypes=AdjArea.somas, inflationRatio=sinflateratio)
        AdjArea.adjust_diameters(
            sectypes=AdjArea.dendrites, inflationRatio=dinflateratio
        )
        AM = AxonMorph(cellID=cname)
        sectypes = []
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            sectypes.append(sectype)
            if sectype not in self.axons:
                continue
            if sectype == 'Axon_Hillock':
                AM.HillockLength.append(sec.L)
                AM.HillockDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.HillockDiam0.append(segs[0].diam)  # end diameters for section
                AM.HillockDiam1.append(segs[-1].diam)
            if sectype in ['Axon_Initial_Segment', 'Unmyelinated_Axon']:
                AM.AISLength.append(sec.L)
                AM.AISDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.AISDiam0.append(segs[0].diam)  # end diameters for section
                AM.AISDiam1.append(segs[-1].diam)
            if sectype == 'Myelinated_Axon':
                AM.MyelLength.append(sec.L)
                AM.MyelDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.MyelDiam0.append(segs[0].diam)  # end diameters for section
                AM.MyelDiam1.append(segs[-1].diam)
            if sectype == 'Axon':
                AM.AxonLength.append(sec.L)
                AM.AxonDiamAvg.append(sec.diam)  # average diameter
                segs = list(sec.allseg())
                AM.AxonDiam0.append(segs[0].diam)  # end diameters for section
                AM.AxonDiam1.append(segs[-1].diam)

                # print('secs: ', secs)
        return AM


if __name__ == '__main__':
    ma = MeasureAxons()
    ma.run()
