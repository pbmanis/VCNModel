from pathlib import Path
import toml
import pprint
import numpy as np
from neuron import h
from cnmodel import cells
from src.vcnmodel.adjust_areas import AdjustAreas
from src.vcnmodel import cell_config as cell_config
from src.vcnmodel import h_reader

pp = pprint.PrettyPrinter(indent=4)

allcells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]

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

    def run(self):
        sl = {}
        for cell in allcells:
            cname = f"VCN_c{cell:02d}"
            cell_hoc = f"{cname}_full.hoc"
            fn = Path(self.baseDirectory, cname, self.morphDirectory, cell_hoc)
            # print(fn.is_file())
            # h.load_file(str(fn))
            # for sec in h.sections:
            #     print(sec, sec.name())
            secs = self.get_axon_lengths(cname, fn)
            sl[cname] = secs
        pp.pprint(sl)
        print(f"Cell       Hillock     AIS      Unmeyl      Myel      Axon  (um)")
        for c in sl:
            print(f"{c:s}   {np.sum(sl[c]['Axon_Hillock']):6.1f}", end="")
            print(f"     {np.sum(sl[c]['Axon_Initial_Segment']):6.1f}", end="")
            print(f"     {np.sum(sl[c]['Unmeylinated_Axon']):6.1f}", end="")
            print(f"     {np.sum(sl[c]['Myelinated_Axon']):6.1f}", end="")
            print(f"     {np.sum(sl[c]['Axon']):6.1f}")
            

    def get_axon_lengths(self, cname, fn):
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
        pt3d = AdjArea.adjust_diameters(sectypes=AdjArea.somas, inflateRatio=sinflateratio)
        pt3d = AdjArea.adjust_diameters(
            sectypes=AdjArea.dendrites, inflateRatio=dinflateratio
        )
        secs = {n:[] for n in self.axons}
        
        for secname, sec in self.HR.sections.items():
            sectype = self.HR.find_sec_group(secname)
            if sectype not in self.axons:
                continue
            secs[sectype] = [sec.L]
        # print('secs: ', secs)

        return secs


if __name__ == '__main__':
    ma = MeasureAxons()
    ma.run()
    