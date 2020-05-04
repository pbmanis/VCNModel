"""
Version 2: read the new format filenames
"""
import sys
import pickle
from pathlib import Path
from dataclasses import dataclass, field
import datetime
import matplotlib
from matplotlib import rc

rc("text", usetex=False)
import matplotlib.pyplot as mpl

from ephys.ephysanalysis import MakeClamps
from ephys.ephysanalysis import RmTauAnalysis
from ephys.ephysanalysis import SpikeAnalysis

AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

import pylibrary.plotting.plothelpers as PH

def grAList():
    return [2, 5, 6, 9, 10, 11, 13, 17, 24, 29, 30]
@dataclass
class PData:
    gradeA: list = field(default_factory = grAList)
    default_modelName:str = "XM13_nacncoop"
    soma_inflate:bool = True
    dend_inflate:bool = True
    basepath:str = "/Users/pbmanis/Desktop/Python/VCN-SBEM-Data"
    
def main():
    PD = PData()
    gbc_names = [f"VCN_c{g:02d}" for g in PD.gradeA]
    (r, c) = PH.getLayoutDimensions(len(gbc_names), pref="height")
    if r * c > len(gbc_names):
        for i in range(len(gbc_names), r * c):
            gbc_names.append("ne")
    P = PH.regular_grid(
        r,
        c,
        figsize=(10, 8),
        order="rowsfirst",
        panel_labels=gbc_names,
        labelposition=(0.05, 0.95),
        margins={
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.12,
            "bottommargin": 0.1,
        },
    )
    chan = "_"  # for certainity in selection, include trailing underscore in this designation


    # default_modelName = 'XM13'
    if len(sys.argv) > 1:
        modelName = sys.argv[1]
    else:
        modelName = PD.default_modelName
    if len(modelName) > 4:
        chan = modelName[4:] + "_"
        modelName = modelName[:4]

    stitle = "unknown scaling"


    for gbc in gbc_names:
        if gbc == "ne":
            continue
        cellpath = f"VCN_Cells/{gbc:2s}"
        cell_ANData = Path(PD.basepath, cellpath, f"Simulations/AN/")
        # print('CellANdata: ', cell_ANData)# ivdatafile = Path(f"VCN_Cells/VCN_c{gbc:2s}/Simulations/IV/VCN_c{gbc:2s}_pulse_XM13{chan:s}_II_PD.soma_inflate=*_monitor.p")
        fng = list(cell_ANData.glob(f"AN_Result_{gbc:s}*.p"))
        # fng = ['VCN_Cells/VCN_c09/Simulations/AN/AN_Result_VCN_c09_2019-10-28-163743-342545.p']
        ivdatafile = None
        # print("\nConditions:  PD.soma_inflate oma_inflate= ", PD.soma_inflate, "  dend=", dend)
        # print('fng: ', fng)
        for fn in fng:
            fnp = Path(fn)
            with (open(fnp, "rb")) as fh:
                d = pickle.load(fh)
            par = d["Params"]
            # print('Params:   ', d["Params"])
            # print('    File: ', fn)
            # print('    Par:  ', par)
            if PD.soma_inflate and PD.dend_inflate:
                if par["soma_inflation"] > 0 and par["dendrite_inflation"] > 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Soma and Dend scaled"
                    print(stitle)
                    break
            elif PD.soma_inflate and not PD.dend_inflate:
                if par["soma_inflation"] > 0 and par["dendrite_inflation"] < 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Soma only scaled"
                    print(stitle)
                    break
            elif PD.dend_inflate and not PD.soma_inflatea:
                if par["soma_inflation"] < 0 and par["dendrite_inflation"] > 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Dend only scaled"
                    print(stitle)
                    break
            elif not par["soma_autoinflate"] and not par["dendrite_autoinflate"]:
                print("\nConditions x: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)
                ivdatafile = Path(fn)
                stitle = "No scaling (S, D)"
                print(stitle)
                break
            else:
                print("continue")
                continue
        if ivdatafile is None:
            print("No simulation found that matched conditions")
            print(fng)
            continue
        print("datafile to read: ", str(ivdatafile))
        if not ivdatafile.is_file():
            print("no file? : ", str(ivdatafile))
            continue
        # ftime = ivdatafile.stat().st_mtime
        # print('  File time: ', datetime.datetime.fromtimestamp(ftime).strftime('%H:%M:%S %d-%b-%Y'))
        print('ivdatafile: ', ivdatafile)
        AR.read_pfile(ivdatafile)
        # for i in range(AR.traces.shape[0]):
        #     mpl.plot(AR.time_base, AR.traces[i])
        # mpl.show()
        bridge_offset = 0.0
        threshold = -32.0
        tgap = 0.0  # gap before fittoign taum
        # print(self.AR.tstart, self.AR.tend)
        RM.setup(AR, SP, bridge_offset=bridge_offset)
        SP.setup(
            clamps=AR,
            threshold=threshold,
            refractory=0.0001,
            peakwidth=0.001,
            interpolate=True,
            verify=False,
            mode="peak",
        )
        SP.analyzeSpikes()
        SP.analyzeSpikeShape()
        # SP.analyzeSpikes_brief(mode='baseline')
        # SP.analyzeSpikes_brief(mode='poststimulus')
        SP.fitOne(function="fitOneOriginal")
        RM.analyze(
            rmpregion=[0.0, AR.tstart - 0.001],
            tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
            to_peak=True,
            tgap=tgap,
        )

        RMA = RM.analysis_summary
        print(RMA)
        with (open(ivdatafile, "rb")) as fh:
            df = pickle.load(fh)
        print("df.keys: ", df.keys())
        r = df["Results"][0]

        for trial in range(len(df["Results"])):
            ds = df["Results"][trial]
            k0 = list(df["Results"][trial].keys())[0]
            dx = ds[k0]["monitor"]
            P.axdict[gbc].plot(dx["time"], dx["postsynapticV"], linewidth=1.0)
            P.axdict[gbc].set_xlim(0.0, 150.0)
            P.axdict[gbc].set_ylim(-200.0, 50.0)
        PH.calbar(
            P.axdict[gbc],
            calbar=[120.0, -95.0, 25.0, 20.0],
            axesoff=True,
            orient="left",
            unitNames={"x": "ms", "y": "mV"},
            font="Arial",
            fontsize=8,
        )
        ftname = str(ivdatafile.name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = f"{ftname:s}"  # "\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']*1e3:.2f} ms"
        P.axdict[gbc].set_title(toptitle, fontsize=5)

    if chan == "_":
        chan = "nav11"
    else:
        chan = chan[:-1]
    P.figure_handle.suptitle(
        f"Model: {modelName:s}  Na Ch: {chan:s} Scaling: {stitle:s} ", fontsize=7
    )
    mpl.show()


if __name__ == "__main__":
    main()
