"""
Version 2: read the new format filenames
"""
import sys
import argparse
import pickle
from pathlib import Path
import datetime
import dataclasses
from dataclasses import dataclass, field
import numpy as np
import matplotlib

import vcnmodel.model_params
from ephys.ephysanalysis import MakeClamps
from ephys.ephysanalysis import RmTauAnalysis
from ephys.ephysanalysis import SpikeAnalysis

AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

from matplotlib import rc

rc("text", usetex=False)
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH

modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
# modetypes = ['find', 'singles', 'IO', 'multi']

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

    parser = argparse.ArgumentParser(description="Plot GBC results")
    parser.add_argument(
        dest="cell",
        action="store",
        nargs="+",
        type=str,
        default=None,
        help="Select the cell(s) or 'A' for all(no default)",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        dest="protocol",
        action="store",
        default="IV",
        help=("Select the protocol (default: IV) from: %s" % runtypes),
    )
    # parser.add_argument('-m', '--modetype', dest='modetype', action='store',
    #                 default='multi', help=('Select the mode type  (default: multi) from: %s' % modetypes))
    parser.add_argument(
        "-M",
        "--modeltype",
        dest="modeltype",
        action="store",
        default="XM13_nacncoop",
        help=("Select the model type (default XM13_nacncoop) from: %s " % modeltypes),
    )

    parser.add_argument(
        "-s",
        "--scaled",
        dest="scaled",
        action="store_true",
        help=("use scaled data or not"),
    )
    args = parser.parse_args()
    args.protocol = args.protocol.upper()
    
    # PD.gradeA = [cn for cn in args.cell]
    print('args.cell: ', args.cell)
    if args.cell[0] in ["A", "a"]:
        pass # (all grade a cells defined in pd dataclass)
    else:
        PD.gradeA = [int(c) for c in args.cell]
    rows, cols = PH.getLayoutDimensions(len(PD.gradeA))
    print(PD.gradeA)
    plabels = [f"VCN_c{g:02d}" for g in PD.gradeA]
    for i in range(rows*cols-len(PD.gradeA)):
        plabels.append(f"empty{i:d}")

    P = PH.regular_grid(
        rows,
        cols,
        order='rowsfirst',
        figsize=(8, 10),
        panel_labels=plabels,
        labelposition=(0.05, 0.95),
        margins={
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.12,
            "bottommargin": 0.1,
        },
    )
    chan = "_"  # for certainity in selection, include trailing underscore in this designation

    modelName = args.modeltype
    print('Scaling: ', args.scaled)
    if args.scaled:
        PD.soma_inflate = True
        PD.dend_inflate = True
    else:
        PD.soma_inflate = False
        PD.dend_inflate = False

    stitle = "unknown scaling"
    # trip filemode based on date of simulatoin
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate,"%Y-%m-%d-%H:%M") 
    changetimestamp = datetime.datetime.timestamp(dts)
    for gbc in PD.gradeA:
        basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
        pgbc = f"VCN_c{gbc:02d}"
        if args.protocol == 'IV':
            name_start = f"VCN_c{gbc:02d}_pulse_*.p"
        elif args.protocol == 'AN':
            name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"
            
        print(f"search for:  {name_start:s}")
        fng = list(Path(basefn).glob(name_start))

        times = np.zeros(len(fng))
        for i, f in enumerate(fng):
            times[i] = f.stat().st_mtime
        # pick most recent file = this should be better managed (master file information)
        ix = np.argmax(times)
        fng = [fng[ix]]
            
        ivdatafile = None
        print("\nConditions: soma= ", PD.soma_inflate, "  dend=",PD.dend_inflate)
        for fn in fng: 
            fnp = Path(fn)
            fns = str(fn)

            if not fnp.is_file():
                print(f"file: {str(fnp):s} NOT FOUND")
                continue
            mtime = fnp.stat().st_mtime
            timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
            print(f"pgbcivr2: Checking file: {fnp.name:s} [{timestamp_str:s}]")
            # print(mtime, changetimestamp)
            if mtime > changetimestamp:
                filemode = 'vcnmodel.v1'
            else:
                filemode = 'vcnmodel.v0'
            print('pgbcivr2: file mode: ', filemode)
            with (open(fnp, "rb")) as fh:
                d = pickle.load(fh)
            # print(d.keys())
            # if "Params" not in list(d.keys()):
  #               print("File missing Params; need to re-run")
  #               continue
            # print(d['Results'][0]['inputSpikeTimes'])
            #
            # exit()
            if filemode in ['vcnmodel.v0']:
                # print(d['runInfo'].keys())
                par = d['runInfo']
                par['soma_inflation'] = False
                par['dendrite_inflation'] = False
                if fns.find('soma=') > -1:
                    par['soma_inflation'] = True
                if fns.find('dend=') > -1:
                    par['dendrite_inflation'] = True
                par["soma_autoinflate"] = False
                par["dendrite_autoinflate"] = False
            elif filemode in ['vcnmodel.v1']:
                try:
                    par = d["Params"]
                except:
                    try:
                        par = d["self.Params"]
                    except:
                        raise ValueError("File missing Params; need to re-run")
                if isinstance(par, vcnmodel.model_params.Params):
                    par = dataclasses.asdict(par)
            # print('pgbcivr2: Params: ', par)
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
            elif PD.dend_inflate and not PD.soma_inflate:
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
        print("\npgbcivr2: datafile to read: ", str(ivdatafile))
        if not ivdatafile.is_file():
            print("no file? : ", str(ivdatafile))
            continue
        
        AR.read_pfile(ivdatafile, filemode=filemode)
        bridge_offset = 0.0
        threshold = -32.0
        tgap = 0.0  # gap before fittoign taum
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
        if args.protocol == 'IV':
            SP.fitOne(function="fitOneOriginal")
            RM.analyze(
                rmpregion=[0.0, AR.tstart - 0.001],
                tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
                to_peak=True,
                tgap=tgap,
            )

            RMA = RM.analysis_summary
            # print("Analysis: RMA: ", RMA)
        v0 = -160.
        trstep = 2.5
        inpstep = 0.5
        for trial in range(len(AR.traces)):
            # ds = df["Results"][trial]
            # k0 = list(df["Results"][trial].keys())[0]
            # dx = ds[k0]["monitor"]
            P.axdict[pgbc].plot(AR.time_base, AR.traces[trial]*1e3, linewidth=0.5)
            if args.protocol == 'AN' and 'inputSpikeTimes' in list(d['Results'][trial].keys()):
                spkt = d['Results'][trial]['inputSpikeTimes']
                # print('input spike trains: ', len(spkt))
                tr_y = trial*(trstep + len(spkt)*inpstep) 
                for ian in range(len(spkt)):
                    vy = v0+tr_y*np.ones(len(spkt[ian]))+inpstep*ian
                    P.axdict[pgbc].scatter(spkt[ian], vy, s=12, marker='|')
                    # print(len(vy), vy)
      #               print(spkt[ian])
                P.axdict[pgbc].set_ylim(-140.0, 40.0)
            else:
                P.axdict[pgbc].set_ylim(-200., 50.)
            P.axdict[pgbc].set_xlim(0.080, np.max(AR.time_base))

        PH.calbar(
            P.axdict[pgbc],
            calbar=[np.max(AR.time_base)-25., 20.0, 25.0, 20.0],
            axesoff=True,
            orient="left",
            unitNames={"x": "ms", "y": "mV"},
            font="Arial",
            fontsize=8,
        )
        ftname = str(ivdatafile.name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = f"{ftname:s}"
        if args.protocol == 'IV':
            toptitle += f"\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']:.2f} ms"
        toptitle += f"\n{timestamp_str:s}"
        P.axdict[pgbc].set_title(toptitle, fontsize=5)

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
