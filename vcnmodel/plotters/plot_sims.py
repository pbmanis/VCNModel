"""
reads the new format filenames

Wrapper for various analysis functions, handles multiple cells.

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
import seaborn
from typing import Union
import vcnmodel.model_params
import vcnmodel.cell_config as cell_config
import vcnmodel.spikestatistics as SPKS
from ephys.ephysanalysis import MakeClamps
from ephys.ephysanalysis import RmTauAnalysis
from ephys.ephysanalysis import SpikeAnalysis
from ephys.ephysanalysis import Utility as EU
# import pylibrary.tools.utility as PU
AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

from matplotlib import rc

rc("text", usetex=False)
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH

modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
experimenttypes = [None, "delays", "largestonly", "removelargest", "mean", "allmean", "twolargest"]
# modetypes = ['find', 'singles', 'IO', 'multi']
analysistypes = ['traces', 'PSTH', 'revcorr', 'SAC']

def grAList():
    return [2, 5, 6, 9, 10, 11, 13, 17, 30]
    
@dataclass
class PData:
    gradeA: list = field(default_factory = grAList)
    default_modelName:str = "XM13_nacncoop"
    soma_inflate:bool = True
    dend_inflate:bool = True
    basepath:str = "/Users/pbmanis/Desktop/Python/VCN-SBEM-Data"

def norm(p, n):
    pmin = np.min(p)
    pmax = np.max(p)
    return (p[n] - pmin) / float(pmax - pmin)

def twinax(fig, ax1, pos=0.):
    ax2 = fig.add_axes(ax1.get_position(True), sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_offset_position('right')
    ax2.tick_params(direction='in', length=5., width=1., labelsize=6)
    ax2.spines['right'].set_position(('data', pos))
    ax2.spines['left'].set_color('none')
    ax2.spines['top'].set_color('none')
    #ax2.set_autoscalex_on(ax1.get_autoscalex_on())
    #ax1.yaxis.tick_left()
    ax2.xaxis.set_visible(False)
    ax2.patch.set_visible(False)
#    PH.adjust_spines(ax2, distance=0.)
    return ax2
                
def clean_spiketimes(spikeTimes, mindT=0.7):
    """
    Clean up spike time array, removing all less than mindT
    spikeTimes is a 1-D list or array
    mindT is difference in time, same units as spikeTimes
    If 1 or 0 spikes in array, just return the array
    """
    if len(spikeTimes) > 1:
        dst = np.diff(spikeTimes)
        st = np.array(spikeTimes[0])  # get first spike
        sok = np.where(dst > mindT)
        st = np.append(st, [spikeTimes[s+1] for s in sok])
        # print st
        spikeTimes = st[~np.isnan(st)]
    return spikeTimes

def get_data_file(fn, changetimestamp, PD) ->Union[None, tuple]:
    """
    Get a data file, and also parse information from the file
    for display
    """
    fnp = Path(fn)
    fns = str(fn)
    ivdatafile = None
    if not fnp.is_file():
      print(f"file: {str(fnp):s} NOT FOUND")
      return None
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
        if par["soma_inflation"] > 0.0 and par["dendrite_inflation"] > 0.0:
            ivdatafile = Path(fn)
            stitle = "Soma and Dend scaled"
            print(stitle)
    elif PD.soma_inflate and not PD.dend_inflate:
      if par["soma_inflation"] > 0.0 and par["dendrite_inflation"] < 0.0:
          ivdatafile = Path(fn)
          stitle = "Soma only scaled"
          print(stitle)
    elif PD.dend_inflate and not PD.soma_inflate:
      if par["soma_inflation"] < 0.0 and par["dendrite_inflation"] > 0.0:
          ivdatafile = Path(fn)
          stitle = "Dend only scaled"
          print(stitle)
    elif not par["soma_autoinflate"] and not par["dendrite_autoinflate"]:
      print("\nConditions x: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)
      ivdatafile = Path(fn)
      stitle = "No scaling (S, D)"
      print(stitle)
    else:
        return None

    if ivdatafile is None:
      print("no file matching conditions : ", str(ivdatafile))
      return None

    if not ivdatafile.is_file():
      print("no file? : ", str(ivdatafile))
      return None
    
    print("\npgbcivr2: datafile to read: ", str(ivdatafile))
    return par, stitle, ivdatafile, filemode, d


def analyze_data(ivdatafile, filemode, protocol): 
    """
    Provide basic spike detection, shape analysis, and
    IV analysis if appropriate
    We use ephys.acq4read.read_pfile to read the pickled data
    file into acq4 format, then we can use the ephys
    analysis tools to analyze the data
    """   
    AR.read_pfile(ivdatafile, filemode=filemode)
    bridge_offset = 0.0
    threshold = -32.0 # mV
    tgap = 0.0  # gap before fittoign taum
    RM.setup(AR, SP, bridge_offset=bridge_offset)
    SP.setup(
      clamps=AR,
      threshold=threshold,
      refractory=0.0001,
      peakwidth=0.002,
      interpolate=True,
      verify=True,
      mode="peak",
    )
    
    SP.set_detector('Kalluri')  # spike detector

    AR.tstart = 1.0
    AR.tend = 0.4*1e3
    AR.tdur = 0.399*1e3
    # print(AR.tstart, AR.tend, AR.tdur, AR.sample_rate)
    SP.analyzeSpikes()
    SP.analyzeSpikeShape()

    RMA = None
    if protocol == 'IV':
      SP.fitOne(function="fitOneOriginal")
      RM.analyze(
          rmpregion=[0.0, AR.tstart - 0.001],
          tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
          to_peak=True,
          tgap=tgap,
      )

      RMA = RM.analysis_summary
      # print("Analysis: RMA: ", RMA)

    return AR, SP, RMA

def plot_traces(P, pgbc, fn, PD, changetimestamp, protocol):
    x = get_data_file(fn, changetimestamp, PD)
    mtime = Path(fn).stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
    if x is None:
      print("No simulation found that matched conditions")
      print(fng)
      return
    # unpack x
    par, stitle, ivdatafile, filemode, d = x 
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
    ntr = len(AR.traces) # number of trials
    v0 = -160.
    trstep = 25.0/ntr
    inpstep = 5.0/ntr
    sz = 50./ntr
    EUD = EU.Utility()
    for trial in range(len(AR.traces)):
        P.axdict[pgbc].plot(AR.time_base, AR.traces[trial]*1e3, linewidth=0.5)
        P.axdict[pgbc].plot(AR.time_base[SP.spikeIndices[trial]], AR.traces[trial][SP.spikeIndices[trial]]*1e3, 'ro', 
            markersize=2.5)
        if protocol == 'AN' and 'inputSpikeTimes' in list(d['Results'][trial].keys()):
            spkt = d['Results'][trial]['inputSpikeTimes']
            # print('input spike trains: ', len(spkt))
            tr_y = trial*(trstep + len(spkt)*inpstep) 
            for ian in range(len(spkt)):
                vy = v0+tr_y*np.ones(len(spkt[ian]))+inpstep*ian
                P.axdict[pgbc].scatter(spkt[ian], vy, s=sz,
                marker='|', linewidths=0.35)
                # print(len(vy), vy)
        #                 print(spkt[ian])
            P.axdict[pgbc].set_ylim(-140.0, 40.0)
        else:
            P.axdict[pgbc].set_ylim(-200., 50.)
    P.axdict[pgbc].set_xlim(0.080, np.max(AR.time_base))

    ftname = str(ivdatafile.name)
    ip = ftname.find("_II_") + 4
    ftname = ftname[:ip] + "...\n" + ftname[ip:]
    toptitle = f"{ftname:s}"
    if protocol == 'IV':
      toptitle += f"\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']:.2f} ms"
    toptitle += f"\n{timestamp_str:s}"
    P.axdict[pgbc].set_title(toptitle, fontsize=5)


def plot_revcorr(P, pgbc, fn, PD, changetimestamp, protocol,
        thr=-20., width=4.0):
    x = get_data_file(fn, changetimestamp, PD)
    mtime = Path(fn).stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
    if x is None:
      print("No simulation found that matched conditions")
      print(fng)
      return
    # unpack x
    par, stitle, ivdatafile, filemode, d = x 
    
    seaborn.set_style('ticks')
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)

    SC = cell_config.CellConfig()
    syninfo = SC.VCN_Inputs[pgbc]

    st = SP.spikeIndices
    trials = list(d['Results'].keys())
    npost = 0  # number of postsynaptic spikes
    npre = 0  # number of presynaptic spikes
    for i, tr in enumerate(trials):
        trd = d['Results'][tr] # trial data
        sv = trd['somaVoltage']
        ti = trd['time']
        dt = ti[1]-ti[0]
        # st[tr] = PU.findspikes(ti, sv, thr, dt=dt, mode='threshold')
        # st[tr] = clean_spiketimes(st[tr])
        npost += len(st[tr])
        for n in range(len(trd['inputSpikeTimes'])):  # for each sgc
            npre += len(trd['inputSpikeTimes'][n])
    ntrials = len(trials)

    print(f'Detected {npost:d} Post spikes')
    print(f'Detected {npre:d} Presynaptic spikes')
    print ("# trials: ", ntrials)
    ninputs = len(syninfo[1])
    sites = np.zeros(ninputs)
    print('ninputs: ', ninputs)
    # clip trace to avoid end effects
    min_time = 10.0 # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time = np.max(ti) - min_time # this window needs to be at least as long as maxwin
    binw = 0.1
    minwin = -5
    maxwin = 2.5
    window = (minwin, maxwin)
    tcwidth = -minwin # msec for total correlation width
    xwidth = -minwin
    tx = np.arange(-xwidth, 0, binw)
    amax = 0.

    for isite in range(ninputs): # precompute areas
        area = syninfo[1][isite][0]
        if area > amax:
            amax = area
        sites[isite] = int(np.around(area*SC.synperum2))
        
    summarySiteTC = {}
    maxtc = 0
    ax = P.axdict[pgbc]
    secax = twinax(P.figure_handle, ax, pos=maxwin)
    sv_sites = []
    C = None
    nspk_plot = 0

    inputlist = range(ninputs)
    spksplotted = False
    for isite in inputlist: # range(ninputs):  # for each ANF input (get from those on first trial)
        # print('isite: ', isite)
        firstwithspikes=False
        sv_avg = np.zeros(1)
        sv_all = np.zeros(1)
        nsp = 0
        n_avg = 0

        for trial in range(ntrials):  # sum across trials
            # print('   trial: ', trial)
            # ax.plot(d['Results'][trial]['time'], d['Results'][trial]['somaVoltage'])
 #            continue
             # stx = st[trial]
            # stx = d['Results'][trial]['spikeTimes']  # get postsynaptic spike train for this trial measured during simulatoin
            stx = AR.time_base[SP.spikeIndices[trial]]
            stx = stx[(stx > min_time) & (stx < max_time)]  # more than minimum delay
            anx = d['Results'][trial]['inputSpikeTimes'][isite]
            anx = anx[(anx > min_time) & (anx < max_time)]
            if len(stx) == 0 or len(anx) == 0:
                print('empty array for stx or anx')
                continue
            andirac = np.zeros(int(200./binw)+1)
            if not firstwithspikes:
                firstwithspikes = True
                C = SPKS.correlogram(stx, anx, width=xwidth, bin=binw, T=None)
                TC = SPKS.total_correlation(anx, stx, width=tcwidth, T=None)
                if np.isnan(TC):
                    TC = 0.
                # definition: spike_triggered_average(spikes,stimulus,max_interval,dt,onset=None,display=False):
                # C = SPKS.spike_triggered_average(st[trial]*ms, andirac*ms, max_interval=width*ms, dt=binw*ms)
            else:
                C = C + SPKS.correlogram(stx, anx, width=xwidth, bin=binw, T=None)
                tct = SPKS.total_correlation(anx, stx, width=tcwidth, T=None)
                if ~np.isnan(tct):
                    TC = TC + tct

            # accumulate postsynaptic spikes
            if nsp == 0:  # first spike in trace
                reltime = np.around(ti,5) - np.around(stx[0], 5)
                areltime = np.argwhere((minwin <= reltime) & (reltime <= maxwin)).squeeze()
                sv_avg = d['Results'][trial]['somaVoltage'][areltime]
                sv_all = sv_avg.copy()
                ti_avg = ti[0:len(areltime)] + minwin
                nsp += 1
    
            else:  # rest of spikes
                for n in range(1,len(stx)):
                    reltime = np.around(ti, 5) - np.around(stx[n], 5)
                    areltime = np.argwhere((minwin <= reltime) & (reltime <= maxwin)).squeeze()
                    if len(areltime) > len(sv_avg):
                        areltime = areltime[0:len(sv_avg)]
                    if len(areltime) < len(sv_avg):
                        nextend = len(sv_avg)-len(areltime)
                        areltime = np.append(areltime, np.arange(areltime[-1]+1, areltime[-1]+nextend+1))
                    if trial == 0:
                        sv_avg = d['Results'][trial]['somaVoltage'][areltime]
                        sv_all = sv_avg.copy()
                        ti_avg = ti[0:len(areltime)] + minwin
                    else:
                        sv_avg += d['Results'][trial]['somaVoltage'][areltime]
                        sv_all = np.vstack((sv_all, d['Results'][trial]['somaVoltage'][areltime]))
                    nsp += 1

            nspk_plot += nsp
            sv_sites.append(sv_all)
        if nsp > 0:
            sv_avg /= nsp
            # print('trial: ', i, sv_avg)
                # C = C + SPKS.spike_triggered_average(st[trial]*ms, andirac*ms, max_interval=width*ms, dt=binw*ms)
        if C is not None:
            nc = int(len(C)/2)
            TC = TC/len(st)
            summarySiteTC[isite] = TC
            color = mpl.cm.viridis(norm(sites, isite))
            ax.set_facecolor((0.7, 0.7, 0.7))
            ax.plot(tx, C[:nc], color=color, label=('Input {0:2d} N={1:3d}'.format(isite, int(sites[isite]))),
                linewidth=0.75)
        if TC > maxtc:
            maxtc = TC
        # only plot 1/10 of input spikes
        stepsize = int(sv_all.shape[0]/20)
        if stepsize > 0:
            sel = list(range(0, sv_all.shape[0], stepsize))
        else:
            sel = list(range(0, sv_all.shape[0], 1))
            
        # if stelen(sel) < 10:
        #     if sv_all.shape[0] > 10:
        #         sel = list(range(10))
        #     else:
        #         sel = list(range(sv_all.shape[0]))
        for t in sel:
            secax.plot(ti_avg, sv_all[t], color='k', linewidth=0.2)
            am = np.argmax(sv_all[t])

        secax.plot(ti_avg, sv_avg, color='r', linewidth=0.75)
                # secax.text(ti_avg[am], np.max(sv_all[t][am])+t*0.5, f"{sv_ntr[t]:.3f}")
            
        # tx2 = np.array([0.2, 0.8])
        # ax2.plot(tx2, TC*np.ones_like(tx2), color=color, linewidth=2)

    print(f"Total spikes plotted: {nsp:d}")
    for trial in range(ntrials):
        stx = st[trial]
    secax.plot([0., 0.], [-120., 10.], 'r', linewidth=0.5)
    print('finished inputs')
    seaborn.despine(ax=ax)
    ax.set_ylabel('Rate of coincidences/bin (Hz)', fontsize=12)
    ax.set_xlabel('T (ms)', fontsize=12)
    ax.set_xlim(window)
    ax.set_ylim(0, 1)
    secax.set_ylim([-120., 10.])
    secax.set_xlim(window)
    secax.set_ylabel('Vm', rotation=-90.)    
    secax.tick_params(direction='in', length=5., width=1., labelsize=12)
    # ax2.set_ylabel('Total Correlation W=%.1f-%0.1f'% (tcwidth[0], tcwidth[1]), fontsize=12)
    # ax2.set_ylim(0, 1.0) # maxtc*1.2)
    # PH.talbotTicks(ax2, axes='xy',
    #                density=(1.0, 1.0), insideMargin=0.05, pointSize=10,
    #                tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1})
                   
    # a = re_self.search(fn[p])
    # b = re_c10.search(fn[p])
    # if a is not None:
    #     inp = 'VCN_c09'
    # elif b is not None:
    #     inp = 'VCN_c10'
    # else:
    #     inp = "Undefined"
    # ax.set_title(f'VCN_{p:s} input={inp:s} [{int(np.min(sites)):d}-{int(np.max(sites)):d}]\nAmax={amax:.1f}', y=0.9, x=0.02,
    #     horizontalalignment='left', fontsize=6)
    return summarySiteTC, sites

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
        choices=runtypes,
        help=("Select the protocol (default: IV) from: %s" % runtypes),
    )
    parser.add_argument(
        "-a",
        "--analysis",
        dest="analysis",
        action="store",
        default=analysistypes[0],
        choices=analysistypes,
        help=(f"Select the analysis type (default: {analysistypes[0]:s}) from: {str(analysistypes):s}" ),
    )
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
    
    parser.add_argument(
        "-e",
        "--experiment",
        dest="experiment",
        action="store",
        default="delays",
        choices=experimenttypes,
        help=("Select the experiment type from: %s " % experimenttypes),
    )
    
    parser.add_argument(
        "-d",
        "--dB",
        dest="dbspl",
        type=float,
        action="store",
        default=None,
        help=("Select the models at specific intensity"),
    )

    parser.add_argument(
        "-r",
        "--nreps",
        dest="nreps",
        type=int,
        action="store",
        default=None,
        help=("Select the models with # reps"),
    )
    
    parser.add_argument(
        "-c",
        "--check",
        dest="check",
        action="store_true",
        help=("Just check selection criteria and return"),
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
        figsize=(14, 10),
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
            name_start = f"IV_Result_VCN_c{gbc:02d}_*.p"
            args.experiment = None
        elif args.protocol == 'AN':
            name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"
            
        print(f"Searching for:  {str(Path(basefn, name_start)):s}")
        fng = list(Path(basefn).glob(name_start))
        
        """ cull list by experiment and db """
        match_index = []
        print(fng)
        print(args.experiment)
        for ix, fn in enumerate(fng):
            ematch = True
            if args.experiment is not None:
                if str(fn).find(args.experiment) < 0:
                    ematch = False

            dmatch = True
            if args.dbspl is not None:
                srch = f"_{int(args.dbspl):03d}dB_"
                if str(fn).find(srch) < 0:
                    dmatch = False

            rmatch = True
            if args.nreps is not None:
                srch = f"_{int(args.nreps):03d}_"
                if str(fn).find(srch) < 0:
                    rmatch = False
         
            if dmatch and ematch and rmatch:
                match_index.append(ix)
        fng = [fng[i] for i in match_index]
        if len(fng) == 0:
            print('No Files matching conditions found')
            continue                
        print('found (filtered):\n', '\n   '.join([str(f) for f in fng]))
        # print(match_index)
        if args.check:
            return

        times = np.zeros(len(fng))
        for i, f in enumerate(fng):
            times[i] = f.stat().st_mtime
        # pick most recent file = this should be better managed (master file information)
        ix = np.argmax(times)
        fng = [fng[ix]]
        
        ivdatafile = None

        print("\nConditions: soma= ", PD.soma_inflate, "  dend=",PD.dend_inflate)
        for fn in fng:
            if args.analysis == "traces":
                plot_traces(P, pgbc, fn, PD, changetimestamp, args.protocol)
            elif args.analysis == "revcorr":
                plot_revcorr(P, pgbc, fn, PD, changetimestamp, args.protocol)

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
