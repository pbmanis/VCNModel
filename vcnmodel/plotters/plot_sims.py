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
import pylibrary.plotting.styler as PLS

modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
experimenttypes = [None, "delays", "largestonly", "removelargest", "mean", "allmean", "twolargest"]
modetypes = ['find', 'singles', 'IO', 'multi']
analysistypes = ['traces', 'PSTH', 'revcorr', 'SAC', 'tuning', 'singles']
dendriteChoices = [
        "normal",
        "passive",
        "active",
    ]

orient_cells = {
        2: [140.,    0.0, -144.0],
        6: [140.,  -59.0,  -12.0],
        5: [140.,  -46.0,  121.0],
        9: [140.,  -74.0,   18.0],
       11: [140.,   -2.0, -181.0],
       10: [140.,    5.0,  -35.0],
       13: [140.,  -22.0,  344.0],
       17: [140., -158.0,   39.0],
       30: [140., -134.0, -181.0],
}
def grAList() -> list:
    """
    Return a list of the 'grade A' cells from the SBEM project
    """
    
    return [2, 5, 6, 9, 10, 11, 13, 17, 30]
    
@dataclass
class PData:
    """
    data class for some parameters that control what we read
    """
    gradeA: list = field(default_factory = grAList)
    default_modelName:str = "XM13_nacncoop"
    soma_inflate:bool = True
    dend_inflate:bool = True
    basepath:str = "/Users/pbmanis/Desktop/Python/VCN-SBEM-Data"
    renderpath:str = "/Users/pbmanis/Desktop/Python/vcnmodel/Renderings"
    thiscell:str=""

def def_empty_np():
    return np.array(0)

def def_empty_list():
    return []

@dataclass
class RevCorrPars:
    ntrials:int =  1
    ninputs:int = 1 

    # clip trace to avoid end effects
    min_time:float = 10.0 # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time:float = 250.0 # this window needs to be at least as long as maxwin
    binw:float = 0.1
    minwin:float = -5
    maxwin:float = 2.5
    amax:float = 0.

@dataclass
class RevCorrData:
    C:list = field(default_factory = def_empty_list)
    TC:list = field(default_factory = def_empty_list)
    st:np.array = field(default_factory = def_empty_np)
    tx:np.array =field(default_factory = def_empty_np)
    ti:np.array =field(default_factory = def_empty_np)
    ti_avg:np.array =field(default_factory = def_empty_np)
    sv_all:np.array = field(default_factory = def_empty_np)
    sv_avg:np.array = field(default_factory = def_empty_np) 
    sites:np.array =  field(default_factory = def_empty_np)
    nsp:int = 0
    max_coin_rate: float=0
    

def norm(p:Union[list, np.ndarray], n:int) -> np.ndarray:
    """
    Simple function to normalize the n'th point of p
    by the min and max
    """
    pmin = np.min(p)
    pmax = np.max(p)
    return (p[n] - pmin) / float(pmax - pmin)

def twinax(fig:object, ax1:object, pos:float=0.) -> object:
    """
    Create a 'twin' axis on the right side of a plot
    Note: pyqtgraph.plotting.styles can also do an inset
    which may be used instead
    """
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

def get_changetimestamp():
    # trip filemode based on date of simulatoin
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate,"%Y-%m-%d-%H:%M") 
    changetimestamp = datetime.datetime.timestamp(dts)
    return(changetimestamp)                
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

def get_data_file(fn:Union[str, Path], changetimestamp:object, PD:dataclass) ->Union[None, tuple]:
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


def analyze_data(ivdatafile:Union[Path, str], filemode, protocol:str) -> tuple: 
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

    # AR.tstart = 1.0
    # AR.tend = 0.4*1e3
    # AR.tdur = 0.399*1e3
    # print(AR.tstart, AR.tend, AR.tdur, AR.sample_rate)
    SP.analyzeSpikes()
    SP.analyzeSpikeShape()

    RMA = None
    if protocol == 'IV':
        # for trial in range(len(AR.traces)):
        #     mpl.plot(AR.time_base, AR.traces[trial]*1e3, linewidth=0.5)
        # mpl.show()
        # exit()
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

def plot_traces(ax:object, fn:Union[Path, str], PD:dataclass, changetimestamp:object, protocol:str) -> None:
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
        ax.plot(AR.time_base, AR.traces[trial]*1e3, linewidth=0.5)
        ax.plot(AR.time_base[SP.spikeIndices[trial]], AR.traces[trial][SP.spikeIndices[trial]]*1e3, 'ro', 
            markersize=2.5)
        if protocol == 'AN' and 'inputSpikeTimes' in list(d['Results'][trial].keys()):
            spkt = d['Results'][trial]['inputSpikeTimes']
            # print('input spike trains: ', len(spkt))
            tr_y = trial*(trstep + len(spkt)*inpstep) 
            for ian in range(len(spkt)):
                vy = v0+tr_y*np.ones(len(spkt[ian]))+inpstep*ian
                ax.scatter(spkt[ian], vy, s=sz,
                marker='|', linewidths=0.35)
                # print(len(vy), vy)
        #                 print(spkt[ian])
            ax.set_ylim(-140.0, 40.0)
        else:
            ax.set_ylim(-200., 50.)
    ax.set_xlim(0.080, np.max(AR.time_base))

    ftname = str(ivdatafile.name)
    ip = ftname.find("_II_") + 4
    ftname = ftname[:ip] + "...\n" + ftname[ip:]
    toptitle = f"{ftname:s}"
    if protocol == 'IV':
        toptitle += f"\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']:.2f} ms"

        secax = PLS.create_inset_axes([0.4, 0, 0.4, 0.4], ax)
        secax.plot(RM.ivss_cmd*1e12, RM.ivss_v*1e3, 'ks-', markersize=3, markerfacecolor='k', zorder=10, clip_on=False)
        ltz = np.where(RM.ivss_cmd <= 0)[0]
        secax.plot(RM.ivss_cmd[ltz]*1e12, RM.ivpk_v[ltz]*1e3, 'ko-', markersize=3, markerfacecolor='w', zorder=10, clip_on=False)
        PH.crossAxes(secax, xyzero=[0., -60.], limits=[-1.0*1e3, -150, 1.0*1e3, -25.])
        PH.talbotTicks(secax, axes='xy',
                density=(1.0, 1.0), insideMargin=0.02, pointSize=6,
                tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
    PH.calbar(ax, calbar=[20., -160., 10., 20.], unitNames={'x': 'ms', 'y': 'mV'}, fontsize=9)
    if RMA is not None:
        PH.referenceline(ax, RMA['RMP'])
        ax.text(-1., RMA['RMP'], f"{RMA['RMP']:.1f}", verticalalignment='center', horizontalalignment='right',
            fontsize=9)
    toptitle += f"\n{timestamp_str:s}"
    ax.set_title(toptitle, fontsize=5)

def plot_revcorr_map(P, pgbc, inputlist, ntrials, C, TC, st, tx, ti_avg, sv_all, sv_avg, sites, nsp, max_coin_rate, window):
    pass

def plot_revcorr2(ax:object, PD:dataclass, RCP:dataclass, RCD:dataclass):
    seaborn.set_style('ticks')
    #secax = twinax(P.figure_handle, ax, pos=maxwin)
    secax = PLS.create_inset_axes([0, 0.5, 1, 0.5], ax)
    PH.noaxes(secax, 'xy')

    secax.set_facecolor((1,1,1,0))
    secax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.set_facecolor((0.7, 0.7, 0.7))
    summarySiteTC = {}
    for isite in range(RCP.ninputs): # range(ninputs):  # for each ANF input (get from those on first trial)
        stepsize = int(RCD.sv_all.shape[0]/20)
        if stepsize > 0:
            sel = list(range(0, RCD.sv_all.shape[0], stepsize))
        else:
            sel = list(range(0, RCD.sv_all.shape[0], 1))
        sel = list(range(0, RCD.sv_all.shape[0], 1))
        # if stelen(sel) < 10:
        #     if sv_all.shape[0] > 10:
        #         sel = list(range(10))
        #     else:
        #         sel = list(range(sv_all.shape[0]))

        if RCD.C[isite] is not None:
            # print('plotting C')
            nc = int(len(RCD.C[isite])/2)
            RCD.TC = RCD.TC/len(RCD.st)
            summarySiteTC[isite] =RCD.TC
            # print(RCD.sites, isite)
            color = mpl.cm.viridis(norm(RCD.sites, isite))
            # print('color', color)
            ax.plot(RCD.tx, RCD.C[isite][:nc], color=color, label=('Input {0:2d} N={1:3d}'.format(isite, int(RCD.sites[isite]))),
                linewidth=1.5, zorder=5)
            RCD.max_coin_rate = np.max((RCD.max_coin_rate, np.max(RCD.C[:nc])))

        if isite == 0:  # only plot the first time through - the APs are the same no matter the trial
            # print('Plotting V')
            for t in sel:
                secax.plot(RCD.ti_avg, RCD.sv_all[t], color="#666666", linewidth=0.2, zorder=1)
                am = np.argmax(RCD.sv_all[t])

            secax.plot(RCD.ti_avg, RCD.sv_avg, color='k', linewidth=0.75, zorder=2)
            PH.calbar(secax, calbar=[1.0, -10, 1.0, 20.0], axesoff=True, orient='right', 
                unitNames={'x': 'ms', 'y': 'mV'} )
            PH.referenceline(secax, -60.)
            seaborn.despine(ax=secax)
    print(f"Total spikes plotted: {RCD.nsp:d}")
    for trial in range(RCP.ntrials):
        stx = RCD.st[trial]
    secax.plot([0., 0.], [-120., 10.], 'r', linewidth=0.5)
    # print('finished inputs')
    seaborn.despine(ax=ax)
    ax.set_ylabel('Rate of coincidences/bin (Hz)', fontsize=10)
    ax.set_xlabel('T (ms)', fontsize=10)
    ax.set_xlim((RCP.minwin, RCP.maxwin))
    if RCD.max_coin_rate > 0.2:
        ax.set_ylim(0, 1)
    else:
        ax.set_ylim(0, 0.25)
    secax.set_ylim([-70., 10.])
    secax.set_xlim((RCP.minwin, RCP.maxwin))
    # secax.set_ylabel('Vm', rotation=-90., fontsize=10)
    secax.tick_params(direction='in', length=5., width=1., labelsize=9)
    ax.tick_params(direction='in', length=5., width=1., labelsize=9)
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
    return summarySiteTC


def get_data(fn:Union[Path, str], PD:dataclass, changetimestamp, protocol):
    
    X = get_data_file(fn, changetimestamp, PD)
    mtime = Path(fn).stat().st_mtime
    timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
    if X is None:
      print("No simulation found that matched conditions")
      return None
    # unpack x
    par, stitle, ivdatafile, filemode, d = X 

    # 2. find spikes
    AR, SP, RMA = analyze_data(ivdatafile, filemode, protocol)
    # set up analysis parameters and result storage
    RCP = RevCorrPars()
    RCD = RevCorrData()    
    RCD.st = SP.spikeIndices
    trials = list(d['Results'].keys())
    RCP.ntrials = len(trials)
    RCD.npost = 0  # number of postsynaptic spikes
    RCD.npre = 0  # number of presynaptic spikes
    for i, tr in enumerate(trials):
        trd = d['Results'][tr] # trial data
        sv = trd['somaVoltage']
        ti = trd['time']
        dt = ti[1]-ti[0]
        # st[tr] = PU.findspikes(ti, sv, thr, dt=dt, mode='threshold')
        # st[tr] = clean_spiketimes(st[tr])
        RCD.npost += len(RCD.st[tr])
        for n in range(len(trd['inputSpikeTimes'])):  # for each sgc
            RCD.npre += len(trd['inputSpikeTimes'][n])
    RCD.ti = ti
    print(f'Detected {RCD.npost:d} Post spikes')
    print(f'Detected {RCD.npre:d} Presynaptic spikes')
    print ("# trials: ", RCP.ntrials)

    # clip trace to avoid end effects
    RCP.max_time = np.max(ti) - RCP.min_time # this window needs to be at least as long as maxwin
    RCD.tx = np.arange(RCP.minwin, 0, RCP.binw)
    summarySiteTC = {}
    return d, AR, SP, RMA, RCP, RCD
    
        
def compute_revcorr(ax, gbc, fn, PD, changetimestamp, protocol,
        thr=-20., width=4.0) -> Union[None, tuple]:
    #
    # 1. Gather data
    #
    SC = cell_config.CellConfig()    
    syninfo = SC.VCN_Inputs[gbc]

    d, AR, SP, RMA, RCP, RCD = get_data(fn, PD, changetimestamp, protocol)
    RCP.ninputs = len(syninfo[1])
    RCD.sites = np.zeros(RCP.ninputs)
    for isite in range(RCP.ninputs): # precompute areas
        area = syninfo[1][isite][0]
        if area > RCP.amax:
            RCP.amax = area
        RCD.sites[isite] = int(np.around(area*SC.synperum2))

    print('ninputs: ', RCP.ninputs)
    maxtc = 0

    RCD.sv_sites = []
    print('inputs:',  RCP.ninputs)
    RCD.C = [None]*RCP.ninputs
    RCD.max_coin_rate = 0.
    nspk_plot = 0
    spksplotted = False
    RCP.min_time = 10. # 220.  # driven window without onset
    RCP.max_time = 300.
    
    for isite in range(RCP.ninputs): # range(ninputs):  # for each ANF input (get from those on first trial)
        # print('isite: ', isite)
        firstwithspikes=False
        RCD.sv_avg = np.zeros(1)
        RCD.sv_all = np.zeros(1)
        RCD.nsp = 0
        n_avg = 0

        for trial in range(RCP.ntrials):  # sum across trials
            stx = AR.time_base[SP.spikeIndices[trial]]
            stx = stx[(stx > RCP.min_time) & (stx < RCP.max_time)]  # more than minimum delay
            anx = d['Results'][trial]['inputSpikeTimes'][isite]
            anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
            if len(stx) == 0 or len(anx) == 0:
               # print('empty array for stx or anx')
                continue
            andirac = np.zeros(int(200./RCP.binw)+1)
            if not firstwithspikes:
                firstwithspikes = True
                RCD.C[isite] = SPKS.correlogram(stx, anx, width=-RCP.minwin, bin=RCP.binw, T=None)
                RCD.TC = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
                if np.isnan(RCD.TC):
                    RCD.TC = 0.
            else:
                RCD.C[isite] += SPKS.correlogram(stx, anx, width=-RCP.minwin, bin=RCP.binw, T=None)
                tct = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
                if ~np.isnan(tct):
                    RCD.TC = RCD.TC + tct

            # accumulate postsynaptic spike waveforms
            if RCD.nsp == 0:  # first spike in trace
                reltime = np.around(RCD.ti,5) - np.around(stx[0], 5)
                areltime = np.argwhere((RCP.minwin <= reltime) & (reltime <=RCP. maxwin)).squeeze()
                RCD.sv_avg = d['Results'][trial]['somaVoltage'][areltime]
                RCD.sv_all = RCD.sv_avg.copy()
                RCD.ti_avg = RCD.ti[0:len(areltime)] + RCP.minwin
                RCD.nsp += 1
    
            else:  # rest of spikes
                for n in range(1,len(stx)):
                    reltime = np.around(RCD.ti, 5) - np.around(stx[n], 5)
                    areltime = np.argwhere((RCP.minwin <= reltime) & (reltime <= RCP.maxwin)).squeeze()
                    if len(areltime) > len(RCD.sv_avg):
                        areltime = areltime[0:len(RCD.sv_avg)]
                    if len(areltime) < len(RCD.sv_avg):
                        nextend = len(RCD.sv_avg)-len(areltime)
                        areltime = np.append(areltime, np.arange(areltime[-1]+1, areltime[-1]+nextend+1))
                    if trial == 0:
                        RCD.sv_avg = d['Results'][trial]['somaVoltage'][areltime]
                        RCD.sv_all = RCD.sv_avg.copy()
                        RCD.ti_avg = RCD.ti[0:len(areltime)] + minwin
                    else:
                        sh = RCD.sv_avg.shape
                        
                        RCD.sv_avg += d['Results'][trial]['somaVoltage'][areltime]
                        RCD.sv_all = np.vstack((RCD.sv_all, d['Results'][trial]['somaVoltage'][areltime]))
                    RCD.nsp += 1

            nspk_plot += RCD.nsp
            RCD.sv_sites.append(RCD.sv_all)
        # RCD.C[isite] /= RCP.ntrials
        if RCD.nsp > 0:
            RCD.sv_avg /= RCD.nsp
            # print('trial: ', i, sv_avg)
                # C = C + SPKS.spike_triggered_average(st[trial]*ms, andirac*ms, max_interval=width*ms, dt=binw*ms)
        if RCD.TC > maxtc:
            maxtc = RCD.TC

    pre_w = [-2.7, -0.7]

    summarySiteTC = plot_revcorr2(ax, PD, RCP, RCD)
    # return summarySiteTC, RCD.sites

################
# now some pariwise, etc. stats on events prior to a spike
################

    tind = np.where((RCD.tx > pre_w[0]) & (RCD.tx < pre_w[1]))[0]
    pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
    participation = np.zeros(RCP.ninputs)
    nperspike = []
    nspikes = 0
    for trial in range(RCP.ntrials):  # accumulate across all trials
        spks = AR.time_base[SP.spikeIndices[trial]]
        for s in spks:  # for each postsynaptic spike
            if s < RCP.min_time or s > RCP.max_time: # trim to those only in a response window
                continue
            # print('pre: ', s)
            nspikes += 1
            nps = 0
            for isite in range(RCP.ninputs):  # test inputs in a window prior
                anxi = d['Results'][trial]['inputSpikeTimes'][isite]  # input spike times for one input
                anxi = anxi[(anxi > RCP.min_time) & (anxi < RCP.max_time)]  # trim to those only in a response window
                ani = anxi - s
                # print('ani: ', ani)
                nevi = len(np.where((ani >= pre_w[0]) & (ani <= pre_w[1]))[0]) # count spikes in ith window
                if nevi > 0:
                    participation[isite] += 1
                    nps += 1
                for jsite in range(isite+1, RCP.ninputs):

                    anxj = d['Results'][trial]['inputSpikeTimes'][jsite]
                    anxj = anxj[(anxj > RCP.min_time) & (anxj < RCP.max_time)]
                    anj = anxj - s
                    nevj = len(np.where((anj >= pre_w[0]) & (anj <= pre_w[1]))[0])
                    if isite != jsite:
                        if nevj > 0 and nevi > 0:
                        # print(nevi, nevi)
                            pairwise[isite, jsite] += 1
                    else:
                        if nevj > 0:
                            pairwise[isite, jsite] += 1
            nperspike.append(nps)

    s_pair = np.sum(pairwise)
    pairwise /= s_pair
    psh = pairwise.shape
    pos = np.zeros((psh[0], psh[1], 2))
    print('# pairwise: ', s_pair, np.sum(pairwise))
    for i in range(RCP.ninputs):
        for j in range(RCP.ninputs):
            # print(f"{pairwise[i,j]:.3f}", end='  ')
            pos[i,j,0] = i+1
            pos[i,j,1] =j+1
        # print()
    # print(nperspike)
    import scipy.stats
    # print(np.unique(nperspike, return_counts=True))
    nperspike = [n for n in nperspike if n != 0]
    nperspike = scipy.stats.itemfreq(nperspike).T
    # print('nperspike counts: ', nperspike)
    # nperspike = np.array(np.unique(nperspike, return_counts=True))/nspikes
    # properly fill out output
    xnspike = np.arange(RCP.ninputs)
    ynspike = np.zeros(RCP.ninputs)
    for j, i in enumerate(nperspike[0]):
        # print(i, j, nperspike[1,j])
        ynspike[i-1] = nperspike[1,j]

    ynspike = np.cumsum(ynspike/nspikes)
        
    # print(RCD.sites)
    # print(pos)
    maxp = np.max(pairwise)
    print(pairwise)
    PSum = PH.regular_grid(rows=2, cols=2, order='rowsfirst', figsize=(6, 6),
         showgrid = False, verticalspacing = 0.08, horizontalspacing = 0.08,
          margins = {'bottommargin': 0.1, 'leftmargin': 0.1, 'rightmargin': 0.1, 'topmargin': 0.1},
          label = ["A", "B", "C", "D"], labelposition = (-0.05, 1.05),)
    sax = PSum.axdict
    # f, sax = mpl.subplots(3,1)
    # f.set_size_inches( w=3.5, h=9)
    sax['A'].plot(np.arange(RCP.ninputs)+1, RCD.sites, 'bo')
    pclip = np.clip(pairwise, np.min(np.where(pairwise > 0)), np.max(pairwise))
    pclip[np.where(pclip == 0)] = np.nan
    pclipped = pclip - np.nanmin(pclip)
    colormap = 'plasma'
    # sax['B'].plot(np.arange(RCP.ninputs)+1, participation/nspikes, 'gx')
    sax['B'].plot(RCD.sites, participation/nspikes, 'gx')
    sax['C'].scatter(pos[:,:,0], pos[:,:,1], s=200*pairwise/maxp, c=pclipped, cmap=colormap)
    cm_sns = mpl.cm.get_cmap(colormap)
    vmax = np.nanmax(pclip)*100
    vmin = np.nanmin(pclip)*100
    print('vmin, vmax: ', vmin, vmax)
    axcbar = PLS.create_inset_axes([0.8, 0.05, 0.05, 0.5], sax['C'])
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    ticks = np.linspace(vmin, vmax, num=4, endpoint=True)
    c2 = matplotlib.colorbar.ColorbarBase(
        axcbar, cmap=cm_sns, ticks=ticks, norm=norm
    )

    # PH.nice_plot(sax['C'], position=-0.2)
    sax['D'].plot(np.arange(RCP.ninputs)+1, ynspike, 'm^-')
    
    sax['A'].set_ylim(bottom=0)
    sax['B'].set_ylim((0, 1.0))
    sax['B'].set_xlim(left=0)
    sax['D'].set_ylim(0, 1.05)
    sax['A'].set_ylabel('# Release Sites')
    sax['A'].set_xlabel('Input #')
    sax['B'].set_xlabel('# Release Sites')
    sax['B'].set_ylabel('Participation')
    sax['C'].set_ylabel('Input #')
    sax['C'].set_xlabel('Input #')
    sax['D'].set_xlabel(f"# Inputs in [{pre_w[0]:.1f} to {pre_w[1]:.1f}] before spike")

    PH.cleanAxes(PSum.axarr.ravel())
    PH.talbotTicks(sax['C'])
    PH.talbotTicks(sax['A'])

    PH.talbotTicks(sax['B'], tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 2})
    PH.talbotTicks(sax['D'], tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 2})
    PH.talbotTicks(axcbar, tickPlacesAdd={'x':0, 'y':2},  floatAdd={'x': 0, 'y': 2}, pointSize=7)
    mpl.show()
    
    return summarySiteTC, RCD.sites


def select_filenames(fng, args)->Union[list, None]:
    # print(fng)
    print('Selct filename for experiment: ', args.experiment)
    match_index = []
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
                
        mmatch = True
        if args.dendriteMode != 'normal':
            srch = f"_{args.dendriteMode:s}_"
            if str(fn).find(srch) < 0:
                mmatch = False

        if dmatch and ematch and rmatch and mmatch:
            match_index.append(ix)
                            
    fng = [fng[i] for i in match_index]
    if len(fng) == 0:
        print('No Files matching conditions found')
        return None             
    # print('found (filtered):\n', '\n   '.join([str(f) for f in fng]))
    return fng


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
        '--dendritemode',
        dest="dendriteMode",
        default="normal",
        choices=dendriteChoices,
        help="Choose dendrite table (normal, active, passive)",
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
    plabels = [f"VCN_c{g:02d}" for g in PD.gradeA]
    for i in range(rows*cols-len(PD.gradeA)):
        plabels.append(f"empty{i:d}")
    
    sizex = cols*3
    sizey = rows*2.5
    if rows * cols < 4:
        sizex *= 2
        sizey *= 2

    chan = "_"  # for certainity in selection, include trailing underscore in this designation

    modelName = args.modeltype
    print('Scaling: ', args.scaled)
    stitle = "unknown scaling"
    if args.scaled:
        PD.soma_inflate = True
        PD.dend_inflate = True
        stitle="scaled"
    else:
        PD.soma_inflate = False
        PD.dend_inflate = False
        stitle= "notScaled"
    # trip filemode based on date of simulatoin
    changedate = "2020-04-29-12:00"
    dts = datetime.datetime.strptime(changedate,"%Y-%m-%d-%H:%M") 
    changetimestamp = datetime.datetime.timestamp(dts)
    for ig, gbc in enumerate(PD.gradeA):
        basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
        pgbc = f"VCN_c{gbc:02d}"
        if args.protocol == 'IV':
            name_start = f"IV_Result_VCN_c{gbc:02d}_*.p"
            args.experiment = None
        elif args.protocol == 'AN':
            name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"
            
        print(f"Searching for:  {str(Path(basefn, name_start)):s}")
        fng = list(Path(basefn).glob(name_start))
        exit()
        """ cull list by experiment and db """


        # print(match_index)
        if args.check:
            return
        ivdatafile = None

        print("\nConditions: soma= ", PD.soma_inflate, "  dend=",PD.dend_inflate)
        if args.analysis in ['traces', 'revcorr']:
            fng = select_filenames(fng, args)
            times = np.zeros(len(fng))
            for i, f in enumerate(fng):
                times[i] = f.stat().st_mtime
            # pick most recent file = this should be better managed (master file information)
      
            ix = np.argmax(times)
            fng = [fng[ix]]
            if ig == 0:
                P = PH.regular_grid(
                    rows,
                    cols,
                    order='rowsfirst',
                    figsize=(sizex, sizey),
                    panel_labels=plabels,
                    labelposition=(0.05, 0.95),
                    margins={
                        "leftmargin": 0.1,
                        "rightmargin": 0.01,
                        "topmargin": 0.15,
                        "bottommargin": 0.15,
                        },
                )
            
            ax = P.axdict[pgbc]
            for fn in fng:
                if args.analysis == "traces":
                    plot_traces(ax, fn, PD, changetimestamp, args.protocol)
                elif args.analysis == "revcorr":
                    compute_revcorr(ax, pgbc, fn, PD, changetimestamp, args.protocol)
        elif args.analysis == 'singles':
            fna = select_filenames(fng, args)
            print(fna)
            
        elif args.analysis == "tuning":
            channel = 'nacncoop'
            cols = 3
            rows = 2
            sizex = cols*3
            sizey = rows*2.5
            P = PH.regular_grid(
                rows,
                cols,
                order='rowsfirst',
                figsize=(sizex, sizey),
                # panel_labels=plabels,
                labelposition=(0.05, 0.95),
                margins={
                    "leftmargin": 0.1,
                    "rightmargin": 0.01,
                    "topmargin": 0.15,
                    "bottommargin": 0.15,
                },
            )
            
            print('...')

            truefile = {'Passive': 'passive', 'Canonical': 'normal', 'Active': 'active'}
            for ic, chdist in enumerate(['Passive', 'Canonical', 'Active']):
                args.dendritemode = truefile[chdist]
                fna = select_filenames(fng, args)
                print('\n plot data from: '.join([str(f) for f in fna]))

                for k, fn in enumerate(fna):
                    plot_traces(P.axarr[0, k], fn, PD, changetimestamp, args.protocol)
                
                P.axarr[0,ic].text(50, 1.0, chdist, horizontalalignment='center')
                basen = fn.parts[-4]
                    
                pngfile = Path(PD.renderpath, f"{basen:s}_{channel:s}_{truefile[chdist]:s}.png")
                print(pngfile)
                imaged = mpl.imread(pngfile)
                P.axarr[1, ic].imshow(imaged)
                

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
