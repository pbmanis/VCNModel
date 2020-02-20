import sys
import pickle
from pathlib import Path
import datetime
import matplotlib

from ephys.ephysanalysis import MakeClamps
from ephys.ephysanalysis import RmTauAnalysis
from ephys.ephysanalysis import SpikeAnalysis
AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

from matplotlib import rc
rc('text', usetex=False)
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH


def main():
    gbc_names = ['08', '09', '11', '14', '16', '17', '18', '19', '20', '21', '22', 'None']
    gbc_names = ['09', '11', '17', '18', 'None']
    P = PH.regular_grid(2, 2, figsize=(6,8), panel_labels=gbc_names, labelposition=(0.05, 0.95),
        margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.12, 'bottommargin': 0.1},)
    chan = '_'  # for certainity in selection, include trailing underscore in this designation

    default_modelName = 'XM13_nacncoop'
    #default_modelName = 'XM13'
    if len(sys.argv) > 1:
        modelName = sys.argv[1]
    else:
        modelName = default_modelName
    if len(modelName) > 4:
        chan = modelName[4:]+'_'
        modelName = modelName[:4]
    soma = True
    dend = True
    stitle="unknown scaling"

    for gbc in gbc_names:
        basefn = f"VCN_Cells/VCN_c{gbc:2s}/Simulations/IV/VCN_c{gbc:2s}_pulse_{modelName:s}{chan:s}"
        # ivdatafile = Path(f"VCN_Cells/VCN_c{gbc:2s}/Simulations/IV/VCN_c{gbc:2s}_pulse_XM13{chan:s}_II_soma=*_monitor.p")
        fng = list(Path('.').glob(basefn+'*.p'))

        ivdatafile = None
        print('\nConditions: soma= ', soma, '  dend=', dend)
        for fnp in fng:
            fn = str(fnp)
            if soma and dend:
                if fn.find('soma=') > 0 and fn.find('dend=') > 0.:
                    ivdatafile = Path(fn)
                    stitle = 'Soma and Dend scaled'
                    print(stitle)
                    break
            elif soma and not dend:
                if fn.find('soma=') > 0 and fn.find('dend=') < 0.:
                    ivdatafile = Path(fn)
                    stitle = 'Soma only scaled'
                    print(stitle)
                    break
            elif dend and not soma:
                if fn.find('soma=') < 0 and fn.find('dend=') > 0.:
                    ivdatafile = Path(fn)
                    stitle = 'Dend only scaled'
                    print(stitle)
                    break
            elif not soma and not dend and fn.find('soma=') == -1 and fn.find('dend=') == -1:
                print('\nConditions x: soma= ', soma, '  dend=', dend)
                ivdatafile = Path(fn)
                stitle = 'No scaling (S, D)'
                print(stitle)
                break
            else:
                print('continue')
                continue
        if ivdatafile is None:
            print('No simulation found that matched conditions')
            print(fng)
            continue
        print('datafile to read: ', str(ivdatafile))
        if not ivdatafile.is_file():
            print('no file? : ', str(ivdatafile))
            continue
        ftime = ivdatafile.stat().st_mtime
        print('  File time: ', datetime.datetime.fromtimestamp(ftime).strftime('%H:%M:%S %d-%b-%Y'))
        AR.read_pfile(ivdatafile)
        # for i in range(AR.traces.shape[0]):
        #     mpl.plot(AR.time_base, AR.traces[i])
        # mpl.show()
        bridge_offset = 0.0
        threshold = -32.
        tgap = 0.  # gap before fittoign taum
            # print(self.AR.tstart, self.AR.tend)
        RM.setup(AR, SP, bridge_offset=bridge_offset)
        SP.setup(clamps=AR, threshold=threshold, 
                refractory=0.0001, peakwidth=0.001, interpolate=True, verify=False, mode='peak')
        SP.analyzeSpikes()
        SP.analyzeSpikeShape()
        # SP.analyzeSpikes_brief(mode='baseline')
        # SP.analyzeSpikes_brief(mode='poststimulus')
        SP.fitOne(function='fitOneOriginal')
        RM.analyze(rmpregion=[0., AR.tstart-0.001],
                    tauregion=[AR.tstart, AR.tstart + (AR.tend-AR.tstart)/5.],
                    to_peak=True, tgap=tgap)

        RMA = RM.analysis_summary
        print(RMA)
        fh = open(ivdatafile, 'rb')
        df = pickle.load(fh)
        r = df['Results'][0]

        for trial in range(len(df['Results'])):
            ds = df['Results'][trial]
            k0 = list(df['Results'][trial].keys())[0]
            dx = ds[k0]['monitor']
            P.axdict[gbc].plot(dx['time'], dx['postsynapticV'], linewidth=1.0)
            P.axdict[gbc].set_xlim(0., 150.)
            P.axdict[gbc].set_ylim(-200., 50.)
        PH.calbar(P.axdict[gbc], calbar=[120., -95., 25., 20.], axesoff=True, orient='left', 
                unitNames={'x': 'ms', 'y': 'mV'}, font='Arial', fontsize=8)
        ftname = str(ivdatafile.name)
        ip = ftname.find('_II_')+4
        ftname = ftname[:ip]+'...\n'+ftname[ip:]
        toptitle = f"{ftname:s}" # "\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']*1e3:.2f} ms"
        P.axdict[gbc].set_title(toptitle, fontsize=5)

    if chan == '_':
        chan = 'nav11'
    else:
        chan = chan[:-1]
    P.figure_handle.suptitle(f"Model: {modelName:s}  Na Ch: {chan:s} Scaling: {stitle:s} ", fontsize=7)
    mpl.show()

if __name__ == '__main__':
    main()
    
