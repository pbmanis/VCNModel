from pathlib import Path
import pickle
import string
import numpy as np
import re
import toml
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')

import pylibrary.plotting.plothelpers as PH

import vcnmodel.cell_config as CellConf

config = toml.load(open('wheres_my_data.toml', 'r'))
dendqual = Path(config['baseDataDirectory'], config['dendriteQualityFile'])

with open(dendqual, 'rb') as fh:
    ASA = pd.read_excel(fh, 'Sheet1')

CC = CellConf.CellConfig() # get configuration methods
    
def main():
    files =["VCN_c09","VCN_c11", "VCN_c17", "VCN_c18"]
    substinput ={'VCN_c17a': 'VCN_c02', 'VCN_c17b': 'VCN_c05', 'VCN_c17c': 'VCN_c10', 'VCN_c17d': 'VCN_c13'}
    reps = 50
    db = 30
    freq = 16000.
    ch='nacncoop'
    srate = 'MS'
    normalPSTHdatasets = {}
    for cellID in files:  # programatically build file names from the configuration for each cell
          somaInflate = CC.get_soma_ratio(cellID)
          dendInflate = CC.get_dendrite_ratio(cellID)
          print(f"Soma, dend ratios: {cellID:s} : {somaInflate:.4f} {dendInflate:4f}")
          normalPSTHdatasets[cellID] = f"VCN_Cells/{cellID:s}/Simulations/AN/AN_Result_{cellID:s}_inp=self_XM13_{ch:s}_II"
          normalPSTHdatasets[cellID] += f"_soma={somaInflate:.3f}_dend={somaInflate:.3f}_ASA={somaInflate:.3f}"
          normalPSTHdatasets[cellID] += f"_delays_multisite_{reps:03d}_tonepip_{db:03d}dB_{freq:.1f}_MS.p"
    print(normalPSTHdatasets)
    # exit()

    """
    Substituions for the input pattern. Using 17 as the target
    """
    # subsPSTHdatasets = {
    # 'VCN_c17a':
    # 'VCN_Cells/VCN_c17/Simulations/AN/AN_Result_VCN_c17_inp=VCN_c02_XM13_{ch:s}_II_soma=1.344_dend=1.344_ASA=1.344_delays_multisite_050_tonepip_040dB_4000.0_MS.p',
    # 'VCN_c17b':
    # 'VCN_Cells/VCN_c17/Simulations/AN/AN_Result_VCN_c17_inp=VCN_c05_XM13_{ch:s}_II_soma=1.344_dend=1.344_ASA=1.344_delays_multisite_050_tonepip_040dB_4000.0_MS.p',
    # 'VCN_c17c':
    # 'VCN_Cells/VCN_c17/Simulations/AN/AN_Result_VCN_c17_inp=VCN_c10_XM13_{ch:s}_II_soma=1.344_dend=1.344_ASA=1.344_delays_multisite_050_tonepip_040dB_4000.0_MS.p',
    # 'VCN_c17d':
    # 'VCN_Cells/VCN_c17/Simulations/AN/AN_Result_VCN_c17_inp=VCN_c13_XM13_{ch:s}_II_soma=1.344_dend=1.344_ASA=1.344_delays_multisite_050_tonepip_040dB_4000.0_MS.p',
    # }
    db = 40
    freq = 4000.0
    subsPSTHdatasets = {}
    for cellIDn in list(substinput.keys()):  # programatically build file names from the configuration for each cell
          print('cellidn: ', cellIDn)
          cellID = cellIDn[:-1] # trim lettter off 
          print('cellid: ', cellID)
          omaInflate = CC.get_soma_ratio(cellID)
          dendInflate = CC.get_dendrite_ratio(cellID)
          # print(f"Soma, dend ratios: {cellID:s} : {somaInflate:.4f} {dendInflate:4f}")
          subsPSTHdatasets[cellID] = f"VCN_Cells/{cellID:7s}/Simulations/AN/AN_Result_{cellID:7s}_inp={substinput[cellIDn]:s}_XM13_{ch:s}_II"
          subsPSTHdatasets[cellID] += f"_soma={somaInflate:.3f}_dend={somaInflate:.3f}_ASA={somaInflate:.3f}"
          subsPSTHdatasets[cellID] += f"_delays_multisite_{reps:03d}_tonepip_{db:03d}dB_{freq:.1f}_MS.p"
    print(subsPSTHdatasets)

    # latest = {}
    # for j, f in enumerate(files):
    #     simdir = Path('VCN_Cells', f, 'Simulations/AN')
    #     # print(simdir)
    #     fs = simdir.glob(f"AN_Result_{f:s}*")
    #     mt = []
    #     fn = []
    #     for fsi in fs:
    #         fn.append(fsi)
    #         mt.append(fsi.stat().st_mtime)
    #
    #     i = np.argmax(mt)
    #     # print('Latest: ', fn[i], mt[i])
    #     latest[files[j]] = fn[i]
    #
    latest = normalPSTHdatasets
    # latest = subsPSTHdatasets
    print(latest)

    re_inp = re.compile('inp=(VCN_c[0-9]{2})')

    # test re
    # m = re.search(re_inp, str(subsPSTHdatasets['VCN_c17a']))
    # if m is not None:
    #     print('m: ', m.group(0))
    # exit()

    # make a plot window

    nrc = PH.getLayoutDimensions(len(latest)*3, pref='height')
    plabels = []
    for i in range(nrc[1]):
        for j in range(nrc[0]):
            plabels.append(f"{string.ascii_uppercase[i]:s}{j+1:d}")
    # print(plabels)
    P = PH.regular_grid(nrc[0], nrc[1], order='columnsfirst', figsize=(10., 6.), showgrid=False,
                    verticalspacing=0.1, horizontalspacing=0.1,
                    margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.08, 'bottommargin': 0.1},
                    labelposition=(-0.15, 1.05), parent_figure=None, panel_labels=plabels)

    n = 0
    names = list(latest.keys())
    for ifx, fx in enumerate(latest):
        filename = latest[fx]
        filename = Path('/Users/pbmanis/Desktop/Python/VCN-SBEM-Data', filename)# print(filename)
        if fx is None:
            continue
        if not Path(filename).is_file():
            print(f'File: {str(filename):s} was not found')
            continue
        with open(filename, 'rb') as fh:
            print('opening specific file: ', filename)
            d = pickle.load(fh)
            # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
        si = d['trials'][0]['stimInfo']
        totaldur = si['pip_start'] + np.max(si['pip_start']) + si['pip_duration'] + si['pip_offduration']
        # print('dur: ', totaldur)
        allt = []
        inp = re.search(re_inp, str(latest[names[ifx]]))
        if inp is not None:
            inp = inp.groups(0)[0]
        else:
            inp = 'input = self'
        P.axdict[plabels[n+0]].set_title(f"{str(filename)[10:17]:s}  input: {inp:s}", fontsize=9)
        for i in range(len(d['trials'])):  # for all trails in the measure.
            trial = d['trials'][i]
            w = trial['stimWaveform']
            stb = trial['stimTimebase']
            edur = np.max(trial['time']/1000.)
            if i < 5:
                P.axdict[plabels[n+0]].plot(trial['time']/1000., trial['somaVoltage'], linewidth=0.5)
            # P.axdict['B'].plot(stb, w, linewidth=0.5)  # stimulus underneath
            P.axdict[plabels[n+1]].plot(trial['spikeTimes']/1000., i*np.ones(len( trial['spikeTimes'])),
                    'o', markersize=1.0, color='k')
            allt.append(trial['spikeTimes'][0]/1000.)
            # inputs = len(trial['inputSpikeTimes'])
            # for k in range(inputs):
            #     tk = trial['inputSpikeTimes'][k]/1000.
            #     y = (i+0.1+k*0.05)*np.ones(len(tk))
            #     P.axdict[plabels[n+2]].plot(tk, y, '|', markersize=2.5, color='r', linewidth=0.5)
        P.axdict[plabels[n+0]].set_xlim(0, edur)
        P.axdict[plabels[n+1]].set_xlim(0, edur)
        P.axdict[plabels[n+2]].set_xlim(0, edur)
        data = d
        allst = []
        for trial in data['trials']:
            # print (trial)
            allst.extend(data['trials'][trial]['spikeTimes']/1000.)
        allst = np.array(allst)
        allst = np.sort(allst)
            # combine all spikes into one array, plot PSTH
            # the histogram of the data
        print(edur)
        P.axdict[plabels[n+2]].hist(allst, int(1000*edur/0.5), range=(0., edur), density=1, facecolor='black', alpha=0.75)
        n += 3


    PH.mpl.show()

        # def show_AN(self, filename=None):
        #     print('Show an')
        #     # bulid plot layout
        #     d = {}
        #     # for p in self.Params['patterns']:  # for each cell/run get the data
        #     #     h = open(self.fn[p], 'rb')
        #     #     d[p] = pickle.load(h)
        #     #     h.close()
        #     if filename is not None:
        #         with open(filename, 'rb') as fh:
        #             print('opening specific file: ', filename)
        #             d['c09'] = pickle.load(fh)
        #         self.Params['patterns'] = ['c09']
        #
        #
        #     # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
        #     sk = sizer.keys()
        #     si = d[self.Params['patterns'][0]]['Params']
        #     totaldur = si['pip_start'] + np.max(si['pip_start']) + si['pip_duration'] + si['pip_offduration']
        #     for j, pattern in enumerate(self.Params['patterns']): # for all cells in the "pattern"
        #
        #         allt = []
        #         for i in range(len(d[pattern]['trials'])):  # for all trails in the measure.
        #             trial = d[pattern]['trials'][i]
        #             w = trial['stimWaveform']
        #             stb = trial['stimTimebase']
        #             # print(trial['somaVoltage'])
        #             P.axdict['A'].plot(trial['time']/1000., trial['somaVoltage'], linewidth=0.5)
        #             print(stb)
        #             print(w)
        #             # P.axdict['B'].plot(stb, w, linewidth=0.5)  # stimulus underneath
        #             P.axdict['C'].plot(trial['spikeTimes']/1000., i*np.ones(len( trial['spikeTimes'])),
        #                     'o', markersize=2.5, color='b')
        #             allt.append(trial['spikeTimes'][0]/1000.)
        #             inputs = len(trial['inputSpikeTimes'])
        #             for k in range(inputs):
        #                 tk = trial['inputSpikeTimes'][k]/1000.
        #                 y = (i+0.1+k*0.05)*np.ones(len(tk))
        #                 P.axdict['E'].plot(tk, y, '|', markersize=2.5, color='r', linewidth=0.5)
        #         # allt = np.array(allt).ravel()
        #         # # print(allt)
        #         # P.axdict['F'].hist(allt, range=(0., totaldur))
        #         # print('plotted hist)')
        #
        #
        #     for a in ['A', 'B', 'C', 'E', 'F']:  # set some common layout scaling
        #         P.axdict[a].set_xlim((0., totaldur))
        #         P.axdict[a].set_xlabel('T (s)')
        #     P.axdict['A'].set_ylabel('mV', fontsize=8)
        #     P.axdict['D'].set_xlabel('Phase', fontsize=8)
        #     P.axdict['C'].set_ylabel('Trial', fontsize=8)
        #     P.axdict['E'].set_ylabel('Trial, ANF', fontsize=8)
        #     P.axdict['B'].set_title('Stimulus', fontsize=9)
        #     P.axdict['E'].set_title('ANF Spike Raster', fontsize=9)
        #     P.axdict['C'].set_title('Bushy Spike Raster', fontsize=9)
        #     P.axdict['F'].set_title('PSTH', fontsize=9)
        #     # combine all spikes into one array, plot PSTH
        #     data = d[pattern]
        #     allst = []
        #     for trial in data['trials']:
        #         # print (trial)
        #         allst.extend(data['trials'][trial]['spikeTimes']/1000.)
        #     allst = np.array(allst)
        #     allst = np.sort(allst)
        #         # combine all spikes into one array, plot PSTH
        #         # the histogram of the data
        #     P.axdict['F'].hist(allst, 100, normed=1, facecolor='blue', alpha=0.75)
        #
        #     if si['soundtype'] == 'SAM':  # calculate vs and plot histogram
        #         # print('allst: ', allst)
        #         phasewin = [data['Params']['pip_start'][0] + 0.25*data['Params']['pip_duration'],
        #             data['Params']['pip_start'][0] + data['Params']['pip_duration']]
        #         # print (phasewin)
        #         spkin = allst[np.where(allst > phasewin[0])]
        #         spikesinwin = spkin[np.where(spkin <= phasewin[1])]
        #
        #         # set freq for VS calculation
        #         if si['soundtype'] == 'tone':
        #             f0 = data['Params']['F0']
        #             print ("Tone: f0=%.3f at %3.1f dbSPL, cell CF=%.3f" %
        #                     (data['Params']['f0'], data['Params']['dB'], data['Params']['F0'] ))
        #         if si['soundtype'] == 'SAM':
        #             f0 = data['Params']['fmod']
        #             tstring = ("SAM Tone: f0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f" %
        #                  (data['Params']['F0'], data['Params']['dB'], data['Params']['fmod'], data['Params']['dmod'], data['Params']['F0']))
        #             print(tstring)
        #             # P.figure_handle.suptitle(tstring + '\n' + self.fn[self.Params['patterns'][0]].replace('_', '\_'), fontsize=9)
        #        # print('spikes: ', spikesinwin)
        #         print((' Sound type: ', si['soundtype']))
        #         if len(spikesinwin) < 10:
        #             vs = 0.
        #             print('AN Vector Strength: Insufficient spike count, n = %d' % (len(spikesinwin)))
        #         else:
        #             vs = vector_strength(spikesinwin*1000., f0)  # vs expects spikes in msec
        #             print('AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (f0, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n']))
        #             vspfile.add_data(self.Params['patterns'][0], self.Params['spirou'], vs['r'])
        #             # fh = open(self.vspfile, 'rb')  # set file handle to write  - open, and append
        #             # vsp = pickle.load(fh)  # get vs data for all cells (outside this file)
        #             # fh.close
        #             # vsp[self.Params['patterns'][0]][self.Params['spirou']] = vs['r']  # update entry for this cell
        #             # fh = open(self.vspfile, 'wb')  # set file handle to write  - open, which will make an empty file
        #             # pickle.dump(vsp, fh)  # save it.
        #             # fh.close()
        #         # print(vs['ph'])
        #             P.axdict['D'].hist(vs['ph'], bins=2*np.pi*np.arange(30)/30.)
        #             P.axdict['D'].set_xlim((0., 2*np.pi))
        #             P.axdict['D'].set_title('Phase (VS = {0:.3f})'.format(vs['r']), fontsize=9, horizontalalignment='center')
        #
        #     # make figure output filename
        #     if self.Params['inflationfactor'] != 1.0:
        #         t = 'Inflated'
        #     else:
        #         t = 'Original'
        #     figname = str(self.fn) # f"{t:s}_{self.Params['patterns'][0]:s}_{self.Params['spirou']:s}.pdf"
        #     # P.figure_handle.suptitle(figname.replace('_', '\_'))
        #     mpl.savefig('figname_1.pdf')  # rasterized to 300 dpi is ok for documentation.
        #     mpl.show()
        
if __name__ == '__main__':
    main()



