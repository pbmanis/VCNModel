from pathlib import Path
import pickle
import string
import numpy as np
import pylibrary.PlotHelpers as PH

files =["VCN_c09","VCN_c11", "VCN_c17", "VCN_c18"]

normalPSTHdatasets = {
'VCN_c09':
PosixPath('VCN_Cells/VCN_c09/Simulations/AN/AN_Result_VCN_c09_inp=self_XM13_II_soma=1.514_dend=1.514_ASA=1.514_delays_multisite_050_tonepip_040dB_4000.0_MS.p'),
'VCN_c11':
PosixPath('VCN_Cells/VCN_c11/Simulations/AN/AN_Result_VCN_c11_inp=self_XM13_II_soma=1.458_dend=1.458_ASA=1.458_delays_multisite_050_tonepip_040dB_4000.0_MS.p'),
'VCN_c17':
PosixPath('VCN_Cells/VCN_c17/Simulations/AN/AN_Result_VCN_c17_inp=self_XM13_II_soma=1.344_dend=1.344_ASA=1.344_delays_multisite_050_tonepip_040dB_4000.0_MS.p'),
'VCN_c18':
PosixPath('VCN_Cells/VCN_c18/Simulations/AN/AN_Result_VCN_c18_inp=self_XM13_II_soma=1.502_dend=1.502_ASA=1.502_delays_multisite_050_tonepip_040dB_4000.0_MS.p')
}

latest = {}
for j, f in enumerate(files):
    simdir = Path('VCN_Cells', f, 'Simulations/AN')
    # print(simdir)
    fs = simdir.glob(f"AN_Result_{f:s}*")
    mt = []
    fn = []
    for fsi in fs:
        fn.append(fsi)
        mt.append(fsi.stat().st_mtime)
    
    i = np.argmax(mt)
    # print('Latest: ', fn[i], mt[i])
    latest[files[j]] = fn[i]

print(latest)
# make a plot window

nrc = PH.getLayoutDimensions(len(latest)*3, pref='height')
plabels = []
for i in range(nrc[0]):
    for j in range(nrc[1]):
        plabels.append(f"{string.ascii_uppercase[j]:s}{i+1:d}")
# print(plabels)
P = PH.regular_grid(nrc[0], nrc[1], order='columnsfirst', figsize=(10., 6.), showgrid=False,
                verticalspacing=0.1, horizontalspacing=0.1,
                margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.08, 'bottommargin': 0.1},
                labelposition=(-0.1, 1.05), parent_figure=None, panel_labels=plabels)

n = 0
for ifx, fx in enumerate(latest):
    filename = latest[fx]
    # print(filename)
    if fx is None:
        continue
    with open(filename, 'rb') as fh:
        print('opening specific file: ', filename)
        d = pickle.load(fh)
        # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
    si = d['trials'][0]['stimInfo']
    totaldur = si['pip_start'] + np.max(si['pip_start']) + si['pip_duration'] + si['pip_offduration']
    # print('dur: ', totaldur)
    allt = []
    P.axdict[plabels[n+0]].set_title(str(filename)[10:17], fontsize=9)
    for i in range(len(d['trials'])):  # for all trails in the measure.
        trial = d['trials'][i]
        w = trial['stimWaveform']
        stb = trial['stimTimebase']
        edur = np.max(trial['time']/1000.)
        if i < 5:
            P.axdict[plabels[n+0]].plot(trial['time']/1000., trial['somaVoltage'], linewidth=0.5)
        # P.axdict['B'].plot(stb, w, linewidth=0.5)  # stimulus underneath
        P.axdict[plabels[n+1]].plot(trial['spikeTimes']/1000., i*np.ones(len( trial['spikeTimes'])),
                'o', markersize=2, color='k')
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
    P.axdict[plabels[n+2]].hist(allst, 100, normed=1, facecolor='black', alpha=0.75)
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
        


