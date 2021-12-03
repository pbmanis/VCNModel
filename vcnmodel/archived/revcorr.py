def plot_revcorr(p, ax, ax2, respike=False, thr=-20., width=4.0):
    d={}
    basepath = 'VCN_Cells/VCN_{0:s}/Simulations/AN'.format(p)
    h = open(os.path.join(basepath, fn[p]), 'rb')
    d[p] = pickle.load(h, )
    h.close()
    seaborn.set_style('ticks')
    syninfo = SC.VCN_Inputs['VCN_{0:s}'.format(p)]
    if not respike:
        st = d[p]['spikeTimes']
    else:
        st = {}
        # print(d.keys())
        # print(d[p].keys())
        # print(d[p]['trials'][0].keys())
        for k in d[p]['trials'].keys():
            trd = d[p]['trials'][k]
            sv = trd['somaVoltage']
            ti = trd['time']
            dt = ti[1]-ti[0]
            st[k] = pu.findspikes(ti, sv, thr, dt=dt, mode='peak')
            st[k] = clean_spiketimes(st[k])
    ndet1 = 0
    for n in st:
        ndet1 = ndet1 + len(st[n])   
    ndet0 = 0
    for n in d[p]['trials'].keys():
        st = d[p]['trials'][n]['spikeTimes']
        ndet0 = ndet0 + len(st)
        
    print(f'Detected {ndet1:d} AN spikes')
    print(f'Detected {ndet0:d} Postsynaptic spikes')
    ntrials = len(d[p]['trials'].keys())
    print ("# trials: ", ntrials)
    ninputs = len(syninfo[1])
    sites = np.zeros(ninputs)
    print('ninputs: ', ninputs)
    binw = 0.1
    tcwidth = width # msec for total correlation width
    xwidth = 5.
    tx = np.arange(-xwidth, 0, binw)
    amax = 0.
    for isite in range(ninputs): # precompute areas
        area = syninfo[1][isite][0]
        if area > amax:
            amax = area
        sites[isite] = int(np.around(area*SC.synperum2))
        
    summarySiteTC = {}
    maxtc = 0
    for isite in range(ninputs):  # for each ANF input (get from those on first trial)
        for trial in range(ntrials):  # sum across trials
            stx = d[p]['trials'][trial]['spikeTimes']  # get postsynaptic spike train for this trial
            anx = d[p]['trials'][trial]['inputSpikeTimes'][isite]

            andirac = np.zeros(int(200./binw)+1)
            if trial == 0:
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
                # C = C + SPKS.spike_triggered_average(st[trial]*ms, andirac*ms, max_interval=width*ms, dt=binw*ms)
        nc = int(len(C)/2)
        TC = TC/len(st)
        summarySiteTC[isite] = TC
        color = plt.cm.viridis(norm(sites, isite))
        ax.set_facecolor((0.7, 0.7, 0.7))
        # print(tx)
        # print('nc: ', nc)
        # print(C[:nc])
        ax.plot(tx, C[:nc], color=color, label=('Input {0:2d} N={1:3d}'.format(isite, int(sites[isite]))),
            linewidth=0.75)
        tx2 = np.array([0.2, 0.8])
        ax2.plot(tx2, TC*np.ones_like(tx2), color=color, linewidth=2)
        if TC > maxtc:
            maxtc = TC

    print('finished inputs')
    seaborn.despine(ax=ax)
    ax.set_ylabel('Rate of coincidences/bin (Hz)', fontsize=12)
    ax.set_xlabel('T (ms)', fontsize=12)
    ax.set_xlim(-5., 1.)
    ax.set_ylim(0, 1)
    

    ax2.set_ylabel('Total Correlation W=%.1f-%0.1f'% (tcwidth[0], tcwidth[1]), fontsize=12)
    ax2.set_ylim(0, 1.0) # maxtc*1.2)
    PH.talbotTicks(ax2, axes='xy',
                   density=(1.0, 1.0), insideMargin=0.05, pointSize=10, 
                   tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1})
                   
    a = re_self.search(fn[p])
    b = re_c10.search(fn[p])
    if a is not None:
        inp = 'VCN_c09'
    elif b is not None:
        inp = 'VCN_c10'
    else:
        inp = "Undefined"
    ax.set_title(f'VCN_{p:s} input={inp:s} [{int(np.min(sites)):d}-{int(np.max(sites)):d}]\nAmax={amax:.1f}', y=0.9, x=0.02,
        horizontalalignment='left', fontsize=6)
    return summarySiteTC, sites

def plot_SAC():
    #fig, ax = plt.subplots(len(pattern_list)+1, 3)
    fig = plt.figure()

    gs_outer = gridspec.GridSpec(4, 1)  # 4 rows, one column
    gs_inner = []
    ax = {}
    for i, outer in enumerate(gs_outer):
        gs_inner.append(gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=outer))
        ax1 = plt.Subplot(fig, gs_inner[-1][:-1, 0])
        fig.add_subplot(ax1)
        ax2 = plt.Subplot(fig, gs_inner[-1][3, 0])
        fig.add_subplot(ax2)
    
        ax3 = plt.Subplot(fig, gs_inner[-1][:, 1])
        fig.add_subplot(ax3)
        ax4 = plt.Subplot(fig, gs_inner[-1][:, 2])
        fig.add_subplot(ax4)
        ax5 = plt.Subplot(fig, gs_inner[-1][:, 3])
        fig.add_subplot(ax5)
        ax[i] = [ax1, ax2, ax3, ax4, ax5]
        for a in ax[i]:
            PH.nice_plot(a)
        #PH.noaxes(ax1, 'x')

    #fig.tight_layout()
    # plt.show()
    # exit()


    # d[cell] has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']
    #print 'data keys: ', d.keys()
    #print len(d['time'])
    fmod = d[patterns[0]]['stimInfo']['mod']  # modulation frequency
    sacht = 2.5
    sacmax = 100
    phaseht = 15

    #print d['spikeTimes']
    for j, pattern in enumerate(patterns):
        w = d[patterns[j]]['stimWaveform'][0][0].generate()
        t = d[patterns[j]]['stimWaveform'][0][0].time
        ax[j][1].plot(t, w)  # stimulus underneath
        for i, st in enumerate(d[pattern]['spikeTimes'].keys()):
            ax[j][0].plot(d[pattern]['spikeTimes'][st]/1000., i*np.ones(len( d[pattern]['spikeTimes'][st])), 'o', markersize=2.5, color='b')
            inputs = len(d[pattern]['inputSpikeTimes'][st])
            # for j in range(inputs):
            #     t = (d['ref']['inputSpikeTimes'][st][j])/1000.
            #     y = (i+0.1+j*0.05)*np.ones(len(t))
            #     ax[0,].plot(t, y, 'o', markersize=2.5, color='r')
    #for i, st in enumerate(d['ref']['somaVoltage'].keys()):
    #    ax[2,].plot(d['ref']['time'], d['ref']['somaVoltage'][st], color='k')
        ax[j][0].set_xlim((0, 0.8))
        ax[j][1].set_xlim((0, 0.8))


    sac = SAC.SAC()
    xf = []
    dt = 0.1
    u = 2.0*np.pi/(1000.0/fmod)
    for j, pattern in enumerate(patterns):
        X = []
        R = {}
        pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                'dur': 1000., 'ddur': 200.}
        pars_x =     pars = {'twin': 500., 'binw': 1., 'ntestrep': 20, 
                'baseper': 1.0, 'stimdur': 800., 'delay': 100.,
                'dur': 1000., 'ddur': 200.}
        tsynch = np.arange(0., 1000., dt)
        x_spl = []
        y_spl = []
        Nspikes = 0
        for i, st in enumerate(d[pattern]['spikeTimes'].keys()):

            X.append(d[pattern]['spikeTimes'][st])
            xw = np.zeros(tsynch.shape[0])
    #        print [int(xx/dt) for xx in X[-1]]
            Nspikes += np.shape(X[-1])[0]
            # calculate vector strength as well.
            x_spl.append(np.cos(np.array(X[-1])*u))
            y_spl.append(np.sin(np.array(X[-1])*u))
        x_spl = np.hstack(x_spl)
        y_spl = np.hstack(y_spl)
        vs = np.sqrt(np.sum(x_spl)**2 + np.sum(y_spl)**2)/Nspikes
        th = np.arctan2(y_spl, x_spl)
        p, z = PCS.rayleigh(th, w=None, d=None, axis=None)
        if j == 0:
            X0 = X  # save the X
        yh, bins = sac.SAC_asm(X, pars)
        print ('mean vs: for %s at fMod = %.2f:  %f angle: %f' % (pattern, fmod, vs, th.mean()))
        print ('rayleigh: p=%f  z=%f (p > 0.05 means data is uniformly distributed)' % (p, z))
        ax[j][2].bar(bins[:-1], yh)
        ax[j][2].set_xlim((-sacmax, sacmax))
        ax[j][2].set_ylim((0, sacht))
        phasebinnum = 101
        phasebins = np.linspace(-0.5, 0.5, num=phasebinnum)
        thhist, thbins = np.histogram(th/np.pi, bins=phasebins, density=False)
        ax[j][3].bar(thbins[:-1]+0.5, thhist, width=1.0/phasebinnum)
        ax[j][3].set_ylim((0, phaseht))

        if j == 0:
            continue
        ycc, ccbins = sac.XAC(X0, X, pars_x, trialbased=True)
        ax[j][4].bar(ccbins[:-1], ycc)
        ax[j][4].set_xlim((-sacmax, sacmax))
    
        # now cross-correlate data...
    

    #PH.nice_plot(ax.flatten().tolist())
    plt.show()