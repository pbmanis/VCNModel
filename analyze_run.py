from __future__ import print_function
__author__ = 'pbmanis'

import os
import pickle
import numpy as np
import pylibrary.Utility as pu
import lmfit
import matplotlib.pylab as PL

verbose = False

class AnalyzeRun():
    def __init__(self, results):
        self.injs = results.keys()
        self.nRun = len(results.keys())
        self.parseResults(results)
        self.parseStim(results)
        self.IVResult={}

    def IV(self):
        """
        Compute the IV of the loaded dataset
        """
        # print 'V shape in IV: ', self.V.shape
        # print 'min v1: ', np.min(self.V[0,:])
        # print 'min v4: ', np.min(self.V[3,:])
        # return
        self.analyzeIV(self.t, self.V, self.I, self.tw, self.thr)

    def parseStim(self, res):
        """
        parse the stimulus information in the results dictionary.
        We only need to look at the first element to get the delay and duration
        """
        try:
            site = res[self.injs[0]].stim
        except:
            site = res[self.injs[0]]['stim']
        self.delay = site['delay']
        self.duration = site['dur']
        self.tw = [self.delay, self.duration, 10.]
        print( 'analyze_run.py: self.tw = ', self.tw)

    def parseResults(self, res):
        """
        Parse the results created by generate_run, building 2d arrays for V and I,
        and a vector for t
        """
        self.somasite=['postsynapticV', 'postsynapticI']
        try:
            msite = res[self.injs[0]].monitor
        except:
            msite = res[self.injs[0]]['monitor']
        vlen = len(msite[self.somasite[1]])
        self.V = np.zeros((self.nRun, vlen))
        self.I = np.zeros((self.nRun, vlen))
        for j,i in enumerate(res.keys()):  # each result key is a current level...
            try:
                msite = res[self.injs[j]].monitor
            except:
                msite = res[self.injs[j]]['monitor']
            if j == 0:
                self.t = msite['time'][0:vlen]
            self.V[j,:] = msite[self.somasite[0]]
            self.I[j,:] = msite[self.somasite[1]]
        self.thr = -30.0  # mV

    def clean_spiketimes(self, spikeTimes, mindT=0.7):
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
            spikeTimes = st
        return spikeTimes

    def analyzeIV(self, t, V, I, tw, thr):
        """ analyze a set of voltage records (IV), with spike threshold
            
            tw is a list of [tdelay, tdur, tssw], where tdelay is the delay to
            the start of the step, tdur is the duration of the step, and tssw is
            the duration of the steady-state window prior to the end of the
            step
            thr is the threshold that will be used for spike detection.
            Returns:
            a dictionary with:
            vmin
            vss
            i for vmin and vss
            spike count
            ispk
            eventually should also include time constant measures,and adaptation ratio
        """
        if verbose:
            print('starting analyzeIV')
        self.thr = thr
        ntraces = np.shape(V)[0]
        vss = np.empty(ntraces)
        vmin = np.zeros(ntraces)
        vrmss = np.zeros(ntraces)
        vm = np.zeros(ntraces)
        ic = np.zeros(ntraces)
        nspikes = np.zeros(ntraces)
        ispikes = np.zeros(ntraces)
        t_minv = np.zeros(ntraces)
        fsl = []
        fisi = []
        spk = {}
        taus = {}
        tauih = {}
        xtfit = {}
        ytfit = {}
        xihfit = {}
        yihfit = {}
        dt = t[1]-t[0]
        tss = [int((tw[1]-tw[2])/dt), int(tw[1]/dt)]
        ts = tw[0]
        te = tw[1]
        td = tw[2]
        for j in range(0, ntraces):
            if verbose:
                print('    analyzing trace: %d' % (j))
            vss[j] = np.mean(V[j,tss[0]:tss[1]])
            ic[j] = np.mean(I[j,tss[0]:tss[1]])
            vm[j] = np.mean(V[j, 0:int((ts-1.0)/dt)])
            mv = np.argmin(V[j,int(ts/dt):int(te/dt)])
            t_minv[j] = t[mv]  # time of minimum
            vmin[j] = V[j, mv]  # value of minimum
            spk[j] = pu.findspikes(t, V[j,:], self.thr, t0=ts, t1=te, dt=1.0, mode='peak')
            spk[j] = self.clean_spiketimes(spk[j])
            nspikes[j] = spk[j].shape[0] # build spike count list
            ispikes[j] = ic[j]  # currents at which spikes were detected
            if nspikes[j] >= 1:  # get FSL
                fsl.append(spk[j][0])
            else:
                fsl.append(None)
            if nspikes[j] >= 2:  # get first interspike interval
                fisi.append(spk[1]-spk[j][0])
            else:
                fisi.append(None)

            # fit the hyperpolarizing responses for Rin and "sag"
            # print 'ssi: ', ssi[0]
            if ic[j] < 0.0 and (t_minv[j]-ts) > 5.*dt: # just for hyperpolarizing pulses...
                if verbose:
                    print('    fitting trace %d' % j)

                    print('t.shape: ', t.shape)
                    print('V.shape: ', V[j,:].shape)
                    print('ts, vmin: ', ts, vmin)
                    print('ic[j]: ', ic[j])

                taus[j], xtfit[j], ytfit[j] = self.single_taufit(t, V[j,:], ts, t_minv[j])
                if verbose:
                    print('     calling fit')
                if (te-t_minv[j]) > 10.*dt:
                    tauih[j], xihfit[j], yihfit[j] = self.single_taufit(t, V[j,:], t_minv[j], te) # fit the end of the trace
                if verbose:
                    print('     completed fit')
            if verbose:
                print('   >>> completed analyzing trace %d' % j)
        if verbose:
            print('done with traces')
        
        RinIVss = (vss - vm)/ic  # measure steady-state input resistance
        RinIVpk = (vmin - vm)/ic  # measure "peak" input resistance
        icn, = np.where(ic < -1e-3)
        # print('ss rin: ', RinIVss[icn])
        # print('ic rin: ', ic[icn])
        if len(RinIVss[icn]) > 0:
            Rinss = np.max(RinIVss[icn])  # extract max
        else:
            Rinss = np.nan
        if len(RinIVpk[icn]) > 0:
            Rinpk = np.max(RinIVpk[icn])  # extract max
        else:
            Rinpk = np.nan
        if verbose:
            print('building IVResult')
        self.IVResult = {'I': ic, 'Vmin': vmin, 'Vss': vss, 'Vrmss': vrmss,
                'Vm': vm, 'Tmin': np.array(t_minv),
                'Ispike': np.array(ispikes), 'Nspike': np.array(nspikes), 'Tspike': np.array(spk),
                'FSL': np.array(fsl), 'FISI': np.array(fisi), 'taus': taus, 'tauih': tauih,
                'Rinss': Rinss, 'Rinpk': Rinpk,
                'taufit': [xtfit, ytfit], 'ihfit': [xihfit, yihfit],
                }

    def saveIVResult(self, name = None):
        """
        Save the result of multiple runs to disk. Results is in a dictionary,
        each element of which is a Param structure, which we then turn into
         a dictionary...
         The file name is a tuple, where the first element is a path and the
         second is a string name that will be decorated with the data type (by appending)
        """
        if name is None:
            return
        fn = os.path.join(name[0], name[1] + '_IVresult.p')
        pfout = open(fn, 'wb')
        pickle.dump({'IVResult': self.IVResult}, pfout)
        pfout.close()

    def expfit(self, p, x, y=None, C=None, sumsq=False, weights=None):
        """
        single exponential time constant with LM algorithm (using lmfit.py)
        'DC', 'a1', 'v1', 'k1', 'a2', 'v2', 'k2'
        """
        yd = p['dc'].value + (p['a'].value * np.exp(-x/p['tau'].value))
        if y is None:
            return yd
        else:
            if sumsq is True:
                return np.sqrt(np.sum((y - yd) ** 2))
            else:
                return y - yd

    def single_taufit(self, x, y, t0, t1, dc = True, fixedDC = -60., verbose=False):
        plot = False  # use to check if this is working.
        p = lmfit.Parameters()
        if dc is True:
            p.add('dc', value = -60., min = -10., max = -150., vary=True)
        else:
            p.add('dc', value = fixedDC, vary=False)

        p.add('a', value = -10., min= -100., max = 100.)
        p.add('tau', value = 25., min = 0.2, max = 100.)

        (cx, cy) = pu.clipdata(y, x, t0, t1, minFlag = False)
        cx -= t0   # fitting is reference to zero time
        try:
            mi = lmfit.minimize(self.expfit, p, args=(cx, cy))
            yfit = self.expfit(mi.params, cx)
            cx += t0  # restore absolute time

            # print 'sorted all data points'
            # for i, v in enumerate(xst):
            #     print '%7.2f\t%7.2f' % (v, yst[i])
            # print '---------------------------'

            if plot:
                PL.figure(314)
                PL.plot(x, y, 'k-')
                PL.plot(cx, yfit, 'r--')
                PL.show()

            fitpars = mi.params

            return fitpars, cx, yfit
        except:
            return None, None, None


