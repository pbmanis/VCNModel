__author__ = 'pbmanis'

import numpy as np
import pylibrary.Utility as pu
import lmfit
import matplotlib.pylab as PL

class AnalyzeRun():
    def __init__(self, results):
        self.injs = results.keys()
        self.nRun = len(results.keys())
        self.parseResults(results)
        self.parseStim(results)


    def IV(self):
        """
        Compute the IV of the loaded dataset
        """
        self.IVresult = self.analyzeIV(self.t, self.V, self.I, self.tw, self.thr)

        # now get input resistance

        print 'Current: ', self.IVresult['I']
        print '# spikes: ', self.IVresult['Nspike']
        return self.IVresult


    def parseStim(self, res):
        """
        parse the stimulus informaiton in the results dictionary.
        We only need to look at the first element to get the delay and duration
        """
        self.delay = res[self.injs[0]].stim['delay']
        self.duration = res[res.keys()[0]].stim['dur']
        self.tw = [self.delay, self.duration, 10.]


    def parseResults(self, res):
        """
        Parse the results created by generate_run, building 2d arrays for V and I,
        and a vector for t
        """
        self.somasite=['postsynapticV', 'postsynapticI']
        vlen = len(res[self.injs[0]].monitor[self.somasite[1]])
        self.V = np.zeros((self.nRun, vlen))
        self.I = np.zeros((self.nRun, vlen))
        for j,i in enumerate(res.keys()):  # each result key is a current level...
            if j == 0:
                self.t = res[self.injs[j]].monitor['time'][0:vlen]
            self.V[j,:] = res[self.injs[j]].monitor[self.somasite[0]]
            self.I[j,:] = res[self.injs[j]].monitor[self.somasite[1]]
        self.thr = -10.0  # mV


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
        ntraces = np.shape(V)[0]
        vss     = []
        vmin    = []
        vrmss    = []
        vm      = []
        ic       = []
        nspikes = []
        ispikes = []
        tmin = []
        fsl = []
        fisi = []
        spk={}
        taus={}
        tauih = {}
        xtfit = {}
        ytfit = {}
        xihfit = {}
        yihfit = {}
        for j in range(0, ntraces):
            ts = tw[0]
            te = tw[1]
            td = tw[2]
            ssv  = pu.measure('mean', t, V[j,:], te-td, te)
            ssi  = pu.measure('mean', t, I[j,:], te-td, te)
            rvm  = pu.measure('mean', t, V[j,:], 0.0, ts-1.0)
            minv = pu.measure('min', t, V[j,:], ts, te)
            spk[j] = pu.findspikes(t, V[j,:], thr, t0=ts, t1=te, dt=1.0, mode='peak')
            nspikes.append(spk[j].shape[0]) # build spike count list
            ispikes.append(ssi[0])  # currents at which spikes were detected
            if nspikes[-1] >= 1:  # get FSL
                fsl.append(spk[j][0])
            else:
                fsl.append(None)
            if nspikes[-1] >= 2:  # get first interspike interval
                fisi.append(spk[1]-spk[j][0])
            else:
                fisi.append(None)
            vm.append(rvm[0])  # rmp

            # fit the hyperpolarizing responses for Rin and "sag"
            if ssi[0] < 0.0: # just for hyperpolarizing pulses...
                taus[j], xtfit[j], ytfit[j] = self.single_taufit(t, V[j,:], ts, minv[1])
                vrmss.append(rvm[0])
                ic.append(ssi[0])
                vss.append(ssv[0]) # get steady state voltage
                vmin.append(minv[0]) # and min voltage
                tmin.append(minv[1]) # and min time
                tauih[j], xihfit[j], yihfit[j] = self.single_taufit(t, V[j,:], minv[1], te) # fit the end of the trace
                print tauih[j]

        vss = np.array(vss)  # steady state during current pulse
        vrmss = np.array(vrmss)  # resting potential for hyperpolarizaing pulses only
        ic = np.array(ic)  # injected current
        vmin = np.array(vmin)  # min voltage, hyperpolarizing only
        tmin = np.array(tmin)  # time of min with hyp (used for Ih fitting)
        RinIVss = (vss - vrmss)/ic  # measure steady-state input resistance
        RinIVpk = (vmin - vrmss)/ic  # measure "peak" input resistance
        Rinss = np.max(RinIVss)  # extract max
        Rinpk = np.max(RinIVpk)  # extract max

        return({'I': ic, 'Vmin': vmin, 'Vss': vss, 'Vrmss': vrmss,
                'Vm': vm, 'Tmin': np.array(tmin),
                'Ispike': np.array(ispikes), 'Nspike': np.array(nspikes), 'Tspike': np.array(spk),
                'FSL': np.array(fsl), 'FISI': np.array(fisi), 'taus': taus, 'tauih': tauih,
                'Rinss': Rinss, 'Rinpk': Rinpk,
                'taufit': [xtfit, ytfit], 'ihfit': [xihfit, yihfit],
                })


    def expfit(self, p, x, y=None, C=None, sumsq=False, weights=None):
        """
        single exponential time constant with LM algorithm (using lmfit.py)
        'DC', 'a1', 'v1', 'k1', 'a2', 'v2', 'k2'
        """
        yd = p['dc'].value + (p['a'].value * np.exp(-x/p['tau'].value))
        if y == None:
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
            p.add('dc', value = -60., min = 0., max = -200., vary=True)
        else:
            p.add('dc', value = fixedDC, vary=False)

        p.add('a', value = -10., min= -100., max = 100.)
        p.add('tau', value = 5., min = 0.5, max = 100.)

        (cx, cy) = pu.clipdata(y, x, t0, t1, minFlag = False)
        cx -= t0   # fitting is reference to zero time
        mi = lmfit.minimize(self.expfit, p, args=(cx, cy))
        mi.leastsq()
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



