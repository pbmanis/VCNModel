from __future__ import print_function
__author__ = 'pbmanis'
"""
Perform analysis of IVs from model_run
Provides basic analysis

"""
import os
import pickle
import numpy as np
import pylibrary.Utility as pu
from lmfit import Model
from lmfit.models import ExponentialModel

verbose = False

class AnalyzeRun():
    def __init__(self, results):
        self.results = results
        self.injs = list(results.keys())
        self.nRun = len(self.injs)
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
        if verbose:
            print( 'analyze_run.py: self.tw = ', self.tw)

    def parseResults(self, res):
        """
        Parse the results created by generate_run, building 2d arrays for V and I,
        and a vector for t
        
        Parameters
        ----------
        res : dict
            The dictionary of the results, as generated by model_run
        
        Returns
        -------
        Nothing
        """

        self.somasite=['postsynapticV', 'postsynapticI']
        inj0 = self.injs[0]
        try:
            msite = self.results[inj0].monitor
        except:
            msite = self.results[inj0]['monitor']
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
        Parameters
        ----------
        spikeTimes : list or numpy array (1-D)
            array of the spike times
        
        mindT : float (default : 0.7)
            minimum time between spikes, in the same units as spikeTimes
            (normally this will be in milliseconds)
        
        Return
        ------
        spikeTimes : list or numpy array (1-D_
            A cleaned list of the spike times where the events are at least
            mindT appart.
            Note: If there is only 1 or 0 spikes in array, just return the array
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
        """
        Analyze a set of voltage records (IV), with spike threshold
            
        Parameters
        ----------
        t : numpy array of floats (1-D)
            time array for the voltage and current
        V : numpy array of floats (2-D)
            Voltage traces to be analyzed. Dimension 0 is trace number,
            dimension 1 corresponds to time evolution of voltage
        I : numpy array of floats (2-D)
            Current traces, corresponding to the voltage traces in V. Should
            be organized in the same way.
        tw : list
            list of [tdelay, tdur, tssw], where:
                tdelay is the delay to the start of the step.
                tdur is the duration of the step
                tssw is the duration of the steady-state window prior
                    to the end of the step
        thr : float
            Voltage threshold that will be used for spike detection.
        
        Returns
        -------
            a dictionary with:
            vmin
            vss
            i for vmin and vss
            spike count
            ispk
            (eventually should also include time constant measures,and adaptation ratio)
        """
        
        if verbose:
            print('starting analyzeIV')
        self.thr = thr
        ntraces = np.shape(V)[0]
        # initialize all result arrays, lists and dicts
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
        # break down the time windows
        tss = [int((tw[1]-tw[2])/dt), int(tw[1]/dt)]
        ts = tw[0]
        te = tw[1]
        td = tw[2]
        for j in range(0, ntraces):
            if verbose:
                print('    analyzing trace: %d' % (j))
            vss[j] = np.mean(V[j,tss[0]:tss[1]])  # steady-state voltage
            ic[j] = np.mean(I[j,tss[0]:tss[1]])  # corresponding currents
            vm[j] = np.mean(V[j, 0:int((ts-1.0)/dt)])  # resting potential - for 1 msec prior to step
            minV = np.argmin(V[j, int(ts/dt):int((ts+te)/dt)])
            t_minv[j] = t[minV+int(ts/dt)]  # time of minimum
            vmin[j] = V[j, minV+int(ts/dt)]  # value of minimum
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

            # fit the hyperpolarizing responses for Tau_m and "sag" due to Ih activation
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
                    tauih[j], xihfit[j], yihfit[j] = self.single_taufit(t, V[j,:], t_minv[j], te+ts) # fit the end of the trace
                if verbose:
                    print('     completed fit')
            if verbose:
                print('   >>> completed analyzing trace %d' % j)
        if verbose:
            print('done with traces')
        
        RinIVss = (vss - vm)/ic  # measure steady-state input resistance
        RinIVpk = (vmin - vm)/ic  # measure "peak" input resistance
        icn, = np.where(ic < -1e-3)
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

    def saveIVResult(self, name):
        """
        Save the result of multiple runs to disk. Results is in a dictionary,
        each element of which is a Param structure, which we then turn into
        a dictionary.
        Parameters
        ----------
        name : tuple of str (default: None)
            The file name is a tuple, where the first element is a path and the
            second is a string name that will be decorated with the data type (by appending)
        """
        fn = os.path.join(name[0], name[1] + '_IVresult.p')
        pfout = open(fn, 'wb')
        pickle.dump({'IVResult': self.IVResult}, pfout)
        pfout.close()

    def single_taufit(self, x, y, t0, t1):
        """
        Perform single exponential fit to voltage traces
        
        Parameters
        ----------
        x : numpy array
            time corresponding to data in y
        y : numpy array
            voltage trace (single trace)
        t0 : float
            start time for fit (ms)
        t1 : float
            end time for fit (ms)
        
        Returns
        -------
        tuple of (best fit parameters as dict, x times, best fit y values)
        If fit fails, return is (None, None, None)
        """
        
        (cx, cy) = pu.clipdata(y, x, t0, t1, minFlag = False)
        expmodel = ExponentialModel() #, prefix='', missing=None, name=None, **kwargs)
        expmodel.set_param_hint('decay', min=0.1, max=50.0)
        cye = np.mean(cy[-5:])
        try:
            result = expmodel.fit(cy-cye, x=cx-t0, amplitude=0., decay=10.)
            if verbose:
                print(result.fit_report())  # print the result
            rbv = result.best_values
            fitr = {'a': 0, 'tau': 5., 'dc': 0.}
            fitr['a'] = rbv['amplitude']
            fitr['tau'] = rbv['decay']
            return fitr, cx, result.best_fit+cye
        except:
            return None, None, None

