"""
GIF fitter wrapper. 
Set up with ln -s giffitterlocation/src gif in the running directory

Basic run will do a test against the datasets provided at the original site.

Fit generalized integrate and fire models to data
Uses Pozzorini et al. PLoS Comp. Bio. 2015 method,
"""
import os
import pickle
import matplotlib.pyplot as mpl
import numpy as np
from gif.Experiment import Experiment
from gif.AEC_Badel import AEC_Badel
from gif.GIF import GIF
from gif.Filter_Rect_LogSpaced import Filter_Rect_LogSpaced
from gif.Filter_Exps import Filter_Exps
import NoiseTrainingGen as NTG
import faulthandler

faulthandler.enable()

class GIFFitter():
    def __init__(self, path, dt=0.1):
        """
        Set up fitter. 
        Parameters
        ----------
        path : str (no default)
            full path to the data to be fit
        dt : float (default: 0.1)
            Time step, in msec
        """
        self.path = path
        self.dt = dt
        self.set_eta_timescales()
        self.set_gamma_timescales()
        self.Exp = Experiment('Experiment 1', self.dt)
        self.GIF = GIF(dt=self.dt)

    def set_templates(self, aec=None, train=None, test=None):
        self.aec_template = aec
        self.train_template=train
        self.test_template = test

    def set_files_test(self, AECtrace=None, trainingset=1, testsets=[], filetype='Igor'):
        if AECtrace is not None:
            self.Exp.setAECTrace(
                os.path.join(self.path, self.aec_template.format(2, AECtrace)), 1.0,
                os.path.join(self.path, self.aec_template.format(3, AECtrace)), 1.0,
                10000.0, FILETYPE=filetype)
            print('Doing AEC')
            self.AEC = AEC_Badel(self.Exp.dt)
            self.AEC.p_expFitRange = [3.0,150.0]
            self.AEC.p_nbRep = 15
            self.Exp.setAEC(self.AEC)
            self.Exp.performAEC()
            # self.AEC.plotKopt()
            # self.AEC.plotKe()
            print('   ...AEC Complete')

        self.Exp.addTrainingSetTrace(
                os.path.join(self.path, self.train_template.format(2, trainingset)), 1.0,
                os.path.join(self.path, self.train_template.format(3, trainingset)), 1.0,
                120000.0, FILETYPE=filetype)

        for ds in testsets:
            self.Exp.addTestSetTrace(
                os.path.join(self.path, self.test_template.format(2, ds)), 1.0,
                os.path.join(self.path, self.test_template.format(3, ds)), 1.0,
                20000.0, FILETYPE=filetype)

    def set_eta_timescales(self, ts=[1.0, 20.]):
        self.eta_timescales = ts

    def set_gamma_timescales(self, ts=[1.0, 20.]):
        self.gamma_timescales = ts

    def fit(self, threshold=0., beforeSpike=5.0, current=None, ax=None, plot=True):
        """
        Perform fit against GIF model. 
        Parameters
        ----------
        threshold : float (default 0)
            Voltage threshold for spike detection, in mV

        beforeSpike : float (default 5.0)
            Time before spike to compute reverse correlation for fitting, in msec
        
        current : array (default None)
            Current injection to use for the fit (in nA)
            If None, this routine generates an arbitrary current.
    
        ax : matplotlib axis object (default None)
            if None, then data is plotted. 
        
        """
        self.Exp.detectSpikes(threshold=threshold, refractory=1.5)
        #self.Exp.plotTrainingSet()

        self.GIF.eta = Filter_Rect_LogSpaced()
        self.GIF.eta.setMetaParameters(length=5000.0, binsize_lb=2.0, binsize_ub=1000.0, slope=4.5)

        self.GIF.gamma = Filter_Rect_LogSpaced()
        self.GIF.gamma.setMetaParameters(length=5000.0, binsize_lb=5.0, binsize_ub=1000.0, slope=5.0)
        self.GIF.eta = Filter_Exps()
        self.GIF.eta.setFilter_Timescales(self.eta_timescales)

        self.GIF.gamma = Filter_Exps()
        self.GIF.gamma.setFilter_Timescales(self.gamma_timescales)

        #To perform the fit using only a specific part of the training set, use the following command before calling self.GIF.fit():
        self.Exp.trainingset_traces[0].setROI([[0,10000.0], [20000.0, 60000.0]])
        
        self.GIF.fit(self.Exp, DT_beforeSpike=beforeSpike, threshold=threshold)
        fittedpars = self.GIF.getParameters()
        self.GIF.plotParameters()
        tsmax = np.max(self.Exp.trainingset_traces[0].getTime())/1000.

        if current is None:
            gen = NTB()
            gen.set_params(skew=8., sigma0=0.2, fmod=0.2, tau=3.0, dt=self.Exp.dt, i0=0, dur=tsmax)
            tb, I = gen.generator()
        else:
            I = current
        V0 = -65.
        (time, V, I_a, V_t, S) = self.GIF.simulate(I, V0, pars=fittedpars)  # simulate response to current trace I with starting voltage V0
#        print('Simulated with fittedpars: \n', fittedpars)
        self.model_traces ={'time': time, 'V': V, 'I_a': I_a, 'I_stim': I, 'V_t': V_t,
            'gifpars': fittedpars._asdict()}
        self.model_traces['ExpData'] = self.Exp.trainingset_traces[0].V
        if plot:
            if ax is None:
                mpl.figure()
                mpl.suptitle('Fitted')
        
                mpl.plot(self.Exp.trainingset_traces[0].getTime(), 
                        self.Exp.trainingset_traces[0].V, 'k-',  linewidth=0.75)
                mpl.plot(time, V, 'r-', linewidth=0.5)

                mpl.show()
            else:
                ax.plot(self.Exp.trainingset_traces[0].getTime(), 
                        self.Exp.trainingset_traces[0].V, 'k-', linewidth=0.75)
                ax.plot(time, V, 'r-', linewidth=0.5)
        self.write_result('GFIT_Original_gifnoise.p')

    def test_simulator(self, current=None, fs=False):
        tsmax = 0.1*np.max(self.Exp.trainingset_traces[0].getTime())/1000.
        if current is None:
            gen = NTB()
            gen.set_params(skew=8., sigma0=0.5, fmod=0.2, tau=3.0, dt=self.Exp.dt, i0=0, dur=tsmax)
            tb, I = gen.generator()

            # tb, I = generator(i0=0, dt=self.Exp.dt, sigma0=0.5, fmod=0.2, tau=3.0,
            #  dur=tsmax)
        else:
            I = current
        V0 = -65.
        if fs is False:
            (time, V, I_a, V_t, S) = self.GIF.simulate(I, V0)  # simulate response to current trace I with starting voltage V0
        else: # "force spikes"
            spks = np.zeros_like(I)
            (time, V, eta_S) = self.GIF.simulateDeterministic_forceSpikes(I, V0, spks)
            I_a = np.zeros_like(V)
            V_t = np.zeros_like(V)
        fittedpars = self.GIF.getParameters()
        self.model_traces ={'time': time, 'V': V, 'I_a': I_a, 'I_stim': I, 'V_t': V_t,
            'gifpars': fittedpars._asdict()}
        self.model_traces['ExpData'] = self.Exp.trainingset_traces[0].V

        mpl.figure()
        mpl.suptitle('Simulated')
    
        mpl.plot(time, V, 'r-', linewidth=0.75)
#        mpl.plot(self.Exp.trainingset_traces[0].getTime(), self.Exp.trainingset_traces[0].V, 'k-', linewidth=0.5)
        mpl.show()

    def write_result(self, fn):
        h = open(fn, 'wb')
        import pickle
        pickle.dump(self.model_traces, h)
        h.close()
            
    def predictor(self):
        print('GIF predictor')
        self.Prediction = self.Exp.compareSpikes(self.GIF, nb_rep=500)
        print('Kistler')
        self.Md = self.Prediction.computeMD_Kistler(4.0, 0.1)    
        self.Prediction.plotRaster(delta=1000.0)

def test_original():
    # templates first argument is the "channel" and the second is the record
    aec_template = 'Cell3_Ger1Elec_ch{0:d}_{1:4d}.ibw'
    train_template = 'Cell3_Ger1Training_ch{0:d}_{1:4d}.ibw'
    test_template = 'Cell3_Ger1Test_ch{0:d}_{1:4d}.ibw'
    basepath = '/Users/pbmanis/Desktop/Python/GIFFittingToolbox/data/fi/'

    GF = GIFFitter(path=basepath)
    GF.GIF.gn = 0.0
    GF.GIF.refract = 8.0
    GF.set_templates(aec=aec_template, train=train_template, test=test_template)
    GF.set_files_test(AECtrace=None, trainingset=1008, testsets = range(1009, 1018), filetype='Igor')
    
    GF.fit(current=GF.Exp.trainingset_traces[0].I, threshold=0.)
    #GIFFitter.GF.plotAverageModel(GF)
    return GF

def fit_cell(cellno=19, gn=0.):
    threshold = -20.
    cell = cellno
    model = 'mGBC'
    cname = 'VCN_c%02d' % cell 
    basepath = '/Users/pbmanis/Desktop/Python/VCNModel'
    cellpath = 'VCN_Cells/%s/Simulations/Noise/' % cname
    
    testset = os.path.join(basepath, cellpath, '%s_%s_gifnoise.p' % (cname, model))
    h = open(testset, 'rb')
    d = pickle.load(h)
    tr = d['Results'][0][0]['monitor']
    V = np.array(tr['postsynapticV'])

    T = np.array(tr['time'])
    I = np.array(d['Results'][0][0]['stim'][1])

    V_units = 1e-3
    I_units = 1e-9

    GF = GIFFitter(path=testset, dt=0.025)
    GF.Exp.addTrainingSetTrace(V, V_units, I, I_units, np.max(T)-0.025*200., FILETYPE='Array')
    GF.Exp.detectSpikes(threshold=threshold, refractory=1.5)
    all_spks_times_trainingset = []

    for tr in GF.Exp.trainingset_traces:
        spks_times = tr.getSpikeTimes()
        all_spks_times_trainingset.append(spks_times)
    mean_firingrate =  np.mean(np.diff(all_spks_times_trainingset))
    GF.GIF.DV = 50.
    GF.GIF.Vt_star = -GF.GIF.DV * np.log(mean_firingrate)
    print GF.GIF.DV, GF.GIF.Vt_star
    # mpl.plot(T, V)
    # mpl.show()
    GF.set_eta_timescales(ts=[0.1, 1., 5., 20.])
    GF.set_gamma_timescales(ts=[0.3, 1.0, 5.0, 20.])
    GF.GIF.C = 0.1
    GF.GIF.gl = 0.05
    GF.GIF.gn = gn
 #   GF.GIF.DV = 4.0
    GF.GIF.lambda0 = 1.
    GF.GIF.dt = 0.1
    GF.GIF.Tref = 2.5
#    GF.test_simulator(current=None, fs=False)
#    exit(1)
    GF.fit(threshold=threshold, current=I, plot=False)
    print 'returned'
    GF.write_result('GFIT_%s_%s_gn=%6.3f_gifnoise.p' % (cname, model, GF.GIF.gn))
       

if __name__ == '__main__':
    # gf = test_original()
    # gf.GIF.plotAverageModel([gf])
    # exit(1)
    for gn in [0., 0.01, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10.]:
        fit_cell(cellno=19, gn=gn)

    