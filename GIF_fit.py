"""
GIF fitter wrapper. 
Set up with ln -s giffitterlocation/src gif in the running directory

Basic run will do a test against the datasets provided at the original site.

"""
import os
from gif.Experiment import *
from gif.AEC_Badel import *
from gif.GIF import *
from gif.Filter_Rect_LogSpaced import *
from gif.Filter_Exps import *
from NoiseTrainingGen import generator

class GIFFitter():
    def __init__(self, path, dt=0.1):
        self.path = path
        self.dt = dt
        self.set_timescales()
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

    def set_timescales(self, ts=[1.0, 5.0, 30.0, 70.0, 100.0, 500.0]):
        self.timescales = ts

    def fit(self, threshold=0., refract=1.0, beforeSpike=4.0):

        self.Exp.detectSpikes(threshold=threshold, ref=refract)
        # self.Exp.plotTrainingSet()
        # self.Exp.plotTestSet()

        self.GIF.Tref = refract    # refractory period

        print('Setting up GIF')
        self.GIF.eta = Filter_Rect_LogSpaced()
        self.GIF.eta.setMetaParameters(length=5000.0, binsize_lb=2.0, binsize_ub=1000.0, slope=4.5)

        self.GIF.gamma = Filter_Rect_LogSpaced()
        self.GIF.gamma.setMetaParameters(length=5000.0, binsize_lb=5.0, binsize_ub=1000.0, slope=5.0)
        self.GIF.eta = Filter_Exps()
        self.GIF.eta.setFilter_Timescales(self.timescales)

        self.GIF.gamma = Filter_Exps()
        self.GIF.gamma.setFilter_Timescales(self.timescales)

        #To perform the fit using only a specific part of the training set, use the following command before calling self.GIF.fit():
        #self.Exp.trainingset_traces[0].setROI([[0,10000.0], [20000.0, 60000.0]])

        self.GIF.fit(self.Exp, DT_beforeSpike=beforeSpike)
        print('GIF fitted')
        # self.GIF.save('./self.GIF.pck')
        # self.GIF_reloaded = GIF.load('./self.GIF.pck')
        # print dir(self.GIF_reloaded)

        self.GIF.printParameters()
        #self.GIF.plotParameters()
        
        # tb, I = generator(i0=0, dt=self.Exp.dt, sigma0=0.2, fmod=0.2, tau=3.0, dur=20.0)
        # V0 = -65
        # (time, V, I_a, V_t, S) = self.GIF.simulate(I, V0)  # simulate response to current trace I with starting voltage V0
        # plt.plot(time, V)
        # plt.show()
        
    def predictor(self):
        print('GIF predictor')
        self.Prediction = self.Exp.compareSpikes(self.GIF, nb_rep=500)
        print('Kistler')
        self.Md = self.Prediction.computeMD_Kistler(4.0, 0.1)    
        self.Prediction.plotRaster(delta=1000.0)

def test_original():
    # templates first argument is the "channel" and teh second is the record
    aec_template = 'Cell3_Ger1Elec_ch{0:d}_{1:4d}.ibw'
    train_template = 'Cell3_Ger1Training_ch{0:d}_{1:4d}.ibw'
    test_template = 'Cell3_Ger1Test_ch{0:d}_{1:4d}.ibw'
    basepath = '/Users/pbmanis/Desktop/Python/GIFFittingToolbox/data/fi/'

    GF = GIFFitter(path=basepath)
    GF.set_templates(aec=aec_template, train=train_template, test=test_template)
    GF.set_files_test(AECtrace=None, trainingset=1008, testsets = range(1009, 1018), filetype='Igor')
    
    GF.fit()    

if __name__ == '__main__':
    test_original()