__author__ = 'pbmanis'

"""
showmap reads the basic file and then shows the maps of the results.

"""
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
import os

simpath = 'Simulations'

fn = 'Stellate_XM13_ihvcn_kht_14.03.02-16.32.11.p'

class ShowMap():
    def __init__(self):

        self.loadfile(fn)

        self.getids()
        self.plotVm()


    def loadfile(self, fn):
        f = open(fn, 'rb')
        self.r = pickle.load(f)
        f.close()

    def getids(self):
        self.ids = []
        for i in self.r['spacemap']:
            self.ids.append(i['id'])


    def plotVm(self):
        res = self.r['result']
        npts = len(self.r['spacemap'])
        idn = np.zeros(npts)
        ih = np.zeros(npts)
        ik = np.zeros(npts)
        vm = np.zeros(npts)
        ns = np.zeros(npts)
        tauih = np.zeros(npts)
        taum = np.zeros(npts)
        bfile = {}
        fx = pl.figure(99)
        afx = fx.add_subplot(111)
        for i in self.r['spacemap']:
            idn = i['id']
            ih[idn]= i['ihvcn']
            ik[idn] = i['kht']
            vm[idn] = res[idn]['Vm']
            ns[idn] = np.max(res[idn]['spikes']['n'])
            tauih[idn] = res[idn]['tauih']
            taum[idn] = res[idn]['taum']
            bfile[idn] = res[idn]['basefile']
            bf = open(bfile[idn]+'_mrun_ID%04d.p' % idn, 'rb')
            bdata = pickle.load(bf)
            #print '\n\n', bdata['Results'][0][0.0]['monitor']['postsynapticV']
            for br in bdata['Results']:
                for inj in br.keys():
                    v = br[inj]['monitor']['postsynapticV']
                    t = br[inj]['monitor']['time']
                    afx.plot(t, v)

        #print 'Vm, ih, ik: ', vm, ih, ik, ns, tauih, taum

            # make grids of the data sets
        nx = np.unique(ih)
        ny = np.unique(ik)

        Gx, Gy = np.meshgrid(nx, ny)
        fig=pl.figure(1)

        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(ik, ih, ns)
        ax.set_xlabel('kht')
        ax.set_ylabel('ihvcn')
        ax.set_zlabel('Spikes')
        pl.show()

    def plotTraces(self):
        pass



if __name__ == "__main__":
    ShowMap()
