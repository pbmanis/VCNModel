__author__ = 'pbmanis'

"""
showmap reads the basic file and then shows the maps of the results.

"""
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
import os
from sets import *

simpath = 'Simulations'
showtraces = False

#fn = 'Stellate_XM13_ihvcn_kht_14.03.02-20.36.31.p'
#fn = 'Stellate_XM13_ihvcn_kht_14.03.03-00.47.17.p'
#fn = 'Stellate_XM13_ihvcn_kht_14.03.06-11.38.48.p'
fn = 'Stellate_XM13_ihvcn_kht_14.03.06-12.05.07.p'

class ShowMap():
    def __init__(self, fname=None):
        if fname == None:
            fname = fn
        print 'loading file: %s' % fname
        self.loadfile(fname)

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
        sptrain = [None]*npts
        iinj = [None]*npts
        tauih = np.zeros(npts)
        taum = np.zeros(npts)
        bfile = {}
        apfig = False
        for i in self.r['spacemap']:
            idn = i['id']
            ih[idn]= i['ihvcn']
            ik[idn] = i['kht']
            vm[idn] = res[idn]['Vm']
            ns[idn] = np.max(res[idn]['spikes']['n'])
            sptrain[idn] = res[idn]['spikes']['n']
            iinj[idn] = res[idn]['spikes']['i']
            tauih[idn] = res[idn]['tauih']
            taum[idn] = res[idn]['taum']
            bfile[idn] = res[idn]['basefile']
            if showtraces:
                bf = open(bfile[idn]+'_mrun_ID%04d.p' % idn, 'rb')
                bdata = pickle.load(bf)
                #print '\n\n', bdata['Results'][0][0.0]['monitor']['postsynapticV']
                if apfig is False:
                    apf = pl.figure(11)
                    apfig = True
                    apfx = apf.add_subplot(111)
                for br in bdata['Results']:
                    for inj in br.keys():
                        v = br[inj]['monitor']['postsynapticV']
                        t = br[inj]['monitor']['time']
                        apfx.plot(t, v)

        #print 'Vm, ih, ik: ', vm, ih, ik, ns, tauih, taum

            # make grids of the data sets
       #ni = np.unique(np.flatten(iinj))
        ihu = np.unique(ih)
        iku = np.unique(ik)
        #Gx, Gy = np.meshgrid(ni, nx)
        vmx = np.zeros((len(ihu), len(iku)))
        taumx = np.zeros((len(ihu), len(iku)))
        fx = pl.figure(99)
        afx1 = fx.add_subplot(221)
        afx2 = fx.add_subplot(222, projection='3d')
        afx3 = fx.add_subplot(223, projection='3d')
        afx4 = fx.add_subplot(224, projection='3d')
        #print iinj
        #print len(iinj)
        for j, iht in enumerate(ihu):
            for i, kht in enumerate(iku): # for each level of kht
                ikindx = [x for x in np.where(ik == kht)[0]]
                ihindx = [y for y in np.where(ih == iht)[0]]
                ind = np.array(list(set(ikindx).intersection(set(ihindx))))
                ix = iinj[ind]
                spx = sptrain[ind]
                iki = np.argsort(ix)
                taumx[j,i] = taum[ind]
                vmx[j,i] = vm[ind]
               # print taumx
                afx1.plot(ix[iki], spx[iki], 'ro-')
        afx2.scatter(ih, ik, taum)
        afx2.set_zlabel('Taum')
        afx2.set_xlabel('x gkht')
        afx2.set_ylabel('x gih')

        afx3.scatter(ih, ik, vm)
        afx3.set_zlabel('Vm')
        afx3.set_xlabel('x gkht')
        afx3.set_ylabel('x gih')

        afx4.scatter(ih, ik, ns)
        afx4.set_zlabel('N Spikes')
        afx4.set_xlabel('x gkht')
        afx4.set_ylabel('x gih')

        Gx, Gy = np.meshgrid(ihu, iku)
        fig=pl.figure(1)

        #ax = fig.add_subplot(221, projection='3d')
        #ax.scatter(ik, ih, ns)
        ax1 = fig.add_subplot(221)
        ax1.contour([ik, ih, ns])
        ax1.set_xlabel('kht')
        ax1.set_ylabel('ihvcn')
        ax1.set_title('N Spikes')
        #ax1.set_zlabel('ns')
        ax2 = fig.add_subplot(222)
        ax2.imshow([ik, ih, taum])
        ax2.set_xlabel('kht')
        ax2.set_ylabel('ihvcn')
        ax2.set_title('Tau m')
        #ax1.set_zlabel('taum')
        ax3 = fig.add_subplot(223)
        ax3.imshow([ik, ih, vm])
        ax3.set_xlabel('kht')
        ax3.set_ylabel('ihvcn')
        ax3.set_title('Vm')
        #ax3.set_zlabel('Vm')
        pl.show()

    def plotTraces(self):
        pass



if __name__ == "__main__":
    ShowMap()
