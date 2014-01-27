__author__ = 'pbmanis'

import matplotlib.pylab as MP
import pickle

def plotResults(res, runInfo):
    fig = MP.figure(100)
    p1 = fig.add_subplot(4,1,1)
    p2 = fig.add_subplot(4,1,2)
    p3 = fig.add_subplot(4,1,3)
    p4 = fig.add_subplot(4,1,4)
    p1.plot(res['vec']['time'], res['vec']['axon'], color=runInfo['clist'][0])
   # p2.plot(self.vec['time'], self.ica['axon'])
    print res.keys()
    for v in res['Voltages']:
        #print v
        p1.plot(res['vec']['time'], res['Voltages'][v], color=runInfo['clist'][v])
    p1.set_ylabel('V')
    for c in res['ICa']:
        p2.plot(res['vec']['time'], res['ICa'][c], color=runInfo['clist'][c])
        p3.plot(res['vec']['time'], res['Cai'][c], color=runInfo['clist'][c])
    p2.set_ylabel('I_{Ca}')
    p3.set_ylabel('[Ca]_i')
    p4.plot(res['vec']['time'], res['vec']['postsynaptic'], color = 'k')
    #p2.set_ylim(-5e-12, 1e-12)

    MP.show()

if __name__ == "__main__":


    fo = open('Canonical/Normal.p')
    d = pickle.load(fo)
    fo.close()
    print d['runInfo'].keys()

    plotResults(d['Results'], d['runInfo'])