#!/usr/bin/env python

# show traces from simulation run (current-voltage - IV)
#
from __future__ import print_function
import os
import sys
import pickle
import matplotlib.pyplot as mpl

def main():
    # modify to point to the file:
    basepath = '/Users/pbmanis/Desktop/Python/VCNModel/VCN_Cells/VCN_c09/Simulations/'

    # fn = os.path.join(basepath, 'IV/VCN_c09_mGBC_monitor.p')  # IV
    # datatype = 'IV'
    fn = os.path.join(basepath, 'AN/AN_Result_VCN_c09_mGBC_delays_N001_035dB_4000.0_HS.p')  # an AN response
    datatype = 'AN'

    fh = open(fn)
    print(' file opened')
    d = pickle.load(fh)


    # FOR IVS:
    # d is a dictionary with a data structure:
    # ['runInfo', 'basename', 'modelPars', 'Results']
    # runInfo is a dict with information about the run (general info)
    # modelpars holds information about the cell used in the run (status from cnmodel.cell)
    # basename is the base file name, referenced to the top directory holding the cells (dir is hierachal)
    # Results holds the results as follows:
    # d['Results'][##] indexes the data for a given current injection, as a dict. The dict has one key,
    # which is the actual current level, in nA, thus,
    # d['Results[0][-1.5] points to the data for the most negative current injection
    # let di = d['Results'][0][-1.5]
    # di is a dictionary with keys: ['stim', 'monitor', 'runInfo', 'distanceMap', 'vec', 'Sections']
    # di['vec'] has the voltage traces for the selected traces, as a dictionary with keys 
    # by the names of the sections. These names match the names used in the hoc file (which is how
    # you would map these back to the original morphology)
    # di['monitor']['time] has the time base for all of the traces in 'vec'
    # Note that there are other keys in the results data, including the always-recorded soma V


    # for AN result data, the structure is different

    # print d.keys()
    # for dx in d.keys():
    #     for k in d[dx]:
    #         if dx not in ['time']:
    #             print dx, k, dir(d[dx][k])
    #             print' '
    #

    if datatype == 'IV' and len(sys.argv) == 1:
        dibase = d['Results'][0]
        k = dibase.keys()[0]
        di = dibase[k]
        ds = di['vec']['sections[0]']  # happens to be the soma section
        dd1 = di['vec']['sections[93]']  # some dendrite sections
        dd2 = di['vec']['sections[211]']
        t = di['monitor']['time']
        mpl.plot(t, ds)
        mpl.plot(t, dd1)
        mpl.plot(t, dd2)
        mpl.show()

    if datatype == 'AN'  and len(sys.argv) == 1:
        di = d  # direct, not a list; reps are stored inside dicts.
        rep = 0  # which repetition to extract
        t = di['time']
        for i in range(0, len(di['allDendriteVoltages'][rep].keys()), 10):  # subset
            ds = di['allDendriteVoltages'][rep]['sections[%i]'%i]  # happens to be the soma section
            mpl.plot(t, ds)
        mpl.show()
    

    if len(sys.argv) > 1:
        fn = sys.argv[1]
        di = d  # direct, not a list; reps are stored inside dicts.
        rep = 0  # which repetition to extract
        t = di['time']
        dout = {'time': t}
        for i in range(0, len(di['allDendriteVoltages'][rep].keys())):  # subset
            ds = di['allDendriteVoltages'][rep]['sections[%i]'%i]  # happens to be the soma section
            dout['sections[%i]'%i] = ds
        with open(fn, 'wb') as f:  # should be platform agnostic
            pickle.dump(dout, f)

if __name__ == '__main__':
    main()



