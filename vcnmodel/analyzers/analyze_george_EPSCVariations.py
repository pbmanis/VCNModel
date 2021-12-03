#!/usr/bin/python
from __future__ import print_function
import os
import subprocess
import logging
import pickle

import vcnmodel.vspfile as vspfile


def main():
    logging.basicConfig(filename='george_run_analyze.log',level=logging.DEBUG, format='%(asctime)s %(message)s', filemode='w')

    cells = [9, 11, 14, 16, 17, 18, 19, 20, 21, 22]
    runtypes = ['all', 'all=mean', 'max=mean']

    #cells = [9]
    runid = ''

    vspfile.init_vsp(cells, runtypes)

    for cellno in cells:
        for runtype in runtypes:
            cmdrun3 = ['python', 'plotters/displayresult_model_run.py', 'VCN_c%02d' % cellno,
             '@spirou-endbulbs.prm', '-P', 'runANPSTH',  '--spirou', '%s'%runtype, '--noplot']

            try:
                res = subprocess.call(cmdrun3)
                logging.debug('RunID = %s, Cell %d (VCN_c%02d) run1 complete' % (runid, cellno, cellno))
            except:
                logging.debug('RunID = %s, Cell %d (VCN_c%02d) run failed' % (runid, cellno, cellno))

        
    print('finis')

    vspfile.print_vsp()

if __name__ == '__main__':
    main()