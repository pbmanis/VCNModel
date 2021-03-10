#!/usr/bin/python
from __future__ import print_function
import os
import subprocess
import logging

logging.basicConfig(filename='george_run.log',level=logging.DEBUG, format='%(asctime)s %(message)s', filemode='w')

cells = [9, 14, 16, 17, 18, 19, 20, 21, 22]
missingcells =[11]
runid = ''
for cellno in cells:
    cmdinit = ['python', 'model_run.py', 'VCN_c%02d' % cellno, '@spirou-endbulbs.prm', '-P', 'initAN']
    cmdrun1 = ['python', 'model_run.py', 'VCN_c%02d' % cellno, '@spirou-endbulbs.prm', '-P', 'runANPSTH', '--spirou', 'all']
    cmdrun2 = ['python', 'model_run.py', 'VCN_c%02d' % cellno, '@spirou-endbulbs.prm', '-P', 'runANPSTH',  '--spirou', 'all=mean']
    cmdrun3 = ['python', 'model_run.py', 'VCN_c%02d' % cellno, '@spirou-endbulbs.prm', '-P', 'runANPSTH',  '--spirou', 'max=mean']
    
    try:
        subprocess.call(cmdinit)
        logging.debug('RunID = %s, Cell %d (VCN_c%02d) init run complete' % (runid, cellno, cellno))
        retflag = subprocess.call(cmdrun1)
        if retflag != 0:
            break
        logging.debug('RunID = %s, Cell %d (VCN_c%02d) run1 complete' % (runid, cellno, cellno))
        subprocess.call(cmdrun2)
        logging.debug('RunID = %s, Cell %d (VCN_c%02d) run2 complete' % (runid, cellno, cellno))
        subprocess.call(cmdrun3)
        logging.debug('RunID = %s, Cell %d (VCN_c%02d) run3 complete' % (runid, cellno, cellno))
    except:
        logging.debug('RunID = %s, Cell %d (VCN_c%02d) run failed' % (runid, cellno, cellno))
        
print('finis')

        