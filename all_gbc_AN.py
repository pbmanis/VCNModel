#!/usr/bin/python
"""
Run one gbc... but know about them all

call: all_gbc_AN.py #

"""
from __future__ import print_function
import os
import sys


forcerun = False
testflag = False
print('sys argv: ', sys.argv)
if len(sys.argv) > 1:
    gbc_no = int(sys.argv[1])
    SR = sys.argv[2]
    if SR not in ['LS', 'MS', 'HS']:
        raise ValueError('SR must be one of [LS, MS or HS]')
    nrep = int(sys.argv[3])
    if nrep > 999 or nrep < 1:
        raise ValueError('nrep must be > 0 and <= 999')
        
    if 'forcerun' in sys.argv[1:]:
        forcerun = True
        
    if 'test' in sys.argv[1:]:
        testflag = True
    
    if gbc_no not in [8, 9, 17, 18, 19, 20, 21, 22]:
        raise ValueError('GBC %d not implemented' % gbc_no)
    gbc_nos = [gbc_no] # [8, 9, 17, 18, 19, 20, 21, 22]
else:  # defaults
    gbc_nos = [8, 9, 17, 18, 19, 20, 21, 22]
    SR = 'MS'
    nrep = 2

print('-'*32)
print('   all_gbc_AN: entry parameters')
print('      gbc#: ', gbc_no)
print('      SR: ', SR)
print('      nrep: ', nrep)
print('      forcerun: ', forcerun)
print('      testflag: ', testflag)
anresultfilenames = []
for i, n in enumerate(gbc_nos):
    anresultfilenames.append('AN_Result_VCN_c{0:02d}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(n, nrep, SR))
print('   anresult files:')
for i in range(len(gbc_nos)):
    print('      ', anresultfilenames[i])
print('-'*32)
if testflag:
    exit()

seeds = [100]*len(gbc_nos)  # use all the same seeds

# create paths to the simulation runs to check for existing IV initialization

for i, n in enumerate(gbc_nos):
    cell = 'VCN_c{0:02d}'.format(n)
    if i == 0:
        import model_run as mrun
    else:
        mrun = reload(mrun)
    M = mrun.ModelRun() # create a new instance for each cell
    ivinitfile = os.path.join(M.initDirectory, 'ANneuronState.dat')
    initf = os.path.join(M.baseDirectory, cell, ivinitfile)
    print ('\nRunning Cell {0:s}: '.format(cell))
    M.Params['cell'] = cell
    M.Params['infile'] = M.Params['cell'] + '.hoc'    
    if not os.path.isfile(initf):
        print('creating new init file for cell, did not find {:s}'.format(initf))
        M.Params['runProtocol'] = 'initAN'
#        try:
        M.run_model(par_map = M.Params)
#        except:
#            print('Failed to run model for {:s}'.format(initf))
    else:
        print('    Initialization file for {:s} exists'.format(cell))
    if not os.path.isfile(initf):
        raise ValueError('Failed to create the AB init state file')
    andatafile = os.path.join(M.baseDirectory, cell, M.simDirectory, 'AN',
             anresultfilenames[i])
    if not os.path.isfile(andatafile) or forcerun:
        print ('Creating AN datafile: {:s}\n'.format(andatafile))
        M = mrun.ModelRun() # create a new instance for each cell
        M.Params['cell'] = cell
        M.Params['infile'] = M.Params['cell'] + '.hoc'    
        M.Params['runProtocol'] = 'runANPSTH'
        M.Params['SGCmodelType'] = 'cochlea'
        M.Params['soundtype'] = 'tonepip'
        M.Params['SR'] = SR
        M.Params['SRType'] = SR
        M.Params['run_duration'] = 0.3
        M.Params['plotFlag'] = False
        M.Params['nReps'] = nrep
        try:
            M.run_model(par_map = M.Params)
        except:  # try to do the initialization again, then run
            M.Params['runProtocol'] = 'initAN'
            M.run_model(par_map = M.Params)
            M.Params['runProtocol'] = 'runANPSTH'
            M.run_model(par_map = M.Params)
            
#        fh = open(M.ANFilename)    
#    else:
#        fh = open(andatafile)


