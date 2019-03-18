#!/usr/bin/python
"""
Run one gbc... but know about them all

call: all_gbc_AN.py #

"""
from __future__ import print_function
import os
import sys
import pickle

class OneANRun(object):
    def __init__(self, gbc_name, modelName, modelType, protocol, SR, nrep, forcerun=False, testing=False, inflateflag=False):
        self.forcerun = forcerun
        self.testflag = testing
        self.inflateflag = inflateflag

        if SR not in ['LS', 'MS', 'HS']:
            raise ValueError('SR must be one of [LS, MS or HS]')
        if nrep > 999 or nrep < 1:
            raise ValueError('nrep must be > 0 and <= 999')
    
        if not isinstance(gbc_name, list):
            gbc_nos = [gbc_name] # [8, 9, 17, 18, 19, 20, 21, 22]
        else:
            gbc_nos = gbc_name

        print('-'*32)
        print('   all_gbc_AN: entry parameters')
        print('      gbc#: ', gbc_nos)
        print('      SR: ', SR)
        print('      nrep: ', nrep)
        print('      forcerun: ', self.forcerun)
        print('      testflag: ', self.testflag)
        # anresultfilenames = {}
        # for i, n in enumerate(gbc_nos):
        #     anresultfilenames[i] = 'AN_Result_VCN_c{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(n, nrep, SR)
        # print('   anresult files:')
        # for i, n in enumerate(gbc_nos):
        #     print('      ', n, anresultfilenames[i])
        # print('-'*32)
        if self.testflag:
            exit()

        seeds = [100]*len(gbc_nos)  # use all the same seeds

        # create paths to the simulation runs to check for existing AN initialization

        for i, n in enumerate(gbc_nos):
            if i == 0:
                import model_run as mrun
            else:
                mrun = reload(mrun)
            M = mrun.ModelRun() # create an instance
            
            cell = f'VCN_c{n:s}'
            cell_ax = f'VCNc{n:s}'
            M.Params['cell'] = cell
            M.Params['Parallel'] = False
            M.Params['modelType'] = modelType
            M.Params['modelName'] = modelName
            M.Params['runProtocol'] = protocol
            M.Params['soma_autoinflate'] = self.inflateflag   
            
            M.Params['SGCmodelType'] = 'cochlea'
            M.Params['soundtype'] = 'tonepip'
            M.Params['SR'] = SR
            M.Params['SRType'] = SR
            M.Params['run_duration'] = 0.3
            M.Params['plotFlag'] = False
            M.Params['nReps'] = nrep
            
            if M.Params['hocfile'] == None: # just use the matching hoc file
                M.Params['hocfile'] = M.Params['cell'] + '.hoc'
            print('hoc file: ', M.Params['hocfile'])
            M.setup_model(par_map=M.Params)
            print('Simulation file: ', M.Params['simulationFilename'])
            print('Initialization ANfile: ', M.Params['initANStateFile'])
            # ivinitfile = os.path.join(M.initDirectory, 'ANneuronState.dat')
            # initf = os.path.join(M.baseDirectory, cell, ivinitfile)
            print ('\nall_gbc_AN Running Cell {0:s}: '.format(cell))
    
            initf = M.Params['initANStateFile']
            if not initf.is_file():
                print(f'creating new init file for cell, did not find {str(initf):s}')
                M.Params['runProtocol'] = 'initAN'
                if not testing:
                    M.run_model(par_map = M.Params)
            elif self.forcerun:
                print(f'creating new init file for cell; forcerun is true {str(initf):s}')
                M.Params['runProtocol'] = 'initAN'
                if not testing:
                    M.run_model(par_map = M.Params)
            else:
                print(f'    Initialization file for {cell:s} exists')
            if not initf.is_file():  # make sure it completed
                raise ValueError(f'Failed to create the initialization state file {str(initf):s}')
        
            ini_timestamp = initf.stat().st_mtime
            simf = M.Params['simulationFilename']
            if simf.is_file():
                if simf.stat().st_mtime < ini_timestamp:
                    forcerun = True
                    print(simf.stat().st_mtime <= ini_timestamp, forcerun)
            # forcerun = True
            if forcerun:
                # try:
                print('Running with protocol: ', protocol)
                M.Params['runProtocol'] = protocol
                M.run_model(par_map = M.Params)
                laststimt = M.Params['simulationFilename']
                with open('lastsim.txt', 'w') as fh:
                    fh.write(str(laststimt))
                    print('wrote simulation file information to ', laststimt)
                # except:  # try to do the initialization again, then run
                #     M.Params['runProtocol'] = 'initAN'
                #     M.run_model(par_map = M.Params)
                #     M.Params['runProtocol'] = 'runANPSTH'
                #     M.run_model(par_map = M.Params)
            
#        fh = open(M.ANFilename)    
#    else:
#        fh = open(andatafile)

if __name__ == '__main__':
    with open('oneanrun.p', 'rb') as fh:
        d = pickle.load(fh)
    
    OneANRun(d['gbc'], d['modelName'], d['modelType'], d['protocol'], d['SR'], d['nrep'], 
        forcerun=d['forcerun'], testing=d['testing'], inflateflag=d['inflateflag'])
    
