#!/usr/bin/python
"""
Run one gbc... but know about them all

call: all_gbc_AN.py #

"""
from __future__ import print_function
import os
import sys
import pickle
from pathlib import Path
import model_run as mrun

class OneANRun(object):
    def __init__(self, gbc_name, modelName, modelType, protocol, SR, nrep, forcerun=False, testing=False, inflateflag=False, initialize=False):
        self.forcerun = forcerun
        self.testflag = testing
        self.inflateflag = inflateflag
        sim_reportfile = Path('lastsim.txt')
        
        if SR not in ['LS', 'MS', 'HS']:
            raise ValueError('SR must be one of [LS, MS or HS]')
        if nrep > 999 or nrep < 1:
            raise ValueError('nrep must be > 0 and <= 999')


        print('-'*32)
        print('   all_gbc_AN: entry parameters')
        print('      gbc#: ', gbc_name)
        print('      SR: ', SR)
        print('      nrep: ', nrep)
        print('      forcerun: ', self.forcerun)
        print('      testflag: ', self.testflag)
        print('      initialize: ', initialize)
        # anresultfilenames = {}
        # for i, n in enumerate(gbc_nos):
        #     anresultfilenames[i] = 'AN_Result_VCN_c{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(n, nrep, SR)
        # print('   anresult files:')
        # for i, n in enumerate(gbc_nos):
        #     print('      ', n, anresultfilenames[i])
        # print('-'*32)
        if self.testflag:
            exit()

        seeds = [100] # use all the same seeds

        # create paths to the simulation runs to check for existing AN initialization


        M = mrun.ModelRun() # create an instance
            
        cell = f'VCN_c{gbc_name:s}'
        cell_ax = f'VCNc{gbc_name:s}'
        M.Params['cell'] = cell
        M.Params['Parallel'] = True
        M.Params['modelType'] = modelType
        M.Params['modelName'] = modelName
        M.Params['runProtocol'] = protocol
        M.Params['soma_autoinflate'] = self.inflateflag   
        M.Params['ANSynapseType'] = 'multisite'
        
        M.Params['SGCmodelType'] = 'cochlea'
        M.Params['soundtype'] = 'tonepip'
        M.Params['SR'] = SR
        M.Params['SRType'] = SR
        M.Params['run_duration'] = 0.3
        M.Params['plotFlag'] = False
        M.Params['nReps'] = nrep
        
        # if M.Params['hocfile'] == None: # just use the matching hoc file
        print('mparams cell: ', M.Params['cell'])
        M.Params['hocfile'] = M.Params['cell'] + '.hoc'
        print(f"{' ':14s}HOC file: {str(M.Params['hocfile']):s}")
        M.setup_model(par_map=M.Params)
        print(f"{' ':14s}Simulation file:        {str(M.Params['simulationFilename']):s}")
        print(f"{' ':14s}Initialization ANfile:  {str(M.Params['initANStateFile']):s}")
        # ivinitfile = os.path.join(M.initDirectory, 'ANneuronState.dat')
        # initf = os.path.join(M.baseDirectory, cell, ivinitfile)
        print (f"\n{' ':14s}all_gbc_AN Running Cell {cell:s}")

        initf = M.Params['initANStateFile']
        if not initf.is_file() or initialize:
            print(f'{" ":14s}Creating new init file for cell {str(initf):s}')
            M.Params['runProtocol'] = 'initAN'
            if not testing:
                M.run_model(par_map = M.Params)
                if sim_reportfile.is_file():
                    sim_reportfile.unlink() # remove the report file, it is not needed here
                return
        elif self.forcerun:
            print(f'{" ":14s}Creating new init file for cell; forcerun is true {str(initf):s}')
            M.Params['runProtocol'] = 'initAN'
            if not testing:
                M.run_model(par_map = M.Params)
                return
        else:
            print(f'{" ":14s}>>>>>Initialization file for {cell:s} exists and will not be updated')
        if not initf.is_file():  # make sure it completed
            raise ValueError(f'{" ":14s}***** Failed to create the initialization state file {str(initf):s} *****')
    
        ini_timestamp = initf.stat().st_mtime
        simf = M.Params['simulationFilename']
        if simf.is_file():
            if simf.stat().st_mtime < ini_timestamp:  # if last run is prior to last initiazization, update run
                forcerun = True
                print(f'{" ":14s}Sim file is older than init; force update')
                print(f"{' ':14s}{str(simf.stat().st_mtime):s} <= {str(ini_timestamp):s}, Forcerun: {str(forcerun):s}")
                simf.unlink()  # remove the file
        # forcerun = True
        if forcerun and not initialize:
            # try:
            print(f'{" ":10s}Running with protocol: {protocol:s}')
            M.Params['runProtocol'] = protocol
            M.run_model(par_map = M.Params)
            sim_reportfile.write_text(M.Params['simulationFilename'])
            # with open('lastsim.txt', 'w') as fh:
            #     fh.write(str(laststimt))
            #     print('wrote simulation file information to ', laststimt)
            # except:  # try to do the initialization again, then run
            #     M.Params['runProtocol'] = 'initAN'
            #     M.run_model(par_map = M.Params)
            #     M.Params['runProtocol'] = 'runANPSTH'
            #     M.run_model(par_map = M.Params)


if __name__ == '__main__':
    with open('oneanrun.p', 'rb') as fh:
        d = pickle.load(fh)
    
    OneANRun(d['gbc'], d['modelName'], d['modelType'], d['protocol'], d['SR'], d['nrep'], 
        forcerun=d['forcerun'], testing=d['testing'], inflateflag=d['inflateflag'], initialize=d['initialize'])
    
