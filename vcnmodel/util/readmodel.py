from pathlib import Path
import numpy as np
import pickle
import matplotlib
from dataclasses import dataclass
import dataclasses
from typing import Union
import datetime

from ephys.ephysanalysis import MakeClamps
import vcnmodel.util.fixpicklemodule as FPM
from vcnmodel.util import Params  # for old simulation results only

class ReadModel():
    def __init__(self):
       self.MC = MakeClamps.MakeClamps()
    
    def read_pfile(self, datasource:Union[str, Path, dict, None]=None, filemode:str='vcnmodel.v0', vscale:float=1e-3, iscale:float=1e-9, plot=False):
        """
        Read a pickled file; optionally plot the data
        Puts the data into a Clamps structure.

        Parameters
        ----------
        filename : str or Path
            The file to be read
        
        mode: str
            version name
        plot: Boolean (default: False)
            Flag to specify plotting the data in a simple mode
        """
        

        """
        The runInfo dictionary holds somethign like this:
        runinfo:  {'folder': PosixPath('VCN_Cells/VCN_c08/Simulations/IV'), 'fileName': 'Normal', 'runName': 'Run', 
        'manipulation': 'Canonical', 'preMode': 'cc', 'postMode': 'cc', 'TargetCellType': 'Bushy', 
        'electrodeSection': 'soma', 'dendriticElectrodeSection': 'dendrite', 
        'dendriticSectionDistance': 100.0, 'celsius': 37, 'nStim': 1, 
        'stimFreq': 200.0, 'stimInj': {'pulse': [-1.0, 2.01, 0.2]}, 
        'stimDur': 100.0, 'stimDelay': 5.0, 'stimPost': 3.0, 
        'vnStim': 1, 'vstimFreq': 200.0, 'vstimInj': 50, 
        'vstimDur': 50.0, 'vstimDelay': 2.0, 'vstimPost': 3.0, 'vstimHolding': -60, 
        'gif_i0': 0.0, 'gif_sigma': 0.5, 'gif_fmod': 0.2, 'gif_tau': 3.0, 
        'gif_dur': 10.0, 'gif_skew': 0.0, 
        'runTime': 'Wed Oct  9 13:05:54 2019', 
        'inFile': None, 'inFileRep': 1, 'spikeTimeList': {}, 
        'v_init': -61.0, 'useSaveState': True, 'tstop': 8.0, 'filename': 'VCN_c08_pulse_'}
        """
        # print('\nmodelPars: ', df['modelPars'])
        """
        The modelPars dict holds the following:
        modelPars:  {'species': 'mouse', 'cellClass': 'bushy', 'modelType': 'II', 
        'modelName': 'mGBC', 'soma': True, 'axon': False, 
        'dendrites': False, 'pumps': False, 'hillock': False, 
        'initialsegment': False, 'myelinatedaxon': False, 
        'unmyelinatedaxon': False, 'na': 'nav11', 'ttx': False, 
        'name': 'bushy', 'morphology': 'VCN_Cells/VCN_c08/Morphology/VCN_c08.hoc', 
        'temperature': 34.0}
        
        Note 10/28/2019 changed structure so that runInfo and modelPars are both 
        subdictionaries of Params (filemode is 'vcnmodel.v0')
        ... and undone later, so that all are top-level (filemode is 'vcnmodel.v1')
        """
        if isinstance(datasource, (str, Path)):
            with open(datasource, 'rb') as fh:
                df = FPM.pickle_load(fh)
                filename = datasource
        else:
            df = datasource
            print("df.Params: ", df.Params)
            raise ValueError()
            
        if filemode in ['vcnmodel.v0']:
            print(f"Reading model file in version v0:, with {len(df['Results']):4d} trials")
        elif filemode in ['vcnmodel.v1']:
            print(f"Reading model file in version v1:, with {len(df['Results']):4d} trials")
        else:
            raise ValueError(f'Unknown file mode: {filemode:s}')
        mtime = Path(filename).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
        if filemode == 'vcnmodel.v0':
            try:
                dinfo = df['Params']['runInfo']
            except:
                try:
                    dinfo = df['runInfo']
                except:
                    raise ValueError ("Cannot read the file in v0 mode?")
            run_protocol = df['Params']['runProtocol']
            if isinstance(dinfo, Params):
                dinfo = dinfo.todict()
            dur = dinfo['stimDur']
            delay = dinfo['stimDelay']
            mode = dinfo['postMode'].upper()
            ntr = len(df['Results'])
            if 'dt' in df['Params'].keys():
                self.rate = df['Params']['dt']
            else:
                self.rate = df['Params'].dt
            V = [[]]*ntr
            I = [[]]*ntr
            for i in range(len(df['Results'])):
                fk = list(df['Results'][i].keys())[0]
                dfx = df['Results'][i][fk]['monitor']
                timebase = dfx['time']
                V[i] = dfx['postsynapticV']*vscale
                I[i] = dfx['i_stim0']*iscale
        else:
            dinfo = df['runInfo']
            x = dir(dinfo)
            if 'stimVC' not in x:  # create fields for missing values from older versions of files.
                dinfo.stimVC = None
            mode = dinfo.postMode.upper()
            dur = dinfo.stimDur
            delay = dinfo.stimDelay
            mode = dinfo.postMode
            try:
                self.rate = df['Params'].dt  # old version, now separated IC and VC
            except:
                if mode == 'VC':
                    self.rate = df['Params'].dtVC
                elif mode == "CC":
                    self.rate = df['Params'].dtIC
                else:
                    raise ValueError("Cannot find rate for data mode: ", mode)

            run_protocol = dinfo.runProtocol
            if dinfo.runProtocol in ['runIV', 'initIV', 'testIV']:
                ntr = len(df['Results'])
                V = [[]]*ntr
                I = [[]]*ntr
                for ii, i in enumerate(df['Results'].keys()):
                    dfx = df['Results'][i]['monitor']
                    timebase = dfx['time']
                    V[ii] = np.array(dfx['postsynapticV'])*vscale
                    I[ii] = np.array(dfx['i_stim0'])*iscale
            elif dinfo.runProtocol in ['runVC', 'initVC', 'testVC']:
                dur = dinfo.vstimDur # msec
                delay = dinfo.vstimDelay # msec
                ntr = len(df['Results'])
                V = [[]]*ntr
                I = [[]]*ntr
                for ii, i in enumerate(df['Results'].keys()):
                    dfx = df['Results'][i]['monitor']
                    timebase = dfx['time']
                    V[ii] = np.array(dfx['postsynapticV'])*vscale
                    I[ii] = np.array(dfx['postsynapticI'])*iscale
                    
            elif dinfo.runProtocol in ['initAN', 'runANPSTH', 'runANIO', 'runANSingles']:

                # two ways data can be organized, so try both
                try:  # cnmodel_models simulations
                    ntr = len(df['Results'])
                    V = [[]]*ntr
                    I = [[]]*ntr
                    for j in list(df['Results'].keys()):
                        dfx = df['Results'][j]
                        timebase = dfx['time']
                        V[j] = np.array(dfx['somaVoltage'])*vscale
                        I[j] = np.zeros_like(V[j])
                except:  # vcnmodel simulatipns
                    ntr = len(df["Results"]['somaVoltage'])
                    V = [[]]*ntr
                    I = [[]]*ntr
                    for j in range(ntr):
                        timebase = df['Results']['time']
                        V[j] = np.array(df['Results']['somaVoltage'][j])*vscale
                        I[j] = np.zeros_like(V[j])

        V = np.array(V)
        I = np.array(I)

        if run_protocol in ['runVC', 'initVC', 'testVC']:
            self.MC.set_clamps(dmode=mode, time=timebase, data=I, cmddata=V, tstart_tdur=[delay, dur])
        else:
            self.MC.set_clamps(dmode=mode, time=timebase, data=V, cmddata=I, tstart_tdur=[delay, dur])
        self.MC.getClampData()
