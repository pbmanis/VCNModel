"""
multi_run does multiple runs with the models
"""
from __future__ import print_function

import os.path
import model_run as MR
import pprint
import json

pp = pprint.PrettyPrinter(indent=4, width=60)

# create an instance of the modelrunner for each cell configuration

# population the configuration parameters

# run the model with the configuration. Data is automatically stored
# set the "auto initialize" flag.

"""

Define run parameters with a dictionary:
common parameters: all keys correspond to the .Params of the model instance.
common = {'FO': F0, 'fmod': fmod, 'dmod': dmod, 'nrep': nrep, 'SRType': SR,
    'run_duration': dur, 'dB': level, 'type': celltype}

For each run, the following values are changed:
{'runID': {'cellID': name, 'inputPattern': inputpattername}}


Parameters:
        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelChoices = ['XM13', 'RM03', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr']
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initAN', 'runANPSTH', 'runANSingles']
        self.soundChoices = ['tone', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        self.speciesChoices = ['mouse', 'guineapig']

        self.srname = ['**', 'LS', 'MS', 'HS']  # runs 1-3, not starting at 0
        self.cellID = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
        self.Params = OrderedDict()

        self.Params['initIVStateFile'] = 'IVneuronState.dat'
        self.Params['initANStateFile'] = 'ANneuronState.dat'
        self.Params['infile ']= None

        self.Params['cellType'] = self.cellChoices[0]
        self.Params['model_type'] = self.modelChoices[0]
        self.Params['species'] = self.speciesChoices[0]
        self.Params['SRType'] = self.SRChoices[2]
        self.Params['SR'] = self.Params['SRType']  # actually used SR this might be cell-defined, rather than command defined
        self.Params['inputPattern'] = None # ID of cellinput pattern (same as cellID): for substitute input patterns.

        self.Params['runProtocol'] = self.protocolChoices[2]  # testIV is default because it is fast and should be run often

        self.Params['run_duration'] = 0.8 # in sec
        self.Params['soundtype'] = 'SAM'  # or 'tonepip'
        self.Params['pip_duration'] = 0.1
        self.Params['pip_start'] = [0.1]
        self.Params['Fs'] = 100e3
        self.Params['F0'] = 4000.
        self.Params['dB'] = 40.
        self.Params['RF'] = 2.5e-3
        self.Params['fmod'] = 20 # hz, modulation if SAM
        self.Params['dmod'] = 0 # percent if SAM
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params['threshold'] = -20
        self.Params['plotFlag'] = False
        self.Params['auto_initialize'] = False

"""
common = {'F0': 4000., 'fmod': 20., 'dmod': 0., 'nReps': 2, 'SRType': 'MS',
    'run_duration': 2.25, 'dB': 40., 'cellType': 'Bushy', 'pip_duration': 2.0}

# For each run, the following values are changed:

runs = {'selfinputs': {'cellID': 'VCN_c18', 'inputPattern': 'VCN_c18'}}


for runid in runs.keys():
    i1 = MR.ModelRun()
    parkeys = i1.Params.keys()
    for p in common.keys():
        if p not in parkeys:
            print('parameter %s not in parkeys list' % p)
            exit()
        else:
            i1.Params[p] = common[p] # set the parameters
        # force a couple of things
    i1.Params['initANStateFile'] = 'initAN_hoc_' + runs[runid]['cellID']+'_input_' + runs[runid]['inputPattern']
    i1.Params['inputPattern'] = runs[runid]['inputPattern']
    i1.Params['cell'] = runs[runid]['cellID']
    i1.Params['infile'] = runs[runid]['cellID'] + '.hoc'

    print(json.dumps(i1.Params, indent=4))

    stateFileExists = i1.check_for_an_statefile()
    
    if not stateFileExists:
        i1.Params['runProtocol'] = 'initAN'
        i1.run_model()
    i1.Params['runProtocol'] = 'runANPSTH'
    i1.run_model()