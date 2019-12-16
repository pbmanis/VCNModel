import sys
import time
import argparse
import numpy as np
from collections import OrderedDict
import cnmodel.cells as cells
from pylibrary.params import Params
import json
import toml


class ModelParams():
    def __init__(self):
        self.parser = None
        
        self.cellID = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
        
        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelNameChoices = ['XM13', 'XM13_nacncoop', 'XM13_nacn', 'XM13_nabu', 'RM03', 'mGBC', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr']
        self.modelTypeChoices = ['II', 'II-I', 'I-II', 'I-c', 'I-t', 'II-o']
        self.SGCmodelChoices = ['Zilany', 'cochlea']  # cochlea is python model of Zilany data, no matlab, JIT computation; Zilany model creates matlab instance for every run.
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initandrunIV',
                                'initAN', 'runANPSTH', 'runANIO', 'runANSingles',
                                'gifnoise']
        self.soundChoices = ['tonepip', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        self.speciesChoices = ['mouse', 'guineapig']
        self.spirouChoices = ['all', 'max=mean', 'all=mean', 'removelargest', 'largestonly']
        self.ANSynapseChoices = ['simple', 'multisite']

        self.cellmap = {'bushy': cells.bushy, 'tstellate': cells.tstellate, 'dstellate': cells.dstellate}
        self.srname = ['LS', 'MS', 'HS']  # runs 0-2, not starting at 0
        
        
        # set up the Params data structure. This should hold everything that might be modified
        # or need to be known when running the model.
        #
        self.Params = OrderedDict()
        self.Params['cell'] = self.cellID
        self.Params['AMPAScale'] = 1.0 # Use the default scale for AMPAR conductances
        self.Params['ANSynapseType'] = 'simple'  # or multisite
        self.Params['ANSynapticDepression'] = 0  # depression calculation is off by default
        self.Params['initIVStateFile'] = None # 'IVneuronState_%s.dat'
        self.Params['initANStateFile'] = None # 'ANneuronState_%s.dat'
        self.Params['simulationFilename'] = None
        self.Params['shortSimulationFilename'] = None
        self.Params['simPath'] = None
        self.Params['hocfile'] = None
        self.Params['usedefaulthoc'] = False
        self.Params['cellType'] = self.cellChoices[0]
        self.Params['modelName'] = self.modelNameChoices[0]
        self.Params['modelType'] = self.modelTypeChoices[0]
        self.Params['SGCmodelType'] = self.SGCmodelChoices[0]
        self.Params['species'] = self.speciesChoices[0]
        self.Params['Ra'] = 150.  # ohm.cm
        self.Params['soma_inflation'] = 1.0  # factor to multiply soma section areas by
        self.Params['soma_autoinflate'] = False  #
        self.Params['dendrite_inflation'] = 1.0
        self.Params['dendrite_autoinflate'] = False
        self.Params['dendrite_fromsoma'] = False
        self.Params['ASA_inflation'] = 1.0
        self.Params['ASA_fromsoma'] = False
        self.Params['lambdaFreq'] = 2000.  # Hz for segment number
        self.Params['sequence'] = '' # sequence for run - may be string [start, stop, step]
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params['SRType'] = self.SRChoices[2]
        # self.Params['SR'] = self.Params['SRType']  # actually used SR: this might be cell-defined, rather than entirely specified from the command line
        self.Params['inputPattern'] = None # ID of cellinput pattern (same as cellID): for substitute input patterns.
        self.Params['synno'] = None  # for selection of single synaptic inputs
        self.Params['lastfile'] = None
        self.Params['Spirou'] = 'all'
        self.Params['runProtocol'] = self.protocolChoices[2]  # testIV is default because it is fast and should be run often
        self.Params['nReps'] = 1
        self.Params['seed'] = 100 # always the same start - avoids lots of recomutation
                                  # in production, vary this or set to None for random values
        self.Params['initialization_time'] = 50.  # nominal time to let system settle, in msec
        self.Params['run_duration'] = 0.25 # in sec
        self.Params['soundtype'] = 'SAM'  # or 'tonepip'
        self.Params['pip_duration'] = 0.1  # duration in seconds for tone pip
        self.Params['pip_start'] = [0.1]  # start (delay to start of pip)
        self.Params['pip_offduration'] = 0.05  # time after pip ends to keep running
        self.Params['Fs'] = 100e3  # cochlea/zilany model rate
        self.Params['F0'] = 16000.  # stimulus frequency
        self.Params['dB'] = 30.  # in SPL
        self.Params['RF'] = 2.5e-3  # rise-fall time
        self.Params['fmod'] = 20 # hz, modulation if SAM
        self.Params['dmod'] = 0 # percent if SAM
        self.Params['threshold'] = -35
        self.Params['signalToMasker'] = 0.
        self.Params['CMMRmode'] = 'CM'
        self.Params['all_modes'] = False
        # parameters for generating a noise signal to generating GIF model of cell
        self.Params['gif_i0'] = 0.  # base current level
        self.Params['gif_sigma'] = 0.5 # std of noise
        self.Params['gif_fmod'] = 0.2 # mod in Hz
        self.Params['gif_tau'] = 3.0 # tau, msec
        self.Params['gif_dur'] = 10. # seconds
        self.Params['gif_skew'] = 0. # a in scipy.skew

        # general control parameters
        self.Params['plotFlag'] = False
        self.Params['auto_initialize'] = False
        self.Params['nWorkers'] = 8
        self.Params['Parallel'] = True
        self.Params['verbose'] = False
        self.Params['save_all_sections'] = False
        self.Params['commandline'] = ''  # store command line on run
        self.Params['checkcommand'] = False
        self.Params['configfile'] = None
        self.Params['tagstring'] = None
        
        
        # set some defaults for modelPar
        self.Params['modelPar'] = {
            'axon': False,
            'cellClass': None,
            'dendrites': False,
            'hillock': False,
            'initialsegment': False,
            'modelName': None,
            'modelType': 'II',
            'morphology': None,
            'myelinatedaxon': False,
            'na': None,
            'name': None,
            'pumps': False,
            'soma': True,
            'species': 'mouse',
            'temperature': 34.0,
            'ttx': False,
            'unmyelinatedaxon': False,
        }

        # runinfo parameters are filled by generate_run at the initialization of a single trace run
        self.Params['runInfo'] = Params(
                              folder="Simulations",
                              fileName='Normal',
                              runName='Run',
                              manipulation="Canonical",
                              preMode="cc",
                              postMode='cc',
                              TargetCellType='', # celltype, # valid are "Bushy", "Stellate", "MNTB"
                              electrodeSection=None, #electrodeSection,
                              dendriticElectrodeSection=None, # dendriticElectrodeSection,
                              dendriticSectionDistance=100., # microns.
                              celsius=37,  # set the temperature.
                              nStim=1,
                              stimFreq=200.,  # hz
                              stimInj={'pulse': np.linspace(-1., 1.00, 11, endpoint=True)}, # iRange,  # nA, a list of test levels for current clamp
                              stimDur=100.0,  # msec
                              stimDelay=5.0,  # msec
                              stimPost=3.0,  # msec
                              vnStim=1,
                              vstimFreq=200.,  # hz
                              vstimInj=50,  # mV amplitude of step (from holding)
                              vstimDur=50.0,  # msec
                              vstimDelay=2.0,  # msec
                              vstimPost=3.0,  # msec
                              vstimHolding=-60, # holding, mV
                              gif_i0 =  self.Params['gif_i0'],
                              gif_sigma =  self.Params['gif_sigma'],
                              gif_fmod =  self.Params['gif_fmod'],
                              gif_tau =  self.Params['gif_tau'],
                              gif_dur =  self.Params['gif_dur'],
                              gif_skew =  self.Params['gif_skew'],
                              runTime = time.asctime(),  # store date and time of run
                              inFile=None,
                              # 'ANFiles/AN10000Hz.txt', # if this is not None, then we will use these spike times...
                              inFileRep=1,  # which rep to use (or array of reps)
                              spikeTimeList={},  # Dictionary of spike times
                              v_init = -61.0,  # from Rothman type II model - not appropriate in all cases
                              useSaveState = True # useSavedState,  # use the saved state.
        )


    def build_parser(self):
        """
        Set up the command line parser for running the model
        This include both verbs for running different types of protocols, and
        parameters that define the cell scaffold (structure, channels) that is being run
        The parameters in this file should also appear in the Params data structure (dict) for the class
        """
        
        parser = argparse.ArgumentParser(description='Simulate activity in a reconstructed model cell',
                        argument_default=argparse.SUPPRESS,
                        fromfile_prefix_chars='@')
        parser.add_argument(dest='cell', action='store',
                       default=None,
                       help='Select the cell (no default)')
        parser.add_argument('--type', '-T', dest='cellType', action='store',
                       default='Bushy', choices=self.cellChoices,
                       help='Define the cell type (default: Bushy)')
        parser.add_argument('--model', '-M', dest='modelName', action='store',
                       default='XM13', choices=self.modelNameChoices,
                       help='Define the model type (default: XM13)')
        parser.add_argument('--modeltype', dest='modelType', action='store',
                       default='II', choices=self.modelTypeChoices,
                       help='Define the model type (default: XM13)')
        parser.add_argument('--sgcmodel', type=str, dest='SGCmodelType', action='store',
                       default='cochlea', choices=self.SGCmodelChoices,
                       help='Define the SGC model type (default: Zilany)')
        parser.add_argument('--protocol', '-P', dest='runProtocol', action='store',
                       default='runIV', choices=self.protocolChoices,
                       help='Protocol to use for simulation (default: IV)')
        parser.add_argument('-H', '--defaulthoc', action='store_true', dest='usedefaulthoc',
                      default=True,
                      help='Use default hoc file for this cell')
        parser.add_argument('--hocfile', dest='hocfile', action='store',
                      default=None,
                      help='hoc file to use for simulation (default is the selected "cell".hoc)')
        parser.add_argument('--inputpattern', '-i', type=str, dest='inputPattern', action='store',
                      default=None,
                      help='cell input pattern to use (substitute) from cell_config.py')
        parser.add_argument('--stimulus', '-s', type=str, dest='soundtype', action='store',
                       default='tonepip', choices=self.soundChoices,
                       help='Define the stimulus type (default: tonepip)')
        parser.add_argument('--check', '-/', action='store_true', default=False, dest='checkcommand',
                       help='Only check command line for valid input; do not run model')
        parser.add_argument('--configfile', type=str, default=None, dest='configfile',
                       help='Read a formatted configuration file (JSON, TOML) for commands')

        # lowercase options are generally parameter settings:
        parser.add_argument('-d', '--dB', type=float, default=30., dest='dB',
                            help='Set sound intensity dB SPL (default 30)')
        parser.add_argument('-f', '--frequency', type=float, default=4000., dest='F0',
                help='Set tone frequency, Hz (default 4000)')
        parser.add_argument('--duration', type=float, default=0.1, dest='pip_duration',
                help='Set sound stimulus duration (sec; default 0.1)')
        parser.add_argument('-r', '--reps', type=int, default=1, dest = 'nReps',
                help='# repetitions')
        parser.add_argument('-S', '--SRType', type=str, default='HS', dest = 'SRType',
                choices=self.SRChoices,
                help=('Specify SR type (from: %s)' % self.SRChoices))
        parser.add_argument('--synapsetype', type=str, default='multisite', dest = 'ANSynapseType',
                choices=self.ANSynapseChoices,
                help=('Specify AN synapse type (from: %s)' % self.ANSynapseChoices))
        parser.add_argument('--depression', type=int, default=0, dest = 'ANSynapticDepression',
                choices=[0, 1],
                help=('Specify AN depression flag for multisite synapses (from: %s)' % str([0, 1])))
        parser.add_argument('--seed', type=int, default=0, dest = 'seed',
                help=('Specify AN seed'))

        parser.add_argument('--fmod', type=float, default=20.0, dest = 'fmod',
                help='Set SAM modulation frequency')
        parser.add_argument('--dmod', type=float, default=100.0, dest = 'dmod',
                help='Set SAM modulation depth (in percent)')
        parser.add_argument('--S2M', type=float, default=0, dest = 'signalToMasker',
                help='Signal to Masker ratio (dB)')
        parser.add_argument('--cmmrmode', type=str, default='CMR', dest = 'CMMRmode',
                choices=self.cmmrModeChoices,
                help=('Specify mode (from: %s)' % self.cmmrModeChoices))

        parser.add_argument('--Spirou', type=str, dest='Spirou', action='store', default='all',
                choices = self.spirouChoices,
                help='Specify Spirou experiment type.... ')
        # Morphology
        parser.add_argument('--soma_inflate', type=float, dest='soma_inflation', action='store', default=1.0,
                help='Specify factor by which to inflate soma AREA')
        parser.add_argument('--soma_autoinflate', action='store_true', dest='soma_autoinflate', default=False,
                help='Automatically inflate soma based on table')
        parser.add_argument('--dendrite_inflate', type=float, dest='dendrite_inflation', action='store', default=1.0,
                help='Specify factor by which to inflate total dendritic AREA')
        parser.add_argument('--dendrite_autoinflate', action='store_true', dest='dendrite_autoinflate', default=False,
                help='Automatically inflate dendrite area based on table')
        parser.add_argument('--dendrite_from_soma', action='store_true', dest='dendrite_fromsoma', default=False,
                help='Automatically inflate dendrite area based on soma inflation')
        parser.add_argument('--ASA_from_soma', action='store_true', dest='ASA_fromsoma', default=False,
                help='Automatically inflate dendrite area based on soma inflation')
            
        parser.add_argument('--tagstring', type=str, default=None, dest='tagstring',
                help="Add a tag string to the output filename to distinguish it")
        parser.add_argument('-a', '--AMPAScale', type=float, default=1.0, dest='AMPAScale',
                help='Set AMPAR conductance scale factor (default 1.0)')
        parser.add_argument('--allmodes', action="store_true", default = False, dest = 'all_modes',
                help=('Force run of all modes (CMR, CMD, REF) for stimulus configuration.'))
        parser.add_argument('--sequence', type=str, default='', dest = 'sequence',
                help=('Specify a sequence for the primary run parameters'))
        parser.add_argument('--plot',  action="store_true", default=False, dest = 'plotFlag',
                help='Plot results as they are generated - requires user intervention... ')
        parser.add_argument('--workers', type=int,  default=4, dest = 'nWorkers',
                help='Number of "workers" for parallel processing (default: 4)')
        parser.add_argument('--noparallel', action='store_false', default=True, dest='Parallel',
                help='Use parallel or not (default: True)')
        parser.add_argument('--auto', action="store_true", default=False, dest='auto_initialize',
                help='Force auto initialization if reading the state fails in initialization')
        parser.add_argument('--saveall', action="store_true", default=False, dest='save_all_sections',
                help='Save data from all sections in model')
        parser.add_argument('--verbose', action="store_true", default=False, dest='verbose',
                help='Print out extra stuff for debugging')

    # Parser arguments for gif noise generator:
    # parameters for generating a noise signal to generating GIF model of cell
        parser.add_argument('--gifi', type=float, default=0.0, dest='gif_i0',
            help='Set Noise for GIF current level (default 0 nA)')
        parser.add_argument('--gifsigma', type=float, default=0.2, dest='gif_sigma',
            help='Set Noise for GIF variance (default 0.2 nA)')
        parser.add_argument('--giffmod', type=float, default=0.2, dest='gif_fmod',
            help='Set Noise for GIF fmod (default 0.2 Hz)')
        parser.add_argument('--giftau', type=float, default=3.0, dest='gif_tau',
            help='Set Noise for GIF tau (default 0.3 ms)')
        parser.add_argument('--gifdur', type=float, default=10., dest='gif_dur',
            help='Set Noise for GIF duration (default 10 s)')
        parser.add_argument('--gifskew', type=float, default=0., dest='gif_skew',
            help='Set Noise for GIF to have skewed distribution (0 = normal)')

        # parser.add_argument('-p', '--print', action="store_true", default=False, dest = 'print_info',
        #     help='Print extra information during analysis')
        #  parser.add_argument('-l', '--list', action="store_true", default=False, dest = 'list_results',
        #     help='List results to screen')
        self.parser = parser

    def parse_parser(self):
        if self.parser is None:
            self.build_parser()
            
        
        args = self.parser.parse_args()

        if args.configfile is not None:
            config = None
            if args.configfile is not None:
                if '.json' in args.configfile:
                    # The escaping of "\t" in the config file is necesarry as
                    # otherwise Python will try to treat is as the string escape
                    # sequence for ASCII Horizontal Tab when it encounters it
                    # during json.load
                    config = json.load(open(args.configfile))
                elif '.toml' in args.configfile:
                    config = toml.load(open(args.configfile))
        
            vargs = vars(args)  # reach into the dict to change values in namespace
            for c in config:
                if c in args:
                    print('c: ', c)
                    vargs[c] = config[c]
        
        self.cmdargs = args  # just save these
        
        # now copy into the model.Params structure
        for key, value in vars(args).items():
            if key in list(self.Params.keys()):
                self.Params[key] = value
            else:
                raise ValueError(f'Parameter {key:s} is not in Params key list')

        # now add the configs if there are any to add.
        # if config is not None:
        #     for k in config.keys():
        #         if k in list(self.Params.keys()):
        #             model.Params[k] = config[k]
        #         else:
        #             raise ValueError(f'Parameter {k:s} is not in Params key list')
        print('Params: ', self.Params)
    
        # print(model.Params['cell'])
        # if model.Params['hocfile'] == None: # just use the matching hoc file
        #     model.Params['hocfile'] = model.Params['cell'] + '.hoc'
        # allargs = ''
        # parsedargs, unparsed = parser.parse_known_args()
        # for parg in vars(parsedargs):
        #     try:
        #         allargs = allargs + ('{0:s}: {1:s}\n'.format(parg, getattr(parsedargs, parg)))
        #     except:
        #         pass
        if self.Params['checkcommand']:
            self.Params['commandline'] = ' '.join(sys.argv)
            print(json.dumps(self.Params, indent=4))  # pprint doesn't work well with ordered dicts
            print(self.Params['commandline'])

    