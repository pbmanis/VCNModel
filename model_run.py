#!/usr/bin/env python
from __future__ import print_function

__author__ = 'pbmanis'
"""
model_run.py

Run a model based on a hoc structure, decorating the structure with ion channels and synapses.

Requires:
Python 2.7.13 (anaconda distribution)
Neuron7.4 (neuron.yale.edu)
pyqtgraph (Campagnola, from github)
neuronvis (Campagnola/Manis, from github)
cnmodel (Campagnola/Manis, from github)
vcnmodel parts:
    generate_run
    cellInitialization
    analysis

cochlea # Rudniki python implementation of Zilany model
thorns  # required for cochlea


Expects the following directory structure:
(example)
VCN_Cells/   # top level for data
    cell_ID/    # VCN_c18, for example (first argument in call should be this directory name)
        MorphologyFiles/  # location for swc and hoc files
            VCN_c18_755V2.hoc (cell body scaled version)
            VCN_c18_755.hoc  (unscaled version)
            misc.swc  (various swc files that were translated to hoc files)
        InitializationFiles/  # location for Init files
            IV/  # initialization for just the IV with no synaptic input
                VCN_c18_755v2.ninit  # different base structures
                VCN_c18755.ninit
            AN/  # initialization for synaptic inputs
                (ditto)  # different base structures of input arrangements
        Simulations/
            IV/  # results from IV simulations
            AN/  # results from AN simulations


usage: model_run.py [-h] [--type {Bushy,TStellate,DStellate}]
                    [--model {XM13,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}]
                    [--sgcmodel {Zilany,cochlea}]
                    [--protocol {initIV,testIV,runIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}]
                    [--hoc HOCFILE] [--inputpattern INPUTPATTERN]
                    [--stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}]
                    [-d DB] [-f F0] [-r NREPS] [-S {LS,MS,HS,fromcell}]
                    [--fmod FMOD] [--depth DMOD] [--S2M SIGNALTOMASKER]
                    [--cmmrmode {CM,CD,REF}] [-a AMPASCALE] [--allmodes]
                    [--sequence SEQUENCE] [--plot] [--workers NWORKERS]
                    [--noparallel] [--auto] [--verbose] [--gifi GIF_I0]
                    [--gifsigma GIF_SIGMA] [--giffmod GIF_FMOD]
                    [--giftau GIF_TAU] [--gifdur GIF_DUR]
                    cell

Simulate activity in a reconstructed model cell

positional arguments:
  cell                  Select the cell (no default)

optional arguments:
  -h, --help            show this help message and exit
  --type {Bushy,TStellate,DStellate}, -T {Bushy,TStellate,DStellate}
                        Define the cell type (default: Bushy)
  --model {XM13,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}, -M {XM13,RM03,mGBC,XM13PasDend,Calyx,MNTB,L23Pyr}
                        Define the model type (default: XM13)
  --sgcmodel {Zilany,cochlea}
                        Define the SGC model type (default: Zilany)
  --protocol {initIV,testIV,runIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}, -P {initIV,testIV,runIV,initAN,runANPSTH,runANIO,runANSingles,runANOmitOne,gifnoise}
                        Protocol to use for simulation (default: IV)
  --hoc HOCFILE, -H HOCFILE
                        hoc file to use for simulation (default is the
                        selected "cell".hoc)
  --inputpattern INPUTPATTERN
                        cell input pattern to use (substitute) from
                        cell_config.py
  --stimulus {tonepip,noise,stationaryNoise,SAM,CMMR}
                        Define the stimulus type (default: tonepip)
  -d DB, --dB DB        Set sound intensity dB SPL (default 30)
  -f F0, --frequency F0
                        Set tone frequency, Hz (default 4000)
  -r NREPS, --reps NREPS
                        # repetitions
  -S {LS,MS,HS,fromcell}, --SRType {LS,MS,HS,fromcell}
                        Specify SR type (from: ['LS', 'MS', 'HS', 'fromcell'])
  --fmod FMOD           Set SAM modulation frequency
  --depth DMOD          Set SAM modulation depth (in percent)
  --S2M SIGNALTOMASKER  Signal to Masker ratio (dB)
  --cmmrmode {CM,CD,REF}
                        Specify mode (from: ['CM', 'CD', 'REF'])
  -a AMPASCALE, --AMPAScale AMPASCALE
                        Set AMPAR conductance scale factor (default 1.0)
  --allmodes            Force run of all modes (CMR, CMD, REF) for stimulus
                        configuration.
  --sequence SEQUENCE   Specify a sequence for the primary run parameters
  --plot                Plot results as they are generated - requires user
                        intervention...
  --workers NWORKERS    Number of "workers" for parallel processing (default:
                        4)
  --noparallel          Use parallel or not (default: True)
  --auto                Force auto initialization if reading the state fails
                        in initialization
  --verbose             Print out extra stuff for debugging
  --gifi GIF_I0         Set Noise for GIF current level (default 0 nA)
  --gifsigma GIF_SIGMA  Set Noise for GIF variance (default 0.2 nA)
  --giffmod GIF_FMOD    Set Noise for GIF fmod (default 0.2 Hz)
  --giftau GIF_TAU      Set Noise for GIF tau (default 0.3 ms)
  --gifdur GIF_DUR      Set Noise for GIF duration (default 10 s)

Example:
Set up initialization:
python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol initIV --model XM13
Then run model:
python model_run.py VCN_c18 --hoc gbc18_w_axon_rescaled.hoc --protocol runIV --model XM13

"""

import sys
from pathlib import Path
import os
import errno
import pickle
import time
import argparse
from collections import OrderedDict
import pprint
import json
import numpy as np
import timeit
import matplotlib
matplotlib.use('Qt4Agg')
# from neuronvis.hoc_viewer import HocViewer

import cell_config
import neuronvis.hoc_graphics as hoc_graphics
from generate_run import GenerateRun
import cellInitialization as cellInit
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.decorator import Decorator
from cnmodel import data as DATA
import pylibrary.Utility as pu  # access to a spike finder routine
import pyqtgraph as pg
import pylibrary.pyqtgraphPlotHelpers as pgh


showCell = True



class ModelRun():
    def __init__(self, args=None):
        
        # use v2 files for model with rescaled soma
        self.cellChoices = ['Bushy', 'TStellate', 'DStellate']
        self.modelNameChoices = ['XM13', 'RM03', 'mGBC', 'XM13PasDend', 'Calyx', 'MNTB', 'L23Pyr']
        self.modelTypeChoices = ['II', 'II-I', 'I-II', 'I-c', 'I-t', 'II-o']
        self.SGCmodelChoices = ['Zilany', 'cochlea']  # cochlea is python model of Zilany data, no matlab, JIT computation; Zilany model creates matlab instance for every run.
        self.cmmrModeChoices = ['CM', 'CD', 'REF']  # comodulated, codeviant, reference
        self.SRChoices = ['LS', 'MS', 'HS', 'fromcell']  # AN SR groups (assigned across all inputs)
        self.protocolChoices = ['initIV', 'testIV', 'runIV', 'initandrunIV', 
                                'initAN', 'runANPSTH', 'runANIO', 'runANSingles', 'runANOmitOne',
                                'gifnoise']
        self.soundChoices = ['tonepip', 'noise', 'stationaryNoise', 'SAM', 'CMMR']
        self.speciesChoices = ['mouse', 'guineapig']
        self.spirouChoices = ['all', 'max=mean', 'all=mean']
        self.ANSynapseChoices = ['simple', 'multisite']

        self.cellmap = {'bushy': cells.bushy, 'tstellate': cells.tstellate, 'dstellate': cells.dstellate}
        self.srname = ['LS', 'MS', 'HS']  # runs 0-2, not starting at 0
        self.cellID = None  # ID of cell (string, corresponds to directory name under VCN_Cells)
        self.setup = False  # require setup ahead of run - but setup can be done separately
        
        # The following are initial values. Some values of some of these are replaced after
        # parsing the command line.
        # See the main code section (at the end of this file).
        self.Params = OrderedDict()
        self.Params['cellID'] = self.cellID
        self.Params['AMPAScale'] = 1.0 # Use the default scale for AMPAR conductances
        self.Params['ANSynapseType'] = 'simple'  # or multisite
        self.Params['SynapticDepression'] = 0  # depression calculation is off by default
        self.Params['initIVStateFile'] = None # 'IVneuronState_%s.dat'
        self.Params['initANStateFile'] = None # 'ANneuronState_%s.dat'
        self.Params['simulationFilename'] = None
        self.Params['hocfile'] = None
        self.Params['usedefaulthoc'] = False
        self.Params['cellType'] = self.cellChoices[0]
        self.Params['modelName'] = self.modelNameChoices[0]
        self.Params['modelType'] = self.modelTypeChoices[0]
        self.Params['SGCmodelType'] = self.SGCmodelChoices[0]
        self.Params['species'] = self.speciesChoices[0]
        self.Params['Ra'] = 150.  # ohm.cm
        self.Params['soma_inflation'] = 1.0
        self.Params['soma_autoinflate'] = False
        self.Params['dendrite_inflation'] = 1.0
        self.Params['dendrite_autoinflate'] = False
        self.Params['lambdaFreq'] = 2000.  # Hz for segment number
        self.Params['sequence'] = '' # sequence for run - may be string [start, stop, step]
        # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        self.Params['SRType'] = self.SRChoices[2]
        self.Params['SR'] = self.Params['SRType']  # actually used SR: this might be cell-defined, rather than entirely specified from the command line
        self.Params['inputPattern'] = None # ID of cellinput pattern (same as cellID): for substitute input patterns.
        self.Params['spirou'] = 'all'
        self.Params['runProtocol'] = self.protocolChoices[2]  # testIV is default because it is fast and should be run often
        self.Params['nReps'] = 1
        self.Params['seed'] = 100 # always the same start - avoids lots of recomutation
                                  # in production, vary this or set to None for random values
        self.Params['run_duration'] = 0.25 # in sec
        self.Params['soundtype'] = 'SAM'  # or 'tonepip'
        self.Params['pip_duration'] = 0.1
        self.Params['pip_start'] = [0.1]
        self.Params['pip_offduration'] = 0.05
        self.Params['Fs'] = 100e3
        self.Params['F0'] = 16000.
        self.Params['dB'] = 30.
        self.Params['RF'] = 2.5e-3
        self.Params['fmod'] = 20 # hz, modulation if SAM
        self.Params['dmod'] = 0 # percent if SAM
        self.Params['threshold'] = -35
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
        self.Params['nWorkers'] = 4
        self.Params['Parallel'] = True
        self.Params['verbose'] = False
        self.Params['save_all_sections'] = False
        self.Params['commandline'] = ''  # store command line on run
        self.baseDirectory = 'VCN_Cells'
        self.morphDirectory = 'Morphology'
        self.initDirectory = 'Initialization'
        self.simDirectory = 'Simulations'


    def print_modelsetup(self):
        """
        Print out all of the parameters in the model
        """
        for p in self.Params.keys():
            print('{:18s} = {:12}'.format(p, self.Params[p]))
        print('-----------')

    def set_celltype(self, cellType):
        """
        Set the cell Type, as requested. The cell type must be in the cell choices
        
        Parameters
        ----------
        cellType : string
            The type of the cell that will be the basis for the model
        
        Returns
        -------
            Nothing
        
        """
        if cellType not in self.cellChoices:
            print('Celltype must be one of: {:s}. Got: {:s}'.format(', '.join(self.cellChoices), cellType))
            exit()
        self.Params['cellType'] = cellType

    def set_model_name(self, model_name):
        """
        Set the model Type, as requested. The model type must be in the model choices
        
        Parameters
        ----------
        model_type : string
            The type of the model that will be used (condutance settings)
        
        Returns
        -------
            Nothing
        
        """
        if model_name not in self.modelNameChoices:
            print('Model type must be one of: {:s}. Got: {:s}'.format(', '.join(self.modelNameChoices), model_type))
            exit()
        self.Params['modelName'] = model_name

    def _make_filenames(self):
        """
        Define program-wide file names (used also in cell_initialization and generate_run) one time for consistencye
        
        This routine generates two names: 
            1. The name of the initizlization file. This name includes the model name, the model type (for that name),
                soma and dendrite inflation factors if relevenat. Other parameters can be added if needed
            2. The name of the simulation data file. This name is similar to the initizlizaiton name, but may include
                information about the type of run, stimuli, etc.
        
        
        """
        self.cellID = Path(self.Params['cell']).stem # os.path.splitext(self.Params['cell'])[0]
        # pick up parameters that should be in both init and run filenames:
        # Run / init type independent parameters:
        namePars = f"{str(self.Params['modelName']):s}_{str(self.Params['modelType']):s}"
        if self.Params['soma_inflation'] != 1.0:
            namePars += f"_soma={self.Params['soma_inflation']:.3f}"
        if self.Params['dendrite_inflation'] != 1.0:
            namePars += f"_dend={self.Params['dendrite_inflation']:.3f}"
        if self.Params['runProtocol'] in ['initIV', 'initandrunIV', 'runIV']:
            if self.Params['initIVStateFile'] is None:
                fn = f"IVneuronState_{namePars:s}.dat"
                ivinitdir = Path(self.baseDirectory, self.cellID,
                                    self.initDirectory)
                self.Params['initIVStateFile'] = Path(ivinitdir, fn)
            print('IV Initialization file: ', self.Params['initIVStateFile'])
            self.mkdir_p(ivinitdir) # confirm existence of that file
        
        if self.Params['runProtocol'] in ['initandrunIV', 'runIV']:
            outPath = Path('VCN_Cells', self.cellID, self.simDirectory, 'IV')
            self.mkdir_p(outPath) # confirm that output path exists
            self.Params['simulationFilename'] = Path(outPath, f"{self.cellID:s}_pulse_{namePars:s}_monitor.p")
            print('Simulation filename: ', self.Params['simulationFilename'])
               
        
        if self.Params['runProtocol'].startswith('runAN'):
            if self.Params['initANStateFile'] is None:
                fn = f"ANneuronState_{namePars:s}_{self.Params['ANSynapseType']:s}.dat"
                aninitdir= Path(self.baseDirectory, self.cellID,
                                    self.initDirectory)
                self.Params['initANStateFile'] = Path(aninitdir, fn)
                self.mkdir_p(aninitdir) # confirm existence of that file  
                print('IV Initialization file: ', self.Params['initANStateFile'])
                   
            outPath = Path('VCN_Cells', self.cellID, self.simDirectory, 'AN')
            self.mkdir_p(outPath) # confirm that output path exists
            ID = self.cellID
            if self.Params['inputPattern'] is not None:
                ID += '_%s' % self.Params['inputPattern']
            addarg = namePars
            if self.Params['spirou'] == 'all':
                addarg += '_all'
            elif self.Params['spirou'] == 'max=mean':
                addarg = '_mean'
            elif self.Params['spirou'] == 'all=mean':
                addarg = '_allmean'
            
            print('soundtype: ', self.Params['soundtype'])
            fn = f"AN_Result_{self.cellID:s}_{addarg:s}_{self.Params['ANSynapseType']:s}"
            fn += f"_{self.Params['nReps']:03d}_{self.Params['soundtype']:s}"
            fn += f"_{int(self.Params['dB']):03d}dB_{self.Params['F0']:06.1f}"
            if self.Params['soundtype'] in ['SAM', 'sam']:
                fn += f"_{self.Params['fmod']:03.1f}_{int(self.Params['dmod']):03d}"
                fn += f"_{self.Params['SR']:2s}.p"
                ofile = Path(outPath, fn)
            else:
                fn += f"_{self.Params['SR']:2s}.p"
                ofile = Path(outPath, fn)
            self.Params['simulationFilename'] = ofile

        
        
    def set_spontaneousrate(self, spont_rate_type):
        """
        Set the SR, overriding SR in the cell_config file. The SR type must be in the SR choices
        
        Parameters
        ----------
        spont_rate_type : string
            The SR type that will be used for AN fibers.
            1 = Low, 2 = medium, 3  = high (following Zilany et al)
        Returns
        -------
            Nothing
        
        """
        if self.spont_rate_type not in self.spont_rate_choices:
            print('SR type must be one of: {:s}. Got: {:s} '.format(', '.join(self.spont_rate_choices), spont_rate_type))
            exit()
        self.Params['SRType'] = spont_rate_type
    
    def set_starttime(self, starttime):
        """
        store the start time... 
        """
        self.Params['StartTime'] = starttime

    def mkdir_p(self, path):
        try:
            Path.mkdir(path, exist_ok=True)
        except:
            raise FileNotFoundError(f"Cannot create path: {str(path):s}")

    def setup_model(self, par_map=None):
        """
        Main entry routine for running all models
        
        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of parameters, passed to models that are run (not used).
        
        Returns
        -------
            Nothing
        
        """
        # if self.Params['ANSynapseType'] == 'simple':
        #     raise ValueError('model_run:: Simple AN synapses are not fully implemented in this version')

        if self.Params['verbose']:
            print('run_model entry')
        if par_map is not None and 'id' in par_map.keys():
            self.idnum = par_map['id']
        else:
            self.idnum = 9999
        self.cellID = Path(self.Params['cell']).stem # os.path.splitext(self.Params['cell'])[0]
        
        print('Morphology directory: ', self.morphDirectory)
        if self.Params['usedefaulthoc']:
            self.Params['hocfile'] = self.Params['cell'] + '.hoc'
        print('Hoc (structure) file: ', self.Params['hocfile'])
        print('Base directory: ', self.baseDirectory)
        print('params: ', self.Params)
        filename = Path(self.baseDirectory, self.cellID, self.morphDirectory, self.Params['hocfile'])
        
        # instantiate cells
        if self.Params['cellType'] in ['Bushy', 'bushy']:
            print('Creating a bushy cell (run_model) ')
            from cnmodel import data

            changes = data.add_table_data('XM13_channels', row_key='field', col_key='model_type', 
                           species='mouse', data=u"""

            This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
            for different cell types based on the Xie and Manis 2013 models for mouse.

            The REFERENCE values are applied to "point" models, and to the soma of
            compartmental models.
            The names of the mechanisms must match a channel mechanism (Neuron .mod files)
            and the following _(gbar, vshift, etc) must match an attribute of that channel
            that can be accessed.

            -----------------------------------------------------------------------------------------------------------------------------------
                           II             II-I           I-c           I-II          I-t       
                                                                                   
            nav11_gbar     1000.  [1]     1000.  [1]     3000.  [1]    1000.  [2]    3000.  [1] 
            kht_gbar       58.0   [3]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
            klt_gbar       80.0   [1]     14.0   [1]     0.0    [1]    20.0   [2]    0.0    [1] 
            ka_gbar        0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0    [1] 
            ihvcn_gbar     30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
            leak_gbar      2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
            leak_erev      -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
            na_type        nav11  [1]     nav11  [1]     nav11  [1]    nav11  [1]    nav11  [1] 
            ih_type        ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
            soma_Cap       26.0   [1]     26.0   [1]     25.0   [1]    26.0   [2]    25.0   [1] 
            nav11_vshift   10.    [1]     4.3    [1]     4.3    [1]    4.3    [1]    4.3    [1]
            e_k            -84    [1]     -84    [1]     -84    [1]    -84    [2]    -84    [1] 
            e_na           50.    [1]     50.    [1]     50.    [1]    50.    [2]    50.    [1] 
            ih_eh          -43    [1]     -43    [1]     -43    [1]    -43    [2]    -43    [1] 

            -----------------------------------------------------------------------------------------------------------------------------------

            [1] Uses channels from Rothman and Manis, 2003
                Conductances are for Mouse bushy cells
                Xie and Manis, 2013
                Age "adult", Temperature=34C
                Units are nS.

            [2] Rothman and Manis, 2003, model I-II
                Some low-voltage K current, based on observations of
                a single spike near threshold and regular firing for higher
                currents (Xie and Manis, 2017)
                
            [3] Increased DR to force AP repolarization to be faster


            """)
            changes = data.add_table_data('XM13_channels_compartments', row_key='parameter', col_key='compartment',
                    species='mouse', model_type='II', data=u"""

            This table describes the ion channel densities relative to somatic densities,
            e.g., relative to REFERENCE densities in the table XM13_channels.
            and voltage shifts, for different compartments of the specified neuron,
            Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse
            (data table: XM13_channels).

            NOTE: unmyelinatedaxon and initialsegment are equivalent in George's models, but only "unmyelinatedaxon" is actually used.
            ------------------------------------------------------------------------------------------------------------------------------------------------------------------
                           axon       unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite

            nav11_gbar     3.0 [1]    15. [1]              0.0 [1]            0.0 [1]           9.0 [1]     2.5 [1]     0.5 [1]          0.25 [1]           0.25 [1]
            kht_gbar       1.0 [1]    1.0 [1]              0.01 [1]           2.0 [1]           1.0 [1]     1.0 [1]     1.0 [1]          0.5 [1]            0.25 [1]
            klt_gbar       1.0 [1]    1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]
            ihvcn_gbar     0.0 [1]    0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]
            leak_gbar      1.0 [1]    0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]
            leak_erev      -65. [1]   -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]
            nav11_vshift   4.3  [1]   4.3 [1]              0.0 [1]            4.3 [1]           4.3 [1]     4.3 [1]     0.0  [1]         0.0  [1]            0.0 [1]
            na_type        nav11      nav11                nav11              nav11             nav11       nav11       nav11            nav11              nav11
            ih_type        ihvcn      ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn
            -------------------------------------------------------------------------------------------------------------------------------------------------------------------

            [1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.


            """)
            data.report_changes(changes)

            self.post_cell = cells.Bushy.create(morphology=str(filename), decorator=Decorator,
                    species=self.Params['species'], 
                    modelName=self.Params['modelName'], modelType=self.Params['modelType'])
        elif self.Params['cellType'] in ['tstellate', 'TStellate']:
            print('Creating a t-stellate cell (run_model) ')
            self.post_cell = cells.TStellate.create(morphology=str(filename), decorator=Decorator,
                    species=self.Params['species'],
                    modelType=self.Params['modelName'], )
        elif self.Params['cellType'] in ['dstellate', 'DStellate']:
            print('Creating a D-stellate cell (run_model)')
            self.post_cell = cells.DStellate.create(morphology=str(filename), decorator=Decorator,
                    species=self.Params['species'],
                    modelType=self.Params['modelName'], )
        else:
            raise ValueError(f"cell type {self.Params['cellType']:s} not implemented")
        
        
        # Set up run parameters
        print('Requested temperature (deg C): ', self.post_cell.status['temperature'])
        self.post_cell.hr.h.celsius = self.post_cell.status['temperature']  # this is set by prepareRun in generateRun. Only place it should be changed
        self.post_cell.hr.h.Ra = self.Params['Ra']
        print('Ra (ohm.cm) = {:8.1f}'.format(self.post_cell.hr.h.Ra))
        print(f'Specified Temperature = {self.post_cell.hr.h.celsius:8.1f} degC ')
        
        if self.Params['soma_autoinflate']: # get values and inflate soma automatically to match mesh
            cconfig = cell_config.CellConfig()
            inflateratio = cconfig.get_soma_ratio(self.cellID)
            self.Params['soma_inflation'] = inflateratio
        
        if self.Params['soma_inflation'] != 1.0:
            print('!!!!!   Inflating soma')
            soma_group = self.post_cell.hr.sec_groups['soma']
            print(' Section in soma group: ', soma_group)
            rtau = self.post_cell.compute_rmrintau(auto_initialize=True, vrange=[-80., -60.])
            print(f"     Original Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}")
            origdiam = {}
            origarea = self.post_cell.somaarea
            self.post_cell.computeAreas()
            a1 = np.sum(list(self.post_cell.areaMap['soma'].values()))
            print('     Original Soma area: ', a1)
            for section in list(soma_group):
                secobj = self.post_cell.hr.get_section(section)
                origdiam[section] = secobj.diam
                print('      diam orig: ', origdiam[section])
                secobj.diam = origdiam[section] * self.Params['soma_inflation']
                print('      diam new:  ', self.post_cell.hr.get_section(section).diam)
                print('      ratio:  ', self.post_cell.hr.get_section(section).diam/origdiam[section])
            self.post_cell.computeAreas()
            a2 = np.sum(list(self.post_cell.areaMap['soma'].values()))
            print('     Revised Soma area: ', a2, '  area ratio: ', a2/a1)

            rtau = self.post_cell.compute_rmrintau(auto_initialize=True, vrange=[-80., -60.])
            print(f"     New Rin: {rtau['Rin']:.2f}, tau: {rtau['tau']*1e3:.2f}, RMP: {rtau['v']:.2f}")

        
        for group in self.post_cell.hr.sec_groups.keys():
            g = self.post_cell.hr.sec_groups[group]
            for section in list(g):
                self.post_cell.hr.get_section(section).Ra = self.post_cell.hr.h.Ra
                if self.Params['verbose']:
                    print('Section: ', section)
                    print('Ra: ', self.post_cell.hr.get_section(section).Ra)
        
        electrode_section = list(self.post_cell.hr.sec_groups['soma'])[0]
        self.electrode_site = self.post_cell.hr.get_section(electrode_section)
        self.electrodeSection = 'soma'
        self.hg = hoc_graphics
        self.get_hoc_file(self.post_cell.hr)

        try:
            dendritic_electrode = list(self.post_cell.hr.sec_groups['dendrite'])[0]
            self.dendriticelectrode_site = self.post_cell.hr.get_section(dendritic_electrode)
            self.dendriticElectrodeSection = 'dendrite'
        except:
            dendritic_electrode = list(self.post_cell.hr.sec_groups['soma'])[0]
            self.dendriticelectrode_site = self.post_cell.hr.get_section(dendritic_electrode)
            self.dendriticElectrodeSection = 'soma'
        
        if self.Params['verbose']:
            if par_map is not None:
                print('Listing par_map (run_model): ', par_map)
            self.post_cell.hr.h.topology()
       
        self.post_cell.set_d_lambda(freq=self.Params['lambdaFreq'])
        
        # handle the following protocols:
        # ['initIV', 'initAN', 'runIV', 'run', 'runANSingles', 'gifnoise']

        self._make_filenames()  # make filenames AFTER all manipulations of the cell
        
        self.setup = True

    def run_model(self, par_map=None):
        if not self.setup:
            self.setup_model(par_mnap=par_map)
        if self.Params['runProtocol'] in ['runANPSTH', 'runANSingles']:
            self.Params['run_duration'] = self.Params['pip_duration'] + np.sum(self.Params['pip_start']) + self.Params['pip_offduration']
            
        if self.Params['runProtocol'] in ['initIV', 'initandrunIV']:
            if self.Params['verbose']:
                print('run_model: protocol is initIV')
            # ivinitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(self.post_cell, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection,
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=self.post_cell.i_test_range,
                             plotting = self.Params['plotFlag'],
                             params = self.Params)
            cellInit.get_initial_condition_state(self.post_cell, tdur=500.,
               filename=self.Params['initIVStateFile'], electrode_site=self.electrode_site)
            print(f'Ran to get initial state for {self.post_cell.hr.h.t:.1f} msec')

        if self.Params['runProtocol'] in ['runIV', 'initandrunIV']:
            if self.Params['verbose']:
                print('Run IV')
            self.iv_run()
        
        if self.Params['runProtocol'] == 'testIV':
            if self.Params['verbose']:
                print( 'test_init')
            # ivinitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initIVStateFile'])
            self.R = GenerateRun(self.post_cell, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection,
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=self.post_cell.irange,
                             plotting = self.Params['plotFlag'], 
                             params = self.params)
            cellInit.test_initial_conditions(self.post_cell, filename=self.Params['initIVStateFile'],
                electrode_site=self.electrode_site)
            self.R.testRun(initfile=self.Params['initIVStateFile'])
            return  # that is ALL, never make testIV/init and then keep running.
        
        if self.Params['runProtocol'] == 'runANPSTH':
            if self.Params['verbose']:
                print('ANPSTH')
            self.an_run(self.post_cell)
        
        if self.Params['runProtocol'] == 'initAN':
            if self.Params['verbose']:
                print('Init AN')
            self.an_run(self.post_cell, make_an_intial_conditions=True)

        if self.Params['runProtocol'] == 'runANIO':
            if self.Params['verbose']:
                print('Run AN IO')
            self.an_run_IO(self.post_cell)

        if self.Params['runProtocol'] == 'runANSingles':
            if self.Params['verbose']:
                print('ANSingles')
            self.an_run_singles(self.post_cell)

        if self.Params['runProtocol'] == 'runANOmitOne':
            if self.Params['verbose']:
                print('ANOmitOne')
            self.an_run_singles(self.post_cell, exclude=True)

        if self.Params['runProtocol'] == 'gifnoise':
            self.noise_run()

        # if showCell:
        #     print(dir(self.post_cell))
        #     self.render = HocViewer(self.post_cell.hr)
        #     cylinder=self.render.draw_cylinders()
        #     cylinder.set_group_colors(self.section_colors, alpha=0.8, mechanism=['nav11', 'gbar'])
        #     pg.show()

    
    def iv_run(self, par_map=None):
        """
        Main entry routine for running all IV (current-voltage relationships with somatic electrode)
        
        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).
        
        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times
        
        """
        print('iv_run: starting')
        start_time = timeit.default_timer()
        if self.Params['verbose']:
            print('iv_run: calling generateRun', self.post_cell.i_test_range)
        # parse i_test_range and pass it here
        if self.Params['sequence'] is '':
            if isinstance(self.post_cell.i_test_range, dict):  # not defined, use default for the cell type
                iinjValues = self.post_cell.i_test_range
            else:  # if it was a list, try to put the list into the pulse range
                iinjValues = {'pulse': self.post_cell.i_test_range}
        else:
            iinjValues = {'pulse': eval(self.Params['sequence'])}
        self.R = GenerateRun(self.post_cell, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection,
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             iRange=iinjValues,
                             plotting = self.Params['plotFlag'],
                             params=self.Params,
                             )
        self.R.runInfo.folder = Path('VCN_Cells', self.cellID, self.simDirectory, 'IV')
        if self.Params['verbose']:
            print('iv_run: calling do_run')
        nworkers = self.Params['nWorkers']
#        print(self.Params['Parallel'])
        if self.Params['Parallel'] == False:
            nworkers = 1
#        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(self.Params['hocfile'], parMap=iinjValues, save='monitor', 
            restore_from_file=True, initfile=self.Params['initIVStateFile'],
            workers=nworkers)
        if self.Params['verbose']:
            print( '   iv_run: do_run completed')
        elapsed = timeit.default_timer() - start_time
        print(f'   iv_rin: Elapsed time: {elapsed:2f} seconds')
        isteps = self.R.IVResult['I']
        if self.Params['verbose']:
            for k, i in enumerate(self.R.IVResult['tauih'].keys()):
                print( '   ih: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['tauih'][i]['tau']))
                print('           dV : %f' % self.R.IVResult['tauih'][i]['a'])
            for k, i in enumerate(self.R.IVResult['taus'].keys()):
                print('   i: %3d (%6.1fnA) tau: %f' % (i, isteps[k], self.R.IVResult['taus'][i]['tau']))
                print( '          dV : %f' % (self.R.IVResult['taus'][i]['a']))
        
        #print('   Nspike, Ispike: ', self.R.IVResult['Nspike'], self.R.IVResult['Ispike'])
        print('   N spikes:   {0:d}'.format(int(np.sum(self.R.IVResult['Nspike']))))
        print('   Rinss:      {0:.1f} Mohm'.format(self.R.IVResult['Rinss']))
        print('   Tau(mean):  {0:.3f} ms'.format(np.mean([self.R.IVResult['taus'][i]['tau'] for i in range(len(self.R.IVResult['taus']))])))
        print('   Vm:         {0:.1f} mV'.format(np.mean(self.R.IVResult['Vm'])))
        if len(self.R.IVResult['taus'].keys()) == 0:
            taum_mean = 0.
            tauih_mean = 0.
        else:
            taum_mean = np.mean([self.R.IVResult['taus'][i]['tau'] for k, i in
                enumerate(self.R.IVResult['taus'].keys())])
            tauih_mean = np.mean([self.R.IVResult['tauih'][i]['tau'] for k, i in
                enumerate(self.R.IVResult['tauih'].keys())])
        # construct dictionary for return results:
        self.IVSummary = {'basefile': self.R.basename,
                          'par_map': par_map, 'ID': self.idnum,
                          'sequence' : iinjValues,
                          'Vm': np.mean(self.R.IVResult['Vm']),
                          'Rin': self.R.IVResult['Rinss'],
                          'taum': taum_mean, 'tauih': tauih_mean,
                          'spikes': {'i': self.R.IVResult['Ispike'], 'n': self.R.IVResult['Nspike']},
                          }
        return self.IVSummary
    
    def noise_run(self, par_map={}):
        """
        Main entry routine for running noise into cell current injection (for generating GIF models)
        
        Parameters
        ----------
        par_map : dict (default: empty)
            A dictionary of paramters, passed to models that are run (not used).
        
        Returns
        -------
            summary : dict
                A summary of the results, including the file, par_map, resting input resistance,
                time constant, and spike times
        
        """
        print('noise_run: starting')
        # parse i_test_range and pass it here
        self.R = GenerateRun(self.post_cell, idnum=self.idnum, celltype=self.Params['cellType'],
                             starttime=None,
                             electrodeSection=self.electrodeSection,
                             dendriticElectrodeSection=self.dendriticElectrodeSection,
                             stimtype='gifnoise',
                             plotting = self.Params['plotFlag'],
                             params=self.Params)
        # ivinitfile = Path(self.baseDirectory, self.cellID,
        #                         self.initDirectory, self.Params['initIVStateFile'])
        self.R.runInfo.folder = Path('VCN_Cells', self.cellID, self.simDirectory, 'Noise')
        if self.Params['verbose']:
            print('noise_run: calling do_run')
        nworkers = self.Params['nWorkers']
#        print(self.Params['Parallel'])
        if self.Params['Parallel'] == False:
            nworkers = 1
#        print('Number of workers available on this machine: ', nworkers)
        self.R.doRun(self.Params['hocfile'], parMap=par_map, save='monitor', restore_from_file=True,
            initfile=self.Params['initIVStateFile'],
            workers=nworkers)
        if self.Params['verbose']:
            print( '   noise_run: do_run completed')


    def check_for_an_statefile(self):
        print('State file name: ', self.Params['initANStateFile'])
        print('Cell: ', self.Params['cell'])
        print('Base Directory: ', self.baseDirectory)
        print('Initialization Directory: ', self.initDirectory)
        statefile = Path(self.baseDirectory, self.Params['cell'],
                            self.initDirectory, self.Params['initANStateFile'])
        return(statefile.is_file())

    def compute_seeds(self, nReps, synapseConfig):
        """
        Generate the seeds for the AN runs
        
        If Params['seed'] is None, we use the randomized seed method - every run is totally
        independent.
        If not, then we generate a consequetive sequence of seeds (so it is controlled and
        reusable) offset by the starting number specified as the seed
        
        Parameters
        ----------
        nReps: int (no default)
            Number of repetions of the stimulus that are needed
        
        synpaseConfig : list (no default)
            The length of the synpaseConfig list indicates how many synapses are on the
            cell, and each synapse is assigned a unique seed.
    
        Returns
        -------
        seeds : numpy array 
            A 2-D numpy array (nReps, len(synapseConfig)) containing the seeds for each run and
            each synapse.
        """
        if self.Params['seed'] is None:  # every run is randomized
            seeds = np.random.randint(32678, size=(nReps, len(synapseConfig)))    
        else:
            startingseed = self.Params['seed']
            seeds = np.arange(0, nReps*len(synapseConfig)).reshape((nReps, len(synapseConfig)))
            seeds = seeds + startingseed
        print('AN Seeds: ', seeds)

        return seeds
    
    def get_synapses(self, synapses):
        gMax = np.zeros(len(synapses))  # total g for each synapse
        nSyn = np.zeros(len(synapses))  # # sites each synapse
        ngMax = np.zeros(len(synapses))
        print('# syn: ', len(synapses))
        if self.Params['ANSynapseType'] == 'simple':
            for i, s in enumerate(synapses):
                gMax[i] = gMax[i] + float(s.psd.terminal.netcon.weight[0])
                nSyn[i] = nSyn[i] + 1

        else:
            for i, s in enumerate(synapses):
                nSyn[i] = len(s.psd.ampa_psd)
                for p in s.psd.ampa_psd:
                    gMax[i] = gMax[i] + p.gmax
                for p in s.psd.nmda_psd:
                    ngMax [i] = ngMax[i] + p.gmax
                print('i nsyn gmax nmdagmax  : ', i, nSyn[i], gMax[i], ngMax[i])
        return (gMax, ngMax, nSyn)
        
    def set_synapse(self, synapse, gampa, gnmda):
        totalg = 0.
        if self.Params['ANSynapseType'] == 'simple':
            synapse.psd.terminal.netcon.weight[0] = gampa
            return
        # multisite:
        for p in synapse.psd.ampa_psd:
            p.gmax = gampa
            totalg = totalg + p.gmax
        totaln = 0.
        for p in synapse.psd.nmda_psd:
            p.gmax = gnmda
            totaln = totaln + p.gmax
        print('total ampa, nmda: ', totalg, totaln)
        
    def an_run(self, post_cell, verify=False, make_an_intial_conditions=False):
        """
        Establish AN inputs to soma, and run the model.
        Requires a synapseConfig list of dicts from cell_config.makeDict()
        each list element represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high], 
                    type (e or i), segment and segment location)
                    Note if segment is None, then the synapse is assigned to the soma.
        
        Parameters
        ----------
        post_cell : Cell object
            Provides access to cell class, as well as neuron and hoc file information
        
        verify : boolean (default: False)
            Flag to control printing of various intermediate results
        
        make_ANIntitialConditions : bool (default : False)
            Flag to control whether the initial conditions need to be recomputed and stored.
        
        Returns
        -------
            Nothing

        """
        print('\n*** an_run\n')
               
        if self.Params['inputPattern'] is not None:
            fromdict = self.Params['inputPattern']
            print('Cell id: %s  using input pattern: %s' % (self.cellID, fromdict))
        else:
            fromdict = self.cellID
        cconfig = cell_config.CellConfig()
        synapseConfig, celltype = cconfig.makeDict(fromdict)
        self.start_time = time.time()
        # compute delays in a simple manner
        # assumption 3 meters/second conduction time
        # delay then is dist/(3 m/s), or 0.001 ms/um of length
        # for i, s in enumerate(synapseConfig):
        #     s[1] = s[3]*0.001/3.0
        #     print 'delay for input %d is %8.4f msec' % (i, s[1])
        
        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV
        
        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(post_cell, synapseConfig, celltype, stimInfo)
        
        # see if we need to save the cell state now.
        if make_an_intial_conditions:
#            print('getting initial conditions for AN')
            # aninitfile = Path(self.baseDirectory, self.cellID,
            #                     self.initDirectory, self.Params['initANStateFile'])
            cellInit.get_initial_condition_state(post_cell, tdur=500.,
                filename=self.Params['initANStateFile'], electrode_site=self.electrode_site, reinit=self.Params['auto_initialize'])
            cellInit.test_initial_conditions(post_cell, filename=self.Params['initANStateFile'],
                electrode_site=self.electrode_site)
            # return
        
        seeds = self.compute_seeds(nReps, synapseConfig)
        stimInfo['seeds'] = seeds  # keep the seed values too.
        # spikeTimes = {}
        # inputSpikeTimes = {}
        # somaVoltage = {}
        # dendriteVoltage = {}
        celltime = []
        # stim = {}
        # stimWaveform = {}
        allDendriteVoltages = {}  # Saving of all dendrite voltages is controlled by --saveall flag
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.

        nWorkers = self.Params['nWorkers']
        TASKS = [s for s in range(nReps)]
        tresults = [None]*len(TASKS)
        result = {}
        
        # manipulate syanptic strengths to test hypotheses here.
        allsf = self.Params['AMPAScale']
        gMax, ngMax, nSyn = self.get_synapses(synapse)
        for i in range(len(nSyn)):
            gMax[i] = gMax[i] * allsf
            ngMax[i] = ngMax[i] * allsf
            self.set_synapse(synapse[i], gMax[i]/nSyn[i], ngMax[i]/nSyn[i])        
        
        if self.Params['spirou'] in ['max=mean']:
            print('setting largest to mean of all inputs')
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print(self.Params['spirou'])
            print('AMPA: mean, gmax, imax: ', meangmax, gMax, imaxgmax)
            print('NMDA: mean, gnmax: ', meangnmax, gnMax)
        
            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            self.set_synapse(synapse[imaxgmax], meangmax/nSyn[imaxgmax], meangnmax/nSyn[imaxgmax])
            print('revised values:')
            self.get_synapses(synapse)
            # for i, s in enumerate(synapse):
            #     if i == imaxgmax:
            #         p.gmax = gMax[i]/nSyn[i]  # except the chosen one
            # for p in s.psd.ampa_psd:
            #     gMax[i] = p.gmax
            #     print('revised gmax i : ', i, gMax[i])
        if self.Params['spirou'] in ['all=mean']:
            print('setting ALL to the mean of all inputs (no variance)')
            gMax, gnMax, nSyn = self.get_synapses(synapse)

            meangmax = np.mean(gMax)  # mean conductance each synapse
            meangnmax = np.mean(gnMax)
            imaxgmax = np.argmax(gMax)  # synapse with largest conductance
            print(self.Params['spirou'])
            print('AMPA: mean, gmax, imax: ', meangmax, gMax, imaxgmax)
            print('NMDA: mean, gnmax: ', meangnmax, gnMax)
        
            gMax[imaxgmax] = meangmax
            gnMax[imaxgmax] = meangnmax
            for i in range(len(synapse)):
                self.set_synapse(synapse[i], meangmax/nSyn[i], meangnmax/nSyn[i])
            print('revised values:')
            self.get_synapses(synapse)
            # for i, s in enumerate(synapse):
            #     if i == imaxgmax:
            #         p.gmax = gMax[i]/nSyn[i]  # except the chosen one
            # for p in s.psd.ampa_psd:
            #     gMax[i] = p.gmax
            #     print('revised gmax i : ', i, gMax[i])

        # run using pyqtgraph's parallel support
        if self.Params['Parallel']:
            with pg.multiprocess.Parallelize(enumerate(TASKS), results=tresults,
                workers=nWorkers) as tasker:
                for j, x in tasker:
                    tresults = self.single_an_run(post_cell, j, synapseConfig,
                        stimInfo, seeds, preCell, self.an_setup_time)
                    tasker.results[j] = tresults
            # retreive the data
            for j, N in enumerate(range(nReps)):

                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'],
                        threshold, t0=0., t1=stimInfo['run_duration']*1000., dt=1.0, mode='peak')
                spikeTimes = self.clean_spiketimes(spikeTimes)
                inputSpikeTimes = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage = np.array(tresults[j]['Vsoma'])
                dendriteVoltage = np.array(tresults[j]['Vdend'])
                stimWaveform = np.array(tresults[j]['stimWaveform'])
                stimTimebase = np.array(tresults[j]['stimTimebase'])
                stim = np.array(tresults[j]['stim']) # save the stimulus
                result[N] = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes,
                    'time': np.array(celltime[0]),'somaVoltage': somaVoltage, 
                    'dendriteVoltage': dendriteVoltage, 'allDendriteVoltages': allDendriteVoltages,
                    'stimWaveform': stimWaveform, 'stimTimebase': stimTimebase,
                     }
                if self.Params['save_all_sections']:  # just save soma sections        for section in list(g):
                    allDendriteVoltages[N] = tresults[j]['allsecVec']

        else:
            # Non parallelized version (with --noparallel flag - useful for debugging):
            for j, N in enumerate(range(nReps)):
                print('Repetition %d' % N)

                tresults[j] = self.single_an_run(post_cell, j, synapseConfig, 
                        stimInfo,  seeds, preCell, self.an_setup_time)

                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'],
                        threshold, t0=0., t1=stimInfo['run_duration']*1000., dt=1.0, mode='peak')
                spikeTimes = self.clean_spiketimes(spikeTimes)
                inputSpikeTimes = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage = np.array(tresults[j]['Vsoma'])
                dendriteVoltage = np.array(tresults[j]['Vdend'])
                stimWaveform = np.array(tresults[j]['stimWaveform'])
                stimTimebase = np.array(tresults[j]['stimTimebase'])
                stim = np.array(tresults[j]['stim']) # save the stimulus
                result[N] = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes,
                    'time': np.array(celltime[0]),'somaVoltage': somaVoltage, 
                    'dendriteVoltage': dendriteVoltage, 'allDendriteVoltages': allDendriteVoltages,
                    'stimWaveform': stimWaveform, 'stimTimebase': stimTimebase,
                     }
                if self.Params['save_all_sections']:  # save data for all sections
                    allDendriteVoltages[N] = tresults[j]['allsecVec']
        
        total_elapsed_time = time.time() - self.start_time
#        total_run_time = time.time() - run_time
        print("Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(total_elapsed_time/60., total_elapsed_time))
        print("Total Setup Time = {:8.2f} min ({:8.0f}s)".format(self.setup_time/60., self.setup_time))
        print("Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(self.an_setup_time/60., self.an_setup_time))
        print("Total Neuron Run Time = %{:8.2f} min ({:8.0f}s)".format(self.nrn_run_time/60., self.nrn_run_time))
        
        if self.Params['save_all_sections']: 
            for n in range(len(allDendriteVoltages)):
                for s in allDendriteVoltages[n].keys():
                    allDendriteVoltages[n][s] = np.array(allDendriteVoltages[n][s])
        self.analysis_filewriter(self.Params['cell'], result, tag='delays')
        if self.Params['plotFlag']:
            print('plotting')
            self.plot_an(celltime, result)


    def an_run_singles(self, post_cell, exclude=False, verify=False):
        """
        Establish AN inputs to soma, and run the model.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=high, 3=high])
        This routine is special - it runs nReps for each synapse, turning off all of the other synapses
        by setting the synaptic conductance to 0 (to avoid upsetting initialization)
        if the "Exclude" flag is set, it turns OFF each synapse, leaving the others running.
        The output file either says "Syn", or "ExclSyn" 
        
        Parameters
        ----------
        post_cell : cnmodel cell object
            Access to neuron and file information
        
        exclude : boolean (default: False)
            Set false to do one synapse at a time.
            Set true to do all BUT one synapse at a time.
        
        verify : boolean (default: False)
            Flag to control printing of various intermediate results
        
        
        Returns
        -------
            Nothing
        
        """
        print('\n***singles\n')        
        self.start_time = time.time()
        synapseConfig, celltype = cell_config.makeDict(self.cellID)
        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV
        
        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(post_cell, synapseConfig, celltype, stimInfo)

        nSyns = len(preCell)
        seeds = self.compute_seeds(nReps, synapseConfig)

        stimInfo['seeds'] = seeds  # keep the seed values too.
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = self.Params['Parallel']
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.
        # get the gMax's
        gMax = np.zeros(len(synapse))
        for i, s in enumerate(synapse):
            for p in s.psd.ampa_psd:
                gMax[i] = p.gmax
        
        for k in range(nSyns):
            # only enable gsyn on the selected input
            if exclude:
                tagname = 'ExcludeSyn%03d'
                for i, s in enumerate(synapse):
                    for p in s.psd.ampa_psd:
                        if i == k:
                            p.gmax = 0.  # disable the one
                        else:
                            p.gmax = gMax[i]  # leaving the others set
            else:
                tagname = 'Syn%03d'
                for i, s in enumerate(synapse):
                    for p in s.psd.ampa_psd:
                        if i != k:
                            p.gmax = 0.  # disable all
                        else:
                            p.gmax = gMax[i]  # except the chosen one
            print('Syn #%d   synapse gMax: %f ' % (k, gMax[k]))
            print('tag: %s' % (tagname % k))
            
            tresults = [None]*nReps
            
            if parallel and self.Params['nWorkers'] > 1:
                nWorkers = self.Params['nWorkers']
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with pg.multiprocess.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
                    for j, x in tasker:
                        tresults = self.single_an_run(post_cell, j, synapseConfig,
                            stimInfo, seeds, preCell, self.an_setup_time)
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    tresults[j] = self.single_an_run(post_cell, j, synapseConfig,
                            stimInfo, seeds, preCell, self.an_setup_time)
            for j, N in enumerate(range(nReps)):
                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'], threshold, t0=0.,
                        t1=stimInfo['run_duration']*1000, dt=1.0, mode='peak')
                spikeTimes[N] = self.clean_spiketimes(spikeTimes[N])
                inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]['Vsoma'])
                dendriteVoltage[N] = np.array(tresults[j]['Vdend'])
                stim[N] = np.array(tresults[j]['stim'])
            stimWaveform = np.array(tresults[0]['stim'])
            
            total_elapsed_time = time.time() - self.start_time
    #        total_run_time = time.time() - run_time
            print("Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(total_elapsed_time/60., total_elapsed_time))
            print("Total Setup Time = {:8.2f}min ({:8.0f}s)".format(self.setup_time/60., self.setup_time))
            print("Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(self.an_setup_time/60., self.an_setup_time))
            print("Total Neuron Run Time = {:8.2f} min ({:8.0f}s)".format(self.nrn_run_time/60., self.nrn_run_time))
            
            result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes,
                'somaVoltage': somaVoltage, 'dendriteVoltage': dendriteVoltage, 'time': np.array(tresults[j]['time']),
                'stimWaveform': stimWaveform}
            
            self.analysis_filewriter(self.Params['cell'], result, tag=tagname % k)
        if self.Params['plotFlag']:
            self.plot_an(celltime, result)

    def an_run_IO(self, post_cell):
        """
        Establish AN inputs to soma, and run the model adjusting gmax over the reps from 0.5 to 4x.
        synapseConfig: list of tuples
            each tuple represents an AN fiber (SGC cell) with:
            (N sites, delay (ms), and spont rate group [1=low, 2=medium, 3=high])
        This routine runs a series of nReps for each synapse, turning off all of the other synapses.
        and varying the synaptic conductance to 0 (to avoid upsetting initialization)
        The output file either says "SynIO" or  
        
        Parameters
        ----------
        post_cell : cnmodel cell object
            Access to neuron and file information
        
        exclude : boolean (default: False)
            Set false to do one synapse at a time.
            Set true to do all BUT one synapse at a time.
        
        verify : boolean (default: False)
            Flag to control printing of various intermediate results
        
        
        Returns
        -------
            Nothing
        
        """
        print('\n*** an_run_IO\n')
        self.start_time = time.time()
        synapseConfig, celltype = cell_config.makeDict(self.cellID)
        nReps = self.Params['nReps']
        threshold = self.Params['threshold'] # spike threshold, mV
        
        stimInfo = self.Params
        preCell, synapse, self.electrode_site = self.configure_cell(post_cell, synapseConfig, celltype, stimInfo)

        nSyns = len(preCell)
        k = 0
        spikeTimes = {}
        inputSpikeTimes = {}
        somaVoltage = {}
        dendriteVoltage = {}
        celltime = []
        parallel = self.Params['Parallel']
        self.setup_time = time.time() - self.start_time
        self.nrn_run_time = 0.0
        self.an_setup_time = 0.
        # get the gMax's
        gMax = np.zeros(len(synapse))
        for i, s in enumerate(synapse):
            for p in s.psd.ampa_psd:
                gMax[i] = p.gmax
        for k in range(nSyns):
            # only enable gsyn on the selected input
            tagname = 'SynIO%03d'
            gmaxs = np.zeros(nReps)
            tresults = [None]*nReps
            
            if parallel and self.Params['nWorkers'] > 1:
                nWorkers = self.Params['nWorkers']
                TASKS = [s for s in range(nReps)]
                # run using pyqtgraph's parallel support
                with pg.multiprocess.Parallelize(enumerate(TASKS), results=tresults, workers=nWorkers) as tasker:
                    for j, x in tasker:
                        for i, s in enumerate(synapse):
                            for p in s.psd.ampa_psd:
                                if i != k:
                                    p.gmax = 0.  # disable all others
                                else:
                                    p.gmax = 4.0*(j+1)*gMax[i]/float(nReps+1)  # except the shosen one
                        tresults = self.single_an_run_fixed(post_cell, j, synapseConfig,
                            stimInfo, preCell, self.an_setup_time)
                        tasker.results[j] = tresults
                # retreive the data
            else:  # easier to debug
                for j, N in enumerate(range(nReps)):
                    for i, s in enumerate(synapse):
                        for p in s.psd.ampa_psd:
                            if i != k:
                                p.gmax = 0.  # disable all others
                            else:
                                p.gmax = 4.0*float(j+1)*gMax[i]/float(nReps+1)  # except the chosen one
                    tresults[j] = self.single_an_run_fixed(post_cell, j, synapseConfig,
                            stimInfo, preCell, self.an_setup_time)
            gmaxs = [4.0*float(j+1)*gMax[0]/float(nReps+1) for j in range(nReps)]
            for j, N in enumerate(range(nReps)):
                celltime.append(tresults[j]['time']) # (self.time)
                spikeTimes[N] = pu.findspikes(tresults[j]['time'], tresults[j]['Vsoma'], threshold, t0=0.,
                        t1=stimInfo['run_duration']*1000, dt=1.0, mode='peak')
                spikeTimes[N] = self.clean_spiketimes(spikeTimes[N])
                inputSpikeTimes[N] = tresults[j]['ANSpikeTimes'] # [tresults[j]['ANSpikeTimes'] for i in range(len(preCell))]
                somaVoltage[N] = np.array(tresults[j]['Vsoma'])
                dendriteVoltage[N] = np.array(tresults[j]['Vdend'])
            stimInfo['gSyns'] = gmaxs
            total_elapsed_time = time.time() - self.start_time
    #        total_run_time = time.time() - run_time
            print("Total Elapsed Time = {:8.2f} min ({:8.0f}s)".format(total_elapsed_time/60., total_elapsed_time))
            print("Total Setup Time = {:8.2f}min ({:8.0f}s)".format(self.setup_time/60., self.setup_time))
            print("Total AN Calculation Time = {:8.2f} min ({:8.0f}s)".format(self.an_setup_time/60., self.an_setup_time))
            print("Total Neuron Run Time = {:8.2f} min ({:8.0f}s)".format(self.nrn_run_time/60., self.nrn_run_time))
            
            result = {'stimInfo': stimInfo, 'spikeTimes': spikeTimes, 'inputSpikeTimes': inputSpikeTimes,
                'somaVoltage': somaVoltage, 'dendriteVoltage': dendriteVoltage, 'time': np.array(tresults[j]['time'])}
            
            self.analysis_filewriter(self.Params['cell'], result, tag=tagname % k)
        if self.Params['plotFlag']:
            self.plot_an(celltime, result)

    def single_an_run_fixed(self, post_cell, j, synapseConfig, stimInfo, preCell, an_setup_time):
        """
        Perform a single run with all AN input on the target cell turned off except for input j.
        
        Parameters
        ----------
        hf : hoc_reader object
            Access to neuron and file information
        
        j : int
            The input that will be active in this run
        
        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.
        
        stimInfo : dict
            A dictionary whose elements include the stimulus delay
        
        preCell : list
            A list of the preCell hoc objects attached to the synapses
        
        post_cell : cells object
            The target cell
        
        Returns
        -------
        anresult : dict
            A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'
        
        """
        print('\n*** single_an_run_fixed\n')
        
        hf = post_cell.hr
        try:
            cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'])
        except:
            self.an_run(post_cell, make_an_intial_conditions=True)
            print('return from an_run, initial conditions')
            try:
                cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'])
            except:
                raise ValueError('Failed initialization for cell: ', self.cellID)

        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.
        stim = {}
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        for i in range(len(preCell)):
            preCell[i].set_spiketrain([10., 20., 30., 40., 50.])
            ANSpikeTimes.append(preCell[i]._spiketrain)
        
        an_setup_time += (time.time() - an0_time)
        nrn_start = time.time()
        if self.Params['save_all_sections']:  # just save soma sections        for section in list(g):
            g = hf.sec_groups[group]
            sec = hf.get_section(section)
            self.allsecVec[sec.name()] = hf.h.Vector()
            self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)  # recording of voltage all set up here
            
        Vsoma = hf.h.Vector()
        Vdend = hf.h.Vector()
        rtime = hf.h.Vector()
        if 'dendrite' in post_cell.all_sections and len(post_cell.all_sections['dendrite']) > 0:
            dendsite = post_cell.all_sections['dendrite'][-1]
            Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        else:
            dendsite = None
            
        Vsoma.record(post_cell.soma(0.5)._ref_v, sec=post_cell.soma)
        rtime.record(hf.h._ref_t)
        hf.h.finitialize()
        hf.h.tstop = 50.
        hf.h.t = 0.
        hf.h.batch_save() # save nothing
        hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
        nrn_run_time += (time.time() - nrn_start)
        if dendsite == None:
            Vdend = np.zeros_like(Vsoma)
        anresult = {'Vsoma': np.array(Vsoma), 'Vdend': np.array(Vdend), 'time': np.array(rtime),
                'ANSpikeTimes': ANSpikeTimes, 'stim': stim}
        
        print('all sec vecs: ', self.allsecVec.keys())
        print(' flag for all sections: ', self.Params['save_all_sections'])
        if self.Params['save_all_sections']: 
            anresult['allsecv'] = self.allsecVec
            
        return anresult


    def configure_cell(self, thisCell, synapseConfig, celltype, stimInfo):
        """
        Configure the cell. This routine builds the cell in Neuron, adds presynaptic inputs
        as described in the synapseConfig, and configures those according to parameters in
        stiminfo.
        
        Parameters
        ----------
        thisCell : hoc_reader object
            Access to neuron and file information
        
        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.
        
        celltype : string
            A string describing the cell type. Determines how the channels are populated
            and how the cell is built, even if it comes from a .hoc structure.
        
        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)
        
        Returns
        -------
            preCell : list
                A list of the preCell hoc objects attached to the synapses
            post_cell : cells object
                The target cell
            synapse : list
                A list of the synapses that were created and connected to this cell
            electrode_site : Neuron hoc object
                The section and location where the electrode is attached to the cell.
        
        
        """
        
        debug = False
        if debug:
            print('hf.sec_groups : ', thisCell.sec_groups.keys())
        if debug:
            print(thisCell.print_all_mechs())
        preCell = []
        synapse = []
        # reconfigure syanpses to set the spont rate group
        stimInfo['SR'] = stimInfo['SRType']
        for i, syn in enumerate(synapseConfig):
            if stimInfo['SRType'] is 'fromcell':  # use the one in the table
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=syn['SR']))
                preCell[-1]._synapsetype = self.Params['ANSynapseType']
                stimInfo['SR'] = self.srname[syn[2]] # use and report value from table
            else:
                try:
                    srindex = self.srname.index(stimInfo['SRType'])
                    print('Retrieved index {:d} with SR type {:s}'.format(srindex, stimInfo['SRType']))
                except:
                    raise ValueError('SR type "%s" not found in Sr type list' % stimInfo['SRType'])
                
                preCell.append(cells.DummySGC(cf=stimInfo['F0'], sr=srindex))  # override
            if self.Params['verbose']:
                print('SRtype, srindex: ', stimInfo['SRType'], srindex)
            for pl in syn['postlocations']:
                postsite = syn['postlocations'][pl]
                # note that we split the number of zones between multiple sites
                if self.Params['ANSynapseType'] == 'simple':
                    print("*******Synapsetype is simple")
                    synapse.append(preCell[-1].connect(thisCell, type='simple',
                        post_opts={'postsec': pl, 'postsite': postsite[0:2]}))
                else:
                    synapse.append(preCell[-1].connect(thisCell, type='multisite',
                        pre_opts={'nzones':int(syn['nSyn']*postsite[2]), 'delay':syn['delay2']},
                        post_opts={'postsec': pl, 'postsite': postsite[0:2]}))
        if self.Params['ANSynapseType'] == 'multisite':
            for i, s in enumerate(synapse):
                s.terminal.relsite.Dep_Flag = self.Params['SynapticDepression']  # turn on or off depression computation

        electrodeSection = list(thisCell.hr.sec_groups['soma'])[0]
        electrode_site = thisCell.hr.get_section(electrodeSection)
        return (preCell, synapse, electrode_site)
    

    def set_dbspl(self, signal, dbspl):
        """Scale the level of `signal` to the given dB_SPL."""
        p0 = 20e-6
        rms = np.sqrt(np.sum(signal**2) / signal.size)
        scaled = signal * 10**(dbspl / 20.0) * p0 / rms
        return scaled


    def single_an_run(self, post_cell, j, synapseConfig, stimInfo, seeds, preCell, an_setup_time):
        """
        Perform a single run with all AN input on the target cell turned off except for input j.
        
        Parameters
        ----------
        hf : hoc_reader object
            Access to neuron and file information
        
        j : int
            The input that will be active in this run
        
        synapseConfig : dict
            A dictionary with information about the synapse configuration to use.
        
        stimInfo : dict
            A dictionary whose elements used include SRType (spont rate type), F0 (stimulus frequency)
        
        seeds : 2d numpy array
            An array listing the seeds for starting each auditory nerve input spike train
        
        preCell : list
            A list of the preCell hoc objects attached to the synapses
        
        post_cell : cells object
            The target cell
        
        an_setup_time : time object
        
        Returns
        -------
        anresult : dict
            A dictionary containing 'Vsoma', 'Vdend', 'time', and the 'ANSpikeTimes'
        
        """
        print('\n*** single_an_run\n')
               
        try:
            cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'],
                autoinit=self.Params['auto_initialize'])
        except:
            print('single_an_run: could not restore initial conditions: try creating again')
            self.an_run(post_cell, make_an_intial_conditions=True)
            print('return from inital run initial conditions #2')
            try:
                cellInit.restore_initial_conditions_state(post_cell, electrode_site=None, filename=self.Params['initANStateFile'])
            except:
                raise ValueError('Failed initialization for cell: ', self.cellID)

        # make independent inputs for each synapse
        ANSpikeTimes = []
        an0_time = time.time()
        nrn_run_time = 0.
        #
        # Generate stimuli - they are always the same for every synaptic input, so just generate once
        #
        if isinstance(stimInfo['pip_start'], list):
            pips = stimInfo['pip_start']
        else:
            pips = [stimInfo['pip_start']]
        if stimInfo['soundtype'] == 'tonepip':
            stim = sound.TonePip(rate=stimInfo['Fs'], duration=stimInfo['run_duration'],
                              f0=stimInfo['F0'], dbspl=stimInfo['dB'],
                              ramp_duration=stimInfo['RF'], pip_duration=stimInfo['pip_duration'],
                              pip_start=pips)
        elif stimInfo['soundtype'] == 'SAM':
            stim = sound.SAMTone(rate=stimInfo['Fs'], duration=stimInfo['run_duration'],
                              f0=stimInfo['F0'], dbspl=stimInfo['dB'],
                              ramp_duration=stimInfo['RF'], 
                              fmod=stimInfo['fmod'], dmod=stimInfo['dmod'],
                              pip_duration=stimInfo['pip_duration'],
                              pip_start=pips)
        else:
            raise ValueError('StimInfo sound type %s not implemented' % stimInfo['soundtype'])
        
        stimWaveform = stim.generate()
        stimTimebase = stim.time
        for i, syn in enumerate(synapseConfig):
            nseed = seeds[j, i]
            if self.Params['SGCmodelType'] in ['Zilany']:
                preCell[i].set_sound_stim(stim, seed=nseed, simulator='matlab')  # generate spike train, connect to terminal
            elif self.Params['SGCmodelType'] in ['cochlea']:
                wf = self.set_dbspl(stim.generate(), stimInfo['dB'])
                stim._sound = wf
                preCell[i].set_sound_stim(stim, seed=nseed, simulator='cochlea')  # generate spike train, connect to terminal
            else:
                raise ValueError('SGC model type type %s not implemented' % self.Params['SGCmodelType'])
            ANSpikeTimes.append(preCell[i]._spiketrain)
        
        an_setup_time += (time.time() - an0_time)
        nrn_start = time.time()
        Vsoma = post_cell.hr.h.Vector()
        Vdend = post_cell.hr.h.Vector()
        rtime = post_cell.hr.h.Vector()
        if 'dendrite' in post_cell.all_sections and len(post_cell.all_sections['dendrite']) > 0:
            dendsite = post_cell.all_sections['dendrite'][-1]
            Vdend.record(dendsite(0.5)._ref_v, sec=dendsite)
        else:
            dendsite = None
            
        if self.Params['save_all_sections']:
            self.allsecVec = OrderedDict()
            for group in post_cell.hr.sec_groups.keys(): # get morphological components
                g = post_cell.hr.sec_groups[group]
                for section in list(g):
                    sec = post_cell.hr.get_section(section)
                    self.allsecVec[sec.name()] = post_cell.hr.h.Vector()
                    self.allsecVec[sec.name()].record(sec(0.5)._ref_v, sec=sec)  # recording of voltage all set up here            
        Vsoma.record(post_cell.soma(0.5)._ref_v, sec=post_cell.soma)
        rtime.record(post_cell.hr.h._ref_t)
        post_cell.hr.h.finitialize()
        post_cell.hr.h.tstop = stimInfo['run_duration']*1000.
        post_cell.hr.h.t = 0.
        post_cell.hr.h.batch_save() # save nothing
        post_cell.hr.h.batch_run(post_cell.hr.h.tstop, post_cell.hr.h.dt, "an.dat")
        nrn_run_time += (time.time() - nrn_start)
        if dendsite == None:
            Vdend = np.zeros_like(Vsoma)
        
        anresult = {'Vsoma': np.array(Vsoma), 'Vdend': np.array(Vdend), 'time': np.array(rtime), 'ANSpikeTimes': ANSpikeTimes,
            'stim': stim, 'stimWaveform': stimWaveform, 'stimTimebase': stimTimebase}
        if self.Params['save_all_sections']:
            anresult['allsecVec'] = self.allsecVec
        return anresult

    def clean_spiketimes(self, spikeTimes, mindT=0.7):
        """
        Clean up spike time array, removing all less than mindT
        spikeTimes is a 1-D list or array
        mindT is difference in time, same units as spikeTimes
        If 1 or 0 spikes in array, just return the array
        """
        if len(spikeTimes) > 1:
            dst = np.diff(spikeTimes)
            st = np.array(spikeTimes[0])  # get first spike
            sok = np.where(dst > mindT)
            st = np.append(st, [spikeTimes[s+1] for s in sok])
            spikeTimes = st
        return spikeTimes
    
    def plot_an(self, celltime, result):
        """
        Plot the cell's voltage reponse to the AN inputs
        
        Parameters
        ----------
        celltime : array (no default)
            time array for the cell voltage data
        somaVoltage : array (no default)
            Voltage recorded at the soma for this cell
        stimInfo : dict (no default)
            stimulus parameters (require 'nReps' and 'threshold' for spikes)
        dendVoltage : array (default: None)
            Voltage recorded at some point in the dendrite for this cell
        
        Returns
        -------
        Nothing
         
        """
        if not self.Params['plotFlag']:
            return
        nReps = self.Params['nReps']
        threshold = self.Params['threshold']
        win = pgh.figure(title='AN Inputs')
        layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
        for j, N in enumerate(range(len(result))):
            layout.plot(0, celltime[N], result[N]['somaVoltage'], pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.plot(0, [np.min(celltime[N]), np.max(celltime[N])], [threshold, threshold],
                 pen=pg.mkPen((0.5, 0.5, 0.5), width=0.5))
            dvdt = np.diff(result[N]['somaVoltage'])/np.diff(celltime[N])  # show in mV/ms
            layout.plot(1, celltime[N][:-1], dvdt, pen=pg.mkPen(pg.intColor(N, nReps)))
            layout.getPlot(0).setXLink(layout.getPlot(1))
        if 'dendVoltage' in result[N].keys() and result[N]['dendVoltage'] is not None:
            for j, N in enumerate(range(len(result))):
                layout.plot(0, celltime[N], result[N]['dendVoltage'], pen=pg.mkPen(pg.intColor(N, nReps)),
                    style=pg.QtCore.Qt.DashLine)
        pgh.show()
    
    def analysis_filewriter(self, filebase, result, tag=''):
        """
        Write the analysis information to a pickled file
        
        Parameters:
        filebase : string (no default)
            base filename - *not used* (value replaced bye cellID)
        result : dict (no default)
            dict hlding results. Must be pickleable
        tag : string (default: '')
            tag to insert in filename string
        """
        k = result[0].keys()
        requiredKeys = ['stimInfo', 'spikeTimes', 'inputSpikeTimes', 'somaVoltage', 'time', 'stimWaveform', 'stimTimebase']
        for rk in requiredKeys:
            assert rk in k
        results = {}
        # results with be a dict with params, stiminfo, and trials as keys
        print('\n*** analysis_filewriter\n')

        results['Params'] = self.Params  # include all the parameters of the run too
        results['trials'] = result
        with open(self.Params['simulationFilename'], 'wb') as fh:
            pickle.dump(results, fh)
        print(f"**** Model output written to: ****\n   {str(self.Params['simulationFilename']):s}")
        self.ANFilename = str(self.Params['simulationFilename'])

    def get_hoc_file(self, hf):
        if hf.file_loaded is False:
            exit()
        self.section_list = hf.get_section_prefixes()
        hf.sec_groups.keys()
        if len(hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.colorMap[i]
#        else: # single section name, assign colors to SectionList types:
#        self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
#            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}
        
        (v, e) = hf.get_geometry()
        self.clist = []
        
        for si in hf.sections: # self.section_list[s]:
            hf.h('access %s' % si)
            sr = hf.h.SectionRef()
            n1 = hf.h.cas().name()
            if sr.has_parent() == 1:
                x=sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])


if __name__ == "__main__":
    curdir = Path.cwd()
    model = ModelRun()  # create instance of the model
    
    parser = argparse.ArgumentParser(description='Simulate activity in a reconstructed model cell',
                    argument_default=argparse.SUPPRESS,
                    fromfile_prefix_chars='@')
    parser.add_argument(dest='cell', action='store',
                   default=None,
                   help='Select the cell (no default)')
    parser.add_argument('--type', '-T', dest='cellType', action='store',
                   default='Bushy', choices=model.cellChoices,
                   help='Define the cell type (default: Bushy)')
    parser.add_argument('--model', '-M', dest='modelName', action='store',
                   default='XM13', choices=model.modelNameChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--modeltype', dest='modelType', action='store',
                   default='II', choices=model.modelTypeChoices,
                   help='Define the model type (default: XM13)')
    parser.add_argument('--sgcmodel', type=str, dest='SGCmodelType', action='store',
                   default='Zilany', choices=model.SGCmodelChoices,
                   help='Define the SGC model type (default: Zilany)')
    parser.add_argument('--protocol', '-P', dest='runProtocol', action='store',
                   default='runIV', choices=model.protocolChoices,
                   help='Protocol to use for simulation (default: IV)')
    parser.add_argument('-H', action='store_true', dest='usedefaulthoc',
                  help='Use default hoc file for this cell')
    parser.add_argument('--hocfile', dest='hocfile', action='store',
                  default=None,
                  help='hoc file to use for simulation (default is the selected "cell".hoc)')
    parser.add_argument('--inputpattern', '-i', type=str, dest='inputPattern', action='store',
                  default=None,
                  help='cell input pattern to use (substitute) from cell_config.py')
    parser.add_argument('--stimulus', '-s', type=str, dest='soundtype', action='store',
                   default='tonepip', choices=model.soundChoices,
                   help='Define the stimulus type (default: tonepip)')
    parser.add_argument('--check', '-/', action='store_true', default=False, dest='checkcommand',
                   help='Only check command line for valid input; do not run model')
                   
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
        choices=model.SRChoices,
        help=('Specify SR type (from: %s)' % model.SRChoices))
    parser.add_argument('--synapsetype', type=str, default='multisite', dest = 'ANSynapseType',
        choices=model.ANSynapseChoices,
        help=('Specify AN synpase type (from: %s)' % model.ANSynapseChoices))
    parser.add_argument('--depression', type=int, default=0, dest = 'ANSynapticDepression',
        choices=[0, 1],
        help=('Specify AN depression flag for multisite synapses (from: %s)' % str([0, 1])))
    
    parser.add_argument('--fmod', type=float, default=20, dest = 'fmod',
        help='Set SAM modulation frequency')
    parser.add_argument('--dmod', type=float, default=100., dest = 'dmod',
        help='Set SAM modulation depth (in percent)')
    parser.add_argument('--S2M', type=float, default=0, dest = 'signalToMasker',
        help='Signal to Masker ratio (dB)')
    parser.add_argument('--cmmrmode', type=str, default='CMR', dest = 'CMMRmode',
        choices=model.cmmrModeChoices,
        help=('Specify mode (from: %s)' % model.cmmrModeChoices))
    
    parser.add_argument('--spirou', type=str, dest='spirou', action='store', default='all',
            choices = model.spirouChoices,
            help='Specify spirou experiment type.... ')
    parser.add_argument('--soma-inflate', type=float, dest='soma_inflation', action='store', default=1.0,
            help='Specify factor by which to inflate soma AREA')
    parser.add_argument('--soma-autoinflate', action='store_true', dest='soma_autoinflate', default=False,
            help='Automatically inflate soma based on table')

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

    args = vars(parser.parse_args())
    
    for k in args.keys():
        model.Params[k] = args[k]
        
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
    if model.Params['checkcommand']:
        model.Params['commandline'] = ' '.join(sys.argv)
        print(json.dumps(model.Params, indent=4))  # pprint doesn't work well with ordered dicts
        print(model.Params['commandline'])

    if not model.Params['checkcommand']:
        model.Params['Icmds'] = '[-2.0, 2.0, 0.5]'
        model.run_model() # then run the model
    
    # if sys.flags.interactive == 0:
    #     pg.QtGui.QApplication.exec_()
    
