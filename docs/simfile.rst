Simulation File Structure
=========================



Simulations are controlled by shell scripts in the `scripts` subdirectory. 
These scripts set up the parameters for the runs, including the *cnmodel* datatables and various stimulus parameters. 
The results are stored on a by-cell basis in the target directory as "VCN-SBEM-Data/VCN_Cells/VCN_Cnn/IV (or AN)" with a long filename. 
The scripts invoke *model_run2.py*, which handles all of the simulations and "experiments", wrapped around cnmodel. 
The data files contain all of the information regarding the model run, as stored in the dataclasses Params and RunInfo 
(there should be no control parameters that are not described in these tables). 
This includes the full text of the "data\_table" that is used ton control cell decoration of channels in cnmodel. 

Analysis is handled through "DataTablesVCN.py", which has a graphical interface. This is described in more detail here.


The simulation result files are Python pickle files that containe both simulation results and the
control parameters. In theory, each result file contains sufficient information that the simulation
could be run (from model_run2.py) without reference to any external files.

The top level dictionary contains the following entries:

    * "basename" - The full path and filename to the file itself. 
    * "Params"  - a dataclass containing the parameters that were passed to model_run2.py. The dataclass 
        contains all
        of the parameters that control a given simulation, including the *cnmodel* datatables (text tables). The Params
        dataclass is defined in vcnmodel/model_params.py. Note that model_params also defines the runInfo class, and 
        contains the command line parser, which fills the dataclass values directly as needed.

    * "runInfo" - a dataclass containing the parmeters for the specific run. This has a different set of details than
        the Params dataclass.

    * "Results" - a dictionary, keyed by trial, containing numpy arrays (vectors) of the simulation results. 
       The simulation results for each run are in a dictionary with keys: 'spikeTimes', 'inputSpikeTimes', 
       'time', 'somaVoltage', 'dendriteVoltage', 'allDendriteVoltages', 'stimWaveform', 'stimTimebase']. `allDendriteVoltages` 
       is only included if Params.save_all_sections is True; Normally, only the `somaVoltage`
       is saved. `spikeTimes` are the cell spike times, `inputSpikeTimes` is a
       list for the spike times of each input (in order). `stimWaveform` and `StimTimebase` are the stimulus. 

    * "modelPars" - a brief structure that holds information also found in Params and runInfo, for example:
    
        {'species': 'mouse', 'cellClass': 'bushy', 'modelType': 'II', 'modelName': 'XM13A_nacncoop',
        'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
        'hillock': False, 'initialsegment': False, 'myelinatedaxon': False,
        'unmyelinatedaxon': False, 'na': 'nacncoop', 'ttx': False,
        'name': 'bushy', 
        'morphology': '/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c17/Morphology/VCN_c17_Full_MeshInflate.hoc',
        'temperature': 34.0}
    
    * "mode" - text string "syn", "reps", etc. Not useful (duplicated information).



Dataclasses Reference (examples)
--------------------------------

Params::

 setup : True
   cellID : VCN_c17
   cell : VCN_c17
   AMPAScale : 1.0
   ANSynapseType : multisite
   ANSynapticDepression : 0
   SynapseConfiguration : None
   synapseConfig : None
   initStateFile : /Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c17/Initialization/AN_Init_VCN_c17_inp=self_XM13A_nacncoop_II_HF=VCN_c17_Full_MeshInflate_normal_all_multisite_HS.p
   simulationFilename : /Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c17/Simulations/AN/runANPSTH-all-2021-11-29.11-43-53/AN_Result_VCN_c17_inp=self_XM13A_nacncoop_II_HF=VCN_c17_Full_MeshInflate_normal_all_multisite_001_tonepip_010dB_16000.0_HS.p
   shortSimulationFilename : None
   simPath : None
   hocfile : VCN_c17_Full_MeshInflate.hoc
   meshInflate : False
   usedefaulthoc : False
   cellType : Bushy
   modelName : XM13A_nacncoop
   dataTable : data_XM13A_nacncoop_normal
   dendriteMode : normal
   dendriteExpt : Full
   axonExpt : default
   modelType : II
   SGCmodelType : cochlea
   species : mouse
   displayStyle : cylinders
   displayMechanism : klt
   displayMode : None
   fullhocfile : False
   dtIC : 0.025
   dtVC : 0.005
   celsius : 37
   Ra : 150.0
   soma_inflation : 1.0
   soma_autoinflate : False
   dendrite_inflation : 1.0
   dendrite_autoinflate : False
   dendrite_fromsoma : False
   ASA_inflation : 1.0
   ASA_fromsoma : False
   lambdaFreq : 2000.0
   area_adjustment_method : pt3d
   SRType : HS
   inputPattern : None
   synno : None
   lastfile : None
   seed : 9
   all_modes : False
   plotFlag : False
   auto_initialize : False
   nWorkers : 4
   Parallel : True
   verbose : False
   save_all_sections : False
   commandline : Namespace(AMPAScale=1.0, ANSynapseType='multisite', ANSynapticDepression=0, ASA_fromsoma=False, 
    CMMRmode='CMR', F0=16000.0, Parallel=True, SGCmodelType='cochlea', SRType='HS', Spirou='all', all_modes=False, 
    auto_initialize=False, axonExpt='default', cell='VCN_c17', cellType='Bushy', checkcommand=False, 
    configfile='xm13a_multisite_parallel.toml', dB=10.0, dataTable='data_XM13A_nacncoop_normal', 
    dendriteExpt='Full', dendriteMode='normal', dendrite_autoinflate=False, 
    dendrite_fromsoma=False, dendrite_inflation=1.0, displayMechanism='klt', displayMode='None', 
    displayStyle='cylinders', displayscale=False, dmod=0.0, fmod=100.0, 
    gif_dur=10.0, gif_fmod=0.2, gif_i0=0.0, gif_sigma=0.2, gif_skew=0.0, gif_tau=3.0, 
    hocfile=None, inputPattern=None, meshInflate=False, modelName='XM13A_nacncoop', 
    modelType='II', nReps=1, nWorkers=4, pip_duration=0.1, pip_offduration=0.1, 
    pip_start=0.2, plotFlag=False, runProtocol='runANPSTH', save_all_sections=False, 
    seed=9, sequence='', signalToMasker=0, soma_autoinflate=False, soma_inflation=1.0, 
    soundtype='tonepip', tagstring=None, testsetup=False, verbose=False, vstimHolding=-80.0)
   commands : ['VCN_c17', '-D', 'Full', '-P', 'runANPSTH', '-r', '1', '--dB', '10', '--Spirou', 'all', 
    '--dendritemode', 'normal', '--configfile', 'xm13a_multisite_parallel.toml',
    '--datatable', 'data_XM13A_nacncoop_normal', '--saveall']
   checkcommand : False
   testsetup : False
   configfile : xm13a_multisite_parallel.toml
   tagstring : None
   initialization_time : 200.0
   SynapseConfig : [
    OrderedDict([('input', 1), ('asa', 278.32), ('synperum2', 0.7686), 
      ('nSyn', 214), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan), 
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'), 
      ('postlocations', {'soma': [461, 0.5, 1.0]})]),
    OrderedDict([('input', 2), ('asa', 261.49), ('synperum2', 0.7686), 
      ('nSyn', 201), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan), 
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'), 
      ('postlocations', {'soma': [461, 0.5, 1.0]})]),
    OrderedDict([('input', 3), ('asa', 105.06), ('synperum2', 0.7686), 
      ('nSyn', 81), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan), 
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'), 
      ('postlocations', {'soma': [461, 0.5, 1.0]})]), 
    OrderedDict([('input', 4), ('asa', 62.36), ('synperum2', 0.7686), 
      ('nSyn', 48), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan),
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'),
      ('postlocations', {'soma': [461, 0.5, 1.0]})]),
    OrderedDict([('input', 5), ('asa', 41.04), ('synperum2', 0.7686),
      ('nSyn', 32), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan), 
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'), 
      ('postlocations', {'soma': [461, 0.5, 1.0]})]),
    OrderedDict([('input', 6), ('asa', 38.19), ('synperum2', 0.7686),
      ('nSyn', 29), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan), 
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'), 
      ('postlocations', {'soma': [461, 0.5, 1.0]})]),
    OrderedDict([('input', 7), ('asa', 36.75), ('synperum2', 0.7686), 
      ('nSyn', 28), ('delay', 0.0), ('SR', 2), ('delay2', 0.0), ('axonLen', nan),
      ('axonR', nan), ('branchLen', nan), ('branchR', nan), ('type', 'AN'),
      ('postlocations', {'soma': [461, 0.5, 1.0]})])]


runInfo::

   folder : Simulations
   fileName : Normal
   runProtocol : runANPSTH
   runName : Run
   manipulation : Canonical
   preMode : CC
   postMode : CC
   TargetCellType : 
   electrodeSection : sections[461]
   electrodeSectionName : soma
   dendriticElectrodeSection : None
   dendriticSectionDistance : 100.0
   ChannelCompartments, 

    This table describes the ion channel densities relative to somatic densities,
    e.g., relative to REFERENCE densities in the table XM13_channels.
    and voltage shifts, for different compartments of the specified neuron,
    Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse
    (data table: XM13_channels).

    NOTE: unmyelinatedaxon and initialsegment are equivalent in George's models, but only "unmyelinatedaxon" is actually used.
    --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                       axon       Unmyelinated_Axon    Myelinated_Axon    Axon_Initial_Segment    Axon_Hillock     soma        Proximal_Dendrite     Distal_Dendrite    Dendritic_Hub     Dendritic_Swelling
                                                                                                                                                      
    nacncoop_gbar      1.0 [1]    100.0 [1]            0.0 [1]            100.0 [1]               5.0 [1]          1.0 [1]     0.5 [1]               0.5 [1]            0.5 [1]           0.5 [1] 
    kht_gbar           1.0 [1]    2.0 [1]              0.01 [1]           2.0 [1]                 1.0 [1]          1.0 [1]     0.5 [1]               0.5 [1]            0.5 [1]           0.5 [1] 
    klt_gbar           1.0 [1]    1.0 [1]              0.01 [1]           2.0 [1]                 1.0 [1]          1.0 [1]     0.5 [1]               0.5 [1]            0.5 [1]           0.5 [1] 
    ihvcn_gbar         0.0 [1]    0.5 [1]              0.0 [1]            0.5 [1]                 0.0 [1]          1.0 [1]     0.5 [1]               0.5 [1]            0.5 [1]           0.5 [1] 
    leak_gbar          1.0 [1]    1.0 [1]              0.25e-3 [1]        1.0 [1]                 1.0 [1]          1.0 [1]     0.5 [1]               0.5 [1]            0.5 [1]           0.5 [1] 
    leak_erev          -65. [1]   -65. [1]             -65. [1]           -65. [1]                -65. [1]         -65. [1]    -65. [1]              -65. [1]           -65. [1]          -65. [1]
    na_type            nacncoop   nacncoop             nacncoop           nacncoop                nacncoop         nacncoop    nacncoop              nacncoop           nacncoop          nacncoop
    nacncoop_vshift    0.  [2]    0. [2]               0. [2]             0. [2]                  0. [2]           0. [2]      0.  [2]               0.  [2]            0.  [2]           0.  [2] 
    ih_type            ihvcn      ihvcn                ihvcn              ihvcn                   ihvcn            ihvcn       ihvcn                 ihvcn              ihvcn             ihvcn   
    --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    [1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.
    [2] Set to 0 (was 4.3 in original model). Matches original Barela et al (2006) scaling.

    
      ChannelData, 

    This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
    for different cell types based on the Xie and Manis 2013 models for mouse.

    The REFERENCE values are applied to "point" models, and to the soma of
    compartmental models.
    The names of the mechanisms must match a channel mechanism (Neuron .mod files)
    and the following _(gbar, vshift, etc) must match an attribute of that channel
    that can be accessed.
    
    -----------------------------------------------------------------------------------------------------------------
                      II                   II-I                I-c                I-II                I-t       
                                                                                                    
    nacncoop_gbar     29.110       [1]     1000.       [1]     3000.       [1]    1000.        [2]    3000.    [1] 
    kht_gbar          1.6884       [3]     58.0        [1]     500.0       [1]    150.0        [2]    500.0    [1] 
    klt_gbar          2.3288       [1]     14.0        [1]     0.0         [1]    20.0         [2]    0.0      [1] 
    ka_gbar           0.0000       [1]     0.0         [1]     0.0         [1]    0.0          [2]    125.0    [1] 
    ihvcn_gbar        0.8733       [1]     30.0        [1]     18.0        [1]    2.0          [2]    18.0     [1] 
    leak_gbar         0.1385       [1]     2.0         [1]     8.0         [1]    2.0          [2]    8.0      [1] 
    leak_erev         -65.0        [1]     -65         [1]     -65         [1]    -65          [2]    -65      [1] 
    na_type           nacncoop     [1]     nacncoop    [1]     nacncoop    [1]    nacncoop     [1]    nacncoop [1] 
    ih_type           ihvcn        [1]     ihvcn       [1]     ihvcn       [1]    ihvcn        [2]    ihvcn    [1] 
    soma_Cap          13.0         [1]     26.0        [1]     25.0        [1]    26.0         [2]    25.0     [1] 
    e_k               -84.0        [1]     -84         [1]     -84         [1]    -84          [2]    -84      [1] 
    e_na              50.0         [1]     50.         [1]     50.         [1]    50.          [2]    50.      [1] 
    ih_eh             -43.0        [1]     -43         [1]     -43         [1]    -43          [2]    -43      [1] 
    nacncoop_vshift   0.           [1]     0.          [1]     0.          [1]    0.           [1]    0.       [1]
    units             mmho/cm2             nS                  nS                 nS                  nS
    -------------------------------------------------------------------------------------------------------------------

    [1] Uses channels from Rothman and Manis, 2003
        Conductances are for Mouse bushy cells
        Xie and Manis, 2013
        Age "adult", Temperature=34C
        Units are nS unless otherwise stated.
        nacn_vshift: was 4.3 in Xie and Manis (2013) for T-stellate cells; 0 for bushy cells
        Here reset to 0 for bushy cells

    [2] Rothman and Manis, 2003, model I-II
        Some low-voltage K current, based on observations of
        a single spike near threshold and regular firing for higher
        currents (Xie and Manis, 2017)
        
    [3] Increased DR to force AP repolarization to be faster


    
   nReps : 1
   seeds : [[ 9 10 11 12 13 14 15]]
   sequence : 
   nStim : 1
   stimFreq : 200.0
   stimInj,  pulse: [-1.  -0.8 -0.6 -0.4 -0.2  0.   0.2  0.4  0.6  0.8  1.   1.2  1.4  1.6 1.8  2. ]
   stimVC,  pulse: [-20. -10.   0.  10.  20.  30.  40.  50.  60.  70.  80.  90. 100. 110. 120.]
   stimDur : 100.0
   stimDelay : 5.0
   stimPost : 3.0
   vnStim : 1
   vstimFreq : 200.0
   vstimInj : 50
   vstimDur : 100.0
   vstimDelay : 5.0
   vstimPost : 25.0
   vstimHolding : -80.0
   initialization_time : 50.0
   run_duration : 0.4
   soundtype : tonepip
   pip_duration : 0.1
   pip_start : 0.2
   pip_offduration : 0.1
   Fs : 100000.0
   F0 : 16000.0
   dB : 10.0
   RF : 0.0025
   fmod : 100.0
   dmod : 0.0
   threshold : -35.0
   signalToMasker : 0
   CMMRmode : CMR
   Spirou : all
   gif_i0 : 0.0
   gif_sigma : 0.2
   gif_fmod : 0.2
   gif_tau : 3.0
   gif_dur : 10.0
   gif_skew : 0.0
   runTime : Mon Nov 29 11:43:52 2021
   inFile : None
   inFileRep : 1
   v_init : -61.0
   useSaveState : True

