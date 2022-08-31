"""make_data_dirs.py
Create a directory structure based on figure_data.py (in the plotters
subdirectory) and populate it with all of the simulation files and intermediate
analyses needed to regenerate the figures using DataTablesVCN.py This also
populates the Morphology directory for each cell with the .hoc and .swc
morphology files.

The directory created by this routine was compressed (zip) and uploaded
to figshare:
Private link: https://figshare.com/s/2f87300f58e1c1e9270c
(DOI: 10.6084/m9.figshare.19520041)

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2014- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
# import multiprocessing
import shutil
import subprocess
from pathlib import Path
from typing import List

import toml
from pylibrary.tools import cprint as CP
from vcnmodel.plotters import \
    figure_data as FD  # table of simulation runs used for plotting figures
import vcnmodel.util.ptree as ptree

cprint = CP.cprint

with open("wheres_my_data.toml", "r") as fh:
    config = toml.load(fh)
sourcepath = Path(config["baseDataDirectory"])

simpath = Path("/Volumes/Pegasus_002/BU_simulation_data")
# simpath = Path("/Volumes/T7SSD/BU_simulation_data")
figpath = Path("/Volumes/Pegasus_002/VCN-SBEM-Data/SBEM-paper Figures")
intermediate = Path(simpath, "IntermediateAnalyses")
intermediate.mkdir(exist_ok=True, parents=True)
intermediate_source = Path(figpath, "IntermediateResults")

simulations = Path(simpath, "Simulations")
simulations.mkdir(exist_ok=True, parents=True)

Zin_dest_dir = Path(simulations, "Impedance_Calculations")
Zin_dest_dir.mkdir(exist_ok=True, parents=True)
Zin_source_dir = Path(sourcepath, "VCN_Cells", "Impedance_Calculations")




readmefilecontents="""This VCN_SBEM_readme.txt file was generated on 2022-08-30
by Paul B. Manis, Ph.D.


GENERAL INFORMATION

1. Title of Dataset: Computational models/results of serial-blockface electron
   microscopic
reconstructions of globular bushy cells from the mouse ventral cochlear nucleus.

2. Author Information
    A. Principal Investigator Contact Information
        Name: Paul B. Manis Institution: University of North Carolina at Chapel
        Hill Address: CB#7070, B027 Marsico Hall, 125 Mason Farm Road, Chapel
        Hill NC 27599-7070. Email: pmanis@med.unc.edu

    B. Principal Investigator Contact Information
        Name: George A. Spirou Institution: Univesity of South Florida Address:
        Tampa FL Email: gspirou@usf.edu

3. Date of data collection (single date, range, approximate date): 2014-01-01
to 2022-05-15

4. Geographic location of data collection : Chapel Hill, NC USA, Morgantown, WV,
    Tampa, FL.

5. Information about funding sources that supported the collection of the data: 
    NIH grants: DC R01 DC015901 (Spirou, Manis, Ellisman), DC R01 DC004551
    (Manis, 2013-2019, Early development) DC R01 DC019053 (Manis, 2020-2025,
    Later development)

SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: MIT License (see LICENSE.txt in the
   github repository)

2. Links to publications that cite or use the data: 

3. Links to other publicly accessible locations of the data: None

4. Links/relationships to ancillary data sets: None

5. Was data derived from another source? No

6. Recommended citation for this dataset: 


DATA & FILE OVERVIEW

1. File List: 
  See the end of this file. 
 
2. Relationship between files, if important: 

3. Additional related data collected that was not included in the current data
   package: 
    A large number of additional simulations were run over the 8 years of this
    project. The deposited simulations are those used in the submitted
    manuscript, which represent the results of several levels of refinement in
    the approach over the years (these include types of annotation of cells, the
    addition of new reconstructed cells, inflation of the cell surface area from
    the reconstructions to match the mesh reconstructions, manipulations of the
    cell structure, exploration of different stimulus types and patterns, a
    change to the file structure for the simulations.) The original volume is
    not included as it is over 1.5 TB in size. 
    

4. Are there multiple versions of the dataset? No
    A. If yes, name of file(s) that was updated: 
        i. Why was the file updated? 
        ii. When was the file updated? 


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
    Reconstructions were made in syGlass (IstoVisio, Morgantown, WV) as SWC
    files (Cannon et al., 1988), and converted to the HOC format for use in
    NEURON for simulations. The SWC files were traced by hand from
    semi-automated mesh and annotated reconstructions of the cell membranes
    (soma, dendrites, axons).

2. Methods for processing the data: 
    The SWC files were converted to the HOC format for use in NEURON for
    simulations. Code for running the simulations was developed
    (github.com/pbmanis/VCN_Model), based on the cnmodel package
    (github.com/cnmodel/cnmodel; Manis and Campagnola, 2018). 

3. Instrument- or software-specific information needed to interpret the data: 
    Python 3.8 or later. The complete set of package requirements and a script
    to build the required environment can be found in requirements.txt at the
    github repository listed in (2).

4. Standards and calibration information, if appropriate: None

5. Environmental/experimental conditions: None

6. Describe any quality-assurance procedures performed on the data: 
    SWC files: visual comparision with mesh files in syGlass. HOC files:
    Integrity of file using the "topology" function in Neuron. Python: Test
    suites in cnmodel; test suite in VCN_Model (for analysis).

7. People involved with sample collection, processing, analysis and/or
   submission: 
    Paul Manis, George Spirou, Mark Ellisman, Matthew Kersting.

DATA-SPECIFIC INFORMATION FOR: [FILENAME] <repeat this section for each dataset,
folder or file, as appropriate>

1. Number of variables: Not applicable.

2. Number of cases/rows: Not applicable.

3. Variable List: Not applicable.

4. Missing data codes: Not applicable

5. Specialized formats or other abbreviations used: 
    The SWC files have an  extended set of descriptors, in addition to the
    standard ones described in Cannon et al. (1988). These are:
        10: "Axon_Hillock", 11: "Dendritic_Swelling", 12: "Dendritic_Hub", 13:
        "Proximal_Dendrite", 14: "Distal_Dendrite", 15: "Axon_Initial_Segment",
        16: "Axon_Heminode", 17: "Axon_Node",
    (in an earlier version the order was different; this is described in
    swc_to_hoc.py in the github repository).

    The simulation result files are held as ".p" files in a date/timestamped
    directory for each simulation. You can run
    vcnmodel/util/inspect_simulation_file.py to get a listing for a specific
    file. The output of that program will be Each file consists of the
    following:
        Parameters (Params) dataclass.  See vcnmodel/model_params.py RunInfo
        (RunInfo) dataclass.   See vcnmodel/model_params.py Results dictionary.
        Generated in vcnmodel/model_run2.py "modelPars" dictionary (somewhat
        redundant with Parameters and RunInfo) basename (name of the file as
        _originally_ stored).
 The data file is a dictionary with several entries:
└── The top level dictionary keys are:  ['basename\n', 'Params\n', 'runInfo\n',
'modelPars\n', 'Results\n', 'mode\n'] └── Key: basename
    └── contains:
    /Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c05/Simulations/AN/runANPSTH-all-2021-08-03.12-53-12/AN_Result_VCN_c05_inp=self_XM13A_nacncoop_II_HF=VCN_c05_Full_MeshInflate_standardized_axon_normal_all_multisite_025_tonepip_030dB_16000.0_HS.p
└── Key: Params
    └── has fields: 
        └── AMPAScale └── ANSynapseType └── ANSynapticDepression └──
        ASA_fromsoma └── ASA_inflation └── Parallel └── Ra └── SGCmodelType └──
        SRType └── SynapseConfig └── SynapseConfiguration └── all_modes └──
        area_adjustment_method └── auto_initialize └── axonExpt └── cell └──
        cellID └── cellType └── celsius └── checkcommand └── commandline └──
        commands └── configfile └── dataTable └── dendriteExpt └── dendriteMode
        └── dendrite_autoinflate └── dendrite_fromsoma └── dendrite_inflation
        └── displayMechanism └── displayMode └── displayStyle └── dtIC └── dtVC
        └── fullhocfile └── hocfile └── hocstring └── initStateFile └──
        initialization_time └── inputPattern └── lambdaFreq └── lastfile └──
        meshInflate └── modelName └── modelType └── nWorkers └── plotFlag └──
        save_all_sections └── seed └── setup └── shortSimulationFilename └──
        simPath └── simulationFilename └── soma_autoinflate └── soma_inflation
        └── species └── srnames └── synapseConfig └── synno └── tagstring └──
        testsetup └── usedefaulthoc └── verbose
└── Key: runInfo
    └── has fields: 
        └── CMMRmode └── F0 └── Fs └── RF └── Spirou └── SpirouSubs └──
        TargetCellType └── clickDuration └── clickRate └── clickStart └──
        clickTrainDuration └── dB └── dendriticElectrodeSection └──
        dendriticSectionDistance └── dmod └── electrodeSection └──
        electrodeSectionName └── fileName └── fmod └── folder └── gif_dur └──
        gif_fmod └── gif_i0 └── gif_sigma └── gif_skew └── gif_tau └── inFile
        └── inFileRep └── initialization_time └── manipulation └── nReps └──
        nStim └── noise_seed └── pip_duration └── pip_offduration └── pip_start
        └── postMode └── preMode └── runName └── runProtocol └── runTime └──
        run_duration └── seeds └── sequence └── signalToMasker └── soundtype └──
        spikeTimeList └── stimDelay └── stimDur └── stimFreq └── stimInj └──
        stimPost └── stimVC └── tableData └── threshold └── useSaveState └──
        v_init └── vnStim └── vstimDelay └── vstimDur └── vstimFreq └──
        vstimHolding └── vstimInj └── vstimPost
└── Key: modelPars
    └── contains:  species = mouse └── contains:  cellClass = bushy └──
    contains:  modelType = II └── contains:  modelName = XM13A_nacncoop
soma, True axon, False dendrites, False pumps, False hillock, False
initialsegment, False myelinatedaxon, False unmyelinatedaxon, False
    └── contains:  na = nacncoop
ttx, False
    └── contains:  name = bushy └── contains:  morphology =
    /Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c05/Morphology/VCN_c05_Full_MeshInflate_standardized_axon.hoc
temperature, 34.0 └── Key: Results
    └── contains: Trials  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24]        └──  One Trial contains: 
            └──  spikeTimes with length of 149 elements └──  inputSpikeTimes
            with length of 7 elements └──  time with length of 50001 elements
            └──  somaVoltage with length of 50001 elements └──  dendriteVoltage
            with length of 50001 elements └──  allDendriteVoltages with length
            of 0 elements └──  stimWaveform with length of 125001 elements └──
            stimTimebase with length of 125001 elements
└── Key: mode
    └── contains: syn

The directory entry for each simulation file has an associated file with the
same name as the directory, but holding some information so that the entire
dataset does not need to be read. These files are generated by DataTablesVCN.py
to speed up building the tables when new simulations are added. The structure is
as follows for an example file: Example file:
/Volumes/Pegasus_002/BU_simulation_data/Simulations/VCN_c05/Simulations/AN/runANPSTH-all-2021-08-03.12-53-12.pkl
    ANSynapticDepression = 0
                  SRType = HS
                axonExpt = standardized cellType = Bushy
            changed_data = None cnmodel_hash =
            b'eee593c673752c19137658d5b9a381ea9ad4580f' command_line = {'cell':
            'VCN_c05', 'cellType': 'Bushy', 'modelName': 'XM13A_nacncoop',
            'modelType': 
                           'II', 'dendriteMode': 'normal', 'hocfile': None,
                           'meshInflate': False, 'dendriteExpt': 'Full',
                           'axonExpt': 'standardized', 'dataTable':
                           'data_XM13A_nacncoop_normal', 'SGCmodelType':
                           'cochlea', 'runProtocol': 'runANPSTH',
                           'displayStyle': 'cylinders', 'displayMechanism':
                           'klt', 'displayMode': 'None', 'displayscale': False,
                           'inputPattern':
                            None, 'soundtype': 'tonepip', 'checkcommand': False, 
                           'testsetup': False, 'configfile':
                           'singles_multisite_parallel.toml', 'dB': 30.0, 'F0':
                           16000.0,
                            'nReps': 25, 'seed': 9, 'SRType': 'HS',
                            'ANSynapseType': 
                           'multisite', 'ANSynapticDepression': 0, 'pip_start':
                           0.2, 'pip_duration': 1.0, 'pip_offduration': 0.05,
                           'fmod': 100.0, 'dmod': 0.0, 'signalToMasker': 0,
                           'CMMRmode': 'CMR', 'Spirou': 'all', 'soma_inflation':
                           1.0, 'soma_autoinflate': False, 'dendrite_inflation':
                           1.0, 'dendrite_autoinflate': False,
                           'dendrite_fromsoma': False, 'ASA_fromsoma': False,
                           'vstimHolding': -80.0, 'tagstring': None,
                           'AMPAScale': 1.0, 'all_modes': False, 'sequence': '',
                           'plotFlag': False, 'nWorkers': 8, 'Parallel': True,
                           'auto_initialize': False, 'save_all_sections': False,
                           'verbose': False, 'gif_i0': 0.0, 'gif_sigma': 0.2,
                           'gif_fmod': 0.2, 'gif_tau': 3.0, 'gif_dur':
                            10.0, 'gif_skew': 0.0}


                   dBspl = 30.0
               dataTable = data_XM13A_nacncoop_normal
                 datestr = 2021-08-03-12:53:12
            dendriteExpt = Full dendriteMode = normal
                      dt = 0.025
                 elapsed = 0.0
                   files =
                   ['/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c05/Simulations/AN/runANPSTH-all-2021-08-03.12-53-12/AN_Result_VCN_c05_inp=self_XM13A_nacncoop_II_HF=VCN_c05_Full_MeshInflate_standardized_axon_normal_delays_multisite_025_tonepip_030dB_16000.0_HS.p']
                filetype = D
                    fmod = 100.0
                 hocfile = VCN_c05_Full_MeshInflate_standardized_axon.hoc
               modelName = XM13A_nacncoop modelType = II
                   nReps = 25
                  pipDur = 1.0
             runProtocol = runANPSTH
         simulation_path = /Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells
               soundType = tonepip
                  stimVC = {'pulse': array([-20., -10.,   0.,  10.,  20.,  30.,
                  40.,  50.,  60.,  70.,  80.,
        90., 100., 110., 120.])}
       synapseExperiment = all
             synapsetype = multisite temperature = 37
               threshold = None
      vcnmodel_code_hash = b'4e287e3ba6ead8037ebfd6fe19fa717be2b54230'
[Note that this file also holds the vcnmodel and cnmodel code hashes.] Generated
by vcnmodel/util/inspect_simulation_file.py



File List: (Appended automatically by the vcnmodel/util/make_data_dirs.py,
write_the_readme function)
"""
class BuildDataSet():
    def __init__(self, testmode:bool=False):
        self.testmode = testmode
        self.ncopy = 0

    def do_copy(self, src, dest):
        if self.testmode:
            cpc = 'm'
            if not dest.is_file() or (src.stat().st_size != dest.stat().st_size) or (src.stat().st_mtime > dest.stat().st_mtime):
                cprint(cpc, f"        Would be Copying:")
                cprint(cpc, f"           from: {str(src.parent):s}")
                cprint(cpc, f"               + {str(src.name):s}" )           
                cprint(cpc, f"             to: {str(dest):s}\n")
            else:
                cpc = "y"   
                cprint(cpc, f"        File exists looks the same, no copy:")
                cprint(cpc, f"           from: {str(src.parent):s}")
                cprint(cpc, f"               + {str(src.name):s}" )     
        else:
            cpc = 'g'
            # first check that file does not exist and if it does,
            # check that the dest date is more recent than dest, otherwise

            if not dest.is_file() or (src.stat().st_size != dest.stat().st_size) or (src.stat().st_mtime > dest.stat().st_mtime):
                shutil.copy2(src, dest)
                cprint(cpc, f"        Copying:")
                cprint(cpc, f"           from: {str(src.parent):s}")
                cprint(cpc, f"               + {str(src.name):s}" )           
                cprint(cpc, f"             to: {str(dest):s}\n")
            else:
                cpc = "y"   
                cprint(cpc, f"        File exists looks the same, no copy:")
                cprint(cpc, f"           from: {str(src.parent):s}")
                cprint(cpc, f"               + {str(src.name):s}" )           

        self.ncopy += 1

    def mk_cellid(self, cellN):   
        if isinstance(cellN, int):
            return f"VCN_c{cellN:02d}", cellN
        elif isinstance(cellN, str):
            return f"VCN_c{int(cellN[0]):02d}", int(cellN[0])
        
    def write_the_readme(self):
        with open(Path(simpath, "README.txt"), "wb") as fh:
            fh.write(bytearray(readmefilecontents.encode('utf-8')))
        # add the full file tree to the readme.
        #r = subprocess.run(["./vcnmodel/util/tree.sh", simpath], capture_output=True)
        r = ptree.ptree(simpath)
        with open(Path(simpath, "README.txt"), "a") as fh:
            #fh.write(r.stdout.decode())
            fh.write(r)

    def copy_morphology(self, allcells, sdirs, cdirs):
        for cn in allcells:
            cell_id = f"VCN_c{cn:02d}"
            source_morph = list(Path(sdirs[cn], "Morphology").glob("*.hoc"))
            source_morph.extend(list(Path(sdirs[cn], "Morphology").glob("*.swc")))
            dest_morph = Path(cdirs[cn], "Morphology")
            for f in source_morph:
                dest_file = Path(dest_morph, f.name)
                print(f"File to copy: {str(dest_file):s}")
                if not dest_file.is_file():  # only copy if not already there
                    print("copying ", f, " to ", dest_morph)
                    self.do_copy(f, dest_morph)
                else:
                    print(f"Morphology file {cell_id:s}  {str(f):s}  already present")
                    print("    dest file was: ", dest_file)

    def copy_simresults(self, FD, simcopy:bool=True, sdirs:List=[], cdirs:List=[]):
        for fig in FD.all_figures:
            figd = FD.all_figures[fig]
            if fig in ["IV_ex", "IV_all", "VC_ex"]:
                for cell_id in list(figd.keys()):
                    filelist = []
                    cell = self.mk_cellid(cell_id)
                    sourcedir = Path(sdirs[cell_id], "Simulations", fig[:2])
                    targetdir = Path(cdirs[cell_id], "Simulations", fig[:2])
                    print("\n", "=" * 80, "\n", fig, " target dir: ", targetdir)
                    for expt in figd[cell_id].keys():  # for each expt in the list
                        sourcefile = Path(sourcedir, figd[cell_id][expt])
                        sourcefilepkl = Path(sourcefile.with_suffix(sourcefile.suffix + ".pkl"))
                        if expt.startswith("Z_"):  # capture impedance runs
                            sourcefile = Path(Zin_source_dir, figd[cell_id][expt])
                            if sourcefile.is_file and sourcefile.suffix == ".pkl":
                                filelist.append({"src": sourcefile, "dest": Zin_dest_dir})
                                # self.do_copy(sourcefile, Zin_dest_dir)

                        elif sourcefile.is_dir() and simcopy:  # make and fill the subdirectory
                            targetsubdir = Path(targetdir, sourcefile.name)
                            targetsubdir.mkdir(exist_ok=True)
                            sourcefiles = sourcefile.glob("*.p")
                            for f in sourcefiles:
                                filelist.append({"src": f, "dest": targetsubdir})
                                # self.do_copy(f, targetsubdir)
                        elif sourcefilepkl.is_file():
                              filelist.append({"src": sourcefilepkl, "dest": targetdir})
                              # self.do_copy(sourcefilepkl, targetdir)
                    # pool = multiprocessing.Pool(processes=16)
                    for file in filelist:
                        self.do_copy(file['src'], file['dest'])
                    #     pool.apply_async(self.do_copy, args=(file['src'], file['dest'],))
                    # pool.close()
                    # pool.join()
            elif fig in [
                "AN_PSTH",
                "AN_revcorr",
                "AN_ex_revcorr",
                "AN_efficacy",
                "AN_SAC",
                "AN_SAM_SAC",
                "AN_BC_09_Pruned",
                # "AN_VS_15", # removed for size
              #  "AN_VS_30",  # removed for size
                "AN_VS_15_BC09",
            ]:
                for cell_idx in list(figd.keys()):
                    cell, cell_id = self.mk_cellid(cell_idx)
                    sourcedir = Path(sdirs[cell_id], "Simulations", fig[:2])
                    targetdir = Path(cdirs[cell_id], "Simulations", fig[:2])
                    print("\n", "=" * 80, "\n", fig, " target dir: ", targetdir)
                    if isinstance(figd[cell_idx], list):
                        expts = figd[cell_idx]
                    elif isinstance(figd[cell_idx], dict):
                        expts = list(figd[cell_idx].keys())
                    filelist = []
                    for expt in expts:  # for each expt in the list
                        print("   Simulation: ", expt)
                        sourcefile = Path(sourcedir, expt)
                        sourcefilepkl = Path(sourcefile.with_suffix(sourcefile.suffix + ".pkl"))
                        if sourcefile.is_dir() and simcopy:  # make and fill the subdirectory
                            targetsubdir = Path(targetdir, sourcefile.name)
                            targetsubdir.mkdir(exist_ok=True)
                            sourcefiles = list(sourcefile.glob("*.p"))
                            for sourcef in sourcefiles:
                                print("        sourcefile: ", sourcef.parent)
                                print("                  + ", sourcef.name)
                            for f in sourcefiles:
                                filelist.append({'src':f, "dest": targetsubdir}) 
                                # self.do_copy(f, targetsubdir)
                        elif sourcefilepkl.is_file():
                            self.do_copy(sourcefilepkl, targetdir)
                    for file in filelist:
                        self.do_copy(file['src'], file['dest'])
                    # pool = multiprocessing.Pool(processes=16)
                    # for file in filelist:
                    #     pool.apply_async(self.do_copy, args=(file['src'], file['dest'],))
                    # pool.close()
                    # pool.join()
    
        # intermediate results:
        # SAM clicks, SAM
        target1 = Path(intermediate, 'SAC_Results_Clicks.pkl')
        source1 = Path(intermediate_source, 'SAC_Results_Clicks.pkl')
        self.do_copy(source1, target1)
        target2 = Path(intermediate, 'SAC_Results_SAM.pkl')
        source2 = Path(intermediate_source, 'SAC_Results_SAM.pkl')
        self.do_copy(source2, target2)
        source_VS15 = Path(config['codeDirectory'], 'VS_data_15dB.py')
        target_VS15 = Path(intermediate, 'VS_data_15dB.py')
        self.do_copy(source_VS15, target_VS15)
        source_VS15_BC09 = Path(config['codeDirectory'], 'VS_data_15dB_BC09.py')
        target_VS15_BC09 = Path(intermediate, 'VS_data_15dB_BC09.py')
        self.do_copy(source_VS15_BC09, target_VS15_BC09)
        source_VS30 = Path(config['codeDirectory'], 'VS_data_30dB.py')
        target_VS30 = Path(intermediate, 'VS_data_30dB.py')
        self.do_copy(source_VS30, target_VS30)

def test_tree():
    r = subprocess.run(["./vcnmodel/util/tree.sh", simpath], capture_output=True)
    x = str(r.stdout).split('\\n')
    for l in x:
        print(l)

def main():

    allcells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
    cdirs = {}
    sdirs = {}
    # make the directory structure for all the cells
    for cn in allcells:
        cell_id = f"VCN_c{cn:02d}"
        celldir = Path(simulations, cell_id)
        cdirs[cn] = Path(simulations, celldir)
        cdirs[cn].mkdir(exist_ok=True)
        for ddir in ["Simulations/AN", "Simulations/IV", "Simulations/VC", "Morphology"]:
            Path(cdirs[cn], ddir).mkdir(exist_ok=True, parents=True)
        sdirs[cn] = Path(sourcepath, "VCN_Cells", cell_id)

    # for p in list(simpath.glob("*/**")):
    #     print(p)
    MD = BuildDataSet(testmode=False)
    MD.copy_morphology(allcells, sdirs=sdirs, cdirs=cdirs)
    MD.copy_simresults(FD, simcopy=True, sdirs=sdirs, cdirs=cdirs)
    # finally, complete the README.txt file
    MD. write_the_readme()
    cprint("g", f"Copied {MD.ncopy:d} files")

if __name__ == "__main__":
    main()

