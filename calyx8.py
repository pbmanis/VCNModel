__author__ = "Paul B. Manis"
# #
# # Analyze voltage distribution in reconstructed MNTB Calyx of Held
# # Paul B. Manis, Ph.D.
# # Dept. Otolaryngology/HNS
# # UNC Chapel Hill
# # pmanis@med.unc.edu
# # July, 2007
# # 
# # October, November 2007 revisions:
# # calyx5.hoc:
# # changed rmp to -80, adjusted Ek to -85. 
# # 11/1/2007 - fixed code for makring swellings.
# # added vax (voltage in axon) as a separate output file... axon[0](0.5) for comparison
# #
# # 11/1/07 - also added GBC model to drive with axon. 
# # 11/5/2007 - separated file into functional groups of code for easier maintenance.
# # 11/6/07 - Added calyx_tune.hoc, a routine to automatically set the densities
# # of the conductances based on some target values from experiments. Does this by
# # simulating a voltage clamp experiment at the base of the calyx under the 
# # current conditions, and running voltage step protocols. The conductance
# # necessary to meet the target currents are then predicted.
# # 11/6/07 - incorporated Ih current, using the Rothman&Manis 2003 descrption
# # 11/10/07 - split parent axon at 30 um from calyx ; distal from calyx is
# # the axon, proximal is a new type "heminode" for insertion of dense Na ch.
# # (Leao et al, 2005; Lindgren and Moore, 79)
# 
# Adapted into python from hoc code, 2013.
# This is a reduced version -
# reads morphology, decorates membrane segments, runs a current or voltage clamp and monitors calcium
# added connections to a postsynaptic neuron, using stochastic synapses, state models of receptors, etc
# as adapted from Xie and Manis (2013).
# Makes use of neuron/hoc, so be careful!

import os, sys, time
import os.path
import pickle
import neuron as h
from neuron import *
import numpy as np
import scipy as sp
import nrnlibrary
from pylibrary.Params import Params
import matplotlib.pylab as MP
import pprint
import re
import calyxPlots_pg
import time
import csv
import hocRender as hr
import pyqtgraph as pg


# GBCFLAG controls whether we use a cut axon or a GBC soma with axon (not actually implemented in this version)
GBCFLAG = 0 # if 0, IC is in cut axon near calyx; otherwise it is in the GBC soma.

topofileList = ["Calyx-S53Acvt3.hoc", "Calyx-68cvt2.hoc"]
monitorSections = {
    "Calyx-S53Acvt3.hoc":  [[0, 149, 149, 5,   88,  130, 117], [0,  0,  1.0, 0.5, 0.5, 0.5, 0.5]],
    "Calyx-68cvt2.hoc": [[0, 85, 85, 5, 26, 40, 67], [0,  0,  1.0, 0.5, 0.5, 0.5, 0.5]]
}
selectedFile = 0 # in reference to topofileList
renderFlag = False
runFlag = True
recordSwellings = True


class calyx8():
    def __init__(self, argsin = None):
        print '%s: init' % (self.__class__.__name__)
        h.load_file("stdrun.hoc")
        h.load_file(os.getcwd()+"/custom_init.hoc") # replace init with one that gets closer to steady state
        #h.load_file(os.getcwd()+"/mesh.hoc")
        # use the Params class to hold program run information and states
        # runInfo holds information that is specific the the run - stimuli, conditions, etc
        # DO NOT add .p to filename...
        self.runInfo = Params(fileName='Normal', runName='CalyxTest', manipulation="Canonical", folder="Canonical",
                         preMode = "cc", postMode = 'vc',
                         pharmManip = {'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False, 'XE': False},
                         TargetCellName='MNTB',
                         TargetCell = None,
                         newCm = 1.0,
                         newRa = 100.0, # // changed 10/20/2007 to center in range')
                         newg_leak = 0.000004935,
                         celsius = 37, #' // set the temperature.
                         eK_def = -85, eNa_def = 50,
                         ca_init = 70e-6, # free calcium in molar
                         v_init = -80, # mV
                         nStim = 1,
                         stimFreq = 200., # hz
                         stimInj = 2.0, # nA
                         stimDur = 0.5, # msec
                         stimDelay = 2.0,# msec
                         stimPost = 3.0, # msec
                         runTime = time.asctime(), # store date and time of run
                         inFile = None, # 'ANFiles/AN10000Hz.txt', # if this is not None, then we will use these spike times...
                         inFileRep = 1, # which rep to use (or array of reps)
                         spikeTimeList = {}, # Dictionary of spike times
                    )

        if self.runInfo.inFile is None:
            self.runInfo.tstop = 1000.0*(self.runInfo.nStim-1)/self.runInfo.stimFreq + self.runInfo.stimPost + self.runInfo.stimDelay
        else:
            maxt = 0.
            with open(self.runInfo.inFile, 'r') as csvfile:
                spks = csv.reader(csvfile, delimiter=',')
                for i, row in enumerate(spks):
                    if i == 0: # ture first line
                        print row
                        maxt = float(row[1])*1000
                        reps = int(row[0])
                        continue
                    if int(row[0]) in self.runInfo.spikeTimeList.keys():
                        self.runInfo.spikeTimeList[int(row[0])].append(float(row[1])*1000.)
                    else:
                        self.runInfo.spikeTimeList[int(row[0])] = [float(row[1])*1000.]
                # for n in self.runInfo.spikeTimeList.keys():
                #     if np.max(self.runInfo.spikeTimeList[n]) > maxt:
                #         maxt = np.max(self.runInfo.spikeTimeList[n])
                # print self.runInfo.spikeTimeList
                self.runInfo.tstop = maxt
        # modelParameters holds information regarding the cell structure itself for this run

        self.modelPars = Params(calyxNames=['parentaxon', 'axon', 'heminode', 'synapse', 'stalk', 'branch', 'neck', 'swelling', 'tip'], # names of calyx morphological parts
                                      calyxColors={'axon': 'red', 'heminode': 'green', 'stalk':'yellow', 'branch': 'green', 'neck': 'blue',
            'swelling': 'magenta', 'tip': 'k', 'parentaxon': 'red', 'synapse': 'cyan'},
                                      mechNames = {},
                                        # variables to control electrode positions
                                        axonselect = 1,
                                        iclocation = 0.5,
                                        vclocation = 1.0,
                                        topofile = topofileList[selectedFile],
                                    # define model parameters. These are straight from EIModel4.py, even though our only target is bushy
                                    # ALL parameters should be defined here that do not depend on calculation or
                                    # which are not modified by input flags

                                    # Auditory nerve synapse information:
                                    # We use cell parameter dictionaries, to define the inputs to primary target cell types:
                                    # sets default values, but these can be overridden

                                        Exc_Thresh = -20,
                                        AN_gEPSC =   {'stellate' : 5, 'bushy': 40, 'test':   2, 'DStellate':   0, 'MNTB': 20},
                                        spike_thr =  {'stellate': -30, 'bushy': -30, 'test': -30, 'DStellate': -20, 'MNTB': -30},
                                        AN_conv =    {'stellate' :  5, 'bushy':   3, 'test':   1, 'DStellate':  25, 'MNTB': 1}, # bushy was 2
                                        AN_zones =   {'stellate' :  5, 'bushy':  90, 'test':   1, 'DStellate':   5, 'MNTB': 1}, # bushy ANZones was 90...
                                        AN_erev =    {'stellate':   7, 'bushy':   7, 'test':   7, 'DStellate':   7, 'MNTB': 7},
                                        AN_gvar =    {'stellate': 0.3, 'bushy': 0.3, 'test': 0.3, 'DStellate':  0.3, 'MNTB': 0.3}, # std of conductance
                                        AN_delay =   {'stellate': 0.0, 'bushy': 0.0, 'test': 0.0, 'DStellate':  0.3, 'MNTB': 0.3}, # netcon delay (fixed)
                                        AN_lat_Flag = {'stellate': 0,    'bushy': 1,    'test': 0,    'DStellate':  0, 'MNTB': 1}, # 1 to allow latency shift during train
                                        AN_lat_A0 =   {'stellate': 0.14, 'bushy': 0.140, 'test': 0.25, 'DStellate':  0.15, 'MNTB': 0.14}, # magnitude of latency change, if allowed
                                        AN_lat_tau =  {'stellate': 21.0, 'bushy': 21.5, 'test': 21.0, 'DStellate':  21.0, 'MNTB': 21.5}, # rate at which latency changes, if allowed
                                        AN_latency =  {'stellate': 0.5,  'bushy': 0.5,  'test': 0.5,  'DStellate':  0.5, 'MNTB': 0.5}, # starting latency
                                        AN_relstd_Flag = {'stellate':   1, 'bushy': 2,    'test': 0,   'DStellate':  0, 'MNTB': 1}, # flag controlling std shift during train (1 to enable)
                                        AN_relstd =      {'stellate': 0.05, 'bushy': 0.05, 'test': 0.05, 'DStellate':  0.05, 'MNTB': 0.05}, # baseline release std
                                        AN_relstd_A0 =   {'stellate': 0.0, 'bushy': 0.05, 'test': 0.0, 'DStellate':  0.0, 'MNTB': 0.05}, # magnitude of std change during train, if allowed
                                        AN_relstd_tau =  {'stellate': 35.0, 'bushy': 35.0, 'test': 35.0, 'DStellate':  35.0, 'MNTB': 35.0}, # time constant of std change during train, if allowed

                                        # define inhibitory input from D-stellate cell onto a target cell
                                        DS_gI = 50.,
                                        Inh_Thresh = -30.0 , # trigger release at this amplitude
                                        INH_conv =   {'stellate':   8, 'bushy':  10, 'test':   1, 'DStellate':    2}, # normally 10
                                        INH_zones =  {'stellate':   5, 'bushy':  10,  'test':   1, 'DStellate':    5}, # normally 10
                                        INH_erev =   {'stellate': -70, 'bushy': -70, 'test': -70., 'DStellate': -70},
                                        INH_gvar =   {'stellate': 0.3, 'bushy': 0.3, 'test': 0.3, 'DStellate':  0.3},
                                        INH_delay =  {'stellate': 0.0, 'bushy': 0.0, 'test': 0.0, 'DStellate':  0.0}, # netcon delay (fixed)
                                        INH_lat_Flag = {'stellate': 0,    'bushy': 0,    'test': 0,    'DStellate':  0}, # 1 to allow latency shift during train
                                        INH_latency = {'stellate': 0.2, 'bushy': 0.2, 'test': 0.2, 'DStellate':  0.2}, # latency for stochastic synapses
                                        INH_relstd_Flag = {'stellate':   1, 'bushy': 2,    'test': 0,   'DStellate':  0}, # flag controlling std shift during train (1 to enable)
                                        INH_relstd = {'stellate': 0.2, 'bushy': 0.2, 'test': 0.3, 'DStellate':  0.2}, # release sigma
                                    )

        # define which sections (and where in the section) to monitor
        # currents, voltages, and ion concentrations can be monitored
        # this should be changed on a per-calyx configuration, but right now we use the tables above.
        #
        self.modelPars.monsecs = monitorSections[self.modelPars.topofile][0]
        self.modelPars.monpos =  monitorSections[self.modelPars.topofile][1]
        self.modelPars.axonnode = self.modelPars.monsecs[1]  #the heminode section (this must be specified in this way in the topofile dict)

        defPars = self.restoreDefaultConductances() # set the default conductance for each section from scratch
        self.modelPars.gPars = self.setDefaultConductances(defPars=defPars) # assign the current "defaults" to the relevant variables...

        self.rundone = 0
        self.thisrep = 1
        self.Mesh_th = 0.1 # minimum dlambda

        self.dia_thresh = 0.25 # minimum diameter of any part... adjustable in menu
        self.vsaveflag = 0
        self.isaveflag = 0
        self.mark_sw = 1
        self.mark_br = 1
        self.mark_st = 1
        self.mark_tp = 1
        self.mark_nk = 1

        # set up arrays for measurement of voltage, calcium, etc.
        self.nswel = 200
        self.nactual_swel = 0
        h.load_file(1, self.modelPars.topofile) # load the morphology file
        # morphology file should include it's own definition call ("celldef")
        self.CalyxStruct={key: None for key in range(self.modelPars.AN_conv[self.runInfo.TargetCellName])} # create a dictionary to access the calyx parts
        self.clist = {}
        for input in range(self.modelPars.AN_conv[self.runInfo.TargetCellName]):
            self.CalyxStruct[input] = {key: None for key in self.modelPars.calyxNames}
            for name in self.modelPars.calyxNames: # for each part
                x=eval("h.%s" % (name)) # find the associated variable
                self.CalyxStruct[input][name]  = list(x) # and populate it with the pointers to the parts
            # h.load_file(1, "calyx_morpho.hoc")
            # h.load_file(1, "calyx_shape.hoc")
            self.biophys(runInfo = self.runInfo, modelPars = self.modelPars,
                     CalyxStruct = self.CalyxStruct[input], createFlag = True)
        # get the map of inputs from the swellings to the axon id's for later.
        self.swellAxonMap = []
        nSwellings = len(self.CalyxStruct[0]['swelling'])
        input = 0
        for swellno in range(nSwellings):
            swelling = self.getAxonSec('swelling', swellno, input)
            self.swellAxonMap.append(swelling) # implicit order

        # now for the postsynaptic cell - just a basic bushy cell, even if it is an "MNTB" model.

        (self.TargetCell, [self.initseg, self.axn, self.internnode]) = nrnlibrary.Cells.bushy(debug=False, ttx=False,
                                        message=None,
                                        nach='jsrnaf',
                                        species='cat', axon=False, dendrite=False,
                                        newModFiles=False, pump=False)

        self.modelPars.stochasticPars = Params(
            LN_Flag = self.modelPars.AN_relstd_Flag[self.runInfo.TargetCellName],
            LN_t0 = 0.,
            LN_A0 = self.modelPars.AN_relstd_A0[self.runInfo.TargetCellName],
            LN_tau = self.modelPars.AN_relstd_tau[self.runInfo.TargetCellName],
            LN_std = self.modelPars.AN_relstd[self.runInfo.TargetCellName],
            Lat_Flag = self.modelPars.AN_lat_Flag[self.runInfo.TargetCellName],
            latency = self.modelPars.AN_latency[self.runInfo.TargetCellName],
            Lat_t0 = 0.,
            Lat_A0 = self.modelPars.AN_lat_A0[self.runInfo.TargetCellName],
            Lat_tau = self.modelPars.AN_lat_tau[self.runInfo.TargetCellName],
            delay=self.modelPars.AN_delay[self.runInfo.TargetCellName])

        self.modelPars.NMDARatio = 0.0
        # make the incoming synapse to the target cell
        # to mimic MNTB cells, we source each presynaptic V from the swellings
        # (implicit assumption is that swellings are synapses).
        # nFibers is set to 1
        # nRZones is set to whatever (let's say, 2 per swelling?)

        results = Params() # put the results in a Param structure as well.

        #print "GBC excitation to target"

        self.synapse = []
        self.coh=[]
        self.psd = []
        self.cleft = []
        self.nc2 = []
        # count swellings first
        nSwellings = len(self.CalyxStruct[0]['swelling']) # note that at present, all inputs have same structure.
        # now modify nr zones if it is an mntb neuron - try to make 400 zones in the number of swelllings
        # magic number of 400 total
        nZonesPerSwelling = int(400./nSwellings)
        print 'Swellings: %d  zones per swelling: %d' % (nSwellings, nZonesPerSwelling)
        #if self.runInfo.TargetCellName == 'MNTB': # only in case of the MMTB
        #    self.modelPars.AN_zones[self.runInfo.TargetCellName] = nZonesPerSwelling
        # Insert synapses in the swellings (but no where else in the calyx)
        for input in range(self.modelPars.AN_conv[self.runInfo.TargetCellName]):
            for swellno, swell in enumerate(self.CalyxStruct[input]['swelling']):
    #        if len(self.CalyxStruct['swelling']) > 1:
                swelling = self.getAxonSec('swelling', swellno, input)
                (calyx, coh, psd, cleft, nc2, par) = nrnlibrary.Synapses.stochastic_synapses(h, parentSection = self.CalyxStruct[input]['axon'][swelling],
                                                                                             targetcell=self.TargetCell,
                                                                                             cellname = self.runInfo.TargetCellName,
                        nFibers = self.modelPars.AN_conv[self.runInfo.TargetCellName], nRZones = self.modelPars.AN_zones[self.runInfo.TargetCellName],
                        message='creating stochastic multisite synapse',
                        thresh = self.modelPars.Exc_Thresh, psdtype = 'ampa', gmax=self.modelPars.AN_gEPSC[self.runInfo.TargetCellName]*1000.,
                        gvar = self.modelPars.AN_gvar[self.runInfo.TargetCellName],
                        eRev = self.modelPars.AN_erev[self.runInfo.TargetCellName], NMDARatio = self.modelPars.NMDARatio,
                        debug = False, Identifier = 1, stochasticPars = self.modelPars.stochasticPars, calciumPars=None)
                self.synapse.append(calyx)
                self.coh.append(coh)
                self.psd.append(psd)
                self.cleft.append(cleft)
                self.nc2.append(nc2)

            #self.modelPars.stochasticPars.show()
        if renderFlag:
            pg.mkQApp()
            pg.dbg()
            render = hr.hocRender(h)
            render.draw_model(modes=['blob'])
            render.getSectionLists(self.modelPars.calyxColors.keys())
            render.paintSectionsByDensity(self.modelPars.calyxColors, self.modelPars.mechNames['CaPCalyx'])
            render.show()

            exit()

        if runFlag:
            self.runModel(runInfo = self.runInfo, modelPars = self.modelPars)

    def restoreDefaultConductances(self, defPars=None):
        """ default conductance values determined by voltage clamp test on model calyx
            canonical settings. This routine restores those values.
            Arguments:
                defPars: if None, we create the result
                 and set tho defaults
                if defPars is not none or empty, se just pass it on.
            Returns:
                defpars, set to their default values
            Side effects: None
        """

        if defPars is None or defPars is {}:
            defPars = Params(gna_def = 0.45,
                             eNa_def = 50.0,
                             eK_def = -85.0,
                             glvk_def = 0.04,
                             ghvkax_def = 0.02,
                             ghvk_def = 0.005,
                             gca_def = 0.003,
                             gh_def = 0.000950,
                             gleak_def = 0.000004935,
                             capump_def = 3e-14, # mol/sec
                             )
            defPars.gPars = self.setDefaultConductances(defPars=defPars)
        return defPars
    
    def setDefaultConductances(self, defPars = None, gPars = {}):
        """ For each type of segment, initialize setting by the default conductances.
            Inputs: defPars or None to use standards
                gPars, a dictionary with keys = parts of calyx, that is filled with the conductances at
                at each location
            Outputs:
                gPars, the conductance parametere list
            Side Effects:
                None
        """
        if defPars is None or defPars == {}:
            raise Exception (' setDefaultConductances requires defPars!')
        for name in self.modelPars.calyxNames:
            gPars[name] = Params(gNabar=defPars.gna_def, ENa=defPars.eNa_def,    gKLbar =defPars.glvk_def, EK = defPars.eK_def,
                gKHbar = defPars.ghvkax_def, gCabar = defPars.gca_def, gHbar = defPars.gh_def, gCaPump = defPars.capump_def).todict()
            # exceptions:
            if name in ['axon', 'heminode']:
                gPars[name]['gCabar'] = 0.0 # no calcium channels in the axon or heminode
            if name in ['axon']:
                gPars[name]['gNabar'] = defPars.gna_def/5.0 # reduce sodium channel conductance in axon
                gPars[name]['gCabar'] = 0.0 # no calcium channels in the axon or heminode
            if name in ['neck', 'branch', 'tip', 'stalk'] :
                gPars[name]['gCabar'] = 0.0
            if name in ['neck', 'branch', 'tip', 'swelling']:
                gPars[name]['gNabar'] = 0.0 # no sodium or lva K channels in the elements of the terminal
                gPars[name]['gKLbar'] = 0.0
        return gPars


    def biophys(self, runInfo=None, modelPars=None, CalyxStruct = None, createFlag = True):
        """
        Inputs: run parameter structure, model parameter structure
        Outputs: None
        Action: Channel insertion into model
        Side Effects:
            Sets conductances in every different kind of section
            Does not update any class variables (via self).

        original hoc code: Paul B. Manis, Ph.D.
        25 Sept. 2007
        Modified to use gca for HH formulation of calcium current
        14 Oct 2007
        converted for Python, 17 Oct 2012 (PB Manis)
        """
        if modelPars is None or runInfo is None:
            raise Exception('calyx::biophys - no parameters or info passed!')
        if createFlag : # first time through insert basics only
            for sec in CalyxStruct['axon']:
                sec.insert('leak')
                sec().g_leak = runInfo.newg_leak
                sec.Ra = runInfo.newRa
                sec.cm = runInfo.newCm
                e_leak = -50
                sec.insert('capmp')
                if 'capmp' not in modelPars.mechNames.keys():
                    modelPars.mechNames['capmp'] = ['capmp', 'pump0']
                pcabar_cachan = 2.5e-5
                sec().cao = 1.0 #insert capump
                sec.cai = 70e-6
                sec.pump0_capmp = 0
 
        for sec in CalyxStruct['axon']: #  always keep these base parameters updated
            sec.g_leak = runInfo.newg_leak
            sec.Ra = runInfo.newRa
            sec.cm = runInfo.newCm
        # for each of the parts of the calyx, insert the appropriate conductances
        # into the model, and set the conductance level and the Nernst potential.
        for name in modelPars.calyxNames:
            gp = modelPars.gPars[name]
            for sec in CalyxStruct[name]:
               # print 'for %s  sec = ' % (name), sec
                sec.insert('na')
                if 'na' not in modelPars.mechNames.keys():
                    modelPars.mechNames['na'] = ['na', 'gnabar']
                if runInfo.pharmManip['TTX']:
                    sec().gnabar_na = 0.0
                else:
                    sec().gnabar_na = gp['gNabar']
                sec().ena = gp['ENa']
                sec.insert('klt')
                if 'klt' not in modelPars.mechNames.keys():
                    modelPars.mechNames['klt'] = ['klt', 'gkltbar']

                if runInfo.pharmManip['DTX']:
                    sec().gkltbar_klt = 0.0
                else:
                    sec().gkltbar_klt = gp['gKLbar']
                sec().ek = gp['EK']
                sec.insert('kht')
                if 'kht' not in modelPars.mechNames.keys():
                    modelPars.mechNames['kht'] = ['kht', 'gkhtbar']

                if runInfo.pharmManip['TEA']:
                    sec().gkhtbar_kht = 0.0
                else:
                    sec().gkhtbar_kht = gp['gKHbar']
                sec.insert('ih')
                if 'ih' not in modelPars.mechNames.keys():
                    modelPars.mechNames['ih'] = ['ih', 'ghbar']

                if runInfo.pharmManip['ZD']:
                    sec().ghbar_ih = 0.0
                else:
                    sec().ghbar_ih = gp['gHbar']
                sec().eh_ih = -43

                if name is 'axon':
                    sec().cm = 0.0001 # mimic myelination
                    sec().pump0_capmp = 0
                else:
                    sec().cm = runInfo.newCm

                if name not in ['axon', 'heminode']:
                    sec.insert('CaPCalyx')
                    if 'CaPCalyx' not in modelPars.mechNames.keys():
                        modelPars.mechNames['CaPCalyx'] = ['CaPCalyx', 'gcapbar']

                    if runInfo.pharmManip['Cd']:
                        sec().gcapbar_CaPCalyx = 0.0
                    else:
                        sec().gcapbar_CaPCalyx = gp['gCabar']
                    sec().eca=43.9


    def calyx_init(self, runInfo=None, modelPars=None):
        """
        Calyx model initialization procedure:
        Set RMP to the resting RMP of the model cell.
        Make sure nseg is large enough in the proximal axon.
        Initialize the leak, and stabilize the pumps, etc.
        Inputs: runInfo and model Param structure.
        Outputs: None.
        Action: Initializes the NEURON state to begin a run.
        """

        if modelPars is None or runInfo is None:
            raise Exception('calyx::calyx_init - no parameters or info passed!')
    # First we set e_leak so that the rmp in each segment is the same
        h.finitialize(runInfo.v_init)
        i = 0
        for input in range(self.modelPars.AN_conv[runInfo.TargetCellName]):
            for sec in self.CalyxStruct[input]['axon']: # set a new value
                if i == 0:
                    i += 1
                e_newleak = sec.v + (sec.ina + sec.ik + +sec.ik_klt + sec.ik_kht + sec.ica + sec.i_ih)/(sec.g_leak*1e3) # hmmm. This is probably NOT right.
                sec().e_leak = e_newleak
        
    # based on those conditions, make sure spatial grid is fine enough for our needs
        #h.Mesh_th = self.Mesh_th
        #h('mesh_init')
        #self.mesh_init()
        h.finitialize(runInfo.v_init)
 #       h('axon[axonnode] ic.loc(iclocation)')
    
        for name in modelPars.calyxNames:
            if name in ['axon', 'heminode']:
                for sec in h.parentaxon:
                    sec.nseg = 11 

    # get ready to run the system to a stable point for the calcium pump
        j = 0
    # save the calcium pump information
        for input in range(self.modelPars.AN_conv[runInfo.TargetCellName]):
            for sec in self.CalyxStruct[input]['axon']:
                savcore = sec().cabulk_capmp
                sec().cabulk_capmp = runInfo.ca_init
                savtau = sec().tau_capmp
                sec().tau_capmp = 1e-6 # make the pump go really fast

    # starting way back in time
        h.t = -1e10
        dtsav = h.dt
        h.dt=1e9 # big time steps for slow process
        temp = h.cvode.active()
        if(temp!= 0):
            h.cvode.active(0) # turn cvode off (note, in this model it will be off because one of the mechanisms is not compatible with cvode at this time
        while(h.t <-1e9): 
            h.fadvance()

    # now restore the pump values and get the system ready for runs 
        if(temp != 0):
            h.cvode.active(1)
        h.dt=dtsav
        h.t = 0
        for input in range(self.modelPars.AN_conv[runInfo.TargetCellName]):
            for sec in self.CalyxStruct[input]['axon']:
                sec().cabulk_capmp = savcore
                sec().tau_capmp = savtau

        if(h.cvode.active()):
            h.cvode.re_init()
        else:
            h.fcurrent()
        h.frecord_init()
        h.finitialize(h.v_init)


    def getMonSec(self, monsec):
        """
        Find the section in the lists that corresponds to the selected monitor section
        Because there can be 2 pointers to the same section (based on the listed elements
        and the original "axon" elements), we have to traverse all of the non-axon elements first
        before assigning any axon elements
        Input: requested monitor section
        Output: None
        Actions/side effects: sets self.clist[monsec] to the appropriate calyx color
        """
        src = re.compile('axon\[(\d*)\]')
        for s in self.modelPars.calyxNames:
            if s in ['axon', 'synapse']:
                continue # skip the catch-all categories
            for element in self.CalyxStruct[0][s]:
                name = element.name()
#                print 's: %s  name: %s' % (s, name)
                g = src.match(name)
                axno = g.groups()[0]
                if (int(axno)) == monsec:
                    self.clist[monsec] = self.modelPars.calyxColors[s]
 #                   print 'match found for %d in %s' % (monsec, s)

    def getAxonSec(self, structure, number, input):
        """
        Get the axon section in the calyx structure at position
        Inputs:
            structure: valid name of dictionary element of calyx tags
            number: the index into the structure
        Returns: the axon number (if found).
        Side effects: None.
        """
        src = re.compile('axon\[(\d*)\]')
        assert structure in self.CalyxStruct[input].keys()
        try:
            element  = self.CalyxStruct[input][structure][number]
        except:
            raise Exception(('Failed to find %d in structure %s' % (number, structure())))
        name = element.name()
        g = src.match(name)
        axno = g.groups()[0]
        return int(axno)


    def calyxrun(self, runInfo = None, modelPars=None):
        """
        Control a single run the calyx model with updated display
        of the voltages, etc.
        Inputs: runInfo and modelPars parameter dictionaries
        Outputs: None
        Actions: displays the results
        Side Effects: A number of class variables are created and modified
            (runInfo and modelPars are not modified).
        """
        if modelPars is None or runInfo is None:
            raise Exception('calyx::calyxrun - no parameters or info passed!')

        # { local i, j, k, b
        print 'calyxrun'
        rundone = 0
        self.vec={}
        for var in ['axon', 'stalk', 'branch', 'neck', 'tip', 'swelling', 'inj', 'time', 'postsynaptic', 'i_stim0']:
            self.vec[var] = h.Vector()

        h.tstop = runInfo.tstop
        h.dt = 0.005 # force small time step. cvode is probably off.
        h.celsius = runInfo.celsius
        npts = 1 + int(h.tstop/h.dt) # number of points in a run
        self.vec2 = {} # voltage in each monitored segment
        self.ica = {} # calcium current in each monitored segment
        self.cai = {} # intracellular calcium (bulk) in each monitored segment

        self.clist = {} # color list
        segarea = {} # np.zeros(len(self.monsecs))
        segvol = {} # np.zeros(len(self.monsecs))
#        pprint.pprint(self.CalyxStruct)
        for m, monsec in enumerate(self.modelPars.monsecs):
            esection=self.CalyxStruct[0]['axon'][monsec]
            self.vec2[monsec] = h.Vector() # set up recording vectors at each site
            self.ica[monsec] = h.Vector()
            self.cai[monsec] = h.Vector()
            self.vec2[monsec].record(esection(self.modelPars.monpos[m])._ref_v, sec=esection)
            self.ica[monsec].record(esection(self.modelPars.monpos[m])._ref_ica, sec=esection)
            self.cai[monsec].record(esection(self.modelPars.monpos[m])._ref_cai, sec=esection)
            segarea[monsec] = esection.L*np.pi*esection.diam * 1e-8 # L and diam in microns, area in cm^2
            segvol[monsec] = esection.L*np.pi*(esection.diam/2.)**2 # leave in microns

            # set color list for plots - each vector is colored by the type of segment it points to
            self.getMonSec(monsec)

        if self.runInfo.postMode in ['vc', 'vclamp']:
            clampV=-60.0
            vcPost = h.SEClamp(0.5, sec=self.TargetCell)
            vcPost.dur1 = runInfo.tstop-3.0
            vcPost.amp1 = clampV
            vcPost.dur2 = 2.0
            vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
            vcPost.dur3 = 1.0
            vcPost.amp3 = clampV
            vcPost.rs = 1e-6
            self.vec['postsynaptic'].record(vcPost._ref_i, sec=self.TargetCell)
        else:
            icPost = h.iStim(0.5, sec=self.TargetCell)
            icPost.delay = 0
            icPost.dur = 1e9 # these actually do not matter...
            icPost.iMax = 0.0
            self.vec['postsynaptic'].record(self.TargetCell()._ref_v, sec=self.TargetCell)

        nactual_swel = 0 # count the number of swellings
        nconv = self.modelPars.AN_conv[runInfo.TargetCellName]
        print 'Number of converging fibers:: ', nconv
        print 'Number of synapses per swelling: ', self.modelPars.AN_zones[runInfo.TargetCellName]
        electrodesite = [None]*nconv
        istim = [None]*nconv
        for i in range(nconv):
            print 'Input 1: # swellings: ', len(self.CalyxStruct[i]['swelling'])
         #   print modelPars.axonnode
         #   print self.CalyxStruct[i]['axon'][modelPars.axonnode]
         #   print len(electrodesite)
            electrodesite[i] = self.CalyxStruct[i]['axon'][modelPars.axonnode]
            istim[i]= h.iStim(0.5, sec=electrodesite[i])
            if self.runInfo.inFile is None:
                stim={}
                stim['NP'] = runInfo.nStim
                stim['Sfreq'] = runInfo.stimFreq # stimulus frequency
                stim['delay'] = runInfo.stimDelay
                stim['dur'] = runInfo.stimDur
                stim['amp'] = runInfo.stimInj
                stim['PT'] = 0.0
                (secmd, maxt, tstims) = nrnlibrary.makestim.makestim(stim, pulsetype='square', dt = h.dt)
                istim[0].delay = 0
                istim[0].dur = 1e9 # these actually do not matter...
                istim[0].iMax = 0.0
            else:
                stim={}
                stim['NP'] = len(runInfo.spikeTimeList[runInfo.inFileRep])
                stim['Sfreq'] = runInfo.stimFreq # stimulus frequency
                stim['delay'] = runInfo.stimDelay
                stim['dur'] = runInfo.stimDur
                stim['amp'] = runInfo.stimInj
                stim['PT'] = 0.0
                #print runInfo.spikeTimeList
                stim['spikeTimes'] = runInfo.spikeTimeList[i+1]
                (secmd, maxt, tstims) = nrnlibrary.makestim.makestim(stim, pulsetype='timedSpikes', dt=h.dt)
                # MP.figure()
                # MP.plot(secmd)
                # MP.show()
                # exit()
                istim[i].delay = 0
                istim[i].dur = 1e9 # these actually do not matter...
                istim[i].iMax = 0.0

           # print 'i: %d  secmd: %d' % (i, len(secmd))
            self.vec['i_stim%d' % i] = h.Vector(secmd)

            self.vec['i_stim%d' % i].play(istim[i]._ref_i, h.dt , 0, sec=electrodesite[i])
            self.vec['axon'].record(electrodesite[i]()._ref_v, sec=electrodesite[i]) # record axon voltage
            self.vec['inj'].record(istim[i]._ref_i,  sec=electrodesite[i])
            self.vec['time'].record(h._ref_t)
       # MP.plot(np.arange(0, maxt-h.dt, h.dt), secmd)
       # MP.show()
        #self.swellAxonMap
        if recordSwellings:
            cinput = 0
            self.vswel = []
            nSwellings = len(self.CalyxStruct[cinput]['swelling']) # record from all swellings...
#            for i, swellno in enumerate(self.swellAxonMap):
            secarray = []
            for swellno in range(len(self.CalyxStruct[0]['axon'])):
                self.vswel.append(h.Vector())
                print swellno
                print self.CalyxStruct[0]['axon'][swellno]
                self.vswel[-1].record(self.CalyxStruct[0]['axon'][swellno]()._ref_v,
                                     sec=self.CalyxStruct[0]['axon'][swellno])
                secarray.append(swellno)

        #exit()
        print 'Running for: ', h.tstop
        h.run()
        npvecs={}
        for k in self.vec.keys():
            npvecs[k] = np.array(self.vec[k])

        npcai={}
        for k in self.cai.keys():
            npcai[k] = np.array(self.cai[k])
        npvec2={}
        for k in self.vec2.keys():
            npvec2[k] = np.array(self.vec2[k])
        npica={}
        for k in self.ica.keys():
            npica[k] = np.array(self.ica[k])
        runInfo.clist = self.clist
        results = Params(Sections=self.modelPars.monsecs, vec= npvecs,
                        Voltages= npvec2, ICa= npica,
                        Cai= npcai,
                        )
        dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        print("CalyxRun done\n")
        if os.path.exists(runInfo.folder) is False:
            os.mkdir(runInfo.folder)
        fn = os.path.join(runInfo.folder, runInfo.fileName+dtime+'.p')
        pfout = open(fn, 'wb')
        pickle.dump({'runInfo': runInfo.todict(),
                     'modelPars': modelPars.todict(),
                     'Results': results.todict()}, pfout)
        pfout.close()

        if recordSwellings:
            fns = os.path.join(runInfo.folder, runInfo.fileName+'_swellings_'+dtime+'.p')
            pfout = open(fns, 'wb')
            pickle.dump({'swellings': secarray,
                         'time': np.array(self.vec['time']),
                         'data': np.array(self.vswel)}, pfout)
            pfout.close()
            plot = pg.plot()
            vs = np.array(self.vswel)
            ts = np.array(self.vec['time'])
            for i in range(len(self.swellAxonMap)):
                plot.plot(ts, vs[i])

        calyxPlots_pg.plotResults(results.todict(), runInfo.todict())



    
    # #compartmentalization function.
    # #To have compartments with a lambda length < $2 
    # #at $1 (Hz) - the value for $2 should be less than 0.1-0.3.
    # #To speedup the demo simulation a value of 0.4@100Hz is used.
    # #Slightly different values for the conductances would be needed
    # #to obtain the same results when using a different compartmentalization.
    # #Contact the author for more information.
    # 
#     def mesh_init(self):
#          h('fast_mesh(500, %f)' % self.Mesh_th)
#          #axon.nseg = 1
#          #print("Number of Compartments: %d\n" % h.nr)
#      
#          # Set up origin in soma
# #        distance(self.axon)
# #         axon
#     # 
    # 
    # 
    def runModel(self, runInfo=None, modelPars=None):
        """
        Main code to be executed when calyxN.hoc is called
        """
        if modelPars is None or runInfo is None:
            raise Exception('calyx::biophys - no parameters or info passed!')
#        h('measure()') # get distances
        defPars = self.restoreDefaultConductances() # canonical conductances...
        self.setDefaultConductances(defPars = defPars) # initialize the conductances
        # print self.CalyxStruct.keys()
        # load anciallary functions that may need all the definitions above
        # h.load_file(1, "calyx_tune.hoc") # requires clamps to be inserted..
        self.calyx_init(runInfo = runInfo, modelPars = modelPars) # this also sets nseg in the axon - so do it before setting up shapes
        # h.shape_graphs() # must happen AFTER calyx_init
        self.calyxrun(runInfo = runInfo, modelPars = modelPars)

