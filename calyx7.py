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
# Adapted into python from hoc code. 
# reduced version - 
# reads morphology, decorates membrane segments, runs a current clamp and monitors calcium 

import os, sys
import neuron as h
from neuron import *
import numpy as np
import scipy as sp
import nrnlibrary
import matplotlib.pylab as MP
import pprint
import re

# GBCFLAG controls whether we use  is used or not:
GBCFLAG = 0 # if 0, IC is in cut axon near calyx; otherwise it is in the GBC soma.

#topofile = "calyx-S53Acvt2.hoc" # this is the membrane topology filename.
topofileList = ["Calyx-S53Acvt3.hoc", "Calyx-68cvt2.hoc"]
monitorSections = {
    "Calyx-S53Acvt3.hoc":  [[0, 149, 149, 5,   88,  130, 117], [0,  0,  1.0, 0.5, 0.5, 0.5, 0.5]],
    "Calyx-68cvt2.hoc": [[0, 85, 85, 5, 26, 40, 67], [0,  0,  1.0, 0.5, 0.5, 0.5, 0.5]]
}

# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class calyx7():
    def __init__(self, argsin = None):
        print 'calyx7: init'
        h.load_file("stdrun.hoc")
        h.load_file(os.getcwd()+"/custom_init.hoc") # replace init with one that gets closer to steady state
        #h.load_file(os.getcwd()+"/mesh.hoc")
        self.clamp_flag = 0 # 0 = current clamp, 1 = voltage clamp
        self.clampstr  = "CC"
        self.calyxNames=['parentaxon', 'axon', 'heminode', 'synapse', 'stalk', 'branch', 'neck', 'swelling', 'tip'] # names of calyx morphological parts
        self.calyxColors={'axon': 'red', 'heminode': 'orange', 'stalk':'yellow', 'branch': 'green', 'neck': 'blue',
            'swelling': 'magenta', 'tip': 'black', 'parentaxon': 'red', 'synapse': 'cyan'}
        self.pharmManip = {'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False, 'XE': False}
        self.eK_def = -85
        self.eNa_def = 50
        self.gPars = {}
        self.restoreDefaultConductances() # set the default condcutances from scratch
        self.setDefaultConductances() # assign the current "defaults" to the relevant variables...

        self.rundone = 0
        self.thisrep = 1
        self.manipulation = "Canonical"
        self.folder = "Canonical"
        self.TargetCellName = 'MNTB'
        self.postMode = 'cc' # current clamp

        self.newCm = 1.0 #// uf/cm2 cap.
        self.newRa = 100.0 # // changed 10/20/2007 to center in range')
        self.newg_leak = 0.000004935 
        self.celsius = 37 #' // set the temperature.')
        self.ca_init = 70e-6
        self.v_init = -80

        self.nStim = 20
        self.stimFreq = 400. # hz
        self.stimInj = 2.0 # nA
        self.stimDur = 1.0 # msec
        self.stimDelay = 2.0 # msecc

        self.Mesh_th = 0.1 # minimum dlambda        
        self.tstop = 1000.0*self.nStim/self.stimFreq + 10. + self.stimDelay

        self.topofile = topofileList[0]
        # define which sections (and where in the section) to monitor
        # currents, voltages, and ion concentrations
        # this should be changed on a per-calyx configuration, but right now
        #
        # sections to monitor. These are selected manually and held in the dictionary above
        self.monsecs = monitorSections[self.topofile][0]
        self.monpos =  monitorSections[self.topofile][1]
        self.axonnode = self.monsecs[1]  #the heminode section (this must be specified in this way in the topofile dict)

        # variables to control electrode positions
        self.axonselect = 1
        self.iclocation = 0.5
        self.vclocation = 1.0
    
        self.dia_thresh = 0.25 # minimum diameter of any part... adjustable in menu
        self.vsaveflag = 0
        self.isaveflag = 0
        self.mark_sw = 1
        self.mark_br = 1
        self.mark_st = 1
        self.mark_tp = 1
        self.mark_nk = 1

    # define model parameters. These are straight from EIModel4.py, even though our only target is bushy
    # ALL parameters should be defined here that do not depend on calculation or
    # which are not modified by input flags

        # Auditory nerve synapse information:
        # We use cell parameter dictionaries, inputs to primary target cell types:
        # sets default values, but these can be overridden by command line parameters, as done below
        # self.AN_gE = 40 # bushy; default for stellate is 5
        # N.B.: self.AN_gEPSC Dstel AN input is 0 for testing TV stim simulation. Normal Value should be 15.

        self.Exc_Thresh = -20
        self.AN_gEPSC =   {'stellate' : 5, 'bushy': 40, 'test':   2, 'DStellate':   0, 'MNTB': 40}
        self.spike_thr =  {'stellate': -30, 'bushy': -30, 'test': -30, 'DStellate': -20, 'MNTB': -30}
        self.AN_conv =    {'stellate' :  5, 'bushy':   2, 'test':   1, 'DStellate':  25, 'MNTB': 1} # bushy was 2
        self.AN_zones =   {'stellate' :  5, 'bushy':  90, 'test':   1, 'DStellate':   5, 'MNTB': 4} # bushy ANZones was 90...
        self.AN_erev =    {'stellate':   7, 'bushy':   7, 'test':   7, 'DStellate':   7, 'MNTB': 7}
        self.AN_gvar =    {'stellate': 0.3, 'bushy': 0.3, 'test': 0.3, 'DStellate':  0.3, 'MNTB': 0.3} # std of conductance
        self.AN_delay =   {'stellate': 0.0, 'bushy': 0.0, 'test': 0.0, 'DStellate':  0.3, 'MNTB': 0.3} # netcon delay (fixed)
        self.AN_lat_Flag = {'stellate': 0,    'bushy': 1,    'test': 0,    'DStellate':  0, 'MNTB': 1} # 1 to allow latency shift during train
        self.AN_lat_A0 =   {'stellate': 0.14, 'bushy': 0.140, 'test': 0.25, 'DStellate':  0.15, 'MNTB': 0.14} # magnitude of latency change, if allowed
        self.AN_lat_tau =  {'stellate': 21.0, 'bushy': 21.5, 'test': 21.0, 'DStellate':  21.0, 'MNTB': 21.5} # rate at which latency changes, if allowed
        self.AN_latency =  {'stellate': 0.5,  'bushy': 0.5,  'test': 0.5,  'DStellate':  0.5, 'MNTB': 0.} # starting latency
        self.AN_relstd_Flag = {'stellate':   1, 'bushy': 2,    'test': 0,   'DStellate':  0, 'MNTB': 2} # flag controlling std shift during train (1 to enable)
        self.AN_relstd =      {'stellate': 0.05, 'bushy': 0.05, 'test': 0.05, 'DStellate':  0.05, 'MNTB': 0.05} # baseline release std
        self.AN_relstd_A0 =   {'stellate': 0.0, 'bushy': 0.05, 'test': 0.0, 'DStellate':  0.0, 'MNTB': 0.05} # magnitude of std change during train, if allowed
        self.AN_relstd_tau =  {'stellate': 35.0, 'bushy': 35.0, 'test': 35.0, 'DStellate':  35.0, 'MNTB': 35.0} # time constant of std change during train, if allowed

        # define inhibitory input from D-stellate cell onto a target cell
        self.DS_gI = 50
        self.Inh_Thresh = -30.0  # trigger release at this amplitude
        self.INH_conv =   {'stellate':   8, 'bushy':  10, 'test':   1, 'DStellate':    2} # normally 10
        self.INH_zones =  {'stellate':   5, 'bushy':  10,  'test':   1, 'DStellate':    5} # normally 10
        self.INH_erev =   {'stellate': -70, 'bushy': -70, 'test': -70., 'DStellate': -70}
        self.INH_gvar =   {'stellate': 0.3, 'bushy': 0.3, 'test': 0.3, 'DStellate':  0.3}
        self.INH_delay =  {'stellate': 0.0, 'bushy': 0.0, 'test': 0.0, 'DStellate':  0.0} # netcon delay (fixed)
        self.INH_lat_Flag = {'stellate': 0,    'bushy': 0,    'test': 0,    'DStellate':  0} # 1 to allow latency shift during train
        self.INH_latency = {'stellate': 0.2, 'bushy': 0.2, 'test': 0.2, 'DStellate':  0.2} # latency for stochastic synapses
        self.INH_relstd_Flag = {'stellate':   1, 'bushy': 2,    'test': 0,   'DStellate':  0} # flag controlling std shift during train (1 to enable)
        self.INH_relstd = {'stellate': 0.2, 'bushy': 0.2, 'test': 0.3, 'DStellate':  0.2} # release sigma

        # set up arrays for measurement of voltage, calcium, etc.
        self.nswel = 200
        self.nactual_swel = 0
        #load_file(1, "calyx_GBC.hoc")
        h.load_file(1, self.topofile) # load the morphology file
        h.celldef()
        # create a dictionary to access the calyx parts
        self.CalyxStruct={}
        self.clist = {}
        for name in self.calyxNames: # for each part
            x=eval("h.%s" % (name)) # find the associated variable
            print 'names: ', name, x
            self.CalyxStruct[name] = list(x) # and populate it with the pointers to the parts
        # h.load_file(1, "calyx_morpho.hoc")
        # h.load_file(1, "calyx_shape.hoc")
        self.biophys(pars = self.gPars, createFlag = True)

        # now for the postsynaptic cell - just a basic bushy cell
        (self.TargetCell, [self.initseg, self.axn, self.internnode]) = nrnlibrary.Cells.bushy(debug=False, ttx=False,
                                        message=None,
                                        nach='jsrnaf',
                                        species='mouse', axon=False, dendrite=False,
                                        newModFiles=False, pump=False)

        sPars = nrnlibrary.Synapses.Params(LN_Flag = self.AN_relstd_Flag[self.TargetCellName],
            LN_t0 = 0.,
            LN_A0 = self.AN_relstd_A0[self.TargetCellName],
            LN_tau = self.AN_relstd_tau[self.TargetCellName],
            LN_std = self.AN_relstd[self.TargetCellName],
            Lat_Flag = self.AN_lat_Flag[self.TargetCellName],
            latency = self.AN_latency[self.TargetCellName],
            Lat_t0 = 0.,
            Lat_A0 = self.AN_lat_A0[self.TargetCellName],
            Lat_tau = self.AN_lat_tau[self.TargetCellName],
            delay=self.AN_delay[self.TargetCellName])

        self.NMDARatio = 0.0
        # make the incoming synapse to the target cell
        # to mimic MNTB cells, we source each presynaptic V from the swellings
        # (implicit assumption is that swellings are synapses).
        # nFibers is set to 1
        # nRZones is set to whatever (let's say, 2 per swelling?)

        print "GBC excitation to target"

        self.calyx = []
        self.coh=[]
        self.psd = []
        self.cleft = []
        self.nc2 = []
        for swellno, swell in enumerate(self.CalyxStruct['swelling']):
#        if len(self.CalyxStruct['swelling']) > 1:
            swelling = self.getAxonSec('swelling', swellno)
            print swelling
            (calyx, coh, psd, cleft, nc2, par) = nrnlibrary.Synapses.stochastic_synapses(h, parentSection = self.CalyxStruct['axon'][swelling],
                                                                                         targetcell=self.TargetCell,
                                                                                         cellname = self.TargetCellName,
                    nFibers = self.AN_conv[self.TargetCellName], nRZones = self.AN_zones[self.TargetCellName],
                    message='creating stochastic multisite synapse',
                    thresh = self.Exc_Thresh, psdtype = 'ampa', gmax=self.AN_gEPSC[self.TargetCellName]*1000.,
                    gvar = self.AN_gvar[self.TargetCellName],
                    eRev = self.AN_erev[self.TargetCellName], NMDARatio = self.NMDARatio,
                    debug = False, Identifier = 1, stochasticPars = sPars, calciumPars=None)
            self.calyx.append(calyx)
            self.coh.append(coh)
            self.psd.append(psd)
            self.cleft.append(cleft)
            self.nc2.append(nc2)
        nrnlibrary.Synapses.printParams(sPars)

        self.runModel()

    def restoreDefaultConductances(self):
        """ default conductance values determined by voltage clamp test on model calyx
            canonical settings. This restores those values.
        """
        self.gna_def = 0.45
        self.glvk_def = 0.04
        self.ghvkax_def = 0.02
        self.ghvk_def = 0.005
        self.gca_def = 0.003
        self.gh_def = 0.000950
        self.gleak_def = 0.000004935
        self.capump_def = 3e-14 # mol/sec
        self.setDefaultConductances()
    
    def setDefaultConductances(self):
        """ For each type of segment, initialize setting by the default conductances """
        for name in self.calyxNames:
            self.gPars[name] = Params(gNabar=self.gna_def, ENa=self.eNa_def,    gKLbar =self.glvk_def, EK = self.eK_def,
                gKHbar = self.ghvkax_def, gCabar = self.gca_def, gHbar = self.gh_def, gCaPump = self.capump_def)
            # exceptions:
            if name in ['axon', 'heminode']:
                self.gPars[name].gCabar = 0.0
            if name in ['axon']:
                self.gPars[name].gNabar = self.gna_def/5.0
            if name in ['neck', 'branch', 'tip', 'swelling']:
                self.gPars[name].gNabar = 0.0
                self.gPars[name].gKLbar = 0.0

# calyx_biophysics.hoc
# The Calyx Biophysics (Channel insertion)
# Every different kind of section has its own conductance levels
# and can have different channels
#
# Paul B. Manis, Ph.D.
# 25 Sept. 2007
# Modified to use gca for HH formulation of calcium current
# 14 Oct 2007
# converted for Python, 17 Oct 2012.
#
# build the biophysics... 

    def biophys(self, pars=None, createFlag = True):
        if createFlag : # first time through insert basics only
            for sec in self.CalyxStruct['axon']:
                sec.insert('leak')
                sec().g_leak = self.newg_leak
                sec.Ra = self.newRa
                sec.cm = self.newCm
                e_leak = -50
                sec.insert('capmp')
                pcabar_cachan = 2.5e-5
                sec().cao = 1.0 #insert capump
                sec.cai = 70e-6
                sec.pump0_capmp = 0
 
        for sec in self.CalyxStruct['axon']: #  always keep these base parameters updated
            sec.g_leak = self.newg_leak
            sec.Ra = self.newRa
            sec.cm = self.newCm
        if pars is None:
            return 
        # for each of the parts of the calyx, insert the appropriate conductances
        # into the model, and set the conductance level and the Nernst potential.
        for name in self.calyxNames:
            gp = self.gPars[name]
            for sec in self.CalyxStruct[name]:
               # print 'for %s  sec = ' % (name), sec
                sec.insert('na')
                if self.pharmManip['TTX']:
                    sec().gnabar_na = 0.0
                else:
                    sec().gnabar_na = gp.gNabar
                sec().ena = gp.ENa
                sec.insert('klt')
                if self.pharmManip['DTX']:
                    sec().gkltbar_klt = 0.0
                else:
                    sec().gkltbar_klt = gp.gKLbar
                sec().ek = gp.EK
                sec.insert('kht')
                if self.pharmManip['TEA']:
                    sec().gkhtbar_kht = 0.0
                else:
                    sec().gkhtbar_kht = gp.gKHbar
                sec.insert('ih')
                if self.pharmManip['ZD']:
                    sec().ghbar_ih = 0.0
                else:
                    sec().ghbar_ih = gp.gHbar
                sec().eh_ih = -43
                if name is 'axon':
                    sec().cm = 0.0001 # mimic myelination
                    sec().pump0_capmp = 0
                else:
                    sec().cm = self.newCm
                if name not in ['axon', 'heminode']:
                    sec.insert('CaPCalyx')
                    if self.pharmManip['Cd']:
                        sec().gcapbar_CaPCalyx = 0.0
                    else:
                        sec().gcapbar_CaPCalyx = gp.gCaPump 
                    sec().eca=43.9


    #*********************CALYX_INIT************************
    # Calyx model initialization procedure:
    # Set RMP to the resting RMP of the model cell.
    # Make sure nseg is large enough in the proximal axon.
    # Initialize the leak, and stabilize the pumps, etc.

    def calyx_init(self):
    # First we set e_leak so that the rmp in each segment is the same
        h.finitialize(self.v_init)
        i = 0
        for sec in self.CalyxStruct['axon']: # set a new value
            #print 'ina: %f  ica: %f  ik: %f  ih: %f g_leak: %f  vinit: %f' % (sec.ina, sec.ica, sec.ik, sec.i_ih, sec.g_leak, h.v_init)
            #print '      ena %f  eca: %f  ek %f  eh %f  eleak %f' % (sec.ena, sec.eca, sec.ek, sec.eh_ih, sec.e_leak)
            if i == 0:
                #print dir(sec().sec)
                i += 1
            e_newleak = sec.v + (sec.ina + sec.ik + +sec.ik_klt + sec.ik_kht + sec.ica + sec.i_ih)/(sec.g_leak*1e3) # hmmm. This is probably NOT right.
            #print 'leak: ', e_newleak
            sec().e_leak = e_newleak
        
    # based on those conditions, make sure spatial grid is fine enough for our needs
        #h.Mesh_th = self.Mesh_th
        #h('mesh_init')
        #self.mesh_init()
        h.finitialize(self.v_init)
 #       h('axon[axonnode] ic.loc(iclocation)')
    
        for name in self.calyxNames:
            if name in ['axon', 'heminode']:
                for sec in h.parentaxon:
                    sec.nseg = 11 

    # get ready to run the system to a stable point for the calcium pump
        j = 0
    # save the calcium pump information
        for sec in self.CalyxStruct['axon']:
            savcore = sec().cabulk_capmp
            sec().cabulk_capmp = self.ca_init
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
        for sec in self.CalyxStruct['axon']:
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
        """
        src = re.compile('axon\[(\d*)\]')
        for s in self.calyxNames:
            if s in ['axon', 'synapse']:
                continue # skip the catch-all categories
            for element in self.CalyxStruct[s]:
                name = element.name()
#                print 's: %s  name: %s' % (s, name)
                g = src.match(name)
                axno = g.groups()[0]
                if (int(axno)) == monsec:
                    self.clist[monsec] = self.calyxColors[s]
                    print 'match found for %d in %s' % (monsec, s)

    def getAxonSec(self, structure, number):
        src = re.compile('axon\[(\d*)\]')
        assert structure in self.CalyxStruct.keys()
        try:
            element  = self.CalyxStruct[structure][number]
        except:
            raise Exception(('Failed to find %d in structure %s' % (number, structure())))

        name = element.name()
#                print 's: %s  name: %s' % (s, name)
        g = src.match(name)
        axno = g.groups()[0]
        return int(axno)

    # ***************CALYXRUN*******************
    # Control a single run the calyx model with updated display
    #  of the voltages, etc. 
    #

    def calyxrun(self):
        # { local i, j, k, b
        print 'calyxrun'
        rundone = 0
        self.vec={}
        for var in ['axon', 'stalk', 'branch', 'neck', 'tip', 'swelling', 'inj', 'time', 'postsynaptic', 'istim']:
            self.vec[var] = h.Vector()

        h.tstop = self.tstop
        h.dt = 0.005 # force small time step. cvode is probably off.
        h.celsius = self.celsius
        npts = 1 + int(h.tstop/h.dt) # number of points in a run
        self.vec2 = {} # voltage in each monitored segment
        self.ica = {} # calcium current in each monitored segment
        self.cai = {} # intracellular calcium (bulk) in each monitored segment
        self.clist = {} # color list
        segarea = {} # np.zeros(len(self.monsecs))
        segvol = {} # np.zeros(len(self.monsecs))
#        pprint.pprint(self.CalyxStruct)
        for m, monsec in enumerate(self.monsecs):
            esection=self.CalyxStruct['axon'][monsec]
            self.vec2[monsec] = h.Vector() # set up recording vectors at each site
            self.ica[monsec] = h.Vector()
            self.cai[monsec] = h.Vector()
            self.vec2[monsec].record(esection(self.monpos[m])._ref_v, sec=esection)
            self.ica[monsec].record(esection(self.monpos[m])._ref_ica, sec=esection)
            self.cai[monsec].record(esection(self.monpos[m])._ref_cai, sec=esection)
            segarea[monsec] = esection.L*np.pi*esection.diam * 1e-8 # L and diam in microns, area in cm^2
            segvol[monsec] = esection.L*np.pi*(esection.diam/2.)**2 # leave in microns

            # set color list for plots - each vector is colored by the type of segment it points to
            self.getMonSec(monsec)

        if self.postMode in ['vc', 'vclamp']:
            clampV=-60.0
            vcPost = h.SEClamp(0.5, sec=self.TargetCell)
            vcPost.dur1 = self.tstop-3.0
            vcPost.amp1 = clampV
            vcPostl.dur2 = 2.0
            vcPost.amp2 = clampV-0.5 # just a tiny step to keep the system honest
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

        electrodesite = self.CalyxStruct['axon'][self.axonnode]
        istim = h.iStim(0.5, sec=electrodesite)
        stim={}
        stim['NP'] = self.nStim
        stim['Sfreq'] = self.stimFreq # stimulus frequency
        stim['delay'] = self.stimDelay
        stim['dur'] = self.stimDur
        stim['amp'] = self.stimInj
        stim['PT'] = 0.0
        (secmd, maxt, tstims) = self.make_pulse(stim, pulsetype='square')
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        self.vec['i_stim'] = h.Vector(secmd)

        self.vec['i_stim'].play(istim._ref_i, h.dt , 0, sec=electrodesite)
        self.vec['axon'].record(electrodesite()._ref_v, sec=electrodesite) # record axon voltage
        self.vec['inj'].record(istim._ref_i,  sec=electrodesite)
        self.vec['time'].record(h._ref_t)
        h.run()
        fig = MP.figure(100)
        p1 = fig.add_subplot(4,1,1)
        p2 = fig.add_subplot(4,1,2)
        p3 = fig.add_subplot(4,1,3)
        p4 = fig.add_subplot(4,1,4)
        p1.plot(self.vec['time'], self.vec['axon'], color=self.clist[0])
       # p2.plot(self.vec['time'], self.ica['axon'])
        for v in self.vec2:
            #print v
            p1.plot(self.vec['time'], self.vec2[v], color=self.clist[v])
        p1.set_ylabel('V')
        for c in self.ica:
            p2.plot(self.vec['time'], self.ica[c], color=self.clist[c])
            p3.plot(self.vec['time'], self.cai[c], color=self.clist[c])
        p2.set_ylabel('I_{Ca}')
        p3.set_ylabel('[Ca]_i')
        p4.plot(self.vec['time'], self.vec['postsynaptic'], color = 'k')
        #p2.set_ylim(-5e-12, 1e-12)


    #     vax[1].record(&axon[axonnode].v(0)) # record axon voltage too
    #     vax[2].record(&axon[axonnode].v(1)) # record axon voltage too
    #     vclampi.record(&vc.i) # record the voltage clamp current
    #     forsec swelling { # save voltage, calcium current and concentration
    #         thissec = secname()
    #         th = j + 1 ncolor = color+j
    #         if(j == mark_sw) { th = 1 ncolor = 2 }
    #         sscanf(thissec, "axon[%d]", &b) # get the axon identification... 
    # #       printf("j = %d,, thissec: %s   b = %d\n", j, thissec, b)
    #         if(monsecs.contains(b)) {
    #             k = monsecs.indwhere("==", b) # is this a selected sectin?
    #             printf("k = %d\n", k)
    #             sprint(scmd, "%s.ica(%.1f)", thissec, monpos.x[k])
    #             sprint(lcmd, "ax %d", monsecs.x[k])
    #             gc.addvar(lcmd, scmd, k+1, th)
    #             sprint(scmd, "%s.cai(%.1f)", thissec, monpos.x[k])
    #             gi.addvar(lcmd, scmd, k+1, th)
    #             k = k + 1 # k counts the monitored sections
    #         }
    #         vswel[j] = new Vector(npts, 0)
    #         caswel[j] = new Vector(npts, 0)
    #         icaswel[j] = new Vector(npts, 0)
    #         sprint(svsw, "vswel[%d].record(&%s.v(0.5))", j, thissec)
    #         sprint(slca, "caswel[%d].record(&%s.cai(0.5))", j, secname())
    #         sprint(sica, "icaswel[%d].record(&%s.ica(0.5))", j, secname())
    #         execute(svsw) # one way to allow us to do this on the fly
    #         execute(slca)
    #         execute(sica)
    #         j = j + 1
    #     }
    # 
    # # note: we do our own run here with fadvance, rather than run()
    #     tdat = new Vector(npts) # tdat keeps the sample times
    #     t=0
    #     calyx_init() # initialize and stabilize up the model
    #     jt = 0
    #     while (t<tstop) { # for the selected time window
    #         fadvance() # step
    #         tdat.x[jt] = t
    #         jt  = jt + 1
    #         gp.plot(t) # plot 
    #         gvc.plot(t)
    #         gc.plot(t)
    #         gi.plot(t)
    #         shv.flush()
    #         shcai.flush()
    #     }
    #     gp.flush()
    #     gvc.flush()
    #     gc.flush()
    #     gi.flush()
    #     shv.flush()
    #     shcai.flush()
    #   doNotify()
        rundone = 1
    #    thisrep = thisrep + 1
    #   vref = axon[0].v(0.5)
    #   for j = 0, nactual_swel-1 { # calculate difference here - 
    #       vmin.x[j] = vswel[j].min - vref
    #   }
        print("CalyxRun done\n")
        MP.show()

    def make_pulse(self, stim, pulsetype = 'square'):
        """Create the stimulus pulse waveform.

            stim is a dictionary with a delay, duration, SFreq,
            and post duration, and number of pulses
            Pulses can be 'square' or 'exponential' in shape.

        """
        delay = int(np.floor(stim['delay']/h.dt))
        ipi = int(np.floor((1000.0/stim['Sfreq'])/h.dt))
        pdur = int(np.floor(stim['dur']/h.dt))
        posttest = int(np.floor(stim['PT']/h.dt))
        NP = int(stim['NP'])
        maxt = h.dt*(stim['delay'] + (ipi*(NP+2)) + posttest + pdur*2)
        w = np.zeros(np.floor(maxt/h.dt))

        #   make pulse
        tstims = [0]*NP
        if pulsetype == 'square':
            for j in range(0, NP):
                t = (delay + j *ipi)*h.dt
                w[delay+ipi*j:delay+(ipi*j)+pdur] = stim['amp']
                tstims[j] = delay+ipi*j
            if stim['PT'] > 0.0:
                send = delay+ipi*j
                for i in range(send+posttest, send+posttest+pdur):
                    w[i] = stim['amp']

        if pulsetype == 'exp':
            for j in range(0, NP):
                for i in range(0, len(w)):
                    if delay+ipi*j+i < len(w):
                        w[delay+ipi*j+i] += stim['amp']*(1.0-np.exp(-i/(pdur/3.0)))*np.exp(-(i-(pdur/3.0))/pdur)
                tstims[j] = delay+ipi*j
            if stim['PT'] > 0.0:
                send = delay+ipi*j
                for i in range(send+posttest, len(w)):
                    w[i] += stim['amp']*(1.0-np.exp(-i/(pdur/3.0)))*np.exp(-(i-(pdur/3.0))/pdur)

        return(w, maxt, tstims)
    # 
    # 
    # 
    # # *************************************************************
    # # Procedure to write data files to disk.
    # # Note:
    # # The "data file" really consists of 3 files. 
    # # 1. A "header" file, in text format, giving information about the
    # # most recent runn
    # # 2. A "data" file, in text format, but just as a matrix of data points
    # # with the first column being time. Only the voltage at the swellings
    # # is in this file at the moment.
    # # 3. A second data file, with the voltage in the axon itself at 3 points
    # # (far end, middle, and junction with the calyx). 
    # # One consideration is whether it might make more sense to write the data
    # # file as one binary file containing the voltages at all the axonal segements.
    # # 
    # # *************************************************************
    # 
    # objref outfile
    # strdef filename, basefilename, expt, today, datafilename, axonfilename
    # 
    # proc write_data() { local i, j
    #     if(rundone == 0) {
    #         return
    #     }
    #     expt = "Calyx5.hoc: voltage, calcium and Ica at swellings"
    #     system("date", today)
    #     basefilename = "C"
    #     maxout = tdat.size
    #     sprint(filename, "%s/%s-%s.txt", folder, basefilename, manipulation)
    #     printf("\nRaw Voltage Data goes into file: %s\n", filename)
    #     sprint(datafilename, "%s/%s-%s.dat", folder, basefilename, manipulation)
    #     sprint(axonfilename, "%s/%s-%s-ax.dat", folder, basefilename, manipulation)
    # # 
    # # the first file is the "header" file 
    # #
    #     outfile = new File()
    #     u = outfile.wopen(filename)
    #     if(u == 0) {
    #         printf("\n UNABLE TO OPEN OUTPUT FILE: %s\n\n", filename)
    #         return # that is an error - all runs are stored to data files
    #     }       
    #     outfile.printf("%s %d %d\n", datafilename, nactual_swel, maxout)
    #     outfile.printf(" Calyx5.hoc (11/2/2007)  Experiment: %s\n", expt)
    #     outfile.printf("Topology File: %s\n", topofile) # indicate topology source
    #     outfile.printf("Data Run: %d Points: %d  Run Executed on: %s\n", thisrep, maxout, today)
    # # write parameters
    #     outfile.printf("Axon:     gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_ax, ghvk_ax, glvk_ax, gca_ax)
    #     outfile.printf("Stalk:    gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_st, ghvk_st, glvk_st, gca_st)
    #     outfile.printf("Swelling: gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_sw, ghvk_sw, glvk_sw, gca_sw)
    #     outfile.printf("Branch:   gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_br, ghvk_br, glvk_br, gca_br)
    #     outfile.printf("Neck:     gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_nk, ghvk_nk, glvk_nk, gca_nk)
    #     outfile.printf("Tip:      gNa  %8.3f   gHVK %8.3f  gLVK %8.3f  gCa %8.3f\n", gna_tp, ghvk_tp, glvk_tp, gca_tp)
    #     outfile.printf("Calcium: ca_init: %f  k1_capmp: %f  pump0: %f\n", ca_init, k1_capmp, pump0_capmp)
    #     outfile.printf("Passive: Ra: %f  g_leak: %f\n", newRa, newg_leak)
    # 
    #     # header line with identification of columns for immediate IGOR import
    # 
    #     outfile.printf("\"tdat\"  ")
    #     for j = 0,(nactual_swel-1){
    #         outfile.printf("\"Vsw%d\" " ,j)
    #     }
    #     for j = 0,(nactual_swel-1){
    #         outfile.printf("\"ICa_sw%d\" " ,j)
    #     }   
    #     for j = 0,(nactual_swel-1){
    #         outfile.printf("\"[Ca]i_sw%d\" " ,j)
    #     }   
    #     outfile.printf("\n")
    #     outfile.printf("Nswel = %d\n", nactual_swel)
    #     measure() # get the swelling map.
    #     for j = 0, (nactual_swel-1) {
    #         outfile.printf("%d %9.3f\n", j, swelldist.x[j]) 
    #     }
    #     outfile.printf("\n")
    #     outfile.close
    #     
    #     # the actual data file (currently voltage at swellings, along with
    #     # calcium channel current and intracellular calcium
    #     #
    #     u = outfile.wopen(datafilename)
    #     if(u == 0) {
    #         printf("\n UNABLE TO OPEN OUTPUT FILE: %s\n\n", filename)
    #         return # that is an error - all runs are stored to data files
    #     }   
    # 
    #     for i = 0, maxout-1 {
    #         ii = i
    #         outfile.printf("%8.3f ", tdat.x[ii])
    #         for j = 0, (nactual_swel-1) {
    #             outfile.printf("%8.6e ", vswel[j].x[ii])
    #         }
    # 
    #         for j = 0, (nactual_swel-1) {
    #             outfile.printf("%8.6e ", icaswel[j].x[ii])
    #         }
    #     
    #         for j = 0, (nactual_swel-1) {
    #             outfile.printf("%8.6e ", caswel[j].x[ii])
    #         }   
    #         outfile.printf("\n")
    #     }
    #     outfile.printf("\n")
    #     outfile.close
    # 
    #     # A second data file with the axon voltage
    #     u = outfile.wopen(axonfilename)
    #     if(u == 0) {
    #         printf("\n UNABLE TO OPEN OUTPUT FILE: %s\n\n", filename)
    #         return # that is an error - all runs are stored to data files
    #     }   
    # 
    #     for i = 0, maxout-1 {
    #         ii = i
    #         outfile.printf("%8.3f ", tdat.x[ii])
    #         for j = 0, 2 {
    #             outfile.printf("%8.6e ", vax[j].x[ii])
    #         }
    #         outfile.printf("\n")
    #     }
    #     outfile.printf("\n")
    #     outfile.close
    # 
    #     printf("Output file write complete\n")
    # }
    # 
    # 
    
    
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
    #*************************************************************************
    # Main code to be executed when calyxN.hoc is called
    # 
    #*************************************************************************
    def runModel(self):
#        h('celldef()') # includes biophysics
    
        # if(GBCFLAG == 1) { defineGBC() } # include the GBC?
#        h('biophys(1)') # make sure we are all set.
#        h('measure()') # get distances
        self.restoreDefaultConductances() # canonical conductances... 
        self.setDefaultConductances() # initialize the conductances
        print("Iclamp in Cut Calyx Axon[1]\n")
        #print self.CalyxStruct.keys()
        #print self.axonnode
        # self.ic = h.IClamp(0.5, sec=self.CalyxStruct['axon'][self.axonnode])
        # self.ic.delay = self.icDelayDef
        # self.ic.dur = self.icDurDef
        # self.ic.amp = self.icAmpDef


        # load anciallary functions that may need all the definitions above
#        h.load_file(1, "calyx_tune.hoc") # requires clamps to be inserted..
        self.calyx_init() # this also sets nseg in the axon - so do it before setting up shapes
#        h.shape_graphs() # must happen AFTER calyx_init
        self.calyxrun()

if __name__ == "__main__":

    print sys.argv[1:]
    calyx7(sys.argv[1:])