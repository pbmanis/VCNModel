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
import os, sys
import neuron as h
from neuron import *
import numpy as np
import scipy as sp
import calyx_biophysics


# # GBCFLAG controls whether we use  is used or not:
h('GBCFLAG = 0') # if 0, IC is in cut axon near calyx; otherwise it is in the GBC soma.
# # note:  gbc model not working correctly 11/4/07 P. Manis
# # cell morph fixed but get a memory error AFTER the run

#strdef ocmd, lfcmd, o2cmd, topofile
topofile = "calyx-S53Acvt2.hoc"
#system("pwd", ocmd) # get the working directory - we move around... 
#sscanf(ocmd, "%s\n", o2cmd) # strip trailing \mn
#sprint(lfcmd, "%s/%s\n", o2cmd, topofile)

# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class calyx6():
    def __init__(self, argsin = None):
        print 'calyx6: init'
        #h.nrn_load_dll(os.getcwd()+'/mechanisms/i386/.libs/libnrnmech.so')
        h.load_file("stdrun.hoc")
        h.load_file(os.getcwd()+"/custom_init.hoc") # replace init with one that gets closer to steady state

        self.clamp_flag = 0 # 0 = current clamp, 1 = voltage clamp
        self.clampstr  = "CC"

        self.rundone = 0
        self.thisrep = 1
        self.manipulation = "Canonical"
        self.folder = "Canonical"

        h('newCm = 1 ') #// uf/cm2 cap.
        h('newRa = 100') # // changed 10/20/2007 to center in range')
        h('newg_leak = 0.000004935') 
        h('axonnode = 149') #'// the heminode section (this may vary with model type)')
        h('celsius = 37') #' // set the temperature.')
        h('ca_init = 70e-6')
        h('v_init = -80')
        h('Ek = -85')
        h('icAmpDef = 5')
        h('icDurDef = 0.25')

        h('Mesh_th = 0.1') # minimum dlambda        
        h('tstop = 10')

        # g_eak set so that input resis (Z_DC) at axon-calyx junction is 
        # 1 Gohm, per Borst and Sakmann, 1978 
        # (Note, did this manually for the passive case, not for the case 
        # with all the conductances in place; default conductance configuration
        # "canonical" then results in an Rin to about 350 Mohm).

    # define which sections (and where in the section) to monitor
    # currents, voltages, and ion concentrations
    # this should be changed on a per-calyx configuration, but right now
    # we nave only one (S53) to work with.
    # sections to monitor. These are selected manually
        self.monsecs = [0, 0, 149, 149, 5, 88, 130, 117]
        self.monpos = [0, 0, 0, 1.0, 0.5, 0.5, 0.5, 0.5]

        # variables to control electrode positions
        h('axonselect = 1')
        h('iclocation = 0.5')
        h('vclocation = 1.0')
    
        # variables for morphology analysis
        # objref shbox, shmorph, slmorph
        #         objref shvbox, shv, slv
        #         objref shcabox, shcai, slcai
        #     
        #         strdef scmd, lcmd, svsw, slca, sica, thissec
        #     
        #         strdef thissecname
        #         objref secid
        #         objref swelldist
        self.dia_thresh = 0.25 # minimum diameter of any part... adjustable in menu
        self.vsaveflag = 0
        self.isaveflag = 0
        self.mark_sw = 1
        self.mark_br = 1
        self.mark_st = 1
        self.mark_tp = 1
        self.mark_nk = 1
    
        # set up arrays for measurement of voltage, calcium, etc.
        h('nswel = 200')
        h('nactual_swel = 0')
        h('objref vswel[nswel], caswel[nswel], icaswel[nswel], vclampi, tdat, vax[3]') 

        # vax[3] 
    
        #load_file(1, "calyx_GBC.hoc")
        ##load_file(1, "calyx-S53Acvt2.hoc") # load the morphology file (old)
        h.load_file(1, "calyx-S53Acvt3.hoc") # load the morphology file
        # important to load these AFTER definitions above.
        # to make the management of the project easier, routines are in separate
        # source files according to general funcitonal groups
        # load these now:
        h.load_file(1, "calyx_biophysics.hoc")
        h.load_file(1, "mesh_init.hoc")
        h.load_file(1, "mesh.hoc")
        h.load_file(1, "calyx_g_control.hoc")
        h.load_file(1, "calyx_experiments.hoc")
        h.load_file(1, "calyx_electrodes.hoc")
        h.load_file(1, "calyx_morpho.hoc")

        h.load_file(1, "calyx_shape.hoc")
        
        self.runModel()

    def reload(self):
      h.load_file(1, "calyx_experiments.hoc")
      h.load_file(1, "calyx_tune.hoc")


# erase the main plots of voltage/current  
# proc ergraphs() {
#   gp.erase_all()
#   gvc.erase_all()
#   gc.erase_all()
#   gi.erase_all()
# }


    #*********************CALYX_INIT************************
    # Calyx model initizlization procedure:
    # Set RMP to the resting RMP of the model cell.
    # Make sure nseg is large enough in the proximal axon.
    # Initialize the leak, and stabilize the pumps, etc.


    def calyx_init(self):
        # local dtsav, temp

    # First we set e_leak so that the rmp in each segment is the same
        h.finitialize(h('v_init'))
        for sec in h.allsec(): # set a new value
            e_newleak = h('v + (ina + ik + ica + i_ih)/g_leak')
    #       printf("%s [%d] -- e_leak: %6.1f  e_newleak: %6.1f, v=%6.2f\n", secname(), j, e_leak, e_newleak, v)
            sec().e_leak = e_newleak
        
    # based on those conditions, make sure spatial grid is fine enough for our needs
        h('mesh_init()')
        h.finitialize(h('v_init')) # recompute
        GBCFLAG = h('GBCFLAG')
        if GBCFLAG == 1:
            h('soma icgbc.loc(0.5)')
        else:
            h('axon[axonnode] ic.loc(iclocation)')
    
        for sec in h.parentaxon:
            sec.nseg = 11 
        for sec in h.heminode:
            sec.nseg = 11
        h('axon[axonnode] vc.loc(vclocation)')

    # get ready to run the system to a stable point for the calcium pump
        j = 0
    # save the calcium pump information
        allsec = h('allsec')
        for sec in allsec:
            savcore = h('cacore_capmp')
            h('cacore_capmp = ca_init')
            savtau = h('tau_capmp')
            h('tau_capmp = 1e-6') # make the pump go really fast

    # starting way back in time
        h.t = -1e10
        dtsav = dt
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
        for sec in h.allsec:
            self.cacore_capmp = savcore
            self.tau_capmp = savtau

        if(cvode.active()):
            h.cvode.re_init()
        else:
            h.fcurrent()
        h.frecord_init()
        h.finitialize(h('v_init'))

    # ***************CALYXRUN*******************
    # Control a single run the calyx model with updated display
    #  of the voltages, etc. 
    #

    def calyxrun(self):
        # { local i, j, k, b
        print 'calyxrun'
        rundone = 0
        h('biophys(0)') # force update of all gmax values 
        h.dt = 0.010 # force small time step. cvode is probably off.
        
    # set up the graphs
        #ergraphs()
        #gp.begin()
        #gvc.begin()
        #gc.begin()
        #gi.begin()
        npts = 1 + int(h.tstop/h.dt) # number of points in a run
        for j in range(self.monsecs.size):
            scmd = "axon[%d].v(%.1f)" % (monsecs.x[j], monpos.x[j])
        #    gp.addexpr(scmd, j+1, 1)
        
    #     if(GBCFLAG == 1) {
    #         gp.addexpr("soma.v(0.5)", 1, 1)
    #     }
    #     gvc.addexpr("vc.i", color, 1) # current (for vclamp)
    # 
    #     nactual_swel = 0 # count the number of swellings
    #     forsec swelling nactual_swel += 1
    #     vmin = new Vector(nactual_swel)
    # 
    #     j = 0 # initialize various counters
    #     k = 0
    #     b = 0
    # 
    #     for i = 0,2 { # save axon voltage near calyx for reference
    #         vax[i] = new Vector(npts, 0)
    #     }
    #     vclampi = new Vector(npts, 0)
    #     vax[0].record(&axon[0].v(0)) # record axon voltage too
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
        thisrep = thisrep + 1
    #   vref = axon[0].v(0.5)
    #   for j = 0, nactual_swel-1 { # calculate difference here - 
    #       vmin.x[j] = vswel[j].min - vref
    #   }
        print("CalyxRun done\n")

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
    # def mesh_init(self):
    #     nra = h.fast_mesh(500, self.Mesh_th)
    #     axon.nseg = 1
    #     print("Number of Compartments: %d\n" % nra)
    # 
    #     # Set up origin in soma
    #     distance(self.axon)
    #     #access axon
    # 
    # 
    # 
    #*************************************************************************
    # Main code to be executed when calyxN.hoc is called
    # 
    #*************************************************************************
    def runModel(self):
        h('celldef()') # includes biophysics
    
        # if(GBCFLAG == 1) { defineGBC() } # include the GBC?
        h('biophys(1)') # make sure we are all set.
        h('measure()') # get distances
        h('RestoreDefaultConductances()') # canonical conductances... 
        h('setdefaultconductances()') # initialize the conductances
        h('insert_clamps()') # electrodes
    
        # load anciallary functions that may need all the definitions above
        h.load_file(1, "calyx_tune.hoc") # requires clamps to be inserted..
        self.calyx_init() # this also sets nseg in the axon - so do it before setting up shapes
#        h.shape_graphs() # must happen AFTER calyx_init
    
        #---------------------------------------------------------
        # Build the main control menu and, separately, the graphs
        #---------------------------------------------------------
        print "build panel\n"
        self.calyxrun()
        # w = new HBox()
        # w.intercept(1)
        # xpanel("Parameters")
        #     
        # #xmenu("Calyx1")
        # #   xbutton("FigureM1", "FigM1()")
        # #   xbutton("Figure2C", "FigM2()")
        # #   xbutton("Figure2D", "FigM3()")
        # #xmenu()
        #     
        # xmenu("Modes")
        # xradiobutton("IClamp", "choose_iclamp()", 1)
        # xradiobutton("Vclamp", "choose_vclamp()", 0)
        # xmenu()
        #     
        # xmenu("Actions", 1)
        # xbutton("Run  ", "calyxrun()")  
        # xbutton("Erase", "ergraphs()")
        # xbutton("Write Data File", "write_data()")
        # xbutton("Mark Swelling", "markshape()")
        # xmenu()
        #     
        # xmenu("Config", 1)
        # xbutton("Canonical", "RestoreDefaultConductances()")
        # xbutton("Passive", "setpassiveconductances()")
        # xbutton("Uniform(from ax)", "setuniformconductances()")
        # xbutton("HighCa", "sethighCaConductances()")
        # xbutton("TipCa", "sethighipca()")
        # xbutton("StalkNa", "setstalkna()")
        # xbutton("TTX", "addttx()")
        # xbutton("DTX", "adddtx()")
        # xbutton("TEA", "addtea()")
        # xbutton("Cd", "addcd()")
        # xbutton("Zd", "addzd()")
        # xbutton("Modify Conductances", "modconductances()")
        # xmenu()
        #     
        # xmenu("Expts", 1)
        # xbutton("Passive", "runpassive()")
        # xbutton("Standard", "runstandard()")
        # xbutton("Uniform", "rununiform()")
        # xbutton("HighCa", "runhighcalcium()")
        # xbutton("HighTipCa", "runhightipca()")
        # xbutton("StalkNa", "runstalkna()")
        # xmenu()
        #     
        # xvalue("IC Loc", "iclocation", 1)
        # xvalue("Iinj (nA)", "ic.amp", 1)
        # xvalue("Idur (ms)", "ic.dur", 1)
        # xvalue("Idel (ms)", "ic.del", 1)
        #     
        # if(GBCFLAG == 1) {
        #     xvalue("Iinjgbc (nA)", "icgbc.amp", 1)
        # xvalue("Idurgbc (ms)", "icgbc.dur", 1)
        # xvalue("Idelgbc (ms)", "icgbc.del", 1)
        # }
        #     
        # xvalue("VC Loc", "vclocation", 1)
        # xvalue("VCHold (mV)", "vc.amp1", 1)
        # xvalue("VCStep (mV)", "vc.amp2", 1)
        # xvalue("VCDur (ms)", "vc.dur2", 1)
        # xvalue("tstop (ms)", "tstop", 1)
        # xvalue("Vrmp (mV)", "v_init", 1)
        #     
        #     
        # xpanel()
        # xpanel("Topology")
        # xvalue("Color","color",1)
        # xvalue("Mark Swelling", "mark_sw", 1, "markswellings()")
        # xvalue("Mark Stalk", "mark_st", 1, "markstalks()")
        # xvalue("Mark Branch", "mark_br", 1, "markbranches()")
        # xvalue("Mark Tip", "mark_tp", 1, "marktips()")
        # xvalue("Mark Neck", "mark_nk", 1, "marknecks()")
        # xvalue("Path from axon", "axonselect", 1, "markpath(0)")
        # xbutton("Set min diam", "markchange()")
        # xvalue("Min. Dia. (um)", "dia_thresh", 1)
        # xbutton("Show Mapping", "showsections()")
        # xbutton("Run  ", "calyxrun()")
        # xbutton("Tune All G and set new", "calyx_tune(0)")
        # xbutton("reload code", "reload()")
        #     
        # xpanel()
        #     
        # w.intercept(0)
        # w.map("MNTB Calyx 1",30, 80, 450, 450)
        #     
        # objref w1, w2
        #     
        # wg = new VBox()
        # wg.intercept(1)
        #     
        # w1 = new HBox()
        # w1.intercept(1)
        # gp = new Graph(0)
        # gp.view(0,-150,170,200,100,0,200,250)
        # gp.size(0, tstop, -100, 20)
        #     
        # gtitle = gp.label(0.2,0.95,"Voltage")
        # gp.exec_menu("Keep Lines")
        #     
        # gvc = new Graph(0)
        # gvc.view(0,-150,170,200,100,0,200,250)
        # gvc.size(0, tstop, -100, 20)
        #     
        # gtitle = gvc.label(0.2,0.95,"Current")
        # gvc.exec_menu("Keep Lines")
        # w1.intercept(0)
        # w1.map()
        #     
        # w2 = new HBox()
        # w2.intercept(1)
        #     
        # gc = new Graph(0)
        # gc.view(0,-150,170,200,100,0,200,250)
        # gc.size(0, tstop, 0.1, -0.5)
        # gctitle = gc.label(0.2, 0.95, "ICa")
        # gc.exec_menu("Keep Lines")
        #     
        # gi = new Graph(0)
        # gi.view(0,-150,170,200,100,0,200,250)
        # gi.size(0, tstop, 0.0, 1.0)
        # gititle = gi.label(0.2, 0.95, "Ca[i]")
        # gi.exec_menu("Keep Lines")
        #     
        # w2.intercept(0)
        # w2.map()
        #     
        # wg.intercept(0)
        # wg.map("MNTB Calyx 1",321, 80, 400,400)
        #     
        #     

if __name__ == "__main__":

    print sys.argv[1:]
    calyx6(sys.argv[1:])