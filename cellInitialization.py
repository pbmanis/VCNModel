"""
cellInitialization provides routines for initializing the membrane potential
of extended cells. 
Once initialized the cells should be run without doing any further
finitialize, and should not be run with h.run(). Use fadvance instead, or 
even better, use the "new" h.batch_run() function in hoc.

"""
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import numpy as np

def initModel(hf, mode='iclamp', vinit=-65., restoreFromFile=False, filename=None, electrodeSite=None):
    """
    Model initialization procedure to set RMP to the resting RMP of the model cell.
    Does not instantiate recording or stimulating.
    
    Params
    ------
    hf : hoc_reader file object
    mode : str
        mode string ('iclamp', 'vc', 'vclamp'). Default: iclamp
    vinit : float
        initial voltage to start initialization, in mV. Default -65 mV
    restoreFromFile : boolean
        flag to cause model to be restored from previously saved state file. 
        The state file must be from the same model construction that exists at the
        time of the call to initModel.
    Returns
    -------
    boolean : Success of initialization. Always True, just to indicate we were called.
    
    """

    if mode in ['vc', 'vclamp']:
        hf.h.finitialize(vinit)
        return True

    # otherwise we are in current clamp
    # Options:
    # 1. adjust e_leak so that the rmp in each segment is the same
    # 2. use ic_constant to inject current in each segment to set rmp
    # 3. allow vm to vary in segments, using existing conductances (may be unstable)
    
    if restoreFromFile:
        restoreInitialConditionsState(hf, electrodeSite=electrodeSite, filename=filename)
        hf.h.frecord_init()
        return True
    
    if hf.h.CVode().active():
        hf.h.CVode().active(0)  # turn cvode off (note, in this model it will be off because one of the mechanisms is not compatible with cvode at this time)

    hf.h.finitialize(vinit)
    hf.h.t = -1e8
    dtsav = hf.h.dt
    hf.h.dt = 1e6  # big time steps for slow process
    n = 0
    while hf.h.t < 0:
        n += 1
        hf.h.fadvance()
    hf.h.dt = dtsav
    hf.h.t = 0
    if hf.h.CVode().active():
        hf.h.CVode().re_init()
    hf.h.fcurrent()
    hf.h.frecord_init()
    
    if electrodeSite is not None:
        vm = electrodeSite.v
    else:
        vm = 0.
    print 'Initialized with finitialize, starting at %8.2f, ending %8.2f ' % (vinit, vm)

    return True

def getInitialConditionsState(hf, tdur=2000., filename=None, electrodeSite=None):
    """
    Run model for a time, and then save the state
    """
    # first to an initialization to get close
    print 'getInitialConditionsState\n'
    print '  starting t = %8.2f' % hf.h.t
    initModel(hf, restoreFromFile=False, electrodeSite=electrodeSite)
    hf.h.tstop = tdur
    print 'running for %8.2f ms' % tdur
    hf.h.run()
    print '  ran until t = %8.2f' % hf.h.t
    if electrodeSite is not None:
        vfinal = electrodeSite.v
    else:
        vfinal = 0.
    print '  V = %8.2f' % vfinal
    state = hf.h.SaveState()
    stateFile = hf.h.File()
    state.save()
    if filename is None:
        filename = 'neuronstate.dat'
    print '  writing state to : %s' % filename
    stateFile.wopen(filename)
    state.fwrite(stateFile)
    stateFile.close()


def restoreInitialConditionsState(hf, electrodeSite=None, filename=None):
    print 'restoring initial conditions'
    hf.h.finitialize()
    stateFile = hf.h.File() # restore state AFTER finitialize
    state = hf.h.SaveState()
    if filename is None:
        filename = 'neuronstate.dat'
    print ' .. restoring file from: %s' % filename
    stateFile.ropen(filename)
    state.fread(stateFile)
    stateFile.close()
    state.restore(1)
    if electrodeSite is not None:
        vm = electrodeSite.v
    else:
        vm = hf.h('v')
#        print 'restored soma v: %8.2f' % vm
#        print 'v_init after restore: %8.2f' % hf.h.v_init
    hf.h.v_init = vm  # note: this leaves a very slight offset...
#        for group in hf.sec_groups.keys():
#            for sec in hf.sec_groups[group]:
#                section = hf.get_section(sec)
#                print 'section: %s  vm=%8.3f' % (section.name(), section(0.5).v)
monitor = {}

def testInitialConditions(hf, electrodeSite=None, filename=None):
    assert electrodeSite is not None
    monitor['time'] = hf.h.Vector()
    monitor['time'].record(hf.h._ref_t)
    monitor['Velectrode'] = hf.h.Vector()
    print 'site: ', electrodeSite
    monitor['Velectrode'].record(electrodeSite(0.5)._ref_v, sec=electrodeSite)

    restoreInitialConditionsState(hf, filename=filename, electrodeSite=electrodeSite)
    hf.h.t = 0.
    hf.h.tstop = 50.
#    hf.h.run()
   # while hf.h.t < hf.h.tstop:
   #     hf.h.fadvance()
    hf.h.batch_save() # save nothing
    hf.h.batch_run(hf.h.tstop, hf.h.dt, "an.dat")
    pg.mkQApp()
    print np.array(monitor['time'])
    print np.array(monitor['Velectrode'])
    pl = pg.plot(np.array(monitor['time']), np.array(monitor['Velectrode']))
    pl.setTitle(filename)
    QtGui.QApplication.instance().exec_()

