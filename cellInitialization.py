#!/usr/bin/python
from __future__ import print_function
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import numpy as np
import cnmodel.util as CU

"""
cellInitialization provides routines for initializing the membrane potential
of extended cells.
Once initialized the cells should be run without any further calls to
finitialize, and should not be run with h.run().
Use fadvance() instead, or
even better, use the "new" h.batch_run() function in hoc, which allow
you to use saved states to continue previous runs.

"""


def init_model(cell, mode='iclamp', vinit=-65., restore_from_file=False, filename=None, electrode_site=None, reinit=False):
    """
    Model initialization procedure to set RMP to the resting RMP of the model cell.
    Does not instantiate recording or stimulating.
    
    Params
    ------
    cell : cell, including hoc_reader file object (as .hr)
    mode : str (default: iclamp)
        mode string ('iclamp', 'vc', 'vclamp'). 
    vinit : float
        initial voltage to start initialization, in mV. Default -65 mV
    restore_from_file : boolean
        flag to cause model to be restored from previously saved state file.
        The state file must be from the same model construction that exists at the
        time of the call to init_model.
    filename : str (default: None)
        Name of a state file to restore the initial conditions from
    electrode_site : NEURON section object (default: None)
        Site for the recording electrode to be located
    reinit : boolean (default: False)
        If restoring from file, setting reinit to True will also force a reinitialization.
    
    Returns
    -------
    boolean : Success of initialization. Always True, just to indicate we were called.
        Errors result in exceptions.
    
    """
    hf = cell.hr
    if mode in ['vc', 'vclamp']:
        hf.h.finitialize(vinit)  # this is sufficient for initialization in voltage clamp
        return True
    if mode not in ['iclamp']:
        raise ValueError('Mode must be "vc", "vclamp" or "iclamp"; got %s' % mode)

    # otherwise we are in current clamp
    # get state if one is specified
    if restore_from_file:
        restore_initial_conditions_state(cell, electrode_site=electrode_site, filename=filename, reinit=reinit)
        try:
            hf.h.frecord_init()  # try an intialization
        except:
            raise ValueError('Unable to restore initial state')
        return True  # much easier here...
    
    # if hf.h.CVode().active():
    #     #hf.cvode.re_init()
    #     hf.h.CVode().active(0)  # turn cvode off (note, in this model it will be off because one of the mechanisms is not compatible with cvode at this time)
    
    CU.custom_init(v_init=vinit)
    # perform initialization by going back in time and adjusting 
    # hf.h.finitialize(vinit)
    # hf.h.t = -1e6
    # dtsav = hf.h.dt
    # hf.h.dt = 1e3  # big time steps for slow process
    # while hf.h.t < 0:
    #     hf.h.fadvance()
    # hf.h.dt = dtsav
    # hf.h.t = 0
    # if hf.h.CVode().active():
    #     hf.h.CVode().re_init()
    # hf.h.fcurrent()
    # hf.h.frecord_init()
    
    if electrode_site is not None:
        vm = electrode_site.v
    else:
        vm = 0.
    print('Initialized with finitialize, starting at %8.2f, ending %8.2f ' % (vinit, vm))
    return True


def printCellInfo(cell):
    """
    Text output of information about cell - just information
    
    Parameters
    ----------
    cell : CNModel cell instance (required)
    
    Returns
    -------
    Nothing
    """
    for st in cell.all_sections.keys():
        for i, s in enumerate(cell.all_sections[st]):
            if i == 0:
                for sg in cell.all_sections[st]:
                    print('nseg: ', sg.nseg)
                print (cell.all_sections[st][0].allseg() ) # iteration object over all segments.
            if i < 20:
                print('%s %d: %d' % (st, i, s.nseg))


def get_initial_condition_state(cell, tdur=2000., filename=None, electrode_site=None, reinit=False, freq=2000.):
    """
    Run model for a time, and then save the state
    
    Parameters
    ----------
    cell : CNModel cell instance (required)
    tdur : float, ms (default: 2000)
        duration to run to obtain an initial, hopefully stable, state
    filename : str (default: None)
        name of the file to save the state into. This file can be
        retrieved and used to restore the initial conditions.
    electrode_site : NEURON section (default: None)
        location for the recording electrode.
    reinit : boolean (default: False)
        Flag to apss to init_model to force reinitialization
    freq : float, Hz (default: 2000.)
        Frequency to use to help set d_lambda for the cell for this run.
        If you change d_lambda from the default after initialization, you should
        consider forcing a reinit.
    
    Returns
    -------
        Nothing
    """
    
#    set_d_lambda(cell, freq=freq)
    cell.set_d_lambda()
    cell.cell_initialize()
    hf = cell.hr
    # first to an initialization to get close
    print('get_initial_condition_state\n')
    print('  starting t = %8.2f' % hf.h.t)
    init_model(cell, restore_from_file=False, electrode_site=electrode_site, reinit=reinit)
    hf.h.tstop = tdur
    print('running for %8.2f ms at %4.1f' % (tdur, hf.h.celsius))
    hf.h.run()
    print('  run completed, t = %8.2f' % hf.h.t)
    if electrode_site is not None:
        vfinal = electrode_site.v
    else:
        vfinal = 0.
    print('  V = %8.2f' % vfinal)
    state = hf.h.SaveState()
    stateFile = hf.h.File()
    state.save()
    if filename is None:
        filename = 'neuronstate.dat'
    print('  writing state to : %s' % filename)
    stateFile.wopen(str(filename))
    state.fwrite(stateFile)
    stateFile.close()


def restore_initial_conditions_state(cell, filename, electrode_site=None, reinit=False):
    """
    Restore initial conditions from a file
    
    Parameters
    ----------
    cell : CNModel cell instance (required)
    filename : str (required)
        name of the file to save the state into. This file can be
        retrieved and used to restore the initial conditions.
    electrode_site : NEURON section (default: None)
        location for the recording electrode.
    reinit : boolean (default: False)
        Flag to apss to init_model to force reinitialization
    
    Returns
    -------
        Nothing
    """
    hf = cell.hr
    hf.h.finitialize()
    stateFile = hf.h.File() # restore state AFTER finitialize
    state = hf.h.SaveState()
    
    print('Restored initial conditions from file: %s' % filename)
    stateFile.ropen(filename)
    state.fread(stateFile)
    stateFile.close()
    state.restore(1)
    if electrode_site is not None:
        vm = electrode_site.v
    else:
        vm = hf.h('v')
#        print 'restored soma v: %8.2f' % vm
#        print 'v_init after restore: %8.2f' % hf.hr.h.v_init
    hf.h.v_init = vm  # note: this leaves a very slight offset...
#        for group in hf.sec_groups.keys():
#            for sec in hf.sec_groups[group]:
#                section = hf.get_section(sec)
#                print 'section: %s  vm=%8.3f' % (section.name(), section(0.5).v)


def test_initial_conditions(cell, electrode_site=None, filename=None):
    """
    Test routine to verify that the initial conditions work.
    
    Parameters
    ----------
   cell : CNModel cell instance (required)
        The cell object to test against
    electrode_site : NEURON section (default: None)
        location for the recording electrode.
    filename : str (default: None; required)
        name of the file to save the state into. This file can be
        retrieved and used to restore the initial conditions.
    
    Returns
    -------
        Nothing
    """
    assert electrode_site is not None  # Make sure that the site has been specified
    monitor = {}
    monitor['time'] = cell.hr.h.Vector()
    monitor['time'].record(cell.hr.h._ref_t)
    monitor['Velectrode'] = cell.hr.h.Vector()
    print('Test Initial Conditions\n   at site: ', electrode_site)
    monitor['Velectrode'].record(electrode_site(0.5)._ref_v, sec=electrode_site)
    
    restore_initial_conditions_state(cell, filename=filename, electrode_site=electrode_site)
    cell.hr.h.t = 0.
    cell.hr.h.tstop = 50.
#    hf.hr.h.run()
   # while hf.hr.h.t < hf.hr.h.tstop:
   #     hf.hr.h.fadvance()
    cell.hr.h.batch_save() # save nothing
    cell.hr.h.batch_run(cell.hr.h.tstop, cell.hr.h.dt, "an.dat")
    print('\ntime: ', np.array(monitor['time']))
    print('\nVelectrode: ', np.array(monitor['Velectrode']))
    # pg.mkQApp()
    # pl = pg.plot(np.array(monitor['time']), np.array(monitor['Velectrode']))
    # pl.setTitle(filename)
    # QtGui.QApplication.instance().exec_()


# def set_d_lambda(cell, freq=100, d_lambda=0.1):
#     """ Sets nseg in each section to an odd value
#         so that its segments are no longer than
#         d_lambda x the AC length constant
#         at frequency freq in that section.
#
#         Be sure to specify your own Ra and cm _before_ calling geom_nseg()
#
#         To understand why this works,
#         and the advantages of using an odd value for nseg,
#         see  Hines, M.L. and Carnevale, N.T.
#         NEURON: a tool for neuroscientists.
#         The Neuroscientist 7:123-135, 2001.
#
#         This routine is a python adaptation of the hoc file.
#         Some of the original code is included and commented
#
#     Parameters
#     ----------
#     cell : CNModel cell instance (required)
#     freq : float, Hz (default: 100.)
#         frequency to use to set d_lambda.
#     d_lambda : float (default: 0.1)
#         initial value for d_lambda
#
#     Returns
#     -------
#         Nothing
#
#     Note that the routine dynamically modifies nseg for every section inspected.
#
#     """
#
# # the defaults are reasonable values for most models
# # freq = 100      # Hz, frequency at which AC length constant will be computed
# # d_lambda = 0.1
#     hf = cell.hr
# #    forall { nseg = int((L/(d_lambda()*lambda_f(freq))+0.9)/2)*2 + 1  }
#     for st in cell.all_sections.keys():
#         for i, section in enumerate(cell.all_sections[st]):
#             nseg  = int((section.L/(d_lambda*lambda_f(hf, freq, section))+0.9)/2)*2 + 1
#             if nseg < 3:
#                 nseg = 3 # ensure at least 3 segments per section...
#             section.nseg = nseg
# #            print '%s nseg=%d' % (section.name(), section.nseg)
#
# def lambda_f(hf, freq, section):
# #     { local i, x1, x2, d1, d2, lam
#     hf.h('access %s' % section.name())
# #    print 'n3d: ', hf.h.n3d()
#     if hf.h.n3d() < 2:
#             return 1e5*np.sqrt(section.diam/(4.0*np.pi*freq*section.Ra*section.cm))
#
#     # above was too inaccurate with large variation in 3d diameter
#     # so now we use all 3-d points to get a better approximate lambda
#     x1 = hf.h.arc3d(0)
#     d1 = hf.h.diam3d(0)
# #    print 'x1, d1: ', x1, d1
#     lam = 0.
#     for i in range(int(hf.h.n3d())-1):
#             x2 = hf.h.arc3d(i)
#             d2 = hf.h.diam3d(i)
# #            print 'x2, d2: ', x2, d2
#             lam = lam + ((x2 - x1)/np.sqrt(d1 + d2))
#             x1 = x2
#             d1 = d2
#     #  length of the section in units of lambda
#     lam = lam * np.sqrt(2.0) * 1e-5*np.sqrt(4.0*np.pi*freq*section.Ra*section.cm)
# #    print lam, section.Ra, section.cm, section.L
#     return section.L/lam





