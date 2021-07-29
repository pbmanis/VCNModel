#!/usr/bin/python
from __future__ import print_function
import numpy as np
import cnmodel.util as CU
from pathlib import Path
from pylibrary.tools.cprint import cprint

"""
cellInitialization provides routines for initializing the membrane potential
of extended cells.
Once initialized the cells should be run without any further calls to
finitialize, and should not be run with h.run().
Use fadvance() instead, or
even better, use the "new" h.batch_run() function in hoc, which allows
you to use saved states to continue previous runs.

"""


def init_model(
    cell: object,
    mode: str = "CC",
    vinit: float = -65.0,
    restore_from_file: bool = False,
    filename: Path = None,
    electrode_site: object = None,
    reinit: bool = False,
) -> bool:
    """
    Model initialization procedure to set RMP to the resting RMP of the model cell.
    Does not instantiate recording or stimulating.

    Params
    ------
    cell : cell, including hoc_reader file object (as .hr)
    mode : str (default: 'CC')
        mode string ('CC', 'VC', 'gifnoise').

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

    if mode not in ['CC', 'VC']:
        raise ValueError('Mode must be "VC" or "CC". Got %s' % mode)
    if mode in ["VC"]:
        cell.hr.h.finitialize(
            vinit
        )  # this is sufficient for initialization in voltage clamp
        cprint("c", "    Initializing for Vclamp")
        return True

    # otherwise we are in current clamp
    # get state if one is specified
    if restore_from_file:
        cprint("c", f"    Restoring initialization from: {str(filename):s}")
        restore_initial_conditions_state(
            cell, electrode_site=electrode_site, filename=filename, reinit=reinit
        )
        try:
            cell.hr.h.frecord_init()  # try an intialization
        except RuntimeError:
            cprint("r", "\nUnable to restore initial state")
            raise RuntimeError("\nUnable to restore initial state")
        cprint("g", "        Initialization was restored OK")
        return True  # much easier here...

    # cprint('g', "    Performing custom initialization")
    CU.custom_init(v_init=vinit)   # this actually will not be run (uses init from below)

    if electrode_site is not None:
        vm = electrode_site.v
    else:
        vm = 0.0
    cprint(
        "g",
        f"    Initialized with finitialize, starting at {vinit:8.2f}, end voltage: {vm:8.2f} mV",
    )
    return True


def printCellInfo(cell: object) -> None:
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
                    print("nseg: ", sg.nseg)
                print(
                    cell.all_sections[st][0].allseg()
                )  # iteration object over all segments.
            if i < 20:
                print("%s %d: %d" % (st, i, s.nseg))


def get_initial_condition_state(
    cell: object,
    filename: Path,
    mode: str = "CC",
    tdur: float = 2000.0,
    electrode_site: object = None,
    reinit: bool = False,
    freq: float = 2000.0
) -> None:
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
    if filename is None:
        raise ValueError("cellInitialization:get_initial_condition_state: Filename must be specified")

    if mode == "VC":
        cell.vm0 = -65.0  #

        cell.i_currents(cell.vm0)
    else:
        cell.cell_initialize()

    # first to an initialization to get close
    cprint("c", f"\nGetting initial_condition_state: file={str(filename):s}")
    print("        starting t = %8.2f" % cell.hr.h.t)
    init_model(
        cell, mode=mode, restore_from_file=False, electrode_site=electrode_site, reinit=reinit
    )
    cell.hr.h.tstop = tdur
    print("        running for %8.2f ms at %4.1f" % (tdur, cell.hr.h.celsius))
    cell.hr.h.run()
    cprint("g", f"  run completed, t = {cell.hr.h.t:8.2f}")
    if electrode_site is not None:
        vfinal = electrode_site.v
    else:
        vfinal = 0.0
        cprint("g", f"Initialized to to voltage: {vfinal:8.2f} mV")
    state = cell.hr.h.SaveState()
    stateFile = cell.hr.h.File()
    state.save()
    cprint("c", f"        writing state to : {str(filename):s}")
    stateFile.wopen(str(filename))
    state.fwrite(stateFile)
    stateFile.close()


def restore_initial_conditions_state(
    cell: object,
    filename: Path,
    electrode_site: object = None,
    reinit: bool = False,
    autoinit: bool = False
) -> None:
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
        Flag to pass on to init_model to force reinitialization

    Returns
    -------
        Nothing
    """
    #    print('restoring from file: {:s}'.format(filename))
    #    cell.set_d_lambda(freq=100, d_lambda=0.1)
    cell.hr.h.finitialize()
    stateFile = cell.hr.h.File()  # restore state AFTER finitialize
    state = cell.hr.h.SaveState()

    stateFile.ropen(str(filename))
    cprint("c", f"Restoring initial conditions from {str(filename):s}")
    try:
        state.fread(stateFile)
    except RuntimeError:
        cprint("r", "stateFile read failed - states do not match with current model")
        raise RuntimeError("stateFile read failed - states do not match")
    stateFile.close()
    cprint("c", "    ... restoring state")
    state.restore(1)
    cprint("g", "    Successfully restored initial conditions.")

    if electrode_site is not None:
        vm = electrode_site.v
    else:
        vm = cell.hr.h("v")
    #        print 'restored soma v: %8.2f' % vm
    #        print 'v_init after restore: %8.2f' % cell.hr.hr.h.v_init
    cell.hr.h.v_init = vm  # note: this leaves a very slight offset...

#        for group in cell.hr.sec_groups.keys():
#            for sec in cell.hr.sec_groups[group]:
#                section = cell.hr.get_section(sec)
#                print 'section: %s  vm=%8.3f' % (section.name(), section(0.5).v)


def test_initial_conditions(
        cell: object,
        electrode_site: object = None,
        filename: Path = None) -> None:
    """
    Test routine to verify that the initial conditions actually are consistent
        with the cell.

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
    monitor["time"] = cell.hr.h.Vector()
    monitor["time"].record(cell.hr.h._ref_t)
    monitor["Velectrode"] = cell.hr.h.Vector()
    cprint(
        "c",
        f"Test Initial Conditions\n   at site: {str(electrode_site):s}\n   with file: {str(filename.name):s}",
    )
    monitor["Velectrode"].record(electrode_site(0.5)._ref_v, sec=electrode_site)

    restore_initial_conditions_state(
        cell, filename=filename, electrode_site=electrode_site
    )
    cell.hr.h.t = 0.0
    cell.hr.h.tstop = 50.0
    #    hf.hr.h.run()
    # while hf.hr.h.t < hf.hr.h.tstop:
    #     hf.hr.h.fadvance()
    cell.hr.h.batch_save()  # save nothing
    cell.hr.h.batch_run(cell.hr.h.tstop, cell.hr.h.dt, "an.dat")
    print(f"    Filename: {str(filename.name):s}")
    print(f"        time: {str(np.array(monitor['time'])):s}")
    print(f"        Velectrode:  {str(np.array(monitor['Velectrode'])):s}")
    # pg.mkQApp()
    # pl = pg.plot(np.array(monitor['time']), np.array(monitor['Velectrode']))
    # pl.setTitle(filename)
    # QtGui.QApplication.instance().exec_()
