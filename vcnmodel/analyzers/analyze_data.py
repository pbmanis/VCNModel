"""
Provides basic spike detection, spike shape analysis, and IV analysis if
appropriate. We use readmodel (ReadModel) to read the pickled data file into
acq4 format, then we can use the ephys analysis tools to analyze the data

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2020 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
from pathlib import Path
from typing import Union
from ephys.ephysanalysis import RmTauAnalysis, SpikeAnalysis
import vcnmodel.util.readmodel as readmodel

SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()


def analyze_data(
    ivdatafile: Union[Path, str],
    filemode: str = "",
    protocol: str = "",
    spike_shape=False,
) -> tuple:
    """
    Perform analysis of spike times and IVs from model data.
    
    Parameters
    ----------
    ivdatafile : path or str path to file
    filemode : str file mode (.v0 or .v1)
    protocol : str
        name of the protocol

    Returns
    -------
    AR : Acq4 Read object (data structure)
    SP : Spike analysis object
    RMA : Analysis summary object
    """
    ReadModel = readmodel.ReadModel()
    ReadModel.read_pfile(ivdatafile, filemode=filemode)

    bridge_offset = 0.0
    threshold = -0.035  # V
    tgap = 0.0  # gap before fittoign taum
    RM.setup(ReadModel.MC, SP, bridge_offset=bridge_offset)
    SP.setup(
        clamps=ReadModel.MC,
        threshold=threshold,  # in units of V
        refractory=0.001,  # In units of sec
        peakwidth=0.001,  # in units of sec
        interpolate=True,
        verify=True,
        mode="peak",
        data_time_units="ms",  # pass units for the DATA
        data_volt_units="V",
    )

    SP.set_detector("Kalluri")  # spike detector: argrelmax, threshold, Kalluri
    SP.analyzeSpikes()
    if spike_shape:
        SP.analyzeSpikeShape()
    SP.analysis_summary["pulseDuration"] = 0.1
    if protocol == "IV":
        print("Read Model MC: ", ReadModel.MC.tstart, ReadModel.MC.tend)
        # SP.fitOne(function="fitOneOriginal")
        RM.analyze(
            rmpregion=[0.0, ReadModel.MC.tstart - 0.001],
            tauregion=[
                ReadModel.MC.tstart,
                ReadModel.MC.tstart
                + (ReadModel.MC.tend - ReadModel.MC.tstart) / 5.0,
            ],
            to_peak=True,
            tgap=tgap,
        )

        # RMA = RM.analysis_summary
    return ReadModel, SP, RM