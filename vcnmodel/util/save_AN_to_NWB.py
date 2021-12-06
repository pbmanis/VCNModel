#!/usr/bin/env python
# coding: utf-8



import matplotlib

matplotlib.use("Qt5Agg")
import datetime as datetime
from pathlib import Path
from typing import Union

import h5py
import numpy as np
import pandas as pd
import pynwb as NWB
from dateutil.tz import tzlocal
from matplotlib import pyplot as mpl

import toml

"""
Save AN data from a model file run to an NWB formatted file
-----------------------------------------------------------
In this script, we read a file generated by model_run with the "--saveall"
flag set so that all section voltages (as taken from the middle of the section) are stored in the file.
The first set of code blocks does this and lets us read and plot the data.

Last update 11/29/2021 pbmanis
Updated for qt5, nwb2.0.0, and matplotlib 3.5
Converted to python 11/30/2021

The top level of the data file contains two specific anciallary data classes: Params and runInfo
------------------------------------------------------------------------------------------------


Plot all locations (for reference) for one stimulus (AN data)
-------------------------------------------------------------
1. Simulation results and input spikes are in the 'Results' list.
2. Each element of the list contains data for a trial.
3. Each trial contains a dict with the following keys:

```
'stimInfo', 'spikeTimes', 'inputSpikeTimes', 'time',
'somaVoltage', 'dendriteVoltage', 'allDendriteVoltages',
'stimWaveform', 'stimTimebase'
```

Depending on the simulation type and run settings, the spikeTimes and allDendriteVoltages entries
may be empty lists.
"""


class SaveNWB(object):
    def __init__(self, basepath:Union[str, Path, None]=None):
        """
        SaveNWB Class
        Provides methods for reading simulations that are output from 
        model_run2.py (vcnmodel), and converting them to NWB format
        Also provides methods for reading the NWB format and plotting
        the data.
        
        Parameters
        ----------
        basepath: Union[str, Path, None] (default: None)
            Sets the base path to the top level of the data set
        
        Returns
        -------
        Nothing
        """
        self.data = None  # storage for file that was read
        if basepath is None:
            self.basepath = Path("/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells")
        else:
            self.basepath = basepath

    def load_simulation(self, fn:Union[str, Path, None]=None):
        """
        Loads a model_run simulation from the specified file/path, relative to the basepath
        Resulting data is in self.data
        
        Parameters
        ----------
        fn : :Union[str, Path, None], default None
        
        Returns
        -------
        Nothing
        Loads from the specified file/path, relative to the basepath
        Resulting data is in self.data
        """

        fnb = Path(self.basepath, fn)
        print("Valid file? ", fnb.is_file())

        self.data = pd.read_pickle(fnb)  # reads the original data file
        print("Keys in data file: ", self.data.keys())
        print(
            "Keys for each trial result; here from the first trial:\n",
            self.data["Results"][0].keys(),
        )

    def get_params_runinfo(self):
        """
        Parse the parameter and run info blocks from the dataset,
        store the values as dictionaries
        
        Parameters
        ----------
        None
        """
        if self.data is None:
            raise ValueError("No data has been read yet")

        params = toml.dumps(self.data["Params"].__dict__)
        # print(params)
        # print("-"*80)
        runinfo = toml.dumps(self.data["runInfo"].__dict__)
        # print(runinfo)
        print(self.data["Params"].cellID)
        # print(d['Params'].hocstring)  # the hoc file for this model run
        # etc.
        print("Parameter names: ", self.data["Params"].__dict__.keys())
        print()
        print("RunInfo names: ", self.data["runInfo"].__dict__.keys())
        self.params = params
        self.runinfo = runinfo

    def getDataArray(self, trial: int = 0):
        """
        Pulls the data for a specified trial into a data array.
        The array is indexed by section number (in spite of the
        name 'allDendriteVoltages', all sections are actually stored).
        
        Parameters
        ----------
        trial: int (default: 0)
            specify with trial from the data set will be loaded
        
        Returns
        -------
        time : np.array
            1d array of time values
        data_array : np.array (2d)
            2d array of data values, sorted by section number
        
        """
        if self.data is None:
            raise ValueError("No data has been read yet")
        dx = self.data["Results"][trial]
        print("trial: ", trial)
        time = dx["time"]
        data_array = None
        secnos = sorted([int(x[9:-1]) for x in list(dx["allDendriteVoltages"].keys())])
        for i, secnum in enumerate(secnos):
            secstr = f"sections[{secnum:d}]"
            y = dx["allDendriteVoltages"][secstr]
            if data_array is None:
                data_array = np.zeros((len(secnos), len(y)))
            data_array[secnum, :] = y  # store data in section order
        return time, data_array

    def plotDendriteVoltages(self, trial: int = 0):
        """
        Simple plot of all the voltages in the sections for a given trial.
        Parameters
        ----------
        trial: int (default: 0)
            specify with trial from the data set will be loaded
        
        Returns
        -------
        nothing
        """
        f, ax = mpl.subplots(1, 1)
        spines_to_remove = ["top", "right"]
        ax = [ax]
        for a in ax:
            for spine in spines_to_remove:
                a.spines[spine].set_visible(False)
        # for trial in d['trials']:
        time, data_array = self.getDataArray(trial=trial)
        for i in range(data_array.shape[0]):
            ax[0].plot(time, data_array[i, :], linewidth=0.25)

    def writeNWB(self, infile, outfile: str = "test", trial: int = 0):
        """
        Save the data in NWB format
        Here we store as much of the ancillary data as we can:
        1. The top level file has the basic information, including a strimg
        version of the params and runinfo dataclasses, concatenated. This preserves key information about 
        the model runs, including the data tables that were used to control the decorations.
        2. AN inputs are stored as a timeseries (these are spike times for each input).
        3. The voltages are stored as "icephys", current clamp traces.
        4. Injected currents are stored as CurrentClampSeries (these may be all 0's if the data is not in IC format)
        
        Assumes (but verifies in getDataArray) that we have already read the data.
        
        Parameters
        ----------
        infile : Path or str
            the input file name
        
        outfile : string
            name of the output file.
        
        trial : int (default = 0)
            The trial to put into the NWB file
        """

        time, data_array = self.getDataArray(trial)
        dx = self.data["Results"][trial]
        """
        Assemble the spike times into a list.
        """
        ANSpikeTimes = []  # store as dict so that input number is an explicit key
        for i in range(len(dx["inputSpikeTimes"])):
            ANSpikeTimes.append(dx["inputSpikeTimes"][i])

        sessionno = 0
        subject = NWB.file.Subject(
            age="0",
            description="vcnmodel",
            genotype="None",
            sex="None",
            species="Computer",
            subject_id="1",
            weight="None",
        )
        self.get_params_runinfo()
        #     print(info)
        nwbfile = NWB.NWBFile(
            "AN Data set",
            str(infile),
            datetime.datetime.now(tzlocal()),
            experimenter="Manis, Paul",
            lab="Manis Lab",
            institution="UNC Chapel Hill",
            experiment_description="Model Output",
            session_id=f"{sessionno:d}",
            notes=self.params + self.runinfo,
            subject=subject,
        )

        device = NWB.device.Device("cnmodel")
        nwbfile.add_device(device)
        # print(an_timestamps)
        ANSpikeOrigin = []
        for i in range(len(ANSpikeTimes)):
            ANSpikeOrigin.append(i * np.ones(len(ANSpikeTimes[i])))
        an_spike_origin = [y for x in ANSpikeOrigin for y in x]
        an_time_stamps = [y for x in ANSpikeTimes for y in x]
        # print(len(an_timestamps))
        # print(an_spike_origin)
        aninputs = NWB.base.TimeSeries(
            name="ANSpikeTimes",
            data=an_spike_origin,
            unit="ms",
            comments="Input AN spike trains",
            description="AN Spike Trains; data is stored in TimeStamps as a dictionary",
            timestamps=an_time_stamps,
        )
        nwbfile.add_stimulus(aninputs)

        elec = NWB.icephys.IntracellularElectrode(
            name="elec0", description="electrodes in middle of section", device=device
        )
        nwbfile.add_ic_electrode(elec)  # not well documented!

        istim = NWB.icephys.CurrentClampStimulusSeries(
            name="Ics",
            data=np.array(dx["stimWaveform"]),
            starting_time=0.0,
            rate=self.data["Params"].__dict__["dtIC"],
            electrode=elec,
            gain=1.0,
            sweep_number=np.uint64(1),
        )
        vdata = NWB.icephys.CurrentClampSeries(
            "CCData",
            data=data_array,
            unit="volts",
            electrode=elec,
            gain=1.0,
            bias_current=0.0,
            bridge_balance=0.0,
            capacitance_compensation=0.0,
            stimulus_description="NA",
            resolution=0.0,
            conversion=1.0,
            timestamps=None,
            starting_time=0.0,
            rate=self.data["Params"].__dict__["dtIC"],
            comments="no comments",
            description="no description",
            control=None,
            control_description=None,
            sweep_number=np.uint64(1),
        )
        nwbfile.add_acquisition(istim)
        nwbfile.add_acquisition(vdata)

        with NWB.NWBHDF5IO(str(Path(outfile)) + ".nwb", "w") as io:
            io.write(nwbfile)
        print(f"Wrote data to NWB file {outfile:s}")

    def readNWB(self, infilename):
        """
        Read and display the data from the specified infilename, and NWB formatted file.
        Provided as an example of how to read the data, and also
        as confirmation that the data contains what is expected.

        Note that the [()] construct is essential for accessing the data
        
        Parameters
        ----------
        infilename: str (no default)
        
        Returns
        -------
        nothing
        """
        
        infile = Path(infilename).with_suffix(".nwb")
        if infile.is_file():
            try:
                io = NWB.NWBHDF5IO(str(infile), "r")
                nwbdata = io.read()

                notes = nwbdata.notes
                vcs = nwbdata.acquisition["CCData"]
                data = vcs.data[()]
                ANSpikes = nwbdata.stimulus["ANSpikeTimes"]
                anspikeorigin = ANSpikes.data[()]
                anspiketimes = ANSpikes.timestamps[()]
            except:
                print("Error reading data file")
            finally:
                io.close()

        print(vcs)

        # print(notes)  # holds most of the run parameters in a text string.

        timebase = np.arange(0, vcs.rate * data.shape[1], vcs.rate)
        self.plotNWB(
            timebase, anspiketimes, anspikeorigin, data, coincidence_window=0.25
        )

    def plotNWB(
        self,
        timebase,
        anspiketimes,
        anspikeorigin,
        data,
        coincidence_window: float = 0.05,
    ):
        """
        Plot the data - used to plot from the NWB file, but is otherwise
        pretty general.
        """
        f, axl = mpl.subplots(2, 1, sharex=True)
        for i in range(data.shape[0]):  # plot traces for every section
            axl[0].plot(timebase, data[i, :], linewidth=0.35)

        ost = np.argsort(anspiketimes)  # ordered sort times
        anspikeorigin = anspikeorigin[ost]
        anspiketimes = anspiketimes[ost]
        ninputs = np.max(anspikeorigin)

        coincident = np.where(np.diff(anspiketimes) < coincidence_window)[0]
        noncoinc = [
            i
            for i, x in enumerate(anspiketimes)
            if x not in coincident and x not in coincident + 1
        ]

        axl[1].scatter(
            anspiketimes[noncoinc],
            anspikeorigin[noncoinc],
            c="k",
            s=25,
            marker="|",
            linewidth=1,
            facecolor="k",
        )

        # plot AN spikes that occur together within a close coincidence window (like Joris and Smith)
        axl[1].scatter(
            anspiketimes[coincident],
            anspikeorigin[coincident],
            c="r",
            s=25,
            marker="|",
            linewidth=2,
            facecolor="r",
        )
        axl[1].scatter(
            anspiketimes[coincident + 1],
            anspikeorigin[coincident + 1],
            c="r",
            s=25,
            marker="|",
            linewidth=2,
            facecolor="r",
        )

        xlim = axl[0].get_xlim()
        axl[1].set_xlim(xlim)
        axl[0].set_ylabel("V(mv)", fontsize=9)
        axl[1].set_ylabel("Input #", fontsize=9)
        axl[1].set_xlabel("T (ms)")
        axl[0].set_title("All section voltages")
        axl[1].set_title(
            f"Input Spikes and Coincidences (dt={coincidence_window:.2f} ms)"
        )

        for ax in axl:
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)


if __name__ == "__main__":
    fn = Path(
        "VCN_c17/Simulations/AN/runANPSTH-all-2021-11-30.12-47-40/AN_Result_VCN_c17_inp=self_XM13A_nacncoop_II_HF=VCN_c17_Full_MeshInflate_normal_delays_multisite_001_tonepip_010dB_16000.0_HS.p"
    )
    SNWB = SaveNWB()
    SNWB.load_simulation(fn)
    SNWB.plotDendriteVoltages()
    outfile = "test"
    SNWB.writeNWB(fn, outfile=outfile)
    SNWB.readNWB(outfile)

    mpl.show()