"""
GenerateRun is a class that sets up for a run after the cell has been decorated with channels.
It requires the celltype, and a section where the electrode will be inserted.
The stimulus in current clamp can consist of single pulses or pulse trains. cd is the
channelDecorator, which is also used to set the current range level for IV's.
Code for reading externally generated spike trains from a file is also included.
Methods:
    
    * do_run(filename) will execute the run. The resulting plot will have the filename text at the top.
    * do_run calls 3 private methods, _prepare_run, _initializeRun, and _execute_run, which respectively,
    * set up the stimuli, initialize the state of the model, and then run the model, generating a plot.

Original version: February 2014, Paul B. Manis UNC Chapel Hill.
    

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2014- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import copy
import csv
import dataclasses
import errno
import multiprocessing as MP
import os
import pickle
import signal
import time
from collections import OrderedDict
from pathlib import Path
from typing import Union

import numpy as np
import pyqtgraph as pg
from cnmodel.util.stim import make_pulse  # makestim
from pylibrary.tools import cprint as CP
from pyqtgraph import multiprocess as mproc
from pyqtgraph.Qt import QtGui

# from vcnmodel import NoiseTrainingGen as NG
from vcnmodel import cellInitialization as cellInit
from vcnmodel.analyzers import analyze_run as ar
from vcnmodel.plotters import IVPlots as IVP

__author__ = "pbmanis"


verbose = False  # use this for testing.
cprint = CP.cprint


class GenerateRun:
    def __init__(
        self,
        Params,
        RunInfo,
        cell,
        idnum=0,
        starttime=None,
    ):
        """Set up for a simulation run

        Parameters
        ----------
        Params : dataclass
            General parameters for the run, defined in model_params.py
        RunInfo : dataclass
            Specific parameters for a simulation run, defined in model_params.py
        cell : object
            cnmodel cell instance that will be used as the scaffold.
        idnum : int, optional
            a run id number, by default 0
        starttime : datetime, optional
            starting time for the run for timing, by default None
        """    
        self.run_initialized = False
        self.Params = Params  # make available to the rest of the classself.RunInfo.
        self.RunInfo = RunInfo
        # print(dir(cell))
        self.cell = cell  # cnmodel cell instance
        self.hf = cell.hr  # get the reader structure and the hoc pointer object locally
        self.basename = "basenamenotset"
        print(" save all sections: ", self.Params.save_all_sections)
        self.idnum = idnum
        self.startTime = starttime

        esecs = list(self.hf.sec_groups[self.RunInfo.electrodeSectionName])
        electrodeSection = esecs[0]
        self.electrode_site = self.hf.get_section(electrodeSection)
        if self.RunInfo.dendriticElectrodeSection is not None:
            print(
                "GenerateRun: Dendrite electrode section: ",
                self.RunInfo.dendriticElectrodeSection,
            )
            dend_sections = list(
                self.hf.sec_groups[self.RunInfo.dendriticElectrodeSection]
            )
            # invert the mappint of sections
            # first, just get the dendrite components, making a dendrite submap
            dendDistMap = {}
            for k in dend_sections:
                dendDistMap[k] = self.hf.distanceMap[k]
            revmap = dict((v, k) for k, v in dendDistMap.items())

            # now find the distlist key that corresponds to the closest value to our desired value
            (num, dist) = min(
                enumerate(revmap),
                key=lambda x: abs(x[1] - self.RunInfo.dendriticSectionDistance),
            )
            print(
                "Monitoring dendrite section number {:d}, at {:6.2f} microns from soma: ".format(
                    num, dist
                )
            )

            dendriticElectrodeSection = list(
                self.hf.sec_groups[self.RunInfo.dendriticElectrodeSection]
            )[num]
            self.dendriticElectrodeSite = self.hf.get_section(dendriticElectrodeSection)

        if self.RunInfo.inFile is None:
            self.RunInfo.tstop = (
                1000.0 * (self.RunInfo.nStim - 1) / self.RunInfo.stimFreq
                + self.RunInfo.stimPost
                + self.RunInfo.stimDelay
            )

        else:
            print("Reading spike train from CSV file")
            maxt = 0.0
            with open(self.RunInfo.inFile, "r") as csvfile:
                spks = csv.reader(csvfile, delimiter=",")
                for i, row in enumerate(spks):
                    if i == 0:  # first line
                        # print(row)
                        maxt = float(row[1]) * 1000
                        continue
                    if int(row[0]) in self.RunInfo.spikeTimeList.keys():
                        self.RunInfo.spikeTimeList[int(row[0])].append(
                            float(row[1]) * 1000.0
                        )
                    else:
                        self.RunInfo.spikeTimeList[int(row[0])] = [
                            float(row[1]) * 1000.0
                        ]
                self.RunInfo.tstop = maxt

        self.monitor = OrderedDict()  # standard monitoring
        self.allsecVec = OrderedDict()  # all section monitoring
        self.mons = OrderedDict()

        if self.Params.verbose:
            print("run-off initialization done")

    def mkdir_p(self, path):

        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _make_filename(self, filename=None, subdir=None):
        """
        Make the filename and a storage directory for the data

        Parameters
        ----------
        filename : str default=None
            The base name of the file to create
        subdir : subdirectory name

        Returns
        -------
        nothing
        """
        self.RunInfo.filename = filename
        if self.startTime is None:
            self.dtime = time.strftime("%y.%m.%d-%H.%M.%S")
        else:  # convert time object
            self.dtime = time.strftime(
                "%y.%m.%d-%H.%M.%S", time.localtime(self.startTime)
            )
        self.mkdir_p(self.RunInfo.folder)
        # f os.path.exists(self.RunInfo.folder) is False:
        #    os.mkdir(self.RunInfo.folder)
        if subdir is not None:
            folder = os.path.join(self.RunInfo.folder, subdir)
        else:
            folder = self.RunInfo.folder
        if os.path.exists(folder) is False:
            os.mkdir(folder)
        # self.basename = os.path.join(folder, filename + self.dtime)
        self.basename = os.path.join(folder, filename)

    def __final_preparation(self, maxt):
        self.hf.h.tstop = maxt + self.RunInfo.stimDelay
        print("PrepareRun, finalprep:\n ")
        print(f"   maxt:     {maxt:8.2f} ms")
        print(f"   delay:    {self.RunInfo.stimDelay:8.2f} ms")
        print(f"   duration: {self.RunInfo.stimDur:8.2f} ms")
        print(f"   tstop:    {self.hf.h.tstop:8.2f} ms")
        print(f"   h.t:      {self.hf.h.t:8.2f} ms")
        print(f"   h.dt:    {self.hf.h.dt:8.3f} ms")
        print("\n----------------\n")
        self.monitor["time"].record(self.hf.h._ref_t)

    def _prepare_run(self, inj=None):
        """
        (private method)
        Control a single run of the model with updated display of the voltages, etc.
        Inputs: inj: override for current injection
        Outputs: None
        Actions: optionally displays the results
        Side Effects: A number of class variables are created and modified, mostly related to the
        generation of stimuli and monitoring of voltages and currents
        """
        cprint("w", f"_prepare_run: {self.RunInfo.postMode:s}, inj={str(inj):s}")
        if self.Params.verbose:
            print("_prepare_run")
        for group in self.hf.sec_groups.keys():  # get morphological components
            if not self.Params.save_all_sections:  # just save soma sections
                if group.rsplit("[")[0] == "soma":
                    self.allsecVec["soma"] = self.hf.h.Vector()
                    section = list(self.hf.sec_groups[group])[0]
                    sec = self.hf.get_section(section)
                    self.allsecVec["soma"].record(sec(0.5)._ref_v, sec=sec)
                    break  # only save the FIRST occurance.
            else:
                g = self.hf.sec_groups[group]
                for section in list(g):
                    sec = self.hf.get_section(section)
                    self.allsecVec[sec.name()] = self.hf.h.Vector()
                    self.allsecVec[sec.name()].record(
                        sec(0.5)._ref_v, sec=sec
                    )  # recording of voltage all set up here

        for var in [
            "time",
            "postsynapticV",
            "dendriteV",
            "postsynapticI",
            "i_stim0",
            "v_stim0",
        ]:  # get standard stuff
            self.monitor[var] = self.hf.h.Vector()

        self.RunInfo.celsius = self.cell.status["temperature"]
        self.hf.h.celsius = self.RunInfo.celsius
        self.clist = {}  # color list
        somasecs = list(self.hf.sec_groups["soma"])
        electrodeSection = somasecs[0]
        self.electrode_site = self.hf.get_section(electrodeSection)

        # electrodeSection = list(self.hf.sec_groups[self.RunInfo.electrodeSection])[0]
        # self.electrode_site = self.hf.get_section(electrodeSection)
        if self.RunInfo.postMode == "VC":
            # Note to self (so to speak): the hoc object returned b
            # by this call must have a life after
            # the routine exits. Thus, it must be "self." Same for the IC stimulus...
            cprint(
                "c", f"Running VClamp experiment, holding={self.RunInfo.vstimHolding}"
            )
            self.hf.h.dt = self.Params.dtVC
            self.vcPost = self.hf.h.SEClamp(
                0.5, sec=self.electrode_site
            )  # self.hf.sections[electrode_site])
            self.vcPost.dur1 = self.RunInfo.vstimDelay
            self.vcPost.amp1 = self.RunInfo.vstimHolding
            self.vcPost.dur2 = 1e9
            self.vcPost.amp2 = (
                self.RunInfo.vstimHolding
            )  # just a tiny step to keep the system honest
            self.vcPost.dur3 = self.RunInfo.vstimPost
            self.vcPost.amp3 = self.RunInfo.vstimHolding
            self.vcPost.rs = 1.0  # value is in meghoms
            stim = {}
            stim["NP"] = self.RunInfo.vnStim
            stim["Sfreq"] = self.RunInfo.vstimFreq  # stimulus frequency
            stim["delay"] = self.RunInfo.vstimDelay
            stim["dur"] = self.RunInfo.vstimDur
            stim["amp"] = self.RunInfo.vstimInj
            stim["PT"] = 0.0
            stim["dt"] = self.hf.h.dt
            if inj is not None:
                stim["amp"] = inj
            else:
                stim["amp"] = self.RunInfo.stimVC["pulse"][0]
            (secmd, maxt, tstims) = make_pulse(stim)
            self.stim = stim
            secmd = secmd + self.RunInfo.vstimHolding  # add holding
            self.monitor["v_stim0"] = self.hf.h.Vector(secmd)
            self.monitor["v_stim0"].play(
                self.vcPost._ref_amp2, self.hf.h.dt, 0, sec=self.electrode_site
            )
            self.monitor["postsynapticV"].record(
                self.electrode_site(0.5)._ref_v, sec=self.electrode_site
            )
            if self.RunInfo.dendriticElectrodeSection is not None:
                self.monitor["dendriteV"].record(
                    self.dendriticElectrodeSite(0.5)._ref_v,
                    sec=self.dendriticElectrodeSite,
                )
            self.monitor["postsynapticI"].record(
                self.vcPost._ref_i, sec=self.electrode_site
            )
            self.mons = ["postsynapticI", "v_stim0"]
            self.__final_preparation(maxt)

        elif self.RunInfo.postMode == "CC":
            cprint("c", "Running IClamp experiment")
            self.hf.h.dt = self.Params.dtIC
            stim = {}
            stim["NP"] = self.RunInfo.nStim
            stim["Sfreq"] = self.RunInfo.stimFreq  # stimulus frequency
            stim["delay"] = self.RunInfo.stimDelay
            stim["dur"] = self.RunInfo.stimDur
            if inj is not None:
                stim["amp"] = inj
            else:
                stim["amp"] = self.RunInfo.stimInj["pulse"][0]
            stim["PT"] = 0.0
            stim["dt"] = self.hf.h.dt
            (secmd, maxt, tstims) = make_pulse(
                stim
            )  # cnmodel.makestim.makestim(stim, pulsetype='square', dt=self.hf.h.dt)
            self.stim = stim
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrode_site)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor["i_stim0"] = self.hf.h.Vector(secmd)
            self.monitor["i_stim0"].play(
                self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrode_site
            )
            self.monitor["postsynapticI"].record(
                self.icPost._ref_i, sec=self.electrode_site
            )
            self.monitor["postsynapticV"].record(
                self.electrode_site(0.5)._ref_v, sec=self.electrode_site
            )
            self.mons = ["postsynapticV", "postsynapticI"]
            if self.RunInfo.dendriticElectrodeSection is not None:
                self.monitor["dendriteV"].record(
                    self.dendriticElectrodeSite(0.5)._ref_v,
                    sec=self.dendriticElectrodeSite,
                )
                self.mons.append("dendriteV")
            self.__final_preparation(maxt)

        elif self.RunInfo.postMode in ["gifnoise"]:
            cprint("c", "Running GIF Noise experiment")
            stim = {}
            self.hf.h.dt = self.Params.dtIC
            self.stim = NG.NoiseGen()
            self.stim.generator(
                dt=self.hf.h.dt,
                i0=self.RunInfo.gif_i0,
                sigma0=self.RunInfo.gif_sigma,
                fmod=self.RunInfo.gif_fmod,
                tau=self.RunInfo.gif_tau,
                dur=self.RunInfo.gif_dur,
                skew=self.RunInfo.gif_skew,
            )
            maxt = 1000.0 * self.RunInfo.gif_dur
            self.icPost = self.hf.h.iStim(0.5, sec=self.electrode_site)
            self.icPost.delay = 2
            self.icPost.dur = 1e9  # these actually do not matter...
            self.icPost.iMax = 1.0
            self.monitor["i_stim0"] = self.hf.h.Vector(self.stim[1])
            self.monitor["i_stim0"].play(
                self.icPost._ref_i, self.hf.h.dt, 0, sec=self.electrode_site
            )
            self.monitor["postsynapticI"].record(
                self.icPost._ref_i, sec=self.electrode_site
            )
            self.monitor["postsynapticV"].record(
                self.electrode_site(0.5)._ref_v, sec=self.electrode_site
            )
            self.mons = ["postsynapticV", "postsynapticI"]
            if self.RunInfo.dendriticElectrodeSection is not None:
                self.monitor["dendriteV"].record(
                    self.dendriticElectrodeSite(0.5)._ref_v,
                    sec=self.dendriticElectrodeSite,
                )
                self.mons.append("dendriteV")
            self.__final_preparation(maxt)
        else:
            print("generate_run.py, mode %s  unknown" % self.RunInfo.postMode)
            return

        # self.hf.h.topology()
        # pg.show()
        # self.hf.h('access %s' % self.hf.get_section(self.electrode_site).name())

    def _clean_neuron_objects(self):
        """
        Remove the NEURON objects from data that will be saved
        This is done by turining them into strings
        """
        if not isinstance(self.RunInfo.electrodeSection, str):
            self.RunInfo.electrodeSection = str(
                self.RunInfo.electrodeSection.name()
            )  # electrodeSection
        # self.RunInfo.electrodeSectionName = str(self.RunInfo.electrodeSection.name())
        if self.RunInfo.dendriticElectrodeSection is not None and not isinstance(
            self.RunInfo.dendriticElectrodeSection, str
        ):
            self.RunInfo.dendriticElectrodeSection = str(
                self.dendriticElectrodeSection.name()
            )  # dendriticElectrodeSection,
        # dendriticSectionDistance = 100.0  # microns.

    def do_run(
        self,
        filename: Union[str, Path, None] = None,
        parMap=None,
        save: Union[str, bool] = False,
        restore_from_file: bool = False,
        initfile: Union[str, Path, None] = None,
    ):
        """
        Perform the run/simulation.
        """
        self._clean_neuron_objects()
        cprint("w", f"do_run: {self.RunInfo.postMode:s}")

        if self.Params.verbose:
            print("generate_run::do_run")
        (p, e) = os.path.splitext(filename)  # make sure filename is clean
        self.RunInfo.filename = (
            p  # change filename in structure, just name, no extension
        )
        if parMap is None or len(parMap) == 0:
            self._make_filename(
                filename=self.RunInfo.filename
            )  # base name pluse underscore
        else:
            mstr = "_"
            for k in parMap.keys():
                if k == "id":
                    continue
                mstr += k + "_"
            # if 'id' in parMap.keys():
            #    mstr += 'ID%04d_' % parMap['id']
            self._make_filename(self.RunInfo.filename + mstr)

        if self.Params.verbose:
            print("genrate_run::do_run: basename is = {:s}".format(self.basename))
        # self.hf.update() # make sure channels are all up to date
        self.results = {}
        CP.cprint("m", self.RunInfo.postMode)
        if self.RunInfo.postMode == "CC":
            ipulses = self.RunInfo.stimInj["pulse"]
            # ipulses = np.arange(s[0], s[1], s[2])
        elif self.RunInfo.postMode == "VC":
            ipulses = self.RunInfo.stimVC["pulse"]
        else:
            ipulses = [0]
        nLevels = len(ipulses)

        if self.Params.Parallel is False:
            nWorkers = 1
        else:
            nWorkers = MP.cpu_count()  # get this automatically
        # workers
        CP.cprint("m", f"Parallel with {nWorkers:d} processes")
        # print(f"do_run: initfile = {str(initfile):s}")
        TASKS = [s for s in range(nLevels)]
        tresults = [None] * len(TASKS)
        signal.signal(
            signal.SIGCHLD, signal.SIG_DFL
        )  # might prevent OSError "no child process"
        # run using pyqtgraph's parallel support
        with mproc.Parallelize(
            enumerate(TASKS), results=tresults, workers=nWorkers
        ) as tasker:
            for i, x in tasker:
                inj = ipulses[i]
                self._prepare_run(inj=inj)  # build the recording arrays
                self.run_initialized = cellInit.init_model(
                    self.cell,
                    vinit=self.RunInfo.vstimHolding,
                    mode=self.RunInfo.postMode,
                    restore_from_file=restore_from_file,
                    filename=initfile,
                )
                tr = {"r": self._execute_run(), "i": inj}  # now you can do the run
                tasker.results[i] = tr
        cprint("c", "Parallel task finished")
        self.results = OrderedDict()

        for i in range(nLevels):
            #  print('level: ', i, '  tresults[i]["i"]: ', tresults[i]['i'])
            self.results[tresults[i]["i"]] = tresults[i]["r"]
        for k, i in enumerate(ipulses):
            if self.Params.plotFlag:
                if k == 0:
                    self.mons = list(self.results[i]["monitor"].keys())
                    self.plotRun(self.results[i], init=True)
                else:
                    self.plotRun(self.results[i], init=False)

        if self.RunInfo.postMode == "CC":
            if self.Params.verbose:
                cprint("r", "do_run, calling IV")
            self.arun = ar.AnalyzeRun(
                self.results
            )  # create an instance of the class with the data
            self.arun.IV()  # compute the IV on the data
            self.IVResult = self.arun.IVResult
            if self.Params.verbose:
                print("do_run, back from IV")
            if save == "monitor":
                self.save_runs(save="monitor")
            if self.Params.verbose:
                print("do_run: ivresult is: {:32}".format(self.IVResult))
            if self.Params.plotFlag:
                self.plotFits("Soma", self.IVResult["taufit"], c="r")
                self.plotFits("Soma", self.IVResult["ihfit"], c="b")
                # print (dir(self.ivplots))
                self.ivplts.show()

        if self.RunInfo.postMode == "VC":
            if self.Params.verbose:
                cprint("r", "do_run, calling VC")
            self.arun = ar.AnalyzeRun(
                self.results
            )  # create an instance of the class with the data
            self.arun.VC()  # compute the IV on the data
            self.VCResult = self.arun.VCResult
            if self.Params.verbose:
                print("do_run, back from VC")
            cprint("m", f"SAVE: {str(save):s}")
            if save == "monitor":
                self.save_runs(save="monitor")
            if self.Params.verbose:
                print("do_run: ivresult is: ", self.VCResult)
            if self.Params.plotFlag:
                # print (dir(self.ivplots))
                self.ivplts.show()
            # mpl.plot(self.monitor['time'], self.monitor['postsyanpticI'])

        if self.RunInfo.postMode in ["gifnoise"]:
            self.IVResult = None
            if save == "monitor":
                self.save_runs("gifnoise")

    def test_run(
        self,
        title: str = "testing...",
        level: float = -0.001,
        dur: float = 50.0,
        initfile=None,
    ):
        """
        Perform a test run
        """
        if initfile is None:
            raise ValueError("generate_run:test_run needs initfile name")
        self._prepare_run(inj=level)
        self.run_initialized = cellInit.init_model(
            self.cell,
            vinit=self.RunInfo.vstimHolding,
            filename=initfile,
            restore_from_file=True,
        )
        self.hf.h.t = 0.0
        self.hf.h.tstop = dur
        self._execute_run(testPlot=True)
        # pg.mkQApp()
        # pl = pg.plot(
        #     np.array(self.monitor["time"]), np.array(self.monitor["postsynapticV"])
        # )
        # pl.setTitle(title)
        # QtGui.QApplication.instance().exec_()

    def _execute_run(self, testPlot=False):
        """
        (private mmethod)
        After prepare run and initialization, this routine actually calls the run method in hoc
        assembles the data, saves it to disk and plots the results.
        Inputs: flag to put up a test plot....
        """
        vfile = "v.dat"
        if self.Params.verbose:
            print("_execute_run")
        assert self.run_initialized is True
        print("Starting Vm at electrode site: {:6.2f}".format(self.electrode_site.v))

        # one way
        self.hf.h.t = 0
        """
        #while (self.hf.h.t < self.hf.h.tstop):
        #    for i=0, tstep/dt {
        #       self.hf.h.fadvance()
        # self.hf.h.run()  # calls finitialize, causes offset
        """
        self.hf.h.batch_save()  # save nothing
        print("Temperature in run at start: {:6.1f}".format(self.hf.h.celsius))
        self.hf.h.batch_run(self.hf.h.tstop, self.hf.h.dt, vfile)
        print("Finishing Vm: {:6.2f}".format(self.electrode_site.v))
        if Path(vfile).is_file():
            Path(vfile).unlink()  # delete the v.dat file once we are done.
        self.monitor["time"] = np.array(self.monitor["time"])
        self.monitor["time"][0] = 0.0
        if self.Params.verbose:
            print("post v: ", self.monitor["postsynapticV"])
        if testPlot:
            pg.plot(
                np.array(self.monitor["time"]), np.array(self.monitor["postsynapticV"])
            )
            QtGui.QApplication.instance().exec_()
            # pg.mkQApp()
            # pl = pg.plot(np.array(self.monitor['time']), np.array(self.monitor['postsynapticV']))
            # if self.Params.cell is not None:
            #     pl.setTitle('%s' % self.Params.cell)
            # else:
            #     pl.setTitle('_execute_run, no filename')
        np_monitor = {}
        for k in list(self.monitor.keys()):
            try:
                np_monitor[k] = np.array(self.monitor[k])
            except: # the monitor dict element may not really have any data, so... skip
                pass

        np_allsecVec = OrderedDict()
        for k in self.allsecVec.keys():
            np_allsecVec[k] = np.array(self.allsecVec[k])
        self.RunInfo.clist = self.clist
        results = {
            "Params": self.Params,
            "RunInfo": self.RunInfo,
            "Sections": list(self.hf.sections.keys()),
            "vec": np_allsecVec,
            "monitor": np_monitor,
            "stim": self.stim,
            "runInfo": self.RunInfo,
            "distanceMap": self.hf.distanceMap,
        }
        if self.Params.verbose:
            print("    _execute_run completed")
        return results

    def save_runs(self, save: Union[str, None] = None):
        """
        Save the result of multiple runs to disk. Results is in a dictionary,
        each element of which is a Param structure, which we then turn into
        a dictionary...
        """
        fn = self.Params.simulationFilename
        cprint("g", f"\nWRITING DATA TO: {str(fn):s}\n")
        pfout = open(fn, "wb")
        mp = copy.deepcopy(self.cell.status)
        del mp["decorator"]
        self._clean_neuron_objects()  # in runinfo
        pickle.dump(
            {
                "basename": self.basename,
                "runInfo": self.RunInfo,
                "modelPars": mp,  # some specific parameters to this run
                "Params": self.Params,  # all the parameters that were passed
                "Results": self.results,  # [{k: x.todict()} for k, x in self.results.items()],
            },
            pfout,
        )

        pfout.close()
        return (
            self.RunInfo.folder,
            self.basename,
        )  # return tuple to assemble name elsewhere

        # if recordSection:
        #     fns = os.path.join(self.RunInfo.folder,
        # self.RunInfo.fileName + '_swellings_' + dtime + '.p')
        #     srsave = sr.SimulationResult()
        #     srsave.save(fns, data=np.array(self.vswel),
        #             time=np.array(self.vec['time']),
        #             hoc_file = 'MorphologyFiles/' + self.modelPars.file,
        #             section_map=secarray,
        #     )
        # pfout = open(fns, 'wb')
        # pickle.dump({'swellings': secarray,
        #              'time': np.array(self.vec['time']),
        #              'data': np.array(self.vswel)}, pfout)
        # pfout.close()
        # plot = pg.plot()
        # vs = np.array(self.vswel)
        # ts = np.array(self.vec['time'])
        # for i in range(len(self.swellAxonMap)):
        #     plot.plot(ts, vs[i])

    def plot_run(self, results, init=True, show=False):
        if init:
            self.ivplts = IVP.IVPlots(title=self.Params.cell, mode="mpl")
        self.ivplts.plotResults(
            results, dataclasses.asdict(self.RunInfo), somasite=self.mons
        )
        if show:
            self.ivplts.show()

    def plot_fits(self, panel, x, c="g"):
        self.ivplts.plotFit(panel, x[0], x[1], c)


if __name__ == "__main__":
    pass
