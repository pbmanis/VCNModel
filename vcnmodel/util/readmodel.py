""" readmodel.py - tools to read model data files.
All reading should be done through this module, as it handles different formats,
and retuns a uniform "ModelData" data structure

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2014-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import dataclasses
import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import numpy as np
import vcnmodel.model_params
import vcnmodel.util.fixpicklemodule as FPM
from ephys.ephysanalysis import MakeClamps
from pylibrary.tools import cprint as CP
from vcnmodel.analyzers import analyze_data
from vcnmodel.analyzers.reverse_correlation import RevCorrData, RevCorrPars
from vcnmodel.util import params  # for old simulation results only
from vcnmodel.util import trace_calls

TRC = trace_calls.TraceCalls
cprint = CP.cprint


def defemptydict():
    return {}


@dataclass
class ModelData:
    success: bool = False
    timestamp: str = ""  # a timestamp string
    data: dict = field(default_factory=defemptydict)  # the basic model data
    SI: object = None  # Params
    RI: object = None  # runInfo
    AR: object = None  # a readmodel instance
    SP: object = None  # a spike analysis instance
    RM: object = None  # a Rm, Tau instance (passive properties)
    RCP: object = None  # RevCorr Pars
    RCD: object = None  # RevCorr Data


@dataclass
class Inflate:
    soma_inflation: float = 0.0
    dendrite_inflation: float = 0.0


class ReadModel:
    """read model data files from two different formats, and put the
    data into a standard format for analysis. The format is very similar to that
    used by ephys ("AR")

    Returns
    -------
    None

    """

    def __init__(self, parent=None, my_parent=None):
        self.MC = MakeClamps.MakeClamps()
        self.parent = parent
        self.my_parent = None
        self.firstline = False

    def set_parent(self, parent=None, my_parent=None):
        self.parent = parent
        self.my_parent = my_parent

    def textclear(self):
        """Clear the text display region

        Raises
        ------
        ValueError
            It is invalid to clear the display if there is
            no specified parent
        """
        if self.parent is None:
            raise ValueError("parent is None")
        else:
            self.parent.textbox.clear()

    def textappend(self, text: str, color="white"):
        """append text to the console or GUI

        Parameters
        ----------
        text : str
            The text to display
        color : str, optional
            color for the text display, by default "white"
        """
        if self.parent is None:
            cprint(color, text)  # just go straight to the terminal
        else:
            self.parent.textbox.setTextColor(self.parent.QColor(color))
            self.parent.textbox.append(text)
            self.parent.textbox.setTextColor(self.parent.QColor("white"))

    def read_pfile(
        self,
        datasource: Union[str, Path, dict, None] = None,
        filemode: str = "vcnmodel.v0",
        vscale: float = 1e-3,
        iscale: float = 1e-9,
    ):
        """
        Read a pickled file; optionally plot the data
        Puts the data into a Clamps structure.

        Parameters
        ----------
        datasource : str/path/dict/None (default None)
            where the data comes from

        filemode: str
            version name for file structure. All contemporary ones are V1.0

        vscale: scaling for voltage
            float, default=1e-3
        iscale: scaling for current
            float, default=1e-9

        """

        """
        The runInfo dictionary holds somethign like this: runinfo:  {'folder':
        PosixPath('VCN_Cells/VCN_c08/Simulations/IV'), 'fileName': 'Normal',
        'runName': 'Run', 'manipulation': 'Canonical', 'preMode': 'cc',
        'postMode': 'cc', 'TargetCellType': 'Bushy', 'electrodeSection': 'soma',
        'dendriticElectrodeSection': 'dendrite', 'dendriticSectionDistance':
        100.0, 'celsius': 37, 'nStim': 1, 'stimFreq': 200.0, 'stimInj':
        {'pulse': [-1.0, 2.01, 0.2]}, 'stimDur': 100.0, 'stimDelay': 5.0,
        'stimPost': 3.0, 'vnStim': 1, 'vstimFreq': 200.0, 'vstimInj': 50,
        'vstimDur': 50.0, 'vstimDelay': 2.0, 'vstimPost': 3.0, 'vstimHolding':
        -60, 'gif_i0': 0.0, 'gif_sigma': 0.5, 'gif_fmod': 0.2, 'gif_tau': 3.0,
        'gif_dur': 10.0, 'gif_skew': 0.0, 'runTime': 'Wed Oct  9 13:05:54 2019',
        'inFile': None, 'inFileRep': 1, 'spikeTimeList': {}, 'v_init': -61.0,
        'useSaveState': True, 'tstop': 8.0, 'filename': 'VCN_c08_pulse_'}
        """
        # print('\nmodelPars: ', df['modelPars'])
        """
        The modelPars dict holds the following: modelPars:  {'species': 'mouse',
        'cellClass': 'bushy', 'modelType': 'II', 'modelName': 'mGBC', 'soma':
        True, 'axon': False, 'dendrites': False, 'pumps': False, 'hillock':
        False, 'initialsegment': False, 'myelinatedaxon': False,
        'unmyelinatedaxon': False, 'na': 'nav11', 'ttx': False, 'name': 'bushy',
        'morphology': 'VCN_Cells/VCN_c08/Morphology/VCN_c08.hoc', 'temperature':
        34.0}
        
        Note 10/28/2019 changed structure so that runInfo and modelPars are both
        subdictionaries of Params (filemode is 'vcnmodel.v0') ... and undone
        later, so that all are top-level (filemode is 'vcnmodel.v1')
        """
        if isinstance(datasource, (str, Path)):
            with open(datasource, "rb") as fh:
                df = FPM.pickle_load(fh)
                filename = datasource
        else:
            df = datasource
            print("df.Params: ", df.Params)
            raise ValueError()

        if filemode in ["vcnmodel.v0"]:
            print(
                f"Reading model file in version v0:, with {len(df['Results']):4d} trials"
            )
        elif filemode in ["vcnmodel.v1"]:
            print(
                f"Reading model file in version v1:, with {len(df['Results']):4d} trials"
            )
        else:
            raise ValueError(f"Unknown file mode: {filemode:s}")
        mtime = Path(filename).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if filemode == "vcnmodel.v0":
            try:
                dinfo = df["Params"]["runInfo"]
            except:
                try:
                    dinfo = df["runInfo"]
                except:
                    raise ValueError("Cannot read the file in v0 mode?")
            run_protocol = df["Params"]["runProtocol"]
            if isinstance(dinfo, params.Params):
                dinfo = dinfo.todict()
            dur = dinfo["stimDur"]
            delay = dinfo["stimDelay"]
            mode = dinfo["postMode"].upper()
            ntr = len(df["Results"])
            if "dt" in df["Params"].keys():
                self.rate = df["Params"]["dt"]
            else:
                self.rate = df["Params"].dt
            V = [[]] * ntr
            I = [[]] * ntr
            for i in range(len(df["Results"])):
                fk = list(df["Results"][i].keys())[0]
                dfx = df["Results"][i][fk]["monitor"]
                timebase = dfx["time"]
                V[i] = dfx["postsynapticV"] * vscale
                I[i] = dfx["i_stim0"] * iscale
        else:
            dinfo = df["runInfo"]
            x = dir(dinfo)
            if (
                "stimVC" not in x
            ):  # create fields for missing values from older versions of files.
                dinfo.stimVC = None
            mode = dinfo.postMode.upper()
            dur = dinfo.stimDur
            delay = dinfo.stimDelay
            mode = dinfo.postMode
            try:
                self.rate = df["Params"].dt  # old version, now separated IC and VC
            except:
                if mode == "VC":
                    self.rate = df["Params"].dtVC
                elif mode == "CC":
                    self.rate = df["Params"].dtIC
                else:
                    raise ValueError("Cannot find rate for data mode: ", mode)

            run_protocol = dinfo.runProtocol
            if dinfo.runProtocol in ["runIV", "initIV", "testIV", "runIVSpikeThreshold"]:
                ntr = len(df["Results"])
                V = [[]] * ntr
                I = [[]] * ntr
                for ii, i in enumerate(df["Results"].keys()):
                    dfx = df["Results"][i]["monitor"]
                    timebase = dfx["time"]
                    V[ii] = np.array(dfx["postsynapticV"]) * vscale
                    I[ii] = np.array(dfx["i_stim0"]) * iscale
            elif dinfo.runProtocol in ["runVC", "initVC", "testVC"]:
                dur = dinfo.vstimDur  # msec
                delay = dinfo.vstimDelay  # msec
                ntr = len(df["Results"])
                V = [[]] * ntr
                I = [[]] * ntr
                for ii, i in enumerate(df["Results"].keys()):
                    dfx = df["Results"][i]["monitor"]
                    timebase = dfx["time"]
                    V[ii] = np.array(dfx["postsynapticV"]) * vscale
                    I[ii] = np.array(dfx["postsynapticI"]) * iscale

            elif dinfo.runProtocol in [
                "initAN",
                "runANPSTH",
                "runANIO",
                "runANSingles",
            ]:

                # two ways data can be organized, so try both
                try:  # cnmodel_models simulations
                    ntr = len(df["Results"])
                    V = [[]] * ntr
                    I = [[]] * ntr
                    for j in list(df["Results"].keys()):
                        dfx = df["Results"][j]
                        timebase = dfx["time"]
                        V[j] = np.array(dfx["somaVoltage"]) * vscale
                        I[j] = np.zeros_like(V[j])
                except:  # vcnmodel simulatipns
                    ntr = len(df["Results"]["somaVoltage"])
                    V = [[]] * ntr
                    I = [[]] * ntr
                    for j in range(ntr):
                        timebase = df["Results"]["time"]
                        V[j] = np.array(df["Results"]["somaVoltage"][j]) * vscale
                        I[j] = np.zeros_like(V[j])

        V = np.array(V)
        I = np.array(I)

        if run_protocol in ["runVC", "initVC", "testVC"]:
            self.MC.set_clamps(
                dmode=mode, time=timebase, data=I, cmddata=V, tstart_tdur=[delay, dur]
            )
        else:
            self.MC.set_clamps(
                dmode=mode, time=timebase, data=V, cmddata=I, tstart_tdur=[delay, dur]
            )
        self.MC.getClampData()
        return self

    def _convert_params(self, data, fns, filemode):
        """
        Capture the parameters in different versions of the files,
        and return as a simple dictionary.

        Parameters
        ----------
        d : dict from the top level of the file
            Must have 'runInfo' and/or 'Params' as entries

        fns: str filename
            used to look for inflation information.

        filemode : str (no default)
            A string corresponding to the file mode for the file
            Must be either "vcnmodel.v0" or "vcnmodel.v1"
        """
        if filemode in ["vcnmodel.v0"]:
            # print(data['runInfo'].keys())
            par = data["runInfo"]
            par["soma_inflation"] = False
            par["dendrite_inflation"] = False
            if fns.find("soma=") > -1:
                par["soma_inflation"] = True
            if fns.find("dend=") > -1:
                par["dendrite_inflation"] = True
            par["soma_autoinflate"] = False
            par["dendrite_autoinflate"] = False
        elif filemode in ["vcnmodel.v1"]:
            try:
                par = data["Params"]
            except ValueError:
                try:
                    par = data["self.Params"]
                except ValueError:
                    raise ValueError("File missing Params; need to re-run")
            if isinstance(par, vcnmodel.model_params.Params):
                par = dataclasses.asdict(par)
        else:
            raise ValueError("File mode must be either vcnmodel.v0 or vcnmodel.v1")

        return par

    def _get_scaling(self, fn, PD, par_in):
        """
        Get the dendrite/soma scaling information from this data file
        and return a string. Also prints out the scaling information as we go.

        Parameters
        ----------
        fn : filename
        PD : Parameter Dataclass (no default)
        par_in : parameter list.
        Returns
        -------
        string with the type of scaling applied to the data.
        """
        stitle = "Bare Scaling"
        if isinstance(par_in, dict):
            par = Inflate(
                par_in["soma_inflation"],
                par_in["dendrite_inflation"],
            )
        if PD.soma_inflate and PD.dend_inflate:
            if par.soma_inflation > 0.0 and par.dendrite_inflation > 0.0:
                ivdatafile = Path(fn)
                stitle = "Soma and Dend scaled"
                print(stitle)
        elif PD.soma_inflate and not PD.dend_inflate:
            if par.soma_inflation > 0.0 and par.dendrite_inflation < 0.0:
                ivdatafile = Path(fn)
                stitle = "Soma only scaled"
                print(stitle)
        elif PD.dend_inflate and not PD.soma_inflate:
            if par.soma_inflation < 0.0 and par.dendrite_inflation > 0.0:
                ivdatafile = Path(fn)
                stitle = "Dend only scaled"
                print(stitle)
        elif not PD.soma_autoinflate and not PD.dendrite_autoinflate:
            print("\nConditions x: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)
            ivdatafile = Path(fn)
            stitle = "No scaling (S, D)"
            print(stitle)
        else:
            ivdatafile = Path(fn)
            self.textappend(
                f"Bare file: no identified soma/dendrite inflation conditions", "red"
            )
            stitle = "Bare Scaling"
        return "      " + stitle

    def _get_changetimestamp(self):
        # trip filemode based on date of simulation
        changedate = "2020-04-29-12:00"
        dts = datetime.datetime.strptime(changedate, "%Y-%m-%d-%H:%M")
        changetimestamp = datetime.datetime.timestamp(dts)
        return changetimestamp

    # @TRC(show=False)
    def get_data_file(
        self,
        fn: Union[str, Path],
        PD: dataclass,
        verbose=False,
    ) -> Union[None, tuple]:
        """
        Get a data file, and also parse information from the file
        for display.
        This routine does not do any analysis.

        Parameters
        ----------
        fn : str or Path
            the file to read
        PD : Parameter Dataclass (no default)
            holds information about scaling, passed to _get_scaling
        verbose: bool
            whether to print stuff to the terminal
        """
        fnp = Path(fn)
        fns = str(fn)
        ivdatafile = Path(fn)
        if self.firstline:
            if not fnp.is_file():
                cprint("r", f"   File: {str(fnp):s} NOT FOUND")
                self.textappend(f"   File: {str(fnp):s} NOT FOUND", color="red")
                return None
            else:
                if verbose:
                    self.textappend(f"   File: {str(fnp):s} OK")
        changetimestamp = self._get_changetimestamp()
        mtime = fnp.stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        if verbose and self.firstline:
            self.textappend(
                f"pgbcivr2: Checking file: {fnp.name:s} [{timestamp_str:s}]"
            )
        if mtime > changetimestamp:
            filemode = "vcnmodel.v1"
        else:
            filemode = "vcnmodel.v0"
        if self.firstline and verbose:
            self.textappend(f"pgbcivr2: file mode: {filemode:s}")
        with (open(fnp, "rb")) as fh:
            data = FPM.pickle_load(fh)
        if verbose:
            self.textappend(f"   ...File read, mode={filemode:s}")
        par = self._convert_params(data, fns, filemode)
        stitle = self._get_scaling(fn, PD, par)

        if verbose:
            print("Read file: ", ivdatafile)
        if ivdatafile is None or not ivdatafile.is_file():
            if self.firstline and verbose:
                self.textappend(f"no file matching conditions : {str(ivdatafile):s}")
            return None

        if self.firstline and verbose:
            self.textappend(f"\npgbcivr2: datafile to read: {str(ivdatafile):s}")
        if isinstance(data["Results"], dict):
            if "time" in list(data["Results"].keys()):
                data["Results"] = self._data_flip(data)
        return par, stitle, ivdatafile, filemode, data

    @TRC(show=False)
    def get_data(
        self,
        fn: Union[Path, str],
        PD: dataclass,
        protocol: str = "",
    ) -> dataclass:
        """get the data from a particular file

        Parameters
        ----------
        fn : Union[Path, str]
            filename
        PD : dataclass
            Parameter Dataclass (no default) holds info about the file
        protocol : str, optional
            Protocol as a string name by default ""

        Returns
        -------
        dataclass
            model_data dataclass structure
            model_data.success is True if the data has been read.
        """
        model_data = ModelData()  # create data structure for results

        X = self.get_data_file(fn, PD=PD)
        if X is None:
            print("No simulation found that matched conditions")
            print("Looking for file: ", fn)
            return model_data  # success flag will be false
        # unpack x
        par, stitle, ivdatafile, filemode, d = X
        if "time" in list(d["Results"].keys()):
            d["Results"] = self._data_flip(d)

        # 2. find spikes
        AR, SP, RM = analyze_data.analyze_data(ivdatafile, filemode, protocol)
        # set up analysis parameters and result storage
        RCP = RevCorrPars()
        RCD = RevCorrData()

        RCD.npost_spikes = 0  # number of postsynaptic spikes
        RCD.npre_spikes = 0  # number of presynaptic spikes

        trials = range(len(d["Results"]))
        RCP.ntrials = len(trials)

        for i, tr in enumerate(list(d["Results"].keys())):
            trd = d["Results"][tr]  # trial data
            time_base = AR.MC.time_base / 1000.0  # convert to seconds
            if "inputSpikeTimes" in list(trd.keys()):
                for n in range(len(trd["inputSpikeTimes"])):  # for each sgc
                    RCD.npre_spikes += len(trd["inputSpikeTimes"][n])
                RCD.npost_spikes += len(trd["spikeTimes"])
                RCD.st = SP.spikeIndices
                # print("TR: ", tr)
                # print("RCD.st: ", RCD.st)
                RCD.npost_spikes += len(RCD.st[tr])
            else:
                RCD.npost_spikes += sum(SP.spikes[i])

        RCD.ti = time_base
        # print(f"Detected {RCD.npost:d} Post spikes")
        # print(f"Detected {RCD.npre:d} Presynaptic spikes")
        # print("# trials: ", RCP.ntrials)

        # clip trace to avoid end effects
        RCP.max_time = (
            np.max(RCD.ti) - RCP.min_time
        )  # this window needs to be at least as long as maxwin
        RCD.tx = np.arange(RCP.minwin, 0, RCP.binw)
        mtime = Path(fn).stat().st_mtime
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
            "%Y-%m-%d-%H:%M"
        )
        model_data.success = True
        model_data.timestamp = timestamp_str
        model_data.SI = d["Params"]
        model_data.RI = d["runInfo"]
        model_data.AR = AR
        model_data.SP = SP
        model_data.RM = RM
        model_data.RCP = RCP
        model_data.RCD = RCD
        model_data.data = d
        return model_data

    def _data_flip(self, d):
        """
        Convert from old data file format to new

        Parameters
        ----------
        d : dict from the top of the data file,
            or from results (d['Results'])

        Returns
        -------
        d['Results'] with the data format adjusted.
        """
        # flip order to put trials first (this was an old format)
        data_res = d["Results"]
        trials = range(d["runInfo"].nReps)

        for tr in trials:
            sv = data_res["somaVoltage"][tr]
            dv = data_res["dendriteVoltage"][tr]
            st = data_res["spikeTimes"][tr]
            isp = data_res["inputSpikeTimes"][tr]
            if len(data_res["stimWaveform"].shape) > 0:
                swv = data_res["stimWaveform"][tr]
            else:
                swv = None
            stb = data_res["stimTimebase"][tr]

            ti = data_res["time"][tr]

            data_res[tr] = {
                "somaVoltage": sv,
                "dendriteVoltage": dv,
                "spikeTimes": st,
                "inputSpikeTimes": isp,
                "stimWaveform": swv,
                "stimTimebase": stb,
                "time": ti,
            }
        delete = [key for key in data_res[0].keys()]  # get from saved keys
        for key in delete:
            del data_res[key]  # but delete top level
        # print("_data_flip: ", data_res.keys())
        return data_res
