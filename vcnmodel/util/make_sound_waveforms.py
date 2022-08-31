 
""" make_sound_waveform.py

 Generate a sound waveform - with different choices of sound. 
 This was taken from model_run2.py, and is a stand-alone version,
 and the parameters are taken from model_params.py

 It is used by vcnmodel.plotters.plot_runs and
 vcnmodel.plotters.plot_sims to show the stimulus waveforms
 associated with a simulation run.

 """
from dataclasses import dataclass, field
from typing import Tuple

import matplotlib.pyplot as mpl
import numpy as np
from cnmodel.util import sound

soundChoices = [
        "tonepip",
        "noise",
        "stationaryNoise",
        "regularClicks",
        "poissonClicks",
        "SAM",
        "CMMR",
    ]

def defstarts():
    return [0.1]

def defmaskerstarts():
    return [0.12]

def defemptydict():
    return {}


def defemptylist():
    return []

@dataclass
class SoundParams:
    # sound parameters
    initialization_time: float = 50.0  # nominal time to let system settle, in msec
    run_duration: float = 0.35  # in sec
    soundtype: str = "tonepip"  # or 'tonepip' or any of the soundChoices
    noise_seed: int = 9  # noise generator seed value
    pip_duration: float = 0.1  # duration in seconds for tone pip
    pip_start: list = field(default_factory=defstarts)  # start (delay to start of pip)
    pip_offduration: float = 0.05  # time after pip ends to keep running

    clickStart: float = 0.1
    clickDuration: float = 1e-4
    clickTrainDuration: float = 0.8
    clickRate: float = 20.0  # Hz  for regular clicks, or mean poisson rate

    masker_start: list=field(default_factory=defmaskerstarts)
    masker_duration: float = 0.1
    Fs: float = 100e3  # cochlea/zilany model rate
    F0: float = 16000.0  # stimulus frequency
    dB: float = 30.0  # in SPL
    RF: float = 2.5e-3  # rise-fall time
    fmod: float = 50  # hz, modulation if SAM
    dmod: float = 100.0  # percent if SAM
 
    signalToMasker: float = 12.0
    CMMRmode: str = "CM"

 
class Make_Sound_Waveform():
    def __init__(self, params=None):
        if params is None:
            self.sound_params = SoundParams()  # use default
        else:
            self.sound_params = params

    def set_sound(self, soundtype:str="tonepip"):
        if soundtype not in soundChoices:
            raise ValueError(f"Sound type {soundtype:s} is not a valid sound choice.")
        self.sound_params.soundtype = soundtype

    def set_dbspl(self, signal, dbspl):
        """Scale the level of `signal` to the given dB_SPL."""
        p0 = 20e-6
        rms = np.sqrt(np.sum(signal ** 2) / signal.size)
        scaled = signal * 10 ** (dbspl / 20.0) * p0 / rms
        return scaled

    def create_sound(self, j:int=0):
        """create the dound

        Parameters
        ----------
        j : int (default = 0)
            index that can be used to control generation of different noise
            samples.

        Raises
        ------
        ValueError
            _description_
        """
        if self.sound_params.soundtype not in soundChoices:
            raise ValueError(f"Sound type {self.sound_params.soundtype:s} is not a valid sound choice.")
        if isinstance(self.sound_params.pip_start, float):
            pips = [self.sound_params.pip_start]
        if self.sound_params.soundtype == "tonepip":
            stim = sound.TonePip(
                rate=self.sound_params.Fs,
                duration=self.sound_params.run_duration,
                f0=self.sound_params.F0,
                dbspl=self.sound_params.dB,
                ramp_duration=self.sound_params.RF,
                pip_duration=self.sound_params.pip_duration,
                pip_start=self.sound_params.pip_start,
            )
        elif self.sound_params.soundtype in ["stationaryNoise", "noise"]:
            if (
                self.sound_params.soundtype == "noise"
            ):  # non-stationary noise generator seed changes on per-run basis
                self.sound_params.noise_seed = self.sound_params.noise_seed + j
            print(
                f" **** Noise type: {self.sound_params.soundtype:s}  seed={self.sound_params.noise_seed}"
            )
            stim = sound.NoisePip(
                rate=self.sound_params.Fs,
                duration=self.sound_params.run_duration,
                dbspl=self.sound_params.dB,
                pip_duration=self.sound_params.pip_duration,
                pip_start=self.sound_params.pip_start,
                ramp_duration=self.sound_params.RF,
                seed=self.sound_params.noise_seed,
            )
        elif self.sound_params.soundtype == "SAM":
            stim = sound.SAMTone(
                rate=self.sound_params.Fs,
                duration=self.sound_params.run_duration,
                f0=self.sound_params.F0,
                dbspl=self.sound_params.dB,
                ramp_duration=self.sound_params.RF,
                fmod=self.sound_params.fmod,
                dmod=self.sound_params.dmod,
                pip_duration=self.sound_params.pip_duration,
                pip_start=self.sound_params.pip_start,
            )
        elif self.sound_params.soundtype in ["regularClicks", "poissonClicks"]:
            if self.sound_params.soundtype == "poissonClicks":
                eventintervals = np.random.exponential(
                    1.0 / self.sound_params.clickRate,
                    int(self.sound_params.clickTrainDuration * self.sound_params.clickRate),
                )
                events = np.cumsum(eventintervals)
                events = events[events < self.sound_params.clickTrainDuration]
            else:
                events = np.linspace(
                    self.sound_params.clickStart,
                    self.sound_params.clickTrainDuration,
                    int(self.sound_params.clickTrainDuration * self.sound_params.clickRate),
                )
            stim = sound.ClickTrain(
                rate=self.sound_params.Fs,
                duration=self.sound_params.clickTrainDuration,
                dbspl=self.sound_params.dB,
                click_duration=self.sound_params.clickDuration,
                click_starts=events,
            )

        elif self.sound_params.soundtype in ["CMMR"]:
            stim=sound.ComodulationMasking(rate=self.sound_params.Fs,
                duration=self.sound_params.run_duration,
                dbspl=self.sound_params.dB,
                pipst=self.sound_params.pip_start,
                pipdu=self.sound_params.pip_duration,
                maskst=self.sound_params.masker_start,
                maskdu=self.sound_params.masker_duration,
                rf=self.sound_params.RF,
                f0=self.sound_params.F0,
                s2n=self.sound_params.signalToMasker,
                fmod=self.sound_params.fmod,
                dmod=self.sound_params.dmod,
                fltype="Tone",
                flspc=1/3.,
                flgap=1,
                flph="Comodulated",
                flspl=25.,
                flN=3,
            )

        else:
            raise ValueError(
                "RunInfo sound type %s not implemented" % self.sound_params.soundtype
            )
        self.stimWaveform = stim.generate()
        self.stimTimebase = stim.time

    def show_waveforms(self, ax:object=None):
        if ax is None:
            fig, ax = mpl.subplots(1, 1)
        ax.plot(self.stimTimebase, self.stimWaveform)
        mpl.show()

    
if __name__ == "__main__":
    MSW = Make_Sound_Waveform()
    for sc in soundChoices:
        MSW.set_sound(sc)
        MSW.create_sound()
        MSW.show_waveforms()
