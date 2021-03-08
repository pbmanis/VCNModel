import dataclasses
from dataclasses import dataclass, field
import numpy as np
import scipy.stats

@dataclass
class VSResult:
    vs : float = np.nan  # vector strength
    n_spikes : int = 0  # number of spikes
    Rayleigh : float = np.nan  # coefficient
    pRayleigh : float = 1.0 # p value for NOT flat
    circ_phase : float = np.nan  # mean phase (circularized phase. radians)
    circ_phaseMean : float=np.nan
    circ_phaseSD: float = np.nan  # SD in phase
    circ_timeSD: float = np.nan # SD in time units
    frequency : float = np.nan  # modulation or VS measure frequency

    # carrier_frequency: float = np.nan  # carrier frequency
    # dmod: float = 0.  # modulation depth
    # dB : float = 0. # stimulus intensity
    # n_inputs : int = 0  # number of inputs to cell
    # an_vs : float = np.nan  # vs for an input (placeholder)
    # an_circ_phaseMean: float=np.nan
    # an_circ_phaseSD: float=np.nan
    # an_circ_timeSD: float = np.nan


class VectorStrength():
    def __init__(self):
        pass
    
    def vector_strength(self, spikes, freq):
        """
        Calculate vector strength and related parameters from a spike train, for the specified frequency
        :param spikes: Spike train, in sec.
            If the data comes from repeated trials, the spike train needs to be flattened into a 1d array before 
            calling. 
        :param freq: Stimulus frequency in Hz
        :return: a dataclass of type VSResult:
        """
        VSR = VSResult()
        amax = 0.
        if spikes is None:
            print('NO spikes')
        
        VSR.frequency = freq
        period = 1.0/freq
        VSR.n_spikes = len(spikes)
        twopi_per  = 2.0*np.pi/period
        phasev = twopi_per*np.fmod(spikes, period)  # convert to time within a cycle period in radians
        VSR.circ_phase = phasev
        sumcos = np.sum(np.cos(phasev))
        sumsin = np.sum(np.sin(phasev))
        mean_phase = np.arctan2(sumsin, sumcos)
 
        # print('circ phase, mean: ', VSR.circ_phase, mean_phase)
        sumc2 = sumcos*sumcos
        sums2 = sumsin*sumsin

        if VSR.n_spikes == 0:
            return VSR
        VSR.vs = (1./VSR.n_spikes)*np.sqrt(sumc2+sums2)  # standard vector strength computation
        # print('VSR nsp @ modf: ', VSR.frequency,  VSR.n_spikes, 'VS: ', VSR.vs, 'sumc2, s2: ', sumc2, sums2)
        VSR.Rayleigh = VSR.n_spikes*VSR.vs*VSR.vs  # Raleigh coefficient (Ashida et al, 2010 and many others: note some papers report 2nvs^2)
        VSR.pRayleigh = np.exp(-VSR.n_spikes*VSR.vs*VSR.vs)  # p value for n > 50 (see Ashida et al. 2010).
        VSR.circ_phaseMean = scipy.stats.circmean(phasev)
        VSR.circ_phaseSD =  scipy.stats.circstd(phasev) # radians
        VSR.circ_timeSD = VSR.circ_phaseSD/(2*np.pi*freq)  # convert to time
        self.spikes = spikes
        # index_row = self.parent.selected_index_rows[0]
        # selected = self.parent.table_manager.get_table_data(
        #     index_row
        # )
        return VSR
        
    def print_VS(self, d):

        self.VS_colnames =  f"Cell,Configuration,carrierfreq,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,AN_phase,AN_phasesd,maxArea,ninputs"
        
        line += f"{d.carrier_frequency:.1f},{d.carrier_frequency:.1f},{d.dmod:.1f},{d.dB:.1f},"
        line += f"{d.vs:.4f},"
        line += f"{d.n_spikes:d},"
        line += f"{d.circ_phaseMean:.4f},"
        line += f"{d.circ_phaseSD:.4f},"
        line += f"{d.Rayleigh:.4f},"
        line += f"{d.pRayleigh:.4e},"
        line += f"{d.an_vs:.4f},"
        line += f"{d.an_circ_phaseMean:.4f},"
        line += f"{d.an_circ_phaseSD:.4f},"
        line += f"{amax:.4f},"
        line += f"{d.n_inputs:d}"
        print(line)


    def vs_test__print(self, d):
        print(f"{'Vector Strength':>24s}: {d['r']:12.4f}")
        print(f"{'Spike Count':>24s}: {d['n']:12d}")
        print(f"{'mean phase (radians)':>24s}: {scipy.stats.circmean(d['ph'])/(2*np.pi):12.4f}")
        print(f"{'SD phase':>24s}: {scipy.stats.circstd(d['ph'])/(2.*np.pi):12.4f}")
        print(f"{'Rayleigh':>24s}: {d['R']:12.4f}")
        print(f"{'p value':>24s}: {d['p']:12.4e}")
        print(f"{'d':>24s}: {d['d']:12.4f}")
    
    