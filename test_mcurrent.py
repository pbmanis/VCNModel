import numpy as np
import matplotlib.pyplot as mpl
import scipy.stats
from gif.Filter_Rect_LogSpaced import *
import timeit

def alphan(v):
    return 0.0001*(v + 45)/(1. - np.exp(-(v + 45.)/9.))

def betan(v):
    return -0.0001*(v + 45)/(1. - np.exp((v + 45.)/9.))

def ngate(v):
    ntau = 1./(3.0*(alphan(v) + betan(v) ))
    ninf = alphan(v)*ntau
    return (ntau, ninf)


T = np.arange(2000)*0.1

en = -80.
dt = 0.1
gl = 0.05
gn = 0.1
C = 0.2
El = -60.
Vr = -70.
Vt_star = -48.0            # mV, steady state voltage threshold VT*
DeltaV      = 0.5  
Trefract = 1.
Trefract_ind = int(Trefract/dt)
lambda0 = 1.0

eta     = Filter_Rect_LogSpaced()    # nA, spike-triggered current (must be instance of class Filter)
gamma   = Filter_Rect_LogSpaced() 
        
T_ind = T.shape[0]

en = -80.;
r = 0.0;
n = 0.0;  # these are all local
n_inf = 0.0;
n_tau  = 1;
dndt = 0.;

print('Refract: ', Trefract_ind)
V = np.zeros_like(T)
V[0] = -65.
       
# initialize n
n_tau, n_inf = ngate(V[0])
n = n_inf
spks_i = np.zeros_like(T)
next_spike = spks_i[0] + Trefract_ind
spks_cnt = 0


(p_eta_support, p_eta) = eta.getInterpolatedFilter(dt)   
p_eta       = p_eta.astype('double')
p_eta_l     = len(p_eta)

(p_gamma_support, p_gamma) = gamma.getInterpolatedFilter(dt)   
p_gamma     = p_gamma.astype('double')
p_gamma_l   = len(p_gamma)

eta_sum = np.zeros(len(T)+2*p_eta_l)
gamma_sum = np.zeros(len(T)+2*p_gamma_l)

spks = np.zeros_like(T)
 
sigma = 20.0
tau=3.0
skew = np.sqrt(2)
#I = 4*(np.random.randn(T.shape[0])+0.5)                                                     
nrand = scipy.stats.skewnorm.rvs(skew, size=int(T.shape[0])) - 0.5
I = nrand*np.sqrt(2.0*(sigma*sigma)*dt/tau)


nt = np.zeros_like(T)
t = 0
n = 0
start_time = timeit.default_timer()

while t < T.shape[0]-1:
    # INTEGRATE VOLTAGE
    V[t+1] = V[t] + dt/C*( -gl*(V[t] - El) - gn*n*(V[t] - en) + I[t] - eta_sum[t] )

    n_tau, n_inf = ngate(V[t])

    # advance n
    dndt = (n_inf - n)/n_tau
    n = n + dndt
    nt[t] = n
    # if t < 100:
    #     print n,
    lambdan = lambda0*np.exp( (V[t+1]-Vt_star-gamma_sum[t])/DeltaV )
    p_dontspike = np.exp(-lambdan*(dt/1000.0))  # since lambda0 is in Hz, dt must also be in Hz (this is why dt/1000.0)
          

    # PRODUCE SPIKE STOCHASTICALLY
    r = np.random.randint(0, 1000)/1000.
    # if t < 100:
    #     print p_dontspike, r

    if (r > p_dontspike):
                        
        if t+1 < T_ind-1:
            spks[t+1] = 1.0 
            nt[t:t+Trefract_ind+1] = nt[t-1]  # save n state as well

        t = t + Trefract_ind   
        
        if t+1 < T_ind-1:
            V[t+1] = Vr
        
        
    #UPDATE ADAPTATION PROCESSES
    for j in range(p_eta_l):
        eta_sum[t+1+j] += p_eta[j]

    for j in range(p_gamma_l):
        gamma_sum[t+1+j] += p_gamma[j]
    
    t = t + 1

    # if t == next_spike:
    #     spks_cnt = spks_cnt + 1
    #     next_spike = spks_i[spks_cnt] + Trefract_ind
    #     V[t-1] = 0.
    #     V[t] = Vr
    #     t = t-1


        
elapsed = timeit.default_timer() - start_time
print ('Elapsed time for simulations: %f' % (elapsed))

Vmmm = -100+np.arange(100)*1.5
m_tau = 1./(3.0*(alphan(Vmmm)+betan(Vmmm)))
m_ninf = alphan(Vmmm)*m_tau

f, ax = mpl.subplots(4,1)
ax =ax.ravel()
ax[0].plot(Vmmm, m_tau, 'r')
ax[1].plot(Vmmm, m_ninf, 'k')
ax[2].plot(T, V, 'g')
ax[3].plot(T, nt, 'c')
mpl.show()