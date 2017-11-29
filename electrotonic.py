#!/usr/bin/env python
from __future__ import print_function
__author__ = 'pbmanis'

"""
Neuron electrotonic structure determination based on charging time
constants in response to a small current injection.

Model of the cell is the "somatic shunt" model (sse Ref 1).

References:
The implementation here is based on these works, but there is a much
larger literature (cited in these works).

1: White JA, Young ED, Manis PB. The electrotonic structure of regular-spiking
neurons in the ventral cochlear nucleus may determine their response properties. 
J Neurophysiol. 1994 May;71(5):1774-86. PubMed PMID: 8064348.

2: White JA, Manis PB, Young ED. The parameter identification problem for the
somatic shunt model. Biol Cybern. 1992;66(4):307-18. PubMed PMID: 1312875.

Requirements:
Beyond a regular scientific Python stack (numpy, scipy, matplotlib):

lmfit 0.96 or later: for curve fitting routines

"""

import os
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.cm as cm
import lmfit
from lmfit import  Model
from cnmodel.util import ExpFitting

# note need to define functions to be fit outside of class - 
# otherwise lmfit gets confused about the arguments (self...)
# 
def doubleexp(x, amp=1., C0=5, C1=1, tau0=20., tau1=0.2):
    """
    Double exponential fit
    amp: I0*Rn
    C0, C1: relative amplitudes of exponential functions
    tau1, tau2 : time constants
    """
    return (amp + C0*np.exp(-x/tau0) + C1*np.exp(-x/tau1))


class Electrotonic():
    
    def __init__(self, I, V, T, tstart):
        """
        Parameters
        ----------
        I : numpy array of floats (1D)
            Current injection waveform, scaled in nA
        V : numpy array of floats
            Voltage waveform recorded at soma, scaled in mV
        T : numpy array of floats
             Time points for the current and voltage arrays
        tstart: float
            Time at which the voltage step starts

        """
        self.V = V
        self.I = I
        self.time = T
        self.stepstart = tstart
        self.itstart = np.argmin(np.fabs(self.time-self.stepstart))
        self.Cm = 1.0
        self.Rn = 50.
        self.taum = 10.
        self.I0 = -0.1
        self.L = None

    def nexp(self, I0, Rn, C=np.zeros(2),
             tau=np.array([10., 0.5]), N=2, noise=0.):
        """
        paremeters:
        N : int (default: 2)
            number of exponential terms
        I0 : current step size
        Rn : input resitance
        
        V = I0 * Rn + sum over i (C[i]*exp(-t/tau[i]))
        
        """
        v = self.V
        t = self.time[self.itstart:] - self.stepstart
        csum = np.zeros_like(t)
        for i in range(N):
            csum = csum + C[i]*np.exp(-t/tau[i])
        v[self.itstart:] = v[self.itstart:] + I0 * Rn + csum
        if noise > 0:
            v = v + np.random.randn(v.shape[0])*noise
        return v

    def fitexps(self, x, y, amp=5):
        dexpmodel = Model(doubleexp, independent_vars=['x'])
        params = dexpmodel.make_params(amp=amp, C0=5, C1=1, tau0=5, tau1=0.5)
        print ('Params: ', [p for p in params.values()])
        self.fitresult = dexpmodel.fit(y, params, x=x)
        
    def coeffs_ratio(self):
         C1 = self.fitresult.params['C1'].value
         C0 = self.fitresult.params['C0'].value
         tau0 = self.fitresult.params['tau0'].value
         tau1 = self.fitresult.params['tau1'].value
         return(C1*tau0/(2.0*C0*tau1))

    def print_pars(self):
        print('---------------------------')
        print('Rn: {0:.3f}'.format(self.Rn))
        print('tau0 (Rn*Cm): {0:.3f}'.format(self.Rn*1e6*self.somaAreaMeasured*self.Cm*1e-8))
        print('alphas: {0:.3f}, {1:.3f}'.format(self.alphas()[0], self.alphas()[1]))
        print('eta: {0:.6f}'.format(self.eta()))
        print('rho: {0:.3f}'.format(self.rho))
        
        print('ASt: {0:.3f}'.format(self.ASt()))
        print('Ct: {0:.3f}, {1:.3f}'.format(self.Ct(0), self.Ct(1)))
        
    def alphas(self):
        """
        eqn 10, White et al. 1992
        """
        a = np.zeros(2)
        for i in range(0, 2):
            a[i] = np.sqrt((self.tau_d/
                        self.fitresult.params['tau%d'%i].value)-1.0)
        return a
        
    def eta(self, alpha=None):
        if alpha is None:
            alpha = self.alphas()
        a = alpha
        at1 = a[1]*np.tan(a[1]*self.L)
        at0 = a[0]*np.tan(a[0]*self.L)
        eta = at1-at0
        etan = (1 + a[0]*a[0])*at1 + (1 + a[1]*a[1])*at0
        eta = eta/etan
        return eta
    
    def ASt(self):
        """
        eq. 12, white et al. 1994
        """
        return self.eta()*self.tau_d/(self.Cm*self.Rn*(self.rho+1.0))

    def betas(self, alpha=None):
        if alpha is None:
            alpha = self.alphas()
        a = alpha
        betas = np.zeros(2)
        for i in range(0,2):
            print(a[i])
            betas[i] = [1.0 - self.eta(alpha=a)*(1.0+a[i]*a[i])]/(a[i]*a[i])
        return betas
    
    def Ct(self, i):
        """ calulate theroetical Ct for a give wet of values
        note: coth = 1/tanh
        """
        a = self.alphas()
        b = self.betas(alpha=a)
        tau0 = self.fitresult.params['tau%d'%0].value
        x1 = -self.I0*self.Rn*(2*(self.rho+1.0)*tau0/self.tau_d)
        x1n = b[i] + 2.0*self.eta()
        x1n = x1n + (np.power(a[i]*b[i]*self.L, 2.0)
                     /(self.rho*self.L/np.tanh(self.L)))
        return x1/x1n

    def Err(self, wa=1.0):
        """
        Eqn. 11, White et al. 1994
        """
        c0m = self.fitresult.params['C0'].value
        c1m = self.fitresult.params['C1'].value
        ASm = self.somaAreaMeasured
        ASt = self.ASt()
        c0t = self.Ct(0)
        c1t = self.Ct(1)
        
        err1 = (c0m - c0t)*(c0m - c0t)
        err2 = (c1m - c1t)*(c1m - c1t)
        err3 = (c0m + c1m)*(c0m + c1m)
        errw = wa*(ASm-ASt)*(ASm-ASt)/(ASm*ASm)
        return errw + (err1 + err2)/err3


if __name__ == '__main__':
    rate = 0.01
    dur = 50.
    npts = int(dur/rate)
    tstart = 5.
    I0 = -0.1  # nA
    Rn = 50.  # Mohm
    I = np.zeros(npts)
    T = np.arange(npts)*rate
    V = np.zeros(npts) - 60.
    C_Ratios = np.array([5., 1.])
    C = -(I0*Rn)*C_Ratios/np.sum(C_Ratios)  # scale ratios to voltage
    print('C: {0:.3f}, {1:.3f}'.format(C[0], C[1]))
    E = Electrotonic(I, V, T, tstart)
    
    ve = E.nexp(I0, Rn, C=C, noise=0.1)
    tf = T[E.itstart:]-tstart
    yf = ve[E.itstart:]
    E.fitexps(tf, yf)
    print ('Fitresult:\n', E.fitresult.fit_report())
    print ('Coeffs Ratio: ', E.coeffs_ratio())
    # mpl.plot(T, ve)
    # mpl.plot(tf+tstart, E.fitresult.best_fit, 'r--')
    # mpl.show()
    
    E.L = 1.0
    E.rho = 1.0
    E.tau_d = 20
    E.somaAreaMeasured = E.taum
    npx = 20
    lr = np.arange(npx)*(3./npx) + 0.01
    tr = np.arange(npx)*(100./npx) + E.fitresult.params['tau%d'%0].value
    esurf = np.zeros((lr.shape[0], tr.shape[0]))
    for i, l in enumerate(lr):
        for j, t in enumerate(tr):
            E.L = l
            E.tau_d = t
            esurf[i,j] = E.Err()
    x, y = np.meshgrid(lr, tr)
    E.print_pars()
    mpl.figure()
    mpl.imshow(esurf, interpolation='bilinear', origin='lower',
                cmap=cm.gray, extent=(np.min(lr), np.max(lr),
                 np.min(tr), np.max(tr)), aspect='auto')
    cs = mpl.contour(x, y, esurf)
    mpl.clabel(cs, inline=1, fontsize=8)
    mpl.show()
    
    
        
        