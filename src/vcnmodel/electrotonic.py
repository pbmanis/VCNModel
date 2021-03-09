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

import argparse
import re
import pandas as pd
from dataclasses import dataclass, field
from typing import Union
import io
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.cm as cm
#import lmfit
from lmfit import  Model
#from cnmodel.util import ExpFitting

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
    return (amp - (C0*np.exp(-x/tau0) + C1*np.exp(-x/tau1)))


class Electrotonic():
    
    def __init__(self, I:Union[float, None]=None,
                       V:Union[float, None]=None, 
                       T:Union[float, None]=None, 
                       tstart:float = 10., tdur:float=5.):
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
            Time at which the voltage step starts (msec)
        tdur: float
            Duration of trace to fit (msec)
        
        """
        self.V = V
        self.doubleexp = doubleexp
        self.I = I
        self.time = T
        self.stepstart = tstart
        self.stepend = tdur + tstart
        self.itstart = np.argmin(np.fabs(self.time-self.stepstart))-1
        self.itend = np.argmin(np.fabs(self.time-self.stepend))
        self.Cm = 1.0
        self.Rn = 50.
        self.taum = 10.
        self.I0 = -0.1
        self.L = None
        self.Vfit = self.V[self.itstart:self.itend]
        self.Tfit = self.time[self.itstart:self.itend]-self.stepstart
    
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
        params = dexpmodel.make_params(amp=amp, min=-100, max=0.)
        params.add('C0', value=0.5, min =-100, max=0.)
        params.add('cratio', value=3.0, min=0, max=100.)
        params.add('C1', expr='C0/cratio')
        params.add('tau0', value=2.0, min=0.2, max=10.0)
        params.add('tratio', value=6.0, min=3.0, max=100.0, vary=True)
        params.add('tau1', expr='tau0/tratio')
        # , C0=5, C1=1, tau0=1, tau1=0.1)
        print ('Initial Params: ', [p for p in params.values()])
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
    
    def alphas_taud(self):
        """
        eqn 10, White et al. 1992
        """
        a = np.zeros(2)
        for i in range(0, 2):
            a[i] = np.sqrt((self.tau_d/
                        self.fitresult.params['tau%d'%i].value)-1.0)
        return a
    
    def alphas_Ctau(self):
        """
        Eq 1/2 Rose and Dagum 1988
        """
        a = np.zeros(2)
        for i in range(0, 2):
            a[i] = (self.fitresult.params['C%d'%i].value/
                        self.fitresult.params['tau%d'%i].value)
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
        eq. 12, White et al. 1994
        """
        return self.eta()*self.tau_d/(self.Cm*self.Rn*(self.rho+1.0))
    
    def betas(self, alpha=None):
        if alpha is None:
            alpha = self.alphas()
        a = alpha
        betas = np.zeros(2)
        for i in range(0,2):
#            print(a[i])
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
    
    def spaceconstant(self):
        """
        Basic space constant from Rall, 1969, for cylinder...
        """
        tau0 = self.fitresult.params['tau%d'%0].value
        tau1 = self.fitresult.params['tau%d'%1].value
        return np.pi/np.sqrt((tau0/tau1)-1)


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
        
# dataclass for data in table 1 - each entry
@dataclass
class CellPars:
    Cell: str=''
    Rn: float=0.  # input resistance
    Vm: float=0.  # resting membrane potential
    Ih: float=0.  # holding current
    C0: float=0.  # amplitude of slow tau
    tau0: float = 0.  # time constant of slow tau
    C0: float=0.  # amplitude of slow tau
    tau0: float = 0.  # time constant of slow tau
    CR: float = 0.  # coefficients ration.

data = """
Cell Rn Vm Ih C0 tau0 C1 tau1 CR
6 137.0 -58.0 0.20 10.2 3.9 2.1 1.2 0.33
11  19.3 -60.0 0.40 1.5 3.9 0.4 1.1 0.47
13  58.4 -56.0 0.3  4.8 9.2 0.7 1.4 0.48
"""
    
    
def WYM94data(data):
    spc = re.compile("[ ;,\t\f\v]+") # format replacing all spaces with tabs
    dataiter = re.finditer(spc, data)
    data = re.sub(spc, ',', data)

    sio = io.StringIO(data)
    df = pd.read_table(sio, sep=',')
    df = df.set_index('Cell')
    # print(df.head())
    dc = {}
    for cell in df.index:
        # print('cell: ', cell)
        # print(df.loc[cell])
        da = CellPars()
        da.Cell = str(cell)
        # print(df.Cell)
        da.Rn = df.loc[cell]['Rn']
        da.Vm = df.loc[cell]['Vm']
        da.Ih = df.loc[cell]['Ih']
        da.C0 = df.loc[cell]['C0']
        da.C1 = df.loc[cell]['C1']
        da.tau0 = df.loc[cell]['tau0']
        da.tau1 = df.loc[cell]['tau1']
        da.CR = df.loc[cell]['CR']
        dc[cell] = da
    # print(dc)
    return dc

def generate_testdata(noise=0.):
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
    
    ve = E.nexp(I0, Rn, C=C, noise=noise)
    tf = T[E.itstart:]-tstart
    yf = ve[E.itstart:]
    return tf, yf, E

def test():
    """
    Generate a test case for measuring electrotonic parameters
    from a current trace.
    """
    
    tf, yf, E = generate_testdata(noise=5.)
    # rate = 0.01
    # dur = 50.
    # npts = int(dur/rate)
    # tstart = 5.
    # I0 = -0.1  # nA
    # Rn = 50.  # Mohm
    # I = np.zeros(npts)
    # T = np.arange(npts)*rate
    # V = np.zeros(npts) - 60.
    # noise = 1.0 # in mV
    # C_Ratios = np.array([5., 1.])
    # C = -(I0*Rn)*C_Ratios/np.sum(C_Ratios)  # scale ratios to voltage
    # print('C: {0:.3f}, {1:.3f}'.format(C[0], C[1]))
    # E = Electrotonic(I, V, T, tstart)
    #
    # ve = E.nexp(I0, Rn, C=C, noise=noise)
    # tf = T[E.itstart:]-tstart
    # yf = ve[E.itstart:]
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
    iax = mpl.imshow(esurf, interpolation='bilinear', origin='lower',
                cmap=cm.gray, extent=(np.min(lr), np.max(lr),
                 np.min(tr), np.max(tr)), aspect='auto')
    ax = iax.axes
    ax.set_title("Error Surface")
    ax.set_xlabel('tau (ms)')
    ax.set_ylabel('Lambda (microns)')
    cs = mpl.contour(x, y, esurf)
    mpl.clabel(cs, inline=1, fontsize=8)
    mpl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Measure electrotonic parameters",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )

    parser.add_argument(
        dest="cell", action="store", default=None, help="Select the cell (no default)"
    )
    args = parser.parse_args()
    
    if args.cell == 'test':
        test()
    
    else:
        WYM94data(data)
    
