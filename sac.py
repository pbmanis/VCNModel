# function [y, yh, hx, mr, X] = sac(varargin)
#    [yan, sachist, sacx, mr] = sac(allsp, width, binw, delay, dur);
#
# Shuffled autocorrelation function
# Based on Louage et al., 2004
# X is an array of N responses x length(t) points in time. The times are in
# msec. Since the number of spikes in X may vary from trial to trial,
# X is assumed to be a cell array of spike times.
# The autocorrelation is calculated for every event in X that is not within
# is not within a single spike train.
# twin is the time window for the autocorrelation to be calculated (in
# msec)
# binw is the bin width for the histogram, in msec.
# delay is the delay to the start of measurement in X; in msec.
# dur is the duration of the measurement window in X; in msec.
# (these last two parameters allow you to pass a response and only analyze
# the portion that is relevant).
#
# y is the output autocorrelation function. It consists of a list of
# the correlation intervals, not corrected at all.
#
# yh is the histogram bin amplitudes, with bin width hx.
# This histogram is normalized.
# mr is the mean rate
#
# 10/5/04 Paul B. Manis, Ph.D.
# 10/5/04 some fixes for an spike trains.
# Note: calling the routine with letters as the first argument will generate the plots
# from Figure 2 of Louage et al.

# adapted from Matlab version.

import numpy as np
import pyqtgraph as pg
import pylibrary.pyqtgraphPlotHelpers as pgh

# global ALLCH DFILE
#
#
# y = []; yh = []; hx = []; mr = 0; X={};
#
# switch nargin
#     case 0
#         return;
#     case 5
#         X = varargin{1};
#         twin = varargin{2};
#         binw = varargin{3};
#         delay = varargin{4};
#         dur = varargin{5};
class SAC(object):

    def __init__(self):
        pass


    def makeTests(self, test):
        baseper = 4./3.
        stimdur = 1000.
        ntestrep = 20
        twin = 5.
        binw = 0.05
        delay = 0.
        dur = stimdur
        ddur = 10.        
        X=[None]*ntestrep

        if test == 'A':
            for i in range(ntestrep):
                X[i] = np.arange(baseper, stimdur, baseper)  # all identical

        elif test == 'B':
            for i in range(ntestrep):
                X[i] = np.arange(baseper, stimdur, baseper)  # all identical
                X[i] = X[i] + 0.08*np.random.random(len(X[i]))

        elif test == 'C':
            for i in range(ntestrep):
                X[i] = np.arange(baseper, stimdur, baseper)  # all identical
                X[i] = X[i] + 0.08*np.random.random(len(X[i]))
                n = len(X[i])
                sig = stimdur*np.sort(np.random.rand(n/4))
                X[i] = np.sort(np.append(X[i], sig))
                
        elif test == 'D':
            for i in range(ntestrep):
                X[i] = np.arange(baseper, stimdur, baseper)  # all identical
                X[i] = X[i] + 0.170*np.random.random(len(X[i]))

        elif test == 'E':
            for i in range(ntestrep):
                bx = np.random.permutation(np.arange(baseper, stimdur, baseper))
                bx = bx[:int(len(bx)/2)]
                X[i] = sorted(bx); # elimnate half by random selection
                X[i] = X[i] + 0.170*np.random.random(len(X[i]))
                
        elif test == 'F':
            binw = 0.15
            x=[None]*ntestrep
            for i in range(ntestrep):
                X[i] = stimdur*np.sort(np.random.random(120)) # corresponds to about 120 s/s
            ddur = 100.

        elif test == 'G':
            binw = 0.15
            sig = stimdur*np.sort(np.random.rand(120))
            for i in range(ntestrep):
                X[i] = sig # same signal every time
            ddur = 100

        pars = {'twin': twin, 'binw': binw, 'ntestrep': ntestrep, 
                'baseper': baseper, 'stimdur': stimdur, 'delay': delay,
                'dur': dur, 'ddur': ddur}
        return(X, pars)   


    def compute(self, X, pars):
        """
        """# now the SAC calculations.
        twin = pars['twin']
        binw = pars['binw']
        delay = pars['delay']
        dur = pars['dur']
        lx = len(X)
        # make sure input bin width and window are integer mulitples.
        if(((twin/binw) - np.floor(twin/binw)) != 0):
            x = np.floor(twin/binw);
            twin = binw * x; # recalculate the window

        maxlx = 100;
        maxn = 100000;
        if (lx > maxlx): # for testing purposes.
            lx = maxlx

        yc = [[]]*lx # Prellocation helps ALOT with speed
        n = 1
        spcount = np.zeros(lx)
        for i in range(lx):
            xi = X[i] # get spike train i
            kxi = np.where((xi >= delay) & (xi < delay+dur)) # window the data to be analyzed
            xi = xi[kxi] # reduce data set
            spcount[i] = len(xi)
        #    fprintf(1, 'i=#d  ns: #d n', i, spcount(i)); # let user know we're working

            yj = np.nan*np.zeros(maxn) # preallocate space for cross correlation
            n = 0
            for j in range(lx):
                if j != i:
                    xj = X[j]; # get spike train j
                    kxj = np.where((xj >= delay) & (xj < delay+dur))
                    xj = xj[kxj[0]] # reduce data set
                    for ti in range(len(xi)): # cross correlate against all spikes in xi
                        iti = xi[ti] # time of ti'th spike in spike train i
                        xjk = np.where((xj >= (iti-twin-binw/2.)) & (xj < (iti+twin+binw/2.))) # find spikes in j, in the window relative to current spike
                        # binw/2 are egde effect corrections for the histogram stuff below.
                        nx = len(xjk[0])
                        if nx > 0 and (n+nx) < maxn:
                            n = n + 1
                            yj[n:n+nx] =  xj[xjk] - iti # compute the difference in time between spike in j and relative to current spike in i
                            n = n + nx-1
            yc[i] = yj
        #    fprintf(1, ' elapsed: #8.3fs  last n = #d\n', toc, n);
        ya = np.array(yc).flatten()
        y = ya[~np.isnan(ya)] # clean it up - and shorten.
#        print 'y: ', y
        if len(y) == 0: # no spikes in this window
            return;

        # now calculate the normalized histogram.
        # normalization is N*(N-1)*r^2*deltat*D
        # where N is the number of presentations (lx), r is the rate (sp/sec),
        # deltat is the bin width for the histogram (sec), and D is the duration of
        # the trace (?).
        # normalization goes to 1 as the time gets larger...
        #
        # rate = sum(spcount)/(lx*dur/1000.); # get average rate
        # print'Mean firing rate: %8.1f' % rate
        # mr = rate
        # nfac = lx*(lx-1)*rate*rate*(binw/1000.)*(dur/1000.); # correction factor
        # yh, bins = np.histogram(y, binw); # generate histogram
        # yh = yh/nfac; # to convert to rate, spikes/second
        #       #m = mean(yh);
        #fprintf(1, 'Mean value of yh: #8.3f, scale factor: #f\n', m, nfac);
        return y

if __name__ == '__main__':

    sac = SAC()
    tests = ['A' ,'B' ,'C', 'D', 'E', 'F', 'G']
    win = pgh.figure(title='AN Inputs')
    layout = pgh.LayoutMaker(cols=2,rows=len(tests), win=win, labelEdges=True, ticks='talbot')
    win.resize(600, 125*len(tests))
    for i, t in enumerate(tests):
        X, p = sac.makeTests(t)
        y = sac.compute(X, p)
        yh, bins = np.histogram(y, bins=np.linspace(-p['ddur'], p['ddur'], int(p['ddur']/0.05)), density=False)
        sach = pg.PlotCurveItem(bins, yh, stepMode=True, fillLevel=0,
            brush=(255, 0, 255, 255), pen=None)
        layout.getPlot((i, 1)).addItem(sach)
        layout.getPlot((i, 1)).setXRange(-5., 5.)
        size = 4
        Y=[[]]*len(X)
        for j in range(len(X)):
            Y[j] = (j+1)*np.ones(len(X[j]))
        raster = pg.ScatterPlotItem(x=np.array(X).flatten(), y=np.array(Y).flatten(), pen='w', brush='b', size=size, pxMode=True)
        layout.getPlot((i, 0)).addItem(raster)
        if t not in ['F', 'G']:
            layout.getPlot((i, 0)).setXRange( 0., 10.)
        else:
            layout.getPlot((i, 0)).setXRange( 0., 100.)
            
    pgh.show()
    