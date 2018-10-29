#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Functions for computing lag matrices following Mitra et al.
"""

import numpy as np

def lag_analysis():
    pass

def corrcoef_lagged(timeseries1, timeseries2, lagmax=5, tmask=None):
    """Compute the lagged correlation between two sets of time series.
   
    
    Parameters
    ----------
    timeseries1
        First set of time series.
    timeseries2
        Second set of time series.
    lagmax
        Maximum lag to consider, in frames.
    tmask
        Temporal mask indicating frames of the input time series that should
        be considered in the lag computation.

    Outputs
    -------
    corr: array
        Lagged correlation matrix with dimensions
        len(timeseries1) x len(timeseries2) x (2*lagmax+1)
    """
    corr = np.zeros([timeseries1.shape[1], timeseries2.shape[1], 2*lagmax+1])

    for k, t in enumerate(np.arange(-lagmax, lagmax+1, 1)):
        tau = abs(t)
        if t > 0:
            ts1_lagged = timeseries1[:-tau, :]
            ts2_lagged = timeseries2[tau:, :]
            observed = (tmask[:-tau]*tmask[tau:])==1
        elif t < 0:
            ts1_lagged = timeseries1[tau:, :]
            ts2_lagged = timeseries2[:-tau, :]
            observed = (tmask[:-tau]*tmask[tau:])==1
        else:
            ts1_lagged = timeseries1
            ts2_lagged = timeseries2
            observed = tmask==1
        ts1_lagged = ts1_lagged[observed,:]
        ts2_lagged = ts2_lagged[observed,:]
        corr[:,:,k] = np.corrcoef(x=ts1_lagged,
                                  y=ts2_lagged,
                                  rowvar=False)[ts1_lagged.shape[1]:,
                                                :ts2_lagged.shape[1]]

    return corr

"""
function  r = lagged_corr2(Avg1,Avg2,L,F)
    % a signal time x ROI
    % L how may lags
    % F logic array of frames passed QC
    L1 = size(Avg1,2);
    L2 = size(Avg2,2);
    r = single(zeros(L1,L2,2*L+1));

    k = 1;
    for i = -L:1:L
        tau = abs(i);
        if i >=0
            Avg1_lagged = Avg1(1:(end-tau),:);
            Avg2_lagged = Avg2((tau+1):end,:);
        else
            Avg1_lagged = Avg1((tau+1):end,:);
            Avg2_lagged = Avg2(1:(end-tau),:); 
        end
        f = (F(1:(end-tau)).*F(tau+1:end))==1;
        Avg1_lagged = Avg1_lagged(f,:);
        Avg2_lagged = Avg2_lagged(f,:);
        r(:,:,k) = correl2(Avg1_lagged, Avg2_lagged);
        k = k+1;
    end

end

function  r = correl2(a,b)
    meana = mean(a,1);
    meanb = mean(b,1);
    a = bsxfun(@minus,a,meana);
    b = bsxfun(@minus,b,meanb);
    maga = sqrt(sum(a.^2,1))';
    magb = sqrt(sum(b.^2,1));
    r = a'*b./(maga * magb);
end
"""

def parabolic_interpolation(timeseries, sample_time, criterion='midpoint'):
    """Use parabolic interpolation to estimate the time at which an input
    time series has a maximum or minimum.

    Parameters
    ----------
    timeseries
        The discretely sampled time series for which an extremum should
        be estimated via parabolic interpolation.
    sample_time
        The duration between samples in the 
    criterion
        The criterion to use in determining whether to compute a maximum
        or minimum. Possible values include:
        * maximum: always compute the maximum.
        * minimum: always compute the minimum.
        * extremum: return whichever has a greater absolute value.
        * midpoint: return the minimum if the time series midpoint is
          negative and the maximum if it is positive.

    Outputs
    -------
    lag
        The temporal lag (in the same units as sample_time) at which the
        input time series has an extremum.
    peak
        The peak value estimated using parabolic interpolation.
    """

    return lag, peak

"""
Siphoned from https://github.com/astronomerdamo/pydcf/blob/master/dcf.py
Come back to this later if there's a sample-related confound.
def sdcf(ts1, ts2, t, dt):

    '''
        Subroutine - sdcf
          DCF algorithm with slot weighting
    '''

    dcf = np.zeros(t.shape[0])
    dcferr = np.zeros(t.shape[0])
    n = np.zeros(t.shape[0])

    dst = np.empty((ts1.shape[0], ts2.shape[0]))
    for i in range(ts1.shape[0]):
        for j in range(ts2.shape[0]):
            dst[i,j] = ts2[j,0] - ts1[i,0]

    for k in range(t.shape[0]):
        tlo = t[k] - dt/2.0
        thi = t[k] + dt/2.0
        ts1idx, ts2idx = np.where((dst < thi) & (dst > tlo))

        mts2 = np.mean(ts2[ts2idx,1])
        mts1 = np.mean(ts1[ts1idx,1])
        n[k] = ts1idx.shape[0]

        dcfdnm = np.sqrt((np.var(ts1[ts1idx,1]) - np.mean(ts1[ts1idx,2])**2) \
                         * (np.var(ts2[ts2idx,1]) - np.mean(ts2[ts2idx,2])**2))

        dcfs = (ts2[ts2idx,1] - mts2) * (ts1[ts1idx,1] - mts1) / dcfdnm

        dcf[k] = np.sum(dcfs) / float(n[k])
        dcferr[k] = np.sqrt(np.sum((dcfs - dcf[k])**2)) / float(n[k] - 1)

    return dcf, dcferr

def gdcf(ts1, ts2, t, dt):

    '''
        Subroutine - gdcf
          DCF algorithm with gaussian weighting
    '''

    h = dt/4.0
    gkrn = lambda x: np.exp(-1.0 * np.abs(x)**2 / (2.0 * h**2)) \
           / np.sqrt(2.0 * np.pi * h)
    cntrbt = gkrn(3.290527*h)

    dcf = np.zeros(t.shape[0])
    dcferr = np.zeros(t.shape[0])
    n = np.zeros(t.shape[0])

    dst = np.empty((ts1.shape[0], ts2.shape[0]))
    for i in range(ts1.shape[0]):
        for j in range(ts2.shape[0]):
            dst[i,j] = ts2[j,0] - ts1[i,0]

    for k in range(t.shape[0]):
        gdst = gkrn(dst - t[k])
        ts1idx, ts2idx = np.where(gdst >= cntrbt)

        mts2 = np.mean(ts2[ts2idx,1])
        mts1 = np.mean(ts1[ts1idx,1])
        n[k] = ts1idx.shape[0]

        dcfdnm = np.sqrt((np.var(ts1[ts1idx,1]) - np.mean(ts1[ts1idx,2])**2) \
                         * (np.var(ts2[ts2idx,1]) - np.mean(ts2[ts2idx,2])**2))

        dcfs = (ts2[ts2idx,1] - mts2) * (ts1[ts1idx,1] - mts1) / dcfdnm
        dcf[k] = np.sum(dcfs) / float(n[k])
        dcferr[k] = np.sqrt(np.sum((dcfs - dcf[k])**2)) / float(n[k] - 1)

    return dcf, dcferr
"""