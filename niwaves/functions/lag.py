#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Functions for computing lag matrices following Mitra et al.
"""

import numpy as np


def lag_analysis(timeseries1, timeseries2=None, tmask=None,
                 lagmax=5, min_block=0, sample_time=1):
    """
    Perform a lag analysis on a single subject by identifying valid epochs
    and computing a correlation weighted by the number of frames in each
    epoch.

    Parameters
    ----------
    timeseries1
        First set of time series.
    timeseries2
        Second set of time series.
    tmask
        Temporal mask indicating frames of the input time series that should
        be considered in the lag computation.
    lagmax
        Maximum lag to consider, in the units of sample_time.
    min_block
        The minimum length of a block, in the same units as sample_time.
        If sample_time is not set explicitly, then this is measured in
        samples. If min_block is not set, then all blocks of valid data
        will be returned.
    sample_time
        The sampling interval of the data.

    Outputs
    -------
    lags
        The temporal lag (in the same units as sample_time) at which the
        input time series has an extremum.
    peaks
        The peak value estimated using parabolic interpolation.
    """
    if timeseries2 is None:
        timeseries2 = timeseries1
    if tmask is None:
        tmask = np.ones(timeseries1.shape[0])
    lags = np.arange(-lagmax, lagmax+1, 1)
    blocks = tmask_blocks(tmask=tmask,
                          min_block=min_block,
                          sample_time=sample_time)
    corr = np.zeros([timeseries1.shape[1], timeseries2.shape[1], 2*lagmax+1])
    for block in blocks:
        tmask_block = np.zeros_like(tmask)
        tmask_block[block] = 1
        corr_lag = corrcoef_lagged(timeseries1=timeseries1,
                                   timeseries2=timeseries2,
                                   tmask=tmask_block,
                                   lagmax=lagmax,
                                   sample_time=sample_time)
        for i, l in enumerate(lags):
            corr_lag[:,:,i] = corr_lag[:,:,i] * (len(block) - abs(l))
        corr = corr + corr_lag

    nframes = sum([len(b) for b in blocks])
    nblocks = len(blocks)
    for i, l in enumerate(lags):
        corr[:,:,i] = corr[:,:,i]/(nframes - abs(l)*nblocks)

    corr = corr.reshape(-1, corr.shape[-1])
    lags, peaks = parabolic_interpolation(timeseries=corr,
                                          sample_time=sample_time,
                                          criterion='midpoint')
    lags = lags.reshape(timeseries1.shape[1]. timeseries2.shape[1])
    peaks = peaks.reshape(timeseries1.shape[1]. timeseries2.shape[1])
    return lags, peaks


def corrcoef_lagged(timeseries1, timeseries2, tmask,
                    lagmax=5, sample_time=1):
    """
    Compute the lagged correlation between two sets of time series.
   
    
    Parameters
    ----------
    timeseries1
        First set of time series.
    timeseries2
        Second set of time series.
    tmask
        Temporal mask indicating frames of the input time series that should
        be considered in the lag computation.
    lagmax
        Maximum lag to consider, in the units of sample_time.
    sample_time
        The sampling interval of the data.

    Outputs
    -------
    corr: array
        Lagged correlation matrix with dimensions
        len(timeseries1) x len(timeseries2) x (2*lagmax+1)
    """
    lagmax = int(np.ceil(lagmax*sample_time))
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
    lag = np.full_like(timeseries[:,0], np.NaN)
    peak = np.full_like(timeseries[:,0], np.NaN)
    if criterion == 'maximum':
        optim = timeseries
    elif criterion == 'minimum':
        optim = -timeseries
    elif criterion == 'extremum':
        optim = abs(timeseries)
    elif criterion == 'midpoint':
        optim = np.expand_dims(np.sign(timeseries[:,timeseries.shape[1]//2])
                ,1)*timeseries
    else:
        raise ValueError('Invalid criterion specified %s' % criterion)
    maxidx = np.argmax(optim, 1)
    use = np.where(np.logical_and(maxidx!=0,
            maxidx!=timeseries.shape[1]-1))[0]
    ts = timeseries[use,:]
    maxidx = maxidx[use]
    shift = maxidx - timeseries.shape[1]//2
    slices = np.take_along_axis(ts,
             np.array([maxidx-1,maxidx,maxidx+1]).T, axis=1)

    # Polynomial coefficients
    a2 = (slices[:,2] + slices[:,0] - 2*slices[:,1])/2
    a1 = (slices[:,2] - slices[:,0])/2
    # axis of symmetry
    lag[use] = -a1 / (2 * a2)
    # vertex
    peak[use] = a2 * lag[use]**2 + a1 * lag[use] + slices[:,1]
    # recentre and convert units
    lag[use] = (lag[use] - shift) * sample_time

    return lag, peak


def tmask_blocks(tmask, min_block=0, sample_time=1):
    """
    Return a list containing the indices of each valid block of a temporal
    mask.

    Parameters
    ----------
    tmask
        The temporal mask, in which 1 denotes valid observations and 0
        denotes invalid or missing observations.
    min_block
        The minimum length of a block, in the same units as sample_time.
        If sample_time is not set explicitly, then this is measured in
        samples. If min_block is not set, then all blocks of valid data
        will be returned.
    sample_time
        The sampling interval of the temporal mask.

    Outputs
    -------
    blocks: list
        A list whose length equals the number of valid blocks in the
        temporal mask and whose entries are lists of the indices in each
        block.
    """
    def _check_block_length(block, blocks, min_block, sample_time):
        if len(block) * sample_time > min_block:
            blocks.append(block)

    blocks = []
    block = []
    for i, t in enumerate(tmask):
        if t:
            block.append(i)
        else:
            _check_block_length(block, blocks, min_block, sample_time)
            block = []
    _check_block_length(block, blocks, min_block, sample_time)
    return blocks


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