#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Functions for visualisation.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def lag_sort(lags, community=None, latency='overall', plot=True):
    """Sort the lag matrix according to community affiliation and latency.
    
    Parameters
    ----------
    lags
        The pairwise matrix of lags. If a community affiliation vector
        if provided, it should partition the first axis of this matrix.
    community
        A numeric community affiliation vector, indicating for each node
        its community membership. If no community affiliation vector is
        provided, only overall latency will be used as a sorting criterion.
    latency
        The latency criterion to be used for sorting lags.
        'community': mean latency within community
        'overall': mean latency within the entire network
    plot
        Indicates that the result should be plotted instead of returned.
    """
    if community is None:
        community = np.ones(lags.shape[0])
    sort_idx = np.empty(len(community), dtype='uint32')
    communities = np.unique(community)
    slc=0
    for c in communities:
        idx = np.where(community==c)[0]
        lat = np.zeros_like(idx, dtype='float')
        for i, node_id in enumerate(idx):
            if latency == 'community':
                node = lags[node_id, idx]
            else:
                node = lags[node_id, :]
            lat[i] = np.nanmean(node)
        sort_idx[slc:slc+len(idx)] = idx[lat.argsort()]
        slc = slc + len(idx)
    lags_sorted = np.take(np.take(lags, sort_idx, axis=0), sort_idx, axis=1)
    if plot:
        plt.subplot(1,2,1)
        plt.imshow(lags_sorted, cmap='jet', vmin=-1.5, vmax=1.5)
        plt.subplot(1,2,2)
        sns.kdeplot(lags.flatten(), shade=True);
    else:
        return lags_sorted
