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


def lag_sort(lags, community=None, latency='overall',
             plot=True, vrange=[-1.5, 1.5]):
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
    vrange
        Indicates the range of disciminable values in the image colourmap.
        Does nothing if `plot` is disabled.
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
        f, (p0, p1) = plt.subplots(1,2,
                                   figsize=(12, 6),
                                   gridspec_kw = {
                                        'width_ratios':[1, 1],
                                        'height_ratios':[1]
                                        })
        p0.imshow(lags_sorted, cmap='jet', vmin=vrange[0], vmax=vrange[1])
        sns.kdeplot(lags.flatten(), shade=True, ax=p1)
        p0.set_xticks([]);
        p0.set_yticks([]);
        p1.set_yticks([]);
        p1.set(xlabel='Time')
        plt.show()
    else:
        return lags_sorted


def simil_plot(simil, within, vrange=[0.01, 0.08]):

    """Prepare a plot of similarity between vectors.
    Intended for making figures similar to those in the Gratton MSC paper.

    Parameters
    ----------
    simil
        A correlation or similarity matrix.
    within
        Number of nodes in a major grouping.
    vrange
        Indicates the range of disciminable values in the image colourmap.
    """

    nodes = simil.shape[0]

    fig_width = 8
    fig_height = fig_width * nodes/(nodes + 2)

    seskey = np.array(list(range(1,within+1))*10, ndmin=2)
    subkey = np.array([i for i in range(1,11) for k in range(within)], ndmin=2)
    keys = np.hstack([subkey.T, (10/within)*seskey.T])
    f, (p0, p1) = plt.subplots(1,2, 
                               figsize=(fig_width, fig_height),
                               gridspec_kw = {'width_ratios':[2, nodes], 'height_ratios':[nodes]})
    p0.imshow(keys, cmap='jet')
    p1.imshow(simil,
              vmin=vrange[0],
              vmax=vrange[1])
    p0.axis('off')

    # Minor ticks
    p1.set_xticks(np.arange(-.5, nodes, 1), minor=True)
    p1.set_yticks(np.arange(-.5, nodes, 1), minor=True)
    p1.grid(color='w', linestyle=':', linewidth=0.5, which='minor')
    p1.set_xticks(np.arange(-.5, nodes, within))
    p1.set_yticks(np.arange(-.5, nodes, within))
    p1.grid(color='w', linestyle='-', linewidth=1)
    for tk in p1.yaxis.get_major_ticks():
        tk.tick1On = False
        tk.tick2On = False
        tk.label1On = False
        tk.label2On = False
    for tk in p1.xaxis.get_major_ticks():
        tk.tick1On = False
        tk.tick2On = False
        tk.label1On = False
        tk.label2On = False

    plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    plt.show()
