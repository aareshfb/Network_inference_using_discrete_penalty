# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 19:39:37 2025

@author: aareshfb
"""

import pandas as pd
import numpy as np
from funcs_cat import *
import time
import multiprocessing as mp

def elem_0_cat(data1,data2,W,nu0,mu,gamma):
    '''
    

    Parameters
    ----------
    data : dict
        A dictionary where each key corresponds to a population index (0, 1, ..., K-1),
        and each value is an n_k x p np.array representing the gene expression data for that category.
    W : np.ndarray size K by K
        A K x K np.array (symmetric) representing the hypergraph structure
    nu0 : float
        Soft-thresholding parameter used in the approximate backward mapping.
    mu : float
        Penalty weight for the sparsity-inducing \(\ell_0\) norm. Larger values encourage sparser solutions.
    gamma : float
        Penalty weight for enforcing similarity between precision matrices of neighboring categories
        as defined by the hypergraph W.

    Returns
    -------
    Theta: np.ndarray of shape (p, p, 3K)
        A stacked set of estimated precision matrices:
        - Theta[:, :, 0:K]   → global (shared) network
        - Theta[:, :, K:2K]  → category 1 networks
        - Theta[:, :, 2K:3K] → category 2 networks
    '''

    S1 = backwardMap(data1, nu0)
    S2 = backwardMap(data2, nu0)
    print('Estimating Theta (Categorial ELEM-0)....')
    Theta = estimateNetwork_Category(S1,S2, W, mu, gamma);
    print('Estimation complete ....')
    return(Theta)

