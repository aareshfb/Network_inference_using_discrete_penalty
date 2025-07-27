# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 19:39:37 2025

@author: aareshfb
"""

import pandas as pd
import numpy as np
from funcs import *
import time
import multiprocessing as mp

def elem_0(data,W,nu0,mu,gamma):
    '''
    

    Parameters
    ----------
    data : dict
        A dictionary where each key corresponds to a population index (0, 1, ..., K-1),
        and each value is an n_k x p np.array representing the gene expression data for that category.
    W : np.ndarray size (K, K)
        A (K, K) np.array (symmetric matrix) representing the hypergraph structure
    nu0 : float
        Soft-thresholding parameter used in the approximate backward mapping.
    mu : float
        Penalty weight for the sparsity-inducing \(\ell_0\) norm. Larger values encourage sparser solutions.
    gamma : float
        Penalty weight for enforcing similarity between precision matrices of neighboring categories
        as defined by the hypergraph W.

    Returns
    -------
    Theta: np.ndarray
        A p x p x K NumPy array containing the estimated precision matrices for all K categories.
        Each slice Theta[:, :, k] corresponds to the precision matrix for category k.

    '''

    S = backwardMap(data, nu0)
    print('Estimating Theta ....')
    Theta,_ = estimateNetwork(S, W, mu, gamma);
    print('Estimation complete ....')
    return(Theta)



    
        
