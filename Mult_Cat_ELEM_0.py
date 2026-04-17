# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 19:39:37 2025

@author: aareshfb
"""

import pandas as pd
import numpy as np
from funcs_cat_Multiple import *
import time
import multiprocessing as mp

def elem_0_catM(data_list,W,nu0,mu,gamma):
    '''
    

    Parameters
    ----------
    data_list : list
        A list of data for different categories. The lenght of the list is the number of categories.
        Each item of data_list (refered to as data_c) must be a dictinary that contains the gene expression data. 
        Additionally, the data_c must be provided as follows
            data_c : dict
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
    Theta: np.ndarray of shape (p, p, (C+1)*K)
        A stacked set of estimated precision matrices:
        - Theta[:, :, 0:K]   → global (shared) network
        - Theta[:, :, K:2K]  → category 1 networks
        - Theta[:, :, 2K:3K] → category 2 networks
        - Theta[:, :, CK:(C+1)K] → category C networks
    '''
    
    S_list = [backwardMap(data_c, nu0) for data_c in data_list]
    print('Estimating Theta (Categorial ELEM-0)....')
    Theta = estimateNetwork_CategoryM(S_list, W, mu, gamma);
    print('Estimation complete ....')
    return(Theta)

