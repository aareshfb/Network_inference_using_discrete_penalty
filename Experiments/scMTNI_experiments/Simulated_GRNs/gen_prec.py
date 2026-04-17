# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 12:58:26 2025

@author: aares
"""

import numpy as np
import pandas as pd


for ct in ["hsc", "cmp", "gmp"]:
    # ---- Step 1: read gene list ----
    with open(ct+"_allGenes.txt") as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(genes)} genes")
    df = pd.read_csv(ct+"_network.tsv", sep="\t", header=None, comment="#")
    df.columns = ["src", "dst"]
    
    n = len(genes)
    idx = {g:i for i,g in enumerate(genes)}
    A = np.zeros((n, n), dtype=int)
    
    for s, t in zip(df["src"], df["dst"]):
        if s+'_'+ct in idx and t+'_'+ct in idx:
            i, j = idx[s+'_'+ct], idx[t+'_'+ct]
            A[i, j] = 1
            A[j, i] = 1
        else:
            print('Error')
    
    A=A+np.eye(n)
    adj_df = pd.DataFrame(A, index=genes, columns=genes)
    adj_df.to_csv(ct+"_precision_matrix.csv")