# -*- coding: utf-8 -*-
"""
Created on Thus Aug 15 09:59:47 2025

@author: aares
"""
from funcs import *
from ELEM_0 import elem_0
import sys 
import pandas as pd
import numpy as np
import time
import os
import re

def read_prec(no_regs,p,t):
    base_path = "../Ground_truth"
    prec = {}

    expected_regs = [f"reg{i}" for i in range(1, no_regs+1)]
    expected_genes = [f"gene{i}" for i in range(1, p-no_regs+1)]

    def natural_key(s):
        m = re.match(r"([A-Za-z]+)(\d+)$", str(s))
        if m:
            return (m.group(1), int(m.group(2)))
        return (str(s), 0)

    ordered_regs = sorted(expected_regs, key=natural_key)
    ordered_cols = ordered_regs + sorted(expected_genes, key=natural_key)

    for i in range(1, t+1):
        path = os.path.join(base_path, f"network-{i}.csv")
        df = pd.read_csv(path, index_col=0)

        row_names = df.index.tolist()
        col_names = df.columns.tolist()

        if row_names != ordered_regs:
            raise ValueError(
                f"{path}: row order mismatch.\n"
                f"Expected rows: {ordered_regs}\n"
                f"Found rows:    {row_names}"
            )

        if col_names != ordered_cols:
            raise ValueError(
                f"{path}: column order mismatch.\n"
                f"Expected cols: {ordered_cols}\n"
                f"Found cols:    {col_names}"
            )

        arr = df.to_numpy(dtype=float)

        if arr.shape != (no_regs, p):
            raise ValueError(
                f"{path}: expected shape (no_regs, p), found {arr.shape}"
            )

        prec[i - 1] = arr
        print(f"Loaded {path}: shape={arr.shape}")

    return prec

def errorMetrics2(Theta1, True_Precission1):


    TP = (Theta1 != 0) * (True_Precission1 != 0)
    FP = (Theta1 != 0) * (True_Precission1 == 0)
    FN = (Theta1 == 0) * (True_Precission1 != 0)

    TP = np.sum(TP)
    FP = np.sum(FP)
    FN = np.sum(FN)

    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    F1 = 2 * Precision * Recall / (Precision + Recall)

    print('Precision:', Precision, '\nRecall:', Recall, "\nF1:", F1, "\n")

    return Precision, Recall, F1



def load_dataset(no_regs,p,t):
    base_path = "../Expression_data"
    rng = np.random.default_rng(0)

    expected_regs = [f"reg{i}" for i in range(1, no_regs+1)]
    expected_genes = [f"gene{i}" for i in range(1, p-no_regs+1)]

    def natural_key(s):
        m = re.match(r"([A-Za-z]+)(\d+)$", str(s))
        if m:
            return (m.group(1), int(m.group(2)))
        return (str(s), 0)

    ordered_regs = sorted(expected_regs, key=natural_key)
    ordered_genes = sorted(expected_genes, key=natural_key)
    ordered_rows = ordered_regs + ordered_genes

    data = {}

    for cluster_idx in range(1, t+1):
        path = os.path.join(
            base_path,
            f"cluster_{cluster_idx}_sparse_expression.csv"
        )

        # First column contains row names; first row contains cell IDs
        df = pd.read_csv(path, index_col=0)

        # Ensure numeric matrix
        df = df.apply(pd.to_numeric, errors="coerce")

        # Add any missing rows with tiny Gaussian noise
        missing_rows = [r for r in ordered_rows if r not in df.index]
        for row_name in missing_rows:
            df.loc[row_name] = rng.normal(
                loc=0.0,
                scale=1e-5,
                size=df.shape[1]
            )

        # Keep only the expected rows, in desired order
        df = df.loc[ordered_rows]

        # Remove row/column labels before storing
        arr = df.to_numpy(dtype=float)

        # Store clusters as 0,...,7 to match your earlier convention
        data[cluster_idx - 1] = arr.T

        print(
            f"Loaded {path}: added missing={len(missing_rows)}, "
            f"final shape={arr.shape}"
        )

    return data

if __name__ == "__main__":

    
    t = 10
    p = 200
    n= 2000
    no_regs = 30
    nu0_list = [float(sys.argv[1])]
    mask_diag = sys.argv[2] == "True"
    # n=200
    # nu0_list=[1e-2]
    mu_list=np.logspace(3, -7, 250)
    gamma=1

    folder = "masked_diag" if mask_diag else "unmasked_diag"
    path = os.path.join("Output", "AUPRC3", folder)
    
    os.makedirs(path, exist_ok=True)

    error_metrics_filename = os.path.join(path, f"Code_L0_n{n}_g{gamma}_nu{nu0_list[0]}.csv")
    if not os.path.isfile(error_metrics_filename):
        with open(error_metrics_filename, 'w') as f:
            f.write('t,p,n,nu_0,mu,gamma,no_edges,bic value,precision,,,,,,,,precision_avg,recall,,,,,,,,recall_avg,F1,,,,,,,,F1_avg,Time(sec),mask_diag')
            f.write('\n')
    
    W = np.eye(10)
    edges = [(0, 5), (0, 8), (5, 7), (1, 5), (1, 6), (6, 9), (4, 8), (3, 4), (2, 8)]
    for i, j in edges:
        W[i, j] = 1
    W=np.triu(W,1)
    print("reading expression data...")
    data = load_dataset(no_regs,p,t)
    
    print('reading true precision matrices')
    True_Precission = read_prec(no_regs,p,t)
    True_Precission=np.stack([True_Precission[i] for i in range(t)],axis=2)
    
    
    for nu0 in nu0_list:
        for mu in mu_list:
        
            start = time.time()
            Theta= elem_0(data,W,nu0,mu,gamma)
            stop = time.time()
            time_network=stop-start
            
            bic=Calculate_Bic(Theta, data, t*[n])
            Theta_eval = Theta.copy()
            True_eval = True_Precission.copy()

            if mask_diag:
                for i in range(Theta_eval.shape[2]):
                    np.fill_diagonal(Theta_eval[:,:,i], 0)
                    np.fill_diagonal(True_eval[:,:,i], 0)
                    
            Precision,Recall,F1=[],[],[]
            for i in range(t):
                # idx=regulator_idx()
                # all regulators are same, taking the index from hsc
                pre,rec,fscore=errorMetrics2(Theta_eval[:no_regs,:,i], True_eval[:,:,i])
                Precision.append(pre); Recall.append(rec); F1.append(fscore)
                
            pre,rec,fscore=errorMetrics2(Theta_eval[:no_regs,:,:], True_eval[:,:,:])
            Precision.append(pre); Recall.append(rec); F1.append(fscore)
    
            #print(f"Precision: {Precision:.4f}")
            #print(f"Recall:    {Recall:.4f}")
            #print(f"F1 Score:  {F1:.4f}")
            
            no_edges=0.5*(np.count_nonzero(Theta)-3*p)
            with open(error_metrics_filename, 'a') as f:
                    #f.write(str(t)+','+str(p)+','+str(n)+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(no_edges)+','+str(bic)+','+str(Precision)+','+str(Recall)+','+str(F1)+','+str(time_network))
                    #f.write('\n')
                    row = [t,p,n,nu0,mu,gamma,no_edges,bic] + Precision + Recall + F1 + [time_network,mask_diag]
                    f.write(','.join(map(str,row))+'\n')
                    