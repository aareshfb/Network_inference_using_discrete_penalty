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

def regulator_idx():
   
    cell_types = ["hsc", "cmp", "gmp"]
    all_indices = {}
    
    for ct in cell_types:
        infile = f"../Simulated_GRNs/{ct}_regulators_idx.txt"
        with open(infile, "r") as f:
            all_indices[ct] = [int(line.strip()) for line in f]
    
    return(all_indices)


def read_prec():
    path="../Simulated_GRNs/hsc_precision_matrix.csv"
    prec={}
    ncols = pd.read_csv(path, nrows=1).shape[1]

    # now load as NumPy
    prec[0] = np.loadtxt("../Simulated_GRNs/hsc_precision_matrix.csv", delimiter=",", skiprows=1, usecols=range(1, ncols))
    prec[1] = np.loadtxt("../Simulated_GRNs/cmp_precision_matrix.csv", delimiter=",", skiprows=1, usecols=range(1, ncols))
    prec[2] = np.loadtxt("../Simulated_GRNs/gmp_precision_matrix.csv", delimiter=",", skiprows=1, usecols=range(1, ncols))
    return(prec)

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


def load_dataset(n, base_path="../"):
    """
    Load simulated datasets depending on n (200, 1000, 2000).
    Returns a dict {0: hsc_df, 1: cmp_df, 2: gmp_df}.
    Each value is a pandas DataFrame (genes x cells).
    """
    if n == 2000:
        files = {
            0: base_path + "hsc.table",
            1: base_path + "cmp.table",
            2: base_path + "gmp.table",
        }
    elif n == 1000:
        files = {
            0: base_path + "hsc_n1000.table",
            1: base_path + "cmp_n1000.table",
            2: base_path + "gmp_n1000.table",
        }
    elif n == 200:
        files = {
            0: base_path + "hsc_n200.table",
            1: base_path + "cmp_n200.table",
            2: base_path + "gmp_n200.table",
        }
    else:
        raise ValueError("n must be one of {200, 1000, 2000}")

    data = {}
    for idx, path in files.items():
        df = pd.read_csv(path, sep="\t", header=None)
        arr = df.iloc[1:, 1:].to_numpy(dtype=float)
        data[idx] = arr.T
        print(f"Loaded {path}: {np.shape(arr)}")

    return data

if __name__ == "__main__":

    
    t = 3
    p = 65
    n= int(sys.argv[1])
    nu0_list = [float(sys.argv[2])]
    mask_diag = sys.argv[3] == "True"
    # n=200
    # nu0_list=[1e-2]
    mu_list=np.logspace(4, 6, 20)
    gamma=1

    folder = "masked_diag" if mask_diag else "unmasked_diag"
    path = os.path.join("Output", "AUPRC", folder)
    
    os.makedirs(path, exist_ok=True)

    error_metrics_filename = os.path.join(path, f"Error_metrics_n{n}.csv")
    if not os.path.isfile(error_metrics_filename):
        with open(error_metrics_filename, 'w') as f:
            f.write('t,p,n,nu_0,mu,gamma,no_edges,bic value,precision,,,precision_avg,recall,,,recall_avg,F1,,,F1_avg,Time(sec),mask_diag')
            f.write('\n')
    
    W = np.eye(3) + np.diag(np.ones(2), k=1) 
    # W=np.triu(W)
    print('reading expression data...')
    data = load_dataset(n)
    
    print('reading true precision matrices')
    True_Precission = read_prec()
    True_Precission=np.stack([x[:,:] for x in True_Precission.values()],2)
    
    
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
                    mask = 1 - np.eye(p)
                    np.fill_diagonal(Theta_eval[:,:,i], 0)
                    np.fill_diagonal(True_eval[:,:,i], 0)
                    
            Precision,Recall,F1=[],[],[]
            for i in range(t):
                idx=regulator_idx()
                # all regulators are same, taking the index from hsc
                pre,rec,fscore=errorMetrics2(Theta_eval[idx['hsc'],:,i], True_eval[idx['hsc'],:,i])
                Precision.append(pre); Recall.append(rec); F1.append(fscore)
                
            pre,rec,fscore=errorMetrics2(Theta_eval[idx['hsc'],:,:], True_eval[idx['hsc'],:,:])
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
                    