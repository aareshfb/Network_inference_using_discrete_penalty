import pandas as pd
import numpy as np
from funcs import *
import time
import multiprocessing as mp
from arboreto.algo import grnboost2
#conda activate conda1-env
import os
import sys
from sklearn.metrics import precision_recall_curve, auc

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

    file_path='../generated_data/'
    file_name = 'simulation_clust_'
    file_ext='.csv'
    n = 2000
    t = 3
    p = 65
    idx=regulator_idx()
    reg_idx=idx['hsc']
    print('generating expression data...')
    data = load_dataset(n)
    
    print('reading true precision matrices')
    True_Precission = read_prec()
    True_Precission=np.stack([x[:,:] for x in True_Precission.values()],2)

    #Create outputfile for BIC
    folder_path = os.path.join('Output/AUPRC')
    os.makedirs(folder_path, exist_ok=True)
    BIC_filename = os.path.join(folder_path, f'GRNBoost_t{t}.csv')
    
    if not os.path.isfile(BIC_filename):
        with open(BIC_filename, 'w') as f:
            f.write('population,t,p,n,precision,recall,F1,Time(sec),AUPRC')
            f.write('\n')
            
    network={}
    TP=0
    FP=0
    FN=0
    time_network=0
    for pop in range(t):
        print("Population:",pop)
        d1=data[pop]
        d1 = pd.DataFrame(data[pop])
        d1.columns = d1.columns.astype(str)
        #d1.columns=d1.columns.astype(str)	
        
        tstart=time.perf_counter()
        network[pop]= grnboost2(d1)
        tend=time.perf_counter()
        
        theta=network[pop].to_numpy()
        theta=theta.astype(float)
	
        theta_mat = np.zeros((p, p))
        theta_mat[theta[:,0].astype(int),theta[:,1].astype(int)] = theta[:,2]
        np.fill_diagonal(theta_mat, 0)
	
        true_pop = True_Precission[:, :, pop].copy()
        np.fill_diagonal(true_pop, 0)
	
        Precision, Recall, F1  = errorMetrics2(theta_mat[reg_idx,:],true_pop[reg_idx,:])
        time_network+=tend-tstart
        
        

        scores = theta_mat[reg_idx, :].flatten()
        truth  = (true_pop[reg_idx, :] != 0).astype(int).flatten()

        precision_curve, recall_curve, _ = precision_recall_curve(truth, scores)
        auprc = auc(recall_curve, precision_curve)
        fname = f"Output/AUPRC/recall_n{n}_pop{pop+1}.npy"
        np.save(fname, recall_curve)

        fname = f"Output/precision_n{n}_pop{pop+1}.npy"
        np.save(fname, precision_curve)
        with open(BIC_filename, 'a') as f:
            f.write(str(pop)+','+str(t)+','+str(p)+','+str(n)+','+str(Precision)+','+str(Recall)+','+str(F1)+','+str(time_network)+','+str(auprc))
            f.write('\n')





    
        
