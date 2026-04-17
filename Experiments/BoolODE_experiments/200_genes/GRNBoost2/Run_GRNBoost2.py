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
import re



def read_prec(no_reg,p,t):
    base_path = "../Ground_truth"
    prec = {}

    expected_regs = [f"reg{i}" for i in range(1, no_reg+1)]
    expected_genes = [f"gene{i}" for i in range(1, p-no_reg+1)]

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

        if arr.shape != (no_reg, p):
            raise ValueError(
                f"{path}: expected shape (no_reg, p), found {arr.shape}"
            )

        prec[i - 1] = arr
        print(f"Loaded {path}: shape={arr.shape}")

    return prec

def load_dataset(no_reg,p,t):
    base_path = "../Expression_data"
    rng = np.random.default_rng(0)

    expected_regs = [f"reg{i}" for i in range(1, no_reg+1)]
    expected_genes = [f"gene{i}" for i in range(1, p-no_reg+1)]

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



if __name__ == "__main__":

    file_path='../generated_data/'
    file_name = 'simulation_clust_'
    file_ext='.csv'
    n = 2000
    t = 10
    p = 200
    no_reg=30
    print('generating expression data...')
    data = load_dataset(no_reg,p,t)
    
    print('reading true precision matrices')
    True_Precission = read_prec(no_reg,p,t)
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
    all_scores = []
    all_truth  = []
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
	
	    # ----- GLOBAL COUNTS -----
        pred = theta_mat[:no_reg, :]
        truth = true_pop[:, :]
    
        pred_bin  = (pred != 0)
        truth_bin = (truth != 0)
    
        TP += np.sum(pred_bin & truth_bin)
        FP += np.sum(pred_bin & ~truth_bin)
        FN += np.sum(~pred_bin & truth_bin)

	
        Precision, Recall, F1  = errorMetrics2(theta_mat[:no_reg,:],true_pop[:,:])
        time_network+=tend-tstart
        
        scores = theta_mat[:no_reg, :].flatten()
        truth  = (true_pop[:, :] != 0).astype(int).flatten()
        
        all_scores.append(scores)
        all_truth.append(truth)

        precision_curve, recall_curve, _ = precision_recall_curve(truth, scores)
        auprc = auc(recall_curve, precision_curve)
        fname = f"Output/AUPRC/recall_n{n}_pop{pop+1}.npy"
        np.save(fname, recall_curve)

        fname = f"Output/precision_n{n}_pop{pop+1}.npy"
        np.save(fname, precision_curve)
        with open(BIC_filename, 'a') as f:
            f.write(str(pop)+','+str(t)+','+str(p)+','+str(n)+','+str(Precision)+','+str(Recall)+','+str(F1)+','+str(time_network)+','+str(auprc))
            f.write('\n')

    # ===== GLOBAL METRICS =====
    Precision_global = TP / (TP + FP)
    Recall_global    = TP / (TP + FN)
    F1_global        = 2 * Precision_global * Recall_global / (Precision_global + Recall_global)
    
    print("\n===== GLOBAL METRICS =====")
    print("Precision:", Precision_global)
    print("Recall:", Recall_global)
    print("F1:", F1_global)
    
    # ===== GLOBAL AUPRC =====
    all_scores = np.concatenate(all_scores)
    all_truth  = np.concatenate(all_truth)
    
    precision_curve_g, recall_curve_g, _ = precision_recall_curve(all_truth, all_scores)
    auprc_global = auc(recall_curve_g, precision_curve_g)
    
    print("Global AUPRC:", auprc_global)
    # Append to file
    
    with open(BIC_filename, 'a') as f:
        f.write("GLOBAL"+','+str(t)+','+str(p)+','+str(n)+','+
                str(Precision_global)+','+str(Recall_global)+','+
                str(F1_global)+','+str(time_network)+','+str(auprc_global))
        f.write('\n')
    
    
    
        
            
