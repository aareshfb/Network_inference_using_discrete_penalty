import os
import numpy as np
import pandas as pd

n_values = [2000]
nu_values = [1e-3,5e-3,1e-4,5e-4,1e-5,5e-5]
method = os.path.basename(os.path.abspath(os.path.dirname(__file__)))

data_dir = "Output/AUPRC/masked_diag"

summary = []

for n in n_values:
    for nu in nu_values:

        file = os.path.join(data_dir, f"{method}_n{n}_g1_nu{nu}.csv")
        if not os.path.isfile(file):
            print(f"Skipping missing file: {file}")
            continue
        df = pd.read_csv(file)
    
        for j in ["avg"]:
    
            recall_col = f"recall_{j}"
            precision_col = f"precision_{j}"
    
            temp = df[[recall_col, precision_col]].copy()
    
            # replace NaN precision with 1
            temp[precision_col] = temp[precision_col].fillna(1)
    
            # drop NaN recall
            temp = temp.dropna(subset=[recall_col])
    
            # sort pairs by recall
            temp = temp.sort_values(by=recall_col)
    
            recall = temp[recall_col].to_numpy()
            precision = temp[precision_col].to_numpy()
    
            # precision envelope
            precision_env = np.maximum.accumulate(precision[::-1])[::-1]
    
    
            # ADDED THIS (0,1) START POINT
            recall = np.insert(recall, 0, 0.0)
            precision_env = np.insert(precision_env, 0, 1.0)
    
            # AUPRC calculation
            if len(recall) <= 1:
                auprc = 0
            else:
                auprc = np.sum(np.diff(recall) * precision_env[1:])
    
            # save arrays
            np.save(os.path.join(data_dir, f"recall_n{n}_j{j}.npy"), recall)
            np.save(os.path.join(data_dir, f"precision_envelope_n{n}_j{j}.npy"), precision_env)
            
            summary.append([method, n, nu, j, round(auprc,2)])

summary_df = pd.DataFrame(summary, columns=["Method","n","nu","Dataset","AUPRC"])
summary_df.to_csv(os.path.join(data_dir, f"{method}_auprc_summary3.csv"), index=False)

print(summary_df)