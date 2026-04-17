import os
import numpy as np
import pandas as pd

n_values = [200, 1000, 2000]
method = os.path.basename(os.path.abspath(os.path.dirname(__file__)))

data_dir = "Output/AUPRC"

summary = []

for n in n_values:

    file = os.path.join(data_dir, f"{method}_{n}.csv")
    df = pd.read_csv(file)

    for j in [1,2,3]:

        recall_col = f"Recall_{j}"
        precision_col = f"Precision_{j}"

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
        
        summary.append([method, n, j, round(auprc,2)])

summary_df = pd.DataFrame(summary, columns=["Method","n","Dataset","AUPRC"])
summary_df.to_csv(os.path.join(data_dir, f"{method}_auprc_summary.csv"), index=False)

print(summary_df)