import os
import numpy as np
import pandas as pd

# method name from parent folder
method = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

# fixed parameters from filename
n_by_p = 20
t = 10

# seeds / experiments
seeds = [1994, 1995, 1996, 1997, 1998]

# base directory
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, "Output", "AUPRC")

file_name = f"{method}_CPU{n_by_p}_L{t}.csv"

summary = []

for seed in seeds:
    seed_dir = os.path.join(data_dir, str(seed))
    file_path = os.path.join(seed_dir, file_name)

    if not os.path.exists(file_path):
        print(f"Missing file: {file_path}")
        continue

    df = pd.read_csv(file_path)

    # expected single precision / recall columns
    recall_col = "Recall"
    precision_col = "Precision"

    # if file uses lowercase instead, handle that too
    if recall_col not in df.columns or precision_col not in df.columns:
        if "recall" in df.columns and "precision" in df.columns:
            recall_col = "recall"
            precision_col = "precision"
        else:
            raise KeyError(
                f"Could not find Recall/Precision columns in {file_path}. "
                f"Columns found: {list(df.columns)}"
            )

    temp = df[[recall_col, precision_col]].copy()

    # replace NaN precision with 1
    temp[precision_col] = temp[precision_col].fillna(1)

    # drop NaN recall
    temp = temp.dropna(subset=[recall_col])

    # sort by recall
    temp = temp.sort_values(by=recall_col)

    recall = temp[recall_col].to_numpy()
    precision = temp[precision_col].to_numpy()

    # precision envelope
    precision_env = np.maximum.accumulate(precision[::-1])[::-1]

    # add (0,1) start point
    recall = np.insert(recall, 0, 0.0)
    precision_env = np.insert(precision_env, 0, 1.0)

    # AUPRC
    if len(recall) <= 1:
        auprc = 0.0
    else:
        auprc = np.sum(np.diff(recall) * precision_env[1:])

    # save curve arrays in the same seed folder
    np.save(os.path.join(seed_dir, f"recall_seed{seed}.npy"), recall)
    np.save(os.path.join(seed_dir, f"precision_envelope_seed{seed}.npy"), precision_env)

    summary.append([method, seed, n_by_p, t, round(auprc, 4)])

summary_df = pd.DataFrame(summary, columns=["Method", "Seed", "CPU", "L", "AUPRC"])

# save per-seed results
summary_df.to_csv(
    os.path.join(data_dir, f"{method}_CPU{n_by_p}_L{t}_auprc_by_seed.csv"),
    index=False
)

# average across seeds
avg_auprc = round(summary_df["AUPRC"].mean(), 4) if not summary_df.empty else np.nan

avg_df = pd.DataFrame(
    [[method, n_by_p, t, avg_auprc]],
    columns=["Method", "CPU", "L", "AUPRC_avg"]
)

avg_df.to_csv(
    os.path.join(data_dir, f"{method}_CPU{n_by_p}_L{t}_auprc_avg.csv"),
    index=False
)

print("Per-seed AUPRC:")
print(summary_df)

print("\nAverage AUPRC across seeds:")
print(avg_df)