import pandas as pd
import os

np_vals = [1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
rd_seeds = [10195, 2338, 339, 40, 51341]

# Output file paths
direct_outfile1 = "Direct_results1.csv"
cat_outfile1 = "Cat_results1.csv"
direct_outfile2 = "Direct_results2.csv"
cat_outfile2 = "Cat_results2.csv"

direct_rows = []
cat_rows = []

for rd_seed in rd_seeds:
    for np in np_vals:
        np_str = str(np).replace('.', 'd')
        file_path = os.path.join(str(rd_seed), f"perturb_np{np_str}.csv")

        if not os.path.exists(file_path):
            print(f"Skipping missing file: {file_path}")
            continue

        df = pd.read_csv(file_path)

        for alg_name, collected_rows in [("Direct_alg", direct_rows), ("Cat_alg", cat_rows)]:
            df_alg = df[df['Alg'] == alg_name].copy()

            # Handle missing/empty BIC_mean
            df_alg = df_alg[df_alg['BIC_mean'].notna()]
            df_alg = df_alg[df_alg['BIC_mean'] != '']
            df_alg['BIC_mean'] = pd.to_numeric(df_alg['BIC_mean'], errors='coerce')
            df_alg = df_alg.dropna(subset=['BIC_mean'])

            if df_alg.empty:
                continue

            # Find row index with lowest BIC_mean
            min_idx = df_alg['BIC_mean'].idxmin()

            # Get that row and the one below (if exists)
            row_indices = [min_idx]
            collected_rows.append(df_alg.loc[row_indices])

# Concatenate and save results
if direct_rows:
    pd.concat(direct_rows).to_csv(direct_outfile1, index=False)
if cat_rows:
    pd.concat(cat_rows).to_csv(cat_outfile1, index=False)
    
for rd_seed in rd_seeds:
    for np in np_vals:
        np_str = str(np).replace('.', 'd')
        file_path = os.path.join(str(rd_seed), f"perturb_np{np_str}.csv")

        if not os.path.exists(file_path):
            print(f"Skipping missing file: {file_path}")
            continue

        df = pd.read_csv(file_path)

        for alg_name, collected_rows in [("Direct_alg", direct_rows), ("Cat_alg", cat_rows)]:
            df_alg = df[df['Alg'] == alg_name].copy()

            # Handle missing/empty BIC_mean
            df_alg = df_alg[df_alg['BIC_mean'].notna()]
            df_alg = df_alg[df_alg['BIC_mean'] != '']
            df_alg['BIC_mean'] = pd.to_numeric(df_alg['BIC_mean'], errors='coerce')
            df_alg = df_alg.dropna(subset=['BIC_mean'])

            if df_alg.empty:
                continue

            # Find row index with lowest BIC_mean
            min_idx = df_alg['BIC_mean'].idxmin()

            # Get that row and the one below (if exists)
            row_indices = [min_idx+1]
            df_alg = df[df['Alg'] == alg_name].copy()
            collected_rows.append(df_alg.loc[row_indices])

# Concatenate and save results
if direct_rows:
    pd.concat(direct_rows).to_csv(direct_outfile2, index=False)
if cat_rows:
    pd.concat(cat_rows).to_csv(cat_outfile2, index=False)

print("Saved results to Direct_results.csv and Cat_results.csv")