import pandas as pd
import os

ptx_vals = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
rd_seeds = [1337, 1338, 1339, 1340, 1341]

# Output files

output_dir = "Output"
os.makedirs(output_dir, exist_ok=True)

# Output file paths
direct_outfile = os.path.join(output_dir, "Direct_results.csv")
cat_outfile = os.path.join(output_dir, "Cat_results.csv")

direct_rows = []
cat_rows = []

for rd_seed in rd_seeds:
    for ptx in ptx_vals:
        ptx_str = str(ptx).replace('.', 'd')
        file_path = os.path.join("Output",str(rd_seed), f"perturb_{ptx_str}.csv")

        if not os.path.exists(file_path):
            print(f"Skipping missing file: {file_path}")
            continue

        df = pd.read_csv(file_path)

        for alg_name, collected_rows in [("Direct_alg", direct_rows), ("Cat_alg", cat_rows)]:
            df_alg = df[df['Alg'] == alg_name].copy()

            # Handle missing/empty bic_mean
            df_alg = df_alg[df_alg['bic_mean'].notna()]
            df_alg = df_alg[df_alg['bic_mean'] != '']
            df_alg['bic_mean'] = pd.to_numeric(df_alg['bic_mean'], errors='coerce')
            df_alg = df_alg.dropna(subset=['bic_mean'])

            if df_alg.empty:
                continue

            # Find row index with lowest bic_mean
            min_idx = df_alg['bic_mean'].idxmin()

            # Get that row and the one below (if exists)
            row_indices = [min_idx]
            if min_idx + 1 in df_alg.index:
                row_indices.append(min_idx + 1)

            collected_rows.append(df_alg.loc[row_indices])

# Concatenate and save results
if direct_rows:
    pd.concat(direct_rows).to_csv(direct_outfile, index=False)
if cat_rows:
    pd.concat(cat_rows).to_csv(cat_outfile, index=False)

print("Saved results to Direct_results.csv and Cat_results.csv")
