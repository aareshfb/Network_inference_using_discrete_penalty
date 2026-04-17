import pandas as pd
import os

# Input lists
seed_list = [1337, 1338, 1339, 1340, 1341]
np_list = [0.25, 0.5, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0 ]

# Output file
results_file = os.path.join("Output", "Vary_np","results.csv")
with open(results_file, 'w') as f:
    f.write('rand_seed,t,p,n,nu_0,mu,gamma,bic value,precision,recall,F1,Time(min)')
    f.write('\n')

# Initialize results DataFrame
results_df = pd.DataFrame()

for seed in seed_list:
    for np_val in np_list:
        np_str = str(np_val).replace('.', 'd')  # e.g., 0.1 → 01
        folder = os.path.join("Output", "Vary_np", str(seed))
        filename = f"ErrorMetrics_L0_np{np_str}.csv"
        file_path = os.path.join(folder, filename)

        try:
            df = pd.read_csv(file_path)
        except FileNotFoundError:
            print(f"File {file_path} not found. Skipping.")
            continue

        if "bic value" not in df.columns:
            print(f"'bic value' column not found in {file_path}. Skipping.")
            continue
        
        
        
        # Get row with smallest BIC
        df["bic value"] = pd.to_numeric(df["bic value"], errors='coerce')
        df = df.dropna(subset=["bic value"])
        
        min_bic_idx = df["bic value"].idxmin()
        best_row = df.loc[min_bic_idx]

        # Tag row with seed and np_val
        best_row["seed"] = seed
        best_row["n_by_p"] = np_val

        results_df = results_df.append(best_row, ignore_index=True)

# Save combined results
results_df.to_csv(results_file, index=False)
print(f"Saved best rows to {results_file}")
