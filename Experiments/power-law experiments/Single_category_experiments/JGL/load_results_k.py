import pandas as pd
import os

# Input lists
seed_list = [1994, 1995, 1996, 1997, 1998]
np_list = [3, 5, 10, 20, 50, 100 ]

# Output file
results_file = os.path.join("Output","results_k.csv")
with open(results_file, 'w') as f:
    f.write('rand_seed,t,p,n,nu_0,mu,gamma,bic value,precision,recall,F1,Time(min)')
    f.write('\n')

# Initialize results DataFrame
results_df = pd.DataFrame()

for seed in seed_list:
    for np_val in np_list:
        np_str = str(np_val).replace('.', 'd')  # e.g., 0.1 → 01
        folder = os.path.join("Output", str(seed))
        filename = f"JGL_CPU20_L{np_str}.csv"
        file_path = os.path.join(folder, filename)

        try:
            df = pd.read_csv(file_path)
        except FileNotFoundError:
            print(f"File {file_path} not found. Skipping.")
            continue

        if " BIC" not in df.columns:
            print(f"'bic value' column not found in {file_path}. Skipping.")
            continue
        
        
        
        # Get row with smallest BIC
        df["bic value"] = pd.to_numeric(df[" BIC"], errors='coerce')
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
