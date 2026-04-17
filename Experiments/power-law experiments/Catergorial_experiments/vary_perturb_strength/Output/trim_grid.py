import os
import pandas as pd

# List of folder names
folders = ['1337', '1338', '1339', '1340', '1341']

# Loop through each folder
for folder in folders:
    if not os.path.isdir(folder):
        print(f"Folder {folder} does not exist. Skipping.")
        continue

    # Walk through files in the folder
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith('.csv'):
                filepath = os.path.join(root, file)

                try:
                    df = pd.read_csv(filepath)
                    original_len = len(df)

                    # Filter out rows where nu_0 == 0.001
                    df = df[df['nu_0'] != 0.005]

                    # Save back to the same file
                    df.to_csv(filepath, index=False)

                    print(f"Processed {filepath}: removed {original_len - len(df)} rows")

                except Exception as e:
                    print(f"Error processing {filepath}: {e}")
