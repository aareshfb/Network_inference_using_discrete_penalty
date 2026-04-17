import os
import pandas as pd
from io import StringIO

# List of folder names
folders = ['1337', '1338', '1339', '1340', '1341']
# List of t values used in file names
t_values = [3, 5, 10, 20, 50, 100]

# Correct column names
column_names = [
    'rd_seed', 'Alg', 't', 'p', 'n', 'Cat', 'nu_0', 'mu', 'gamma',
    'precision', 'recall', 'F1', 'Time', 'BIC', 'BIC_mean'
]

for folder in folders:
    folder_path = os.path.join('.', folder)

    for t in t_values:
        file_name = f'perturb_T{t}.csv'
        file_path = os.path.join(folder_path, file_name)

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        with open(file_path, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            print(f"Too few lines to process: {file_path}")
            continue

        # Skip the first line (malformed header)
        data = ''.join(lines[1:])
        df = pd.read_csv(StringIO(data), header=None)

        # Assign correct column names
        df.columns = column_names

        # Save the corrected file
        df.to_csv(file_path, index=False)
        print(f"Fixed and saved: {file_path}")
