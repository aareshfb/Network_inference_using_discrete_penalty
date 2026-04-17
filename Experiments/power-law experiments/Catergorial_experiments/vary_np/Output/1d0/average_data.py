import pandas as pd

# List of files and output names
files = [
    ('Cat_results1.csv', 'Cat_results1_averaged.csv'),
    ('Cat_results2.csv', 'Cat_results2_averaged.csv'),
    ('Direct_results1.csv', 'Direct_results1_averaged.csv'),
    ('Direct_results2.csv', 'Direct_results2_averaged.csv')
]

for input_file, output_file in files:
    df = pd.read_csv(input_file)

    # Group by 'n' and average numeric columns
    # This will average over all 'rd_seed' values for each 'n'
    avg_df = df.groupby('n', as_index=False).mean()

    # Save averaged dataframe to CSV
    avg_df.to_csv(output_file, index=False)

    print(f"Averaged file saved: {output_file}")
