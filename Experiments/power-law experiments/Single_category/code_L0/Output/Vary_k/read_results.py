import pandas as pd

# Load your CSV file
df = pd.read_csv("results_k.csv")  # Replace with your actual file name

# Group by 'T' and compute mean, min, max for the metrics
metrics = ["precision", "recall", "F1"]
agg_funcs = ['mean', 'min', 'max']

# Aggregate
results = df.groupby("t")[metrics].agg(agg_funcs)

# Flatten the multi-level column index
results.columns = ['_'.join(col).strip() for col in results.columns.values]
results.reset_index(inplace=True)

# Save to CSV
results.to_csv("plot_results.csv", index=False)
print("Saved results to plot_results.csv")
