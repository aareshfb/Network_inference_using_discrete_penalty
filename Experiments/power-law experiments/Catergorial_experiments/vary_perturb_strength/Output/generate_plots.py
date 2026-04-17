import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 15})

# Create 'Plots' directory if it doesn't exist
os.makedirs("Plots", exist_ok=True)

# Load data
cat1 = pd.read_csv('Cat_results1_averaged.csv')
cat2 = pd.read_csv('Cat_results2_averaged.csv')
direct1 = pd.read_csv('Direct_results1_averaged.csv')
direct2 = pd.read_csv('Direct_results2_averaged.csv')

# Scale specified columns
# for df in [cat1, cat2, direct1, direct2]:
#     df['perturb_strength'] *= 100
#     df['precision'] *= 100

# Compute averages


# --- Plot 1: F1 vs perturb_strength ---
plt.figure()
x = cat1['perturb_strength']
x = x / (1 + x)
plt.plot(x*100, cat1['F1'], 'bo', markersize=10)#,label='Cat 1'
plt.plot(x*100, cat2['F1'], 'b^', markersize=10)#,label='Cat 2
plt.plot(x*100, direct1['F1'], 'ro', markersize=10)
plt.plot(x*100, direct2['F1'], 'r^', markersize=10)
plt.plot(x*100, 0.5*(cat1['F1'].values+cat2['F1'].values), 'b--', label='categorical ELEM-0')
plt.plot(x*100, 0.5*(direct1['F1']+direct2['F1']), 'r--', label='ELEM-0')
plt.xlabel(r'Local Edge Ratio ($\delta$ %)')
plt.ylabel('F1 score')
# plt.title('F1 vs Perturb Strength')
plt.legend(loc='lower right')
plt.ylim(0., 1.005)
# plt.grid(True)
plt.savefig("Plots/f1.png")
plt.savefig("Plots/f1.tiff")
plt.close()

# --- Plot 2: F1 vs precision ---
plt.figure()
plt.plot(x*100, cat1['precision'], 'bo', markersize=10)
plt.plot(x*100, cat2['precision'], 'b^', markersize=10)
plt.plot(x*100, direct1['precision'], 'ro', markersize=10)
plt.plot(x*100, direct2['precision'], 'r^', markersize=10)
plt.plot(x*100, 0.5*(cat1['precision']+cat2['precision']), 'b--', label='categorical ELEM-0')
plt.plot(x*100, 0.5*(direct1['precision']+direct2['precision']), 'r--', label='ELEM-0')
plt.xlabel(r'Local Edge Ratio ($\delta$ %)')
plt.ylabel('precision')
# plt.title('F1 vs Precision')
plt.legend(loc='lower right')
plt.ylim(0.0, 1.005)
# plt.grid(True)
plt.savefig("Plots/precision.png")
plt.savefig("Plots/precision.tiff")
plt.close()

# --- Plot 3: Recall vs perturb_strength ---
plt.figure()
plt.plot(x*100, cat1['recall'], 'bo', markersize=10)
plt.plot(x*100, cat2['recall'], 'b^', markersize=10)
plt.plot(x*100, direct1['recall'], 'ro', markersize=10)
plt.plot(x*100, direct2['recall'], 'r^', markersize=10)
plt.plot(x*100, 0.5*(cat1['recall']+cat2['recall']), 'b--', label='categorical ELEM-0')
plt.plot(x*100, 0.5*(direct1['recall']+direct2['recall']), 'r--', label='ELEM-0')
plt.xlabel(r'Local Edge Ratio ($\delta$ %)')
plt.ylabel('recall')
# plt.title('Recall vs Perturb Strength')
plt.legend(loc='lower right')
plt.ylim(0.0, 1.005)
# plt.grid(True)
plt.savefig("Plots/recall.png")
plt.savefig("Plots/recall.tiff")
plt.close()

print('Plots saved.')