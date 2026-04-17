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
#     df['t'] *= 100
#     df['precision'] *= 100

# Compute averages


# --- Plot 1: F1 vs t ---
plt.figure()
plt.plot(cat1['t'], cat1['F1'], 'bo', markersize=10 )#,label='Cat 1'
plt.plot(cat2['t'], cat2['F1'], 'b^', markersize=10 )#,label='Cat 2
plt.plot(direct1['t'], direct1['F1'], 'ro', markersize=10)
plt.plot(direct2['t'], direct2['F1'], 'r^', markersize=10)
plt.plot(cat1['t'], 0.5*(cat1['F1'].values+cat2['F1'].values), 'b--', label='categorical ELEM-0')
plt.plot(direct1['t'], 0.5*(direct1['F1']+direct2['F1']), 'r--', label='ELEM-0')
plt.xlabel(r'K')
plt.ylabel('F1 score')
# plt.title('F1 vs K')
plt.legend(loc='upper right')
plt.ylim(0., 1.005)
# plt.grid(True)
plt.savefig("Plots/f1_t.png", bbox_inches='tight')
plt.savefig("Plots/f1_t.tiff")
plt.close()

# --- Plot 2: F1 vs t ---
plt.figure()
plt.plot(cat1['t'], cat1['precision'], 'bo', markersize=10)
plt.plot(cat2['t'], cat2['precision'], 'b^', markersize=10)
plt.plot(direct1['t'], direct1['precision'], 'ro', markersize=10)
plt.plot(direct2['t'], direct2['precision'], 'r^', markersize=10)
plt.plot(cat1['t'], 0.5*(cat1['precision']+cat2['precision']), 'b--', label='categorical ELEM-0')
plt.plot(direct1['t'], 0.5*(direct1['precision']+direct2['precision']), 'r--', label='ELEM-0')
plt.xlabel(r'K')
plt.ylabel('precision')
# plt.title(' Precision vs K')
plt.legend(loc='lower right')
plt.ylim(0., 1.005)
# plt.grid(True)
plt.savefig("Plots/precision_t.png", bbox_inches='tight')
plt.savefig("Plots/precision_t.tiff")
plt.close()

# --- Plot 3: Recall vs t ---
plt.figure()
plt.plot(cat1['t'], cat1['recall'], 'bo', markersize=10)
plt.plot(cat2['t'], cat2['recall'], 'b^', markersize=10)
plt.plot(direct1['t'], direct1['recall'], 'ro', markersize=10)
plt.plot(direct2['t'], direct2['recall'], 'r^', markersize=10)
plt.plot(cat1['t'], 0.5*(cat1['recall']+cat2['recall']), 'b--', label='categorical ELEM-0')
plt.plot(direct1['t'], 0.5*(direct1['recall']+direct2['recall']), 'r--', label='ELEM-0')
plt.xlabel(r'K')
plt.ylabel('recall')
# plt.title('Recall vs K')
plt.legend(loc='upper right')
plt.ylim(0., 1.005)
# plt.grid(True)
plt.savefig("Plots/recall_t.png", bbox_inches='tight')
plt.savefig("Plots/recall_t.tiff")
plt.close()

print('Plots saved.')