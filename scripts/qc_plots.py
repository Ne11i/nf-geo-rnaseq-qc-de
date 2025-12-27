import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text

counts_tsv = sys.argv[1]
cpm_tsv = sys.argv[2]
tpm_tsv = sys.argv[3]
meta_tsv = sys.argv[4]

counts_df = pd.read_csv(counts_tsv, sep='\t', index_col=0)
cpm_df = pd.read_csv(cpm_tsv, sep='\t', index_col=0)
meta_df = pd.read_csv(meta_tsv, sep='\t', index_col=0)

print(f"[INFO] Counts shape: {counts_df.shape}")
print(f"[INFO] CPM shape: {cpm_df.shape}")

samples = cpm_df.columns.tolist()

# 1. Total counts per sample
plt.figure(figsize=(6,4))
palette = sns.color_palette("viridis", n_colors=len(counts_df.columns))
counts_df.sum(axis=0).plot(kind='bar', color=palette)
plt.title('Total counts per sample')
plt.ylabel('Counts')
plt.xlabel('Samples')
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('total_counts.png', dpi=300, bbox_inches='tight')
plt.savefig('total_counts.pdf', bbox_inches='tight')
plt.close()

# 2. log2CPM distributions
plt.figure(figsize=(8,6))
for sample in samples:
    sns.kdeplot(np.log2(cpm_df[sample] + 1), label=sample)
plt.xlabel('log2(CPM)')
plt.ylabel('Frequency')
plt.title('Distribution of gene expression across samples (log2CPM)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('log2cpm_distributions.png', dpi=300, bbox_inches='tight')
plt.savefig('log2cpm_distributions.pdf', bbox_inches='tight')
plt.close()

# 3. PCA plot
cpm_df_log = np.log2(cpm_df + 1)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(cpm_df_log.T)
pca = PCA(n_components=2)
pcs = pca.fit_transform(scaled_data)

plt.figure(figsize=(9,6))

foam_samples = meta_df[meta_df['condition'] == 'foam'].index
non_foam_samples = meta_df[meta_df['condition'] == 'non_foam'].index

foam_idx = [samples.index(s) for s in foam_samples]
non_foam_idx = [samples.index(s) for s in non_foam_samples]

plt.scatter(pcs[foam_idx,0], pcs[foam_idx,1], 
           c='#498a52', edgecolor='black', linewidth=0.3, s=50, label='foam')
plt.scatter(pcs[non_foam_idx,0], pcs[non_foam_idx,1], 
           c='#81498a', edgecolor='black', linewidth=0.3, s=50, label='non_foam')

texts = []
for i, name in enumerate(samples):
    texts.append(plt.text(pcs[i,0]+0.02, pcs[i,1]+0.02, name, fontsize=9))

adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.title('PCA on log2CPM')
plt.legend(title='Condition', frameon=True, edgecolor='gray')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pca_log2cpm.png', dpi=300, bbox_inches='tight')
plt.savefig('pca_log2cpm.pdf', bbox_inches='tight')
plt.close()

print("[INFO] QC plots saved: total_counts.*, log2cpm_distributions.*, pca_log2cpm.*")
