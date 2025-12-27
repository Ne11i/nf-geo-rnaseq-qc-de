import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import mygene
import warnings
warnings.filterwarnings('ignore')

counts_tsv = sys.argv[1]
meta_tsv = sys.argv[2]

counts_df = pd.read_csv(counts_tsv, sep='\t',index_col=0)
meta_df = pd.read_csv(meta_tsv, sep='\t')
meta_df = meta_df.set_index('sample')

print(f"[INFO] Counts shape: {counts_df.shape}")
print(f"[INFO] Meta shape: {meta_df.shape}")

counts_filtered = counts_df[(counts_df > 10).sum(axis=1) >= 3]
print(f"[INFO] Filtered counts shape: {counts_filtered.shape}")

deseq_data = counts_filtered.T
deseq_meta = meta_df

dds = DeseqDataSet(counts=deseq_data.astype(int), 
                   metadata=deseq_meta, 
                   design_factors="condition", 
                   refit_cooks=True)
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition", "foam", "non_foam"])
stat_res.summary()
results_df = stat_res.results_df

results_df.to_csv("deseq_results.csv")
print("[INFO] DESeq2 results saved")

mg = mygene.MyGeneInfo()
ens_ids = results_df.index.tolist()
query_res = mg.querymany(ens_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
id_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_res}
results_df.index = results_df.index.map(id_to_symbol)

results_filtered = results_df.dropna(subset=['padj', 'log2FoldChange']).copy()
results_filtered['-log10(padj)'] = -np.log10(results_filtered['padj'])
results_filtered['change_status'] = 'ns'




results_filtered.loc[(results_filtered["padj"] < 0.05) & 
                    (results_filtered["log2FoldChange"] > 1), "change_status"] = "Up"
results_filtered.loc[(results_filtered["padj"] < 0.05) & 
                    (results_filtered["log2FoldChange"] < -1), "change_status"] = "Down"

top_genes = results_filtered.nsmallest(10, "padj")
top_genes.to_csv("top_genes.csv")
print("[INFO] Top genes saved")

plt.figure(figsize=(9,6))
sns.scatterplot(data=results_filtered, 
                x="log2FoldChange", 
                y="-log10(padj)", 
                hue="change_status",
                palette={"Up": "#ad455c", "Down": "#4561ad", "ns": "grey"}, 
                s=30, alpha=0.6)

for _, row in top_genes.iterrows():
    plt.text(row["log2FoldChange"] + 0.15, 
             row["-log10(padj)"], 
             row.name, fontsize=8, ha='left')

plt.axvline(x=1, color='grey', linestyle='--', alpha=0.7)
plt.axvline(x=-1, color='grey', linestyle='--', alpha=0.7)
plt.axhline(y=-np.log10(0.05), color='grey', linestyle='--', alpha=0.7)
plt.title("Volcano plot: foam vs non_foam")
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(FDR adjusted p-value)")
plt.legend(title="Expression change")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
plt.savefig('volcano_plot.pdf', bbox_inches='tight')
plt.close()

print("[INFO] Volcano plot saved")
print(f"[INFO] Found {len(results_filtered[results_filtered['change_status'] != 'ns'])} significant genes")
print("[INFO] Differential expression analysis complete!")
