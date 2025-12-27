import sys
import pandas as pd

expression_pkl = sys.argv[1]
samples_tsv = sys.argv[2]  

expression_data = pd.read_pickle(expression_pkl)
samples_df = pd.read_csv(samples_tsv, sep="\t")

counts_temp = []
cpm_temp = []
rpkm_temp = []
meta = []

for sample_name, gene_df in expression_data.items():
    
    counts_col = gene_df['Count'].rename(sample_name)
    cpm_col = gene_df['CPM'].rename(sample_name)
    rpkm_col = gene_df['RPKM'].rename(sample_name)

    counts_temp.append(counts_col)
    cpm_temp.append(cpm_col)
    rpkm_temp.append(rpkm_col)
    
    sample_info = samples_df[samples_df['sample'] == sample_name].iloc[0]
    meta.append({
        'sample': sample_name,
        'condition': sample_info['condition'],
        'replicate': sample_info['replicate']
    })

counts_df = pd.concat(counts_temp, axis=1, join='inner')
cpm_df = pd.concat(cpm_temp, axis=1, join='inner')
tpm_df = pd.concat(rpkm_temp, axis=1, join='inner')
tpm_df = tpm_df.div(tpm_df.sum(axis=0), axis=1) * 1e6

meta_df = pd.DataFrame(meta).set_index('sample')

counts_df.to_csv("counts.tsv", sep="\t")
cpm_df.to_csv("cpm.tsv", sep="\t")
tpm_df.to_csv("tpm.tsv", sep="\t")
meta_df.to_csv("meta.tsv", sep="\t")

print("[INFO] Expression matrices generated: counts, CPM, TPM, meta")
