#!/usr/bin/python3
import pandas as pd

# Load dfs
utr_len = pd.read_csv(snakemake.input.utr_len, sep='\t')
utr_len = utr_len[['gene_id', 'gene_name', 'associated_protein', 'UTR_length']]

half_life = pd.read_csv(snakemake.input.half_life, sep='\t')
half_life = half_life[['gene_id', 'half-life_PC1']]

uorf_nb = pd.read_csv(snakemake.input.uORF_number, sep='\t')
uorf_nb = uorf_nb[['gene_id', 'uORF_number']]

# Merge dfs
print(utr_len)
print(half_life)
print(uorf_nb)

df = utr_len.merge(half_life, how='left', on='gene_id')
df = df.merge(uorf_nb, how='left', on='gene_id')
df = df.sort_values(by='associated_protein')
print(df)
df.to_csv(snakemake.output.df, index=False, sep='\t')