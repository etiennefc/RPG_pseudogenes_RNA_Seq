#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import functions as ft
import collections as coll
from scipy.stats import mannwhitneyu


df = pd.read_csv(snakemake.input.df, sep='\t')
print(df)
gtf = read_gtf(snakemake.input.gtf)

# Select only 5'UTR features

gtf = gtf[gtf['feature'] == 'five_prime_utr']

# Merge 5'UTR info to df
df = df.merge(gtf[['gene_id', 'start', 'end']], how='left', on='gene_id')
df['UTR_length'] = df['end'] - df['start'] + 1

# Keep longest 5'UTR per gene
df = df.sort_values(by=['gene_id', 'UTR_length'], ascending=[True, False])
df = df.drop_duplicates(subset='gene_id')

df.to_csv(snakemake.output.df, sep='\t', index=False)
print(df)
print(coll.Counter(df[df['UTR_length'].isna()]['associated_protein']))

df = df[~df['UTR_length'].isna()]

# Create density
df_list = [df[df['associated_protein'] == 'RPS27A']['UTR_length'], 
        df[df['associated_protein'] == 'RPS27AP5']['UTR_length']]

ft.density_x(df_list, "5'UTR length", "Density", 'linear', "", ['tab:blue', 'tab:orange'], [f'RPS27A-associated ({len(df_list[0])})', 
                f'RPS27AP5-associated ({len(df_list[1])})'], snakemake.output.density)

U, p = mannwhitneyu(list(df_list[0]), list(df_list[1]))
print(p)
