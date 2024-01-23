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

half_life = pd.read_csv(snakemake.input.half_life, sep='\t')
print(list(half_life.gene_id))
print([i for i in half_life.gene_id if i in df.gene_id])

# Merge dfs
df = df.merge(half_life[['gene_id', 'half-life_PC1']], how='left', on='gene_id')
df['half-life_PC1'] = df['half-life_PC1'].astype(float)
df.to_csv(snakemake.output.df, sep='\t', index=False)
print(df)


print(coll.Counter(df[df['half-life_PC1'].isna()]['associated_protein']))
df = df[~df['half-life_PC1'].isna()]

# Create density

df_list = [df[df['associated_protein'] == 'RPS27A']['half-life_PC1'], 
        df[df['associated_protein'] == 'RPS27AP5']['half-life_PC1']]

ft.density_x(df_list, "Transcript half-life (PC1)", "Density", 'linear', "", ['tab:blue', 'tab:orange'], [f'RPS27A-associated ({len(df_list[0])})', 
                f'RPS27AP5-associated ({len(df_list[1])})'], snakemake.output.density)

U, p = mannwhitneyu(list(df_list[0]), list(df_list[1]))
print(p)

