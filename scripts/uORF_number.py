#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import functions as ft
import collections as coll
from scipy.stats import mannwhitneyu
from pybedtools import BedTool
import subprocess as sp 


df = pd.read_csv(snakemake.input.df, sep='\t')
gtf = read_gtf(snakemake.input.gtf)

# Select only gene features

gtf = gtf[gtf['feature'] == 'gene']

# Create bed of all these genes of interest
df = df.merge(gtf[['gene_id', 'seqname', 'start', 'end', 'strand']], how='left', on='gene_id')
df['score'] = '.'
cols = ['seqname', 'start', 'end', 'gene_id', 'score', 'strand']
df = df[cols + ['associated_protein']]
df.to_csv('temp_bed_genes.bed', sep='\t', index=False, header=False)
bed_genes = BedTool('temp_bed_genes.bed')

# Create bed of all the uORFs (remove duplicates and filter for uORF strength)
temp_uorf = pd.read_csv(snakemake.input.uORF, sep='\t', names=cols + ['uORF_kozak_strength'])
temp_uorf = temp_uorf.drop_duplicates(subset=['seqname', 'start', 'end', 'strand'])
temp_uorf = temp_uorf[temp_uorf['uORF_kozak_strength'] != 'weak']  # romove uORF with weak Kozak sequence
temp_uorf.to_csv('temp_uorf.bed', sep='\t', index=False, header=False)
temp_uorf_bed = BedTool('temp_uorf.bed')

# Return intersect between uORF and the genes
pos = bed_genes.intersect(temp_uorf_bed, c=True, s=True, F=1).to_dataframe()
pos = pos.rename(columns={'thickStart': 'associated_protein', 'thickEnd': 'uORF_number', 'name': 'gene_id'})



pos[['gene_id', 'uORF_number', 'associated_protein']].to_csv(snakemake.output.df, sep='\t', index=False)
print(pos)

# Create density
df_list = [pos[pos['associated_protein'] == 'RPS27A']['uORF_number'], 
        pos[pos['associated_protein'] == 'RPS27AP5']['uORF_number']]

ft.density_x(df_list, "uORF number per gene", "Density", 'linear', "", ['tab:blue', 'tab:orange'], [f'RPS27A-associated ({len(df_list[0])})', 
                f'RPS27AP5-associated ({len(df_list[1])})'], snakemake.output.density)

U, p = mannwhitneyu(list(df_list[0]), list(df_list[1]))
print(p)

sp.call('rm temp_bed_genes.bed temp_uorf.bed', shell=True)