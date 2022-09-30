#!/usr/bin/python3
import pandas as pd

comparisons = str(snakemake.wildcards.comparisons)
c1, c2 = comparisons.split('-')
deseq = pd.read_csv(snakemake.input.deseq_output, 
                names=['gene_id', 'baseMean', 'log2FoldChange', 
                'lfcSE', 'stat', 'pval', 'padj'])
gtf = pd.read_csv(snakemake.input.gtf, header=4, sep='\t', 
            names=['chr', 'source', 'feature', 'start', 
                    'end', 'score', 'strand', 'frame', 
                    'attributes'])
tpm = pd.read_csv(snakemake.input.tpm, sep='\t')

# Get gene_name and biotype per gene from gtf
gtf = gtf[gtf['feature'] == 'gene']
attributes = list(gtf.attributes)

gene_name_d, biotype_d = {}, {}
for attri in attributes:
    att = attri.split(';')
    for a in att:
        if 'gene_id' in a:
            gene_id = a.split('"')[-2]
        elif 'gene_name' in a:
            gene_name = a.split('"')[-2]
            gene_name_d[gene_id] = gene_name
        elif 'gene_biotype' in a:
            biotype = a.split('"')[-2]
            biotype_d[gene_id] = biotype

# Add gene_name and gene_biotype columns to deseq df
deseq['gene_name'] = deseq['gene_id'].map(gene_name_d)
deseq['gene_name'] = deseq['gene_name'].fillna(deseq['gene_id'])
deseq['gene_biotype'] = deseq['gene_id'].map(biotype_d)

# Merge dfs
df = deseq.merge(tpm, how='left', left_on='gene_id', right_on='gene')
df = df.drop(columns=['gene', 'baseMean', 'lfcSE', 'stat'])

# Sort and filter df
df = df[['gene_id', 'gene_name', 'gene_biotype', 'log2FoldChange', 
        'pval', 'padj', f'{c1}_1', f'{c1}_2', f'{c1}_3', f'{c2}_1', 
        f'{c2}_2', f'{c2}_3']]
df = df.sort_values(by=['padj', 'pval'])
df = df.dropna()

print(df)
df.to_csv(snakemake.output.df, sep='\t', index=False)



