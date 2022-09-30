#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Configure comparison, pval_threshold and colors
comparisons = str(snakemake.wildcards.comparisons)
comp1, comp2 = comparisons.split('-')
DE_tool, quantifier = 'DESeq2', 'Kallisto'
pval_threshold = 0.05
colors = {'|log2FC| > 1 & padj<'+str(pval_threshold): 'lightgreen',
            'Not significative': 'grey'}

# Load DE df
df = pd.read_csv(snakemake.input.deseq, sep='\t')

# Drop genes/transcripts with NaN in log2FoldChange, pvalue and/or padj
df = df.dropna(subset=['log2FoldChange', 'pval', 'padj'])

# Create -log10padj column
df['-log10padj'] = -np.log10(df['padj'])

# Create hue column for significative points (|log2FC| > 1 & padj<0.05)
df.loc[(df['padj'] < pval_threshold) & (np.abs(df['log2FoldChange']) > 1),
        'Statistical significance'] = '|log2FC| > 1 & padj<'+str(pval_threshold)
df['Statistical significance'] = df['Statistical significance'].fillna('Not significative')

# Create volcano function
def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            pval_threshold, **kwargs):
    """
    Creates a violin plot (using a x, y and hue column).
    """
    sns.set_theme()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["legend.loc"] = 'center'

    plt.suptitle(title, fontsize=20)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    # Add threshold lines (padj)
    plt.axhline(y=-np.log10(pval_threshold), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(2), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(0.5), color='black', ls='--', lw=0.5)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.savefig(path, bbox_inches='tight', dpi=600)

# Create volcano
volcano(df, 'log2FoldChange', '-log10padj', 'Statistical significance',
        'log2(Fold change)', '-log10(FDR-adjusted p-value)',
        f'Comparison between {comp1} and {comp2}\nusing {DE_tool} after {quantifier}',
        colors, snakemake.output.volcano, pval_threshold)