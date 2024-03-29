import os
from pathlib import Path

configfile: "config.json"

original_id = list(config['dataset'].values())
simple_id = list(config['dataset'].keys())

# Get DeSeq comparisons (c1 vs c2, fold change with regards to c1)
comparisons, comparisons_full = [], []
with open('data/comparisons.tsv', 'r') as f:
    for line in f:
        if 'cdn1' not in line:  # skip header
            c1, c2 = line.strip('\n').split(' ')  # get the 2 conditions in comparison
            comparisons.append(f'{c1}-{c2}')
            comparisons_full.append(f'results/DESeq2_tximport/{c1}-{c2}.csv')


include: "rules/downloads.smk"
include: "rules/qc_trimming.smk"
include: "rules/kallisto.smk"
include: "rules/deseq2.smk"
include: "rules/final_formatting_and_figures.smk"
include: "rules/transcript_UTR_analysis.smk"

rule all:
    input:
        df = 'results/transcript_analyses.tsv'
        qc_before_trim = expand("data/FastQC/Before_trim/{sample_id}_R1_fastqc.html",
            sample_id=simple_id),
        qc_after_trim = expand("data/FastQC/After_trim/{sample_id}_R1_fastqc.html", 
            sample_id=simple_id),
        pca_plot = 'results/figures/pca_abundance.svg',
        volcano = expand('results/figures/volcano_{comparisons}.svg', comparisons=comparisons)


rule all_downloads:
    input:
        reference_genome = config['path']['reference_genome'],
        gtf = config['path']['reference_annotation']



