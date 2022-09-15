import os
from pathlib import Path

configfile: "config.json"

original_id = list(config['dataset'].values())
simple_id = list(config['dataset'].keys())

include: "rules/downloads.smk"
include: "rules/qc_trimming.smk"
include: "rules/kallisto.smk"
include: "rules/deseq2.smk"

rule all:
    input:
        qc_before_trim = expand("data/FastQC/Before_trim/{id}_R1_fastqc.html",
            id=simple_id),
        qc_after_trim = expand("data/FastQC/After_trim/{id}_R1_fastqc.html", 
            id=simple_id),
        kallisto_quant = expand("results/kallisto/{id}/abundance.tsv", id=simple_id),
        deseq = "results/DESeq2/genes"


rule all_downloads:
    input:
        reference_genome = config['path']['reference_genome'],
        gtf = config['path']['reference_annotation']



