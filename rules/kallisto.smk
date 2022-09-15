rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = rules.download_genome.output.genome,
        gtf = rules.download_annotation.output.gtf
    output:
        transcriptome = 'data/references/transcriptome.fa'
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.transcriptome}"

rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        transcriptome = rules.create_transcriptome.output.transcriptome
    output:
        idx = 'data/references/kallisto.idx'
    params:
        kmer = "31"
    log:
        "logs/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"

rule kallisto_quant:
    """ Generates counts per transcripts using Kallisto pseudo-alignment."""
    input:
        idx = rules.kallisto_index.output.idx,
        trimmed_fq_R1 = rules.fastp.output.fastq1,
        trimmed_fq_R2 = rules.fastp.output.fastq2
    output:
        quant = "results/kallisto/{sample_id}/abundance.tsv"
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{sample_id}"
    log:
        "logs/kallisto/{sample_id}.log"
    threads:
        1
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.trimmed_fq_R1} {input.trimmed_fq_R2} "
        "&> {log}"

rule generate_transcriptID_geneID:
    """
    Generating a two-column text file containing the gene -> transcript
    relationship
    """
    input:
        gtf = rules.download_annotation.output.gtf
    output:
        map = "data/references/transcriptID_geneID.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_transcriptID_geneID.py"

rule kallisto_combine_quantification:
    """
    Custom Python script to collect and format Kallisto results for further
    processing.
    """
    input:
        datasets = expand("results/kallisto/{sample_id}/abundance.tsv",
                           sample_id=simple_id),
        map = rules.generate_transcriptID_geneID.output.map
    output:
        tpm = "results/kallisto_combined/tpm.tsv",
        est_counts = "results/kallisto_combined/est_counts.tsv",
        transcript_tpm = "results/kallisto_combined/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto_combined/transcript_est_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/kallisto-salmon_combine_quantification.py"