rule format_output_df:
    """ Add gene_name, gene_biotype and TPM 
        to DESeq2 output. Also, filter DESeq2 
        output."""
    input:
        deseq_output = "results/DESeq2_tximport/{comparisons}.csv",
        gtf = 'data/references/human_ensembl_107.gtf',
        tpm = "results/kallisto_combined/tpm.tsv"
    output:
        df = 'results/tables/final_df_{comparisons}.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_output_df.py"

rule pca_abundance:
    """ Create pca plot from abundance (in TPM) values 
        of all genes in all 6 samples."""
    input:
        tpm = "results/kallisto_combined/tpm.tsv"
    output:
        pca_plot = 'results/figures/pca_abundance.svg'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/pca_abundance.py"

rule volcano:
    """ Create volcano plot from comparison between RPS"""
    input:
        deseq = rules.format_output_df.output.df
    output:
        volcano = 'results/figures/volcano_{comparisons}.svg'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/volcano.py"