rule DESeq2_genes:
    """ Differential expression per gene for the different conditions. The
        samples in design.tsv must match exactly the values (not the keys) of
        the 'datasets' dictionary located in the config.json file."""
    input:
        counts = rules.kallisto_combine_quantification.output.est_counts,
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/genes"),
    params:
        dir = "results/DESeq2/genes"
    log:
        "logs/DESeq2/genes.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/DESeq2_genes.R"  ##TO modify