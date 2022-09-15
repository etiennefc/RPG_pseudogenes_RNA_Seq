rule DESeq2:
    """ Differential expression for the different conditions """
    input:
        quant = expand("results/kallisto/{sample_id}/abundance.tsv",
                        sample_id=simple_id),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = rules.generate_transcriptID_geneID.output.map
    output:
        results = directory("results/DESeq2_tximport/"),
        out_files = comparisons_full
    params:
        kallisto_dir = "results/kallisto"
    log:
        "logs/DESeq2.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/DESeq2_kallisto_tximport.R"
