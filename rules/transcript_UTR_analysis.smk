rule UTR_length:
    """ Compare UTR length of genes differentially associated with RPS27A vs RPS27AP5. """
    input:
        gtf = 'data/references/human_ensembl_107.gtf',
        df = "results/differentially_associated_genes.tsv"
    output:
        density = 'results/figures/density/UTR_length.svg',
        df = 'results/UTR_length.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/UTR_length.py"

rule half_life:
    """ Compare half-life of genes differentially associated with RPS27A vs RPS27AP5. """
    input:
        half_life = 'data/references/half_life_Agarwal_et_al_2022.tsv',
        df = "results/differentially_associated_genes.tsv"
    output:
        density = 'results/figures/density/half_life.svg',
        df = 'results/half_life.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/half_life.py"

rule uORF_number:
    """ Compare the number of uORF present in genes differentially associated with RPS27A vs RPS27AP5. 
        The uORF number comes from the uORFdb database (downloaded on Jan 23th 2024) and was prefiltered 
        on the cluster to keep only human uORFs."""
    input:
        uORF = 'data/references/uORFdb_homo_sapiens_filtered.bed',
        gtf = 'data/references/human_ensembl_107.gtf',
        df = "results/differentially_associated_genes.tsv"
    output:
        density = 'results/figures/density/uORF_number.svg',
        df = 'results/uORF.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/uORF_number.py"

rule merge_dfs:
    """ Merge the 3 previous dfs together"""
    input:
        utr_len = rules.UTR_length.output.df,
        half_life = rules.half_life.output.df,
        uORF_number = rules.uORF_number.output.df 
    output:
        df = 'results/transcript_analyses.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_dfs.py"