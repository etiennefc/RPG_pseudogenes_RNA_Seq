rule download_genome:
    """Download the reference genome (fasta file) used for this analysis from
        ENSEMBL ftp servers."""
    output:
        genome = config['path']['reference_genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gunzip {output.genome}.gz"

rule download_annotation:
    """Download the annotation (gtf file) used for this analysis."""
    output:
        gtf = config['path']['reference_annotation'],
    params:
        link_annotation = config['download']['annotation'],
    shell:
        "wget --quiet -O {output.gtf}.gz {params.link_annotation} && "
        "gunzip {output.gtf}.gz"

