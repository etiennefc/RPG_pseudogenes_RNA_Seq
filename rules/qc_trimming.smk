import os
rule rename_samples:
    """Rename samples with a nicer understandable name"""
    output:
        renamed_fastq = expand(
            "data/references/fastq/{sample_id}_R{pair}.fastq.gz",
            sample_id=simple_id, pair=[1, 2])
    run:
        for new_name, old_name in config['dataset'].items():
            for num in [1, 2]:
                old = "data/references/fastq/{}_R{}.fastq.gz".format(old_name, num)
                new = "data/references/fastq/{}_R{}.fastq.gz".format(new_name, num)
                print(old, new)
                os.rename(old, new)

rule qc_before_trim:
    "Assess fastq quality before trimming reads"
    input:
        fastq1 = "data/references/fastq/{sample_id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{sample_id}_R2.fastq.gz"
    output:
        qc_report1 = "data/FastQC/Before_trim/{sample_id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/Before_trim/{sample_id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/Before_trim"
    log:
        "logs/fastqc/before_trim/{sample_id}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"

rule fastp:
    """ Trim reads from fastq files using fastp."""
    input:
        fastq1 = "data/references/fastq/{sample_id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{sample_id}_R2.fastq.gz"
    output:
        fastq1 = "results/fastp/{sample_id}_R1.fq",
        fastq2 = "results/fastp/{sample_id}_R2.fq",
        unpaired_fastq1 = "results/fastp/{sample_id}_R1.unpaired.fq",
        unpaired_fastq2 = "results/fastp/{sample_id}_R2.unpaired.fq",
        html_report = "results/fastp/{sample_id}_fastp.html",
        json_report = "results/fastp/{sample_id}_fastp.json"
    threads:
        8
    params:
        options = ["--qualified_quality_phred 30", "--length_required 20", 
                "--cut_window_size 1", "--cut_mean_quality 30", "--cut_front",
                "--cut_tail"]
    log:
        "logs/fastp/{sample_id}.log"
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.fastq1} -I {input.fastq2} "
        "-o {output.fastq1} -O {output.fastq2} "
        "--unpaired1 {output.unpaired_fastq1} "
        "--unpaired2 {output.unpaired_fastq2} "
        "--thread {threads} "
        "-h {output.html_report} "
        "-j {output.json_report} "
        "{params.options} "
        "&> {log}"

rule qc_after_trim:
    "Assess fastq quality after trimming reads"
    input:
        fastq1 = rules.fastp.output.fastq1,
        fastq2 = rules.fastp.output.fastq2
    output:
        qc_report1 = "data/FastQC/After_trim/{sample_id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/After_trim/{sample_id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/After_trim"
    log:
        "logs/fastqc/after_trim/{sample_id}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"
