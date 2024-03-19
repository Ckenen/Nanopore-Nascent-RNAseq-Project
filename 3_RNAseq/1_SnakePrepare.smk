#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "data/datasets"
outdir = "results/prepare"

rule all:
    input:
        # expand(outdir + "/cutadapt/{sample}_1.fastq.gz", sample=samples),
        expand(outdir + "/bowtie2/{sample}.bam", sample=samples),

rule cutadapt:
    input:
        fq1 = indir + "/{sample}_1.fastq.gz",
        fq2 = indir + "/{sample}_2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/{sample}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/{sample}_2.fastq.gz"
    log:
        outdir + "/cutadapt/{sample}.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2: # remove rRNA
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = config["bt2_rrna_idx"]
    output:
        bam = outdir + "/bowtie2/{sample}.bam"
    log:
        outdir + "/bowtie2/{sample}.log"
    threads:
        8
    shell:
        """
        ./scripts/prepare/remove_rrna.pe.sh {input} {threads} {output} &> {log}
        """