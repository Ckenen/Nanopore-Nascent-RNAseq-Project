#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/rmdup"
OUTDIR = "results/expression"

rule all:
    input:
        expand(OUTDIR + "/fpkm/{sample}.{species}.tsv", sample=SAMPLES, species=SPECIES)

rule calculate_fpkm:
    input:
        bam = INDIR + "/{sample}.{species}.bam",
        bed = lambda wildcards: get_transcript_bed(wildcards.species),
        tsv = lambda wildcards: get_annotation_tsv(wildcards.species)
    output:
        tsv = OUTDIR + "/fpkm/{sample}.{species}.tsv"
    log:
        OUTDIR + "/fpkm/{sample}.{species}.log"
    threads:
        THREADS
    shell:
        """
        nasctools CalculateFPKM -t {threads} -s R -l PE -a {input.tsv} \
            {input.bam} {input.bed} {output.tsv} &> {log}
        """
