#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
TYPES = ["all", "rmdup"]
OUTDIR = "results/expression"

rule all:
    input:
        expand(OUTDIR + "/fpkm/{sample}.{species}.{t}.tsv", sample=SAMPLES, species=SPECIES, t=TYPES),
        expand(OUTDIR + "/feature_count/{sample}.{species}.{t}.tsv", sample=SAMPLES, species=SPECIES, t=TYPES),

def get_input_bam(wildcards):
    if wildcards.t == "all":
        return "results/mapping/filtered/%s.%s.bam" % (wildcards.sample, wildcards.species)
    elif wildcards.t == "rmdup":
        return "results/mapping/rmdup/%s.%s.bam" % (wildcards.sample, wildcards.species)
    assert False

rule calculate_fpkm:
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
        bed = lambda wildcards: get_transcript_bed(wildcards.species),
        tsv = lambda wildcards: get_annotation_tsv(wildcards.species)
    output:
        tsv = OUTDIR + "/fpkm/{sample}.{species}.{t}.tsv"
    log:
        OUTDIR + "/fpkm/{sample}.{species}.{t}.log"
    threads:
        THREADS
    shell:
        """
        nasctools CalculateFPKM \
            --threads {threads} \
            --strand R \
            --layout PE \
            --annotation {input.tsv} \
            {input.bam} {input.bed} {output.tsv} &> {log}
        """

rule feature_count:
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
        gtf = lambda wildcards: get_gtf(wildcards.species)
    output:
        tsv = OUTDIR + "/feature_count/{sample}.{species}.{t}.tsv"
    log:
        OUTDIR + "/feature_count/{sample}.{species}.{t}.log"
    conda:
        "subread"
    threads:
        THREADS
    shell:
        """
        featureCounts -T {threads} \
            -s 2 -p -B \
            -a {input.gtf} \
            -o {output.tsv} \
            {input.bam} &> {log}
        """
