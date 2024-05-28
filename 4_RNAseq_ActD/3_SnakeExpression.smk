#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
TYPES = ["all", "rmdup"]
OUTDIR = "results/expression"

rule all:
    input:
        expand(OUTDIR + "/fpkm/{sample}.{species}.{t}.tsv", sample=SAMPLES, species=SPECIES, t=TYPES),

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
        nasctools CalculateFPKM -t {threads} -s R -l PE -a {input.tsv} \
            {input.bam} {input.bed} {output.tsv} &> {log}
        """
