#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
types = ["all", "rmdup"]
outdir = "results/expression"

rule all:
    input:
        expand(outdir + "/fpkm/{sample}.{species}.{t}.tsv", sample=samples, species=species_list, t=types),
        expand(outdir + "/feature_count/{sample}.{species}.{t}.tsv", sample=samples, species=species_list, t=types),

def get_bam(wildcards):
    if wildcards.t == "all":
        return "results/mapping/filtered/%s.%s.bam" % (wildcards.sample, wildcards.species)
    elif wildcards.t == "rmdup":
        return "results/mapping/rmdup/%s.%s.bam" % (wildcards.sample, wildcards.species)
    assert False

rule calculate_fpkm:
    input:
        bam = lambda wildcards: get_bam(wildcards),
        bed = lambda wildcards: FILES[wildcards.species]["TRANSCRIPT_BED_GZ"],
        tsv = lambda wildcards: FILES[wildcards.species]["ANNOTATION_TSV"]
    output:
        tsv = outdir + "/fpkm/{sample}.{species}.{t}.tsv"
    log:
        outdir + "/fpkm/{sample}.{species}.{t}.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM --threads {threads} --strand R --layout PE --annotation {input.tsv} {input.bam} {input.bed} {output.tsv} &> {log}
        """

rule feature_count:
    input:
        bam = lambda wildcards: get_bam(wildcards),
        gtf = lambda wildcards: FILES[wildcards.species]["ANNOTATION_GTF"]
    output:
        tsv = outdir + "/feature_count/{sample}.{species}.{t}.tsv"
    log:
        outdir + "/feature_count/{sample}.{species}.{t}.log"
    threads:
        8
    shell:
        """
        featureCounts -T {threads} -s 2 -p -B -a {input.gtf} -o {output.tsv} {input.bam} &> {log}
        """