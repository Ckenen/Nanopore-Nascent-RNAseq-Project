#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/rmdup"
OUTDIR = "results/introns"

rule all:
    input:
        expand(OUTDIR + "/counts/{sample}.tsv", sample=SAMPLES),

rule quantify_introns:
    input:
        bam = INDIR + "/{sample}.human.bam",
        fa = get_fasta("human"),
        bed = get_transcript_bed("human")
    output:
        tsv = OUTDIR + "/counts/{sample}.tsv"
    log:
        OUTDIR + "/counts/{sample}.log"
    threads:
        THREADS
    shell:
        """
        ./scripts/quantify_introns.py -t {threads} -d R \
            -f {input.fa} -b {input.bed} \
            {input.bam} {output.tsv} &> {log}
        """