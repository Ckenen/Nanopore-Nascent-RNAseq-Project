#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/rmdup"
OUTDIR = "results/introns"

rule all:
    input:
        expand(OUTDIR + "/counts/{sample}.tsv", sample=SAMPLES),
        expand(OUTDIR + "/novel_introns/{sample}.bed", sample=SAMPLES),

rule quantify_introns:
    input:
        bed = get_transcript_bed("human"),
        bam = INDIR + "/{sample}.human.bam"
    output:
        txt = OUTDIR + "/counts/{sample}.tsv"
    log:
        OUTDIR + "/counts/{sample}.log"
    threads:
        THREADS
    shell:
        """
        ./scripts/quantify_introns.py {input.bed} \
            {input.bam} {threads} {output.txt} &> {log}
        """

rule quantify_introns_novel:
    input:
        bam = INDIR + "/{sample}.human.bam",
        fasta = get_fasta("human")
    output:
        bed = OUTDIR + "/novel_introns/{sample}.bed"
    log:
        OUTDIR + "/novel_introns/{sample}.log"
    threads:
        THREADS
    shell:
        """
        ./scripts/quantify_introns.novel.py \
            -t {threads} \
            -d R \
            {input.bam} {input.fasta} {output.bed} &> {log}
        """