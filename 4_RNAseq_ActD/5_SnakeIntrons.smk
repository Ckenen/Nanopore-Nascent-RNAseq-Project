#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/rmdup"
outdir = "results/introns"

rule all:
    input:
        expand(outdir + "/counts/{sample}.tsv", sample=samples),
        expand(outdir + "/novel_introns/{sample}.bed", sample=samples),

rule quantify_introns:
    input:
        bed = outdir + "/introns.bed.gz",
        bam = indir + "/{sample}.human.bam"
    output:
        txt = outdir + "/counts/{sample}.tsv"
    log:
        outdir + "/counts/{sample}.log"
    threads:
        8
    shell:
        """
        ./scripts/quantify_introns.py {input.bed} {input.bam} {threads} {output.txt} &> {log}
        """

rule quantify_introns_novel:
    input:
        bam = indir + "/{sample}.human.bam",
        fasta = FILES["human"]["GENOME_FASTA"]
    output:
        bed = outdir + "/novel_introns/{sample}.bed"
    log:
        outdir + "/novel_introns/{sample}.log"
    shell:
        """
        ./scripts/quantify_introns.novel.py {input.bam} {input.fasta} R {output.bed} &> {log}
        """