#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/marked_strand"
outdir = "results/expression"

rule all:
    input:
        expand(outdir + "/fpkm/{run_cell}.tsv", run_cell=run_cells),

rule calculate_fpkm:
    input:
        bam = indir + "/{run}/{cell}.bam",
        bed = TRANSCRIPT_BED_GZ,
        anno = ANNOTATION_TSV
    output:
        txt = outdir + "/fpkm/{run}/{cell}.tsv"
    log:
        outdir + "/fpkm/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM --threads {threads} --layout PE --strand TAG --strand-tag ST \
            --annotation {input.anno} {input.bam} {input.bed} {output.txt} &> {log}
        """