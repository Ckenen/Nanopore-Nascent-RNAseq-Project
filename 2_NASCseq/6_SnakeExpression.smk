#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mismatch/marked_nascent"
OUTDIR = "results/expression"
# RUN_CELLS = list(filter(lambda x: "_SE" not in x, RUN_CELLS))

rule all:
    input:
        expand(OUTDIR + "/fpkm/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/genes/{run_cell}.tsv", run_cell=RUN_CELLS),

rule calculate_fpkm:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        bed = config["TRANSCRIPT_BED_GZ"],
        txt = config["ANNOTATION_TSV"]
    output:
        txt = OUTDIR + "/fpkm/{run}/{cell}.tsv"
    log:
        OUTDIR + "/fpkm/{run}/{cell}.log"
    params:
        layout = lambda wildcards: get_layout(wildcards.cell)
    threads:
        THREADS
    shell:
        """
        nasctools CalculateFPKM \
            --threads {threads} \
            --layout {params.layout} \
            --strand TAG \
            --strand-tag ST \
            --nascent \
            --annotation {input.txt} \
            {input.bam} {input.bed} {output.txt} &> {log}
        """

rule report_gene_number:
    input:
        rules.calculate_fpkm.output.txt
    output:
        tsv = OUTDIR + "/genes/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/report_gene_number.py {input} > {output}
        """