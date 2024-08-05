#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mismatch/ratio_consensus"
OUTDIR = "results/signal2noise"
RUN_CELLS = RUN_CELLS_CELLLINE

rule all:
    input:
        expand(OUTDIR + "/pc/{run_cell}.tsv", run_cell=RUN_CELLS),


rule estimate_pc:
    input:
        tsv1 = lambda wildcards: get_estimate_pe_model(wildcards.cell),
        tsv2 = INDIR + "/{run}/{cell}.tsv",
        tsv3 = INDIR + "/{run}/{cell}.events.tsv"
    output:
        tsv = OUTDIR + "/pc/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/signal2noise/estimate_pc.py -m long -e {input.tsv3} {input.tsv1} {input.tsv2} > {output}
        """
