#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mismatch/ratio_consensus"
outdir = "results/signal2noise"
run_cells = run_cells_cellline

rule all:
    input:
        expand(outdir + "/pc/{run_cell}.tsv", run_cell=run_cells),


rule estimate_pc:
    input:
        tsv1 = indir + "/{run}/{cell}.tsv",
        tsv2 = lambda wildcards: get_estimate_pe_model(wildcards.cell),
        tsv3 = indir + "/{run}/{cell}.events.tsv",
    output:
        tsv = outdir + "/pc/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/signal2noise/estimate_pc.py {input} > {output}
        """
