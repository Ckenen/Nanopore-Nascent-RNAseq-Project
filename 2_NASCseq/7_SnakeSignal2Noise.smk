#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mismatch/marked_nascent"
OUTDIR = "results/signal2noise"

rule all:
    input:
        expand(OUTDIR + "/pc/{run_cell}.tsv", run_cell=RUN_CELLS),

rule estimate_pc:
    input:
        tsv1 = "reports/Estimate.Pe.Model.K562.PE.tsv",
        tsv2 = "results/mismatch/ratio/{run}/{cell}.tsv",
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        txt = OUTDIR + "/pc/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/estimate_pc.py -m short -b {input.bam} \
            {input.tsv1} {input.tsv2} > {output}
        """
