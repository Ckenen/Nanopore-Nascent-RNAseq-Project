#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/mark_duplicate"
OUTDIR = "results/consensus"
RUN_CELLS = list(DAT[DAT["Run"] == "20220719_K562R3"]["RunCell"])

rule all:
    input:
        expand(OUTDIR + "/consensus_accuracy/{run_cell}", run_cell=RUN_CELLS),

rule estimate_accuracy:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        mmi = lambda wildcards: get_genome_splice_mmi(wildcards.cell),
        bed = lambda wildcards: get_transcript_bed(wildcards.cell)
    output:
        out = directory(OUTDIR + "/consensus_accuracy/{run}/{cell}")
    log:
        OUTDIR + "/consensus_accuracy/{run}/{cell}.log"
    threads:
        24
    shell:
        """
        ./scripts/consensus/estimate_accuracy.py {input.bam} {input.mmi} {input.bed} {threads} {output} &> {log}
        """

