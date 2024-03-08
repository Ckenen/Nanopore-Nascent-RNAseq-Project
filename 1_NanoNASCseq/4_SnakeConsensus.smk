#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_duplicate"
outdir = "results/consensus"
d = dat[dat["Run"] == "20220719_K562R3"]
run_cells = list(d["RunCell"])

rule all:
    input:
        expand(outdir + "/consensus_accuracy/{run_cell}", run_cell=run_cells),

rule estimate_accuracy:
    input:
        bam = indir + "/{run}/{cell}.bam",
        mmi = lambda wildcards: FILES[get_species(wildcards.cell)]["GENOME_SPLICE_MMI"],
        bed = lambda wildcards: FILES[get_species(wildcards.cell)]["TRANSCRIPT_BED"]
    output:
        out = directory(outdir + "/consensus_accuracy/{run}/{cell}")
    log:
        outdir + "/consensus_accuracy/{run}/{cell}.log"
    threads:
        24
    shell:
        """
        ./scripts/estimate_accuracy.py {input.bam} {input.mmi} {input.bed} {threads} {output} &> {log}
        """

