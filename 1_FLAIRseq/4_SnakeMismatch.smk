#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/mark_duplicate"
OUTDIR = "results/mismatch"

rule all:
    input:
        expand(OUTDIR + "/events/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/ratio_all/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/ratio_rmdup/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/ratio_consensus/{run_cell}.tsv", run_cell=RUN_CELLS),

rule get_events:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        bed = lambda wildcards: get_snp_bed(wildcards.cell)
    output:
        bam = OUTDIR + "/events/{run}/{cell}.bam"
    log:
        OUTDIR + "/events/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools GetEvent -t {threads} -s {input.bed} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule report_mismatch_all:
    input:
        bam = rules.get_events.output.bam
    output:
        tsv = OUTDIR + "/ratio_all/{run}/{cell}.tsv"
    log:
        OUTDIR + "/ratio_all/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch -t {threads} -s F {input.bam} {output.tsv} &> {log}
        """

rule report_mismatch_uniq:
    input:
        bam = rules.get_events.output.bam
    output:
        tsv = OUTDIR + "/ratio_rmdup/{run}/{cell}.tsv"
    log:
        OUTDIR + "/ratio_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch -t {threads} -s F --discard-duplicates {input.bam} {output.tsv} &> {log}
        """

rule report_mismatch_consensus:
    input:
        bam = rules.get_events.output.bam,
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        tsv = OUTDIR + "/ratio_consensus/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/ratio_consensus/{run}/{cell}.events.tsv"
    log:
        OUTDIR + "/ratio_consensus/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_consensus_mismatch.py {input} {output} &> {log}
        """

