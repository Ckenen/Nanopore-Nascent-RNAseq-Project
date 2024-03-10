#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_duplicate"
outdir = "results/mismatch"
# run_cells = run_cells[:2]

rule all:
    input:
        expand(outdir + "/events/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/ratio_all/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/ratio_rmdup/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/ratio_consensus/{run_cell}.tsv", run_cell=run_cells),

rule get_events:
    input:
        bam = indir + "/{run}/{cell}.bam",
        bed = lambda wildcards: get_snp_bed(wildcards.cell)
    output:
        bam = outdir + "/events/{run}/{cell}.bam"
    log:
        outdir + "/events/{run}/{cell}.log"
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
        tsv = outdir + "/ratio_all/{run}/{cell}.tsv"
    log:
        outdir + "/ratio_all/{run}/{cell}.log"
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
        tsv = outdir + "/ratio_rmdup/{run}/{cell}.tsv"
    log:
        outdir + "/ratio_rmdup/{run}/{cell}.log"
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
        tsv = outdir + "/ratio_consensus/{run}/{cell}.tsv",
        tsv2 = outdir + "/ratio_consensus/{run}/{cell}.events.tsv"
    log:
        outdir + "/ratio_consensus/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_consensus_mismatch.py {input} {output} &> {log}
        """

