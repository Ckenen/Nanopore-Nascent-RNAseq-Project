#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/marked_strand"
OUTDIR = "results/mismatch"

rule all:
    input:
        expand(OUTDIR + "/events/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/marked_nascent/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/ratio/{run_cell}.tsv", run_cell=RUN_CELLS),

rule get_events:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        bed = config["SNP_BED_GZ"]
    output:
        bam = OUTDIR + "/events/{run}/{cell}.bam"
    log:
        OUTDIR + "/events/{run}/{cell}.log"
    threads:
        THREADS
    shell:
        """
        nasctools GetEvent \
            --threads {threads} \
            --snp {input.bed} \
            {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_nascent:
    input:
        bam = rules.get_events.output.bam
    output:
        bam = OUTDIR + "/marked_nascent/{run}/{cell}.bam"
    log:
        OUTDIR + "/marked_nascent/{run}/{cell}.log"
    params:
        layout = lambda wildcards: get_layout(wildcards.cell)
    shell:
        """
        nasctools MarkNascent \
            --platform NGS \
            --layout {params.layout} \
            {input.bam} {output.bam} &> {log}
        samtools index {output.bam}
        """

rule report_ratio:
    input:
        bam = rules.mark_nascent.output.bam
    output:
        tsv = OUTDIR + "/ratio/{run}/{cell}.tsv"
    log:
        OUTDIR + "/ratio/{run}/{cell}.log"
    threads:
        THREADS
    shell:
        """
        nasctools ReportMismatch \
            --threads {threads} \
            --strand TAG \
            --strand-tag ST \
            {input.bam} {output.tsv} &> {log}
        """
