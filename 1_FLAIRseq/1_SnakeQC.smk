#!/usr/bin/env RUNSnakemake
include: "0_SnakeCommon.smk"
INDIR = "data/datasets"
OUTDIR = "results/qc"

rule all:
    input:
        expand(OUTDIR + "/check_barcode_list/{run}.txt", run=RUNS),
        expand(OUTDIR + "/read_length/{run}.tsv", run=RUNS),
        expand(OUTDIR + "/read_length/{run}.pdf", run=RUNS),
        expand(OUTDIR + "/library_summary/{run}.tsv", run=RUNS),

# Ensure the barcodes in the reads are consistent with what is expected.

def get_run_barcode_list(run):
    tmp = DAT[DAT["Run"] == run]
    assert len(tmp) > 0
    return ["Bar%d" % x for x in tmp["Barcode"]]

rule check_barcode_list:
    input:
        fq = INDIR + "/{run}.fastq.gz",
        fa = config["BARCODES"]
    output:
        txt = OUTDIR + "/check_barcode_list/{run}.txt"
    log:
        OUTDIR + "/check_barcode_list/{run}.log"
    params:
        barcodes = lambda wildcards: ",".join(get_run_barcode_list(wildcards.run))
    threads:
        4
    shell:
        """
        ./scripts/qc/check_barcode_list.py {input.fq} {input.fa} \
            {threads} {params.barcodes} {output.txt} &> {log}
        """

# Analysis of the length of raw reads, including cDNA sequences and exogenous sequences.

rule stat_read_length:
    input:
        fq = INDIR + "/{run}.fastq.gz"
    output:
        txt = OUTDIR + "/read_length/{run}.tsv"
    shell:
        """
        ./scripts/qc/stat_read_length.sh {input.fq} > {output.txt}
        """

rule plot_read_length:
    input:
        txt = rules.stat_read_length.output.txt
    output:
        pdf = OUTDIR + "/read_length/{run}.pdf"
    shell:
        """
        ./scripts/qc/plot_read_length.py {input.txt} {output.pdf}
        """

# Report several length metrics of reads.

rule report_library_summary:
    input:
        txt = rules.stat_read_length.output.txt
    output:
        txt = OUTDIR + "/library_summary/{run}.tsv"
    shell:
        """
        ./scripts/qc/report_library_summary.py {input.txt} > {output.txt}
        """
