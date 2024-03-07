#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/qc"

rule all:
    input:
        expand(outdir + "/check_barcode_list/{run}.txt", run=runs),
        expand(outdir + "/read_length/{run}.tsv", run=runs),
        expand(outdir + "/read_length/{run}.pdf", run=runs),
        outdir + "/read_length.all_runs.pdf",
        expand(outdir + "/library_summary/{run}.tsv", run=runs),
        outdir + "/library_summary.all_runs.tsv",

rule check_barcode_list:
    input:
        fq = indir + "/{run}.fastq.gz",
        fa = config["barcodes"]
    output:
        txt = outdir + "/check_barcode_list/{run}.txt"
    params:
        barcodes = lambda wildcards: ",".join(get_barcode_list(wildcards.run))
    threads:
        8
    shell:
        """
        ./scripts/qc/check_barcode_list.py {input.fq} {input.fa} \
            {threads} {params.barcodes} {output.txt}
        """

rule stat_read_length:
    input:
        fq = indir + "/{run}.fastq.gz"
    output:
        txt = outdir + "/read_length/{run}.tsv"
    shell:
        """
        ./scripts/qc/stat_read_length.sh {input.fq} > {output.txt}
        """

rule plot_read_length:
    input:
        txt = rules.stat_read_length.output.txt
    output:
        pdf = outdir + "/read_length/{run}.pdf"
    shell:
        """
        ./scripts/qc/plot_read_length.py {input.txt} {output.pdf}
        """

rule merge_pdf:
    input:
        pdfs = expand(outdir + "/read_length/{run}.pdf", run=runs)
    output:
        pdf = outdir + "/read_length.all_runs.pdf"
    shell:
        """
        merge_pdf.py `dirname {input.pdfs[0]}`/*.pdf {output.pdf}
        """

rule report_library_summary:
    input:
        txt = rules.stat_read_length.output.txt
    output:
        txt = outdir + "/library_summary/{run}.tsv"
    shell:
        """
        ./scripts/qc/report_library_summary.py {input.txt} > {output.txt}
        """

rule merge_library_summary:
    input:
        txts = expand(rules.report_library_summary.output.txt, run=runs)
    output:
        txt = outdir + "/library_summary.all_runs.tsv"
    shell:
        """
        head -n 1 {input.txts[0]} > {output.txt}
        for f in {input.txts}; do tail -n 1 $f >> {output.txt}; done
        """