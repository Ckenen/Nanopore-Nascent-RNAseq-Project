#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/demux"

rule all:
    input:
        expand(outdir + "/barcodes/{run}.fa", run=runs),
        expand(outdir + "/fbilr/{run}.tsv.gz", run=runs),
        expand(outdir + "/fbilr/{run}.stats.tsv.gz", run=runs),
        expand(outdir + "/splitted/{run}", run=runs),
        expand(outdir + "/splitted/{run}.pdf", run=runs),
        outdir + "/splitted_cell_reads.all_runs.pdf",
        expand(outdir + "/trimmed/{run_cell}.fastq.gz", run_cell=run_cells),

rule get_barcodes:
    input:
        fa = config["barcodes"]
    output:
        fa = outdir + "/barcodes/{run}.fa",
        tsv = outdir + "/barcodes/{run}.tsv"
    run:
        import subprocess
        d = dat[dat["Run"] == wildcards.run]
        bcs = ["Bar%d" % bc for bc in sorted(set(d["Barcode"]))]
        cmd = "samtools faidx %s %s > %s" % (input.fa, " ".join(bcs), output.fa)
        subprocess.check_call(cmd, shell=True)
        with open(output.tsv, "w+") as fw:
            for cell, bc in d[["Cell", "Barcode"]].values:
                fw.write("%s\tBar%s\n" % (cell, bc))

rule fbilr:
    input:
        fq = indir + "/{run}.fastq.gz",
        fa = rules.get_barcodes.output.fa
    output:
        txt = outdir + "/fbilr/{run}.tsv.gz"
    log:
        outdir + "/fbilr/{run}.log"
    threads:
        24
    shell:
        """(
        fbilr -t {threads} -m PE -b {input.fa} {input.fq} \
            | gzip -c > {output.txt} ) &> {log}
        """

rule stat_matrix2:
    input:
        txt = rules.fbilr.output.txt
    output:
        txt = outdir + "/fbilr/{run}.stats.tsv.gz"
    shell:
        """
        ./scripts/demux/stat_matrix2.sh {input.txt} | gzip -c > {output.txt}
        """

rule split_reads:
    input:
        fq = indir + "/{run}.fastq.gz",
        mtx = rules.fbilr.output.txt,
        txt = rules.get_barcodes.output.tsv
    output:
        out = directory(outdir + "/splitted/{run}")
    log:
        outdir + "/splitted/{run}.log"
    threads:
        12
    shell:
        """
        ./scripts/demux/split_reads.py {input} {output} &> {log}
        pigz -p {threads} {output}/*/*.fastq
        """

rule plot_cell_reads: # optional
    input:
        txtdir = rules.split_reads.output.out
    output:
        pdf = outdir + "/splitted/{run}.pdf"
    shell:
        """
        ./scripts/demux/plot_cell_reads.py {input.txtdir}/reads.tsv {output.pdf}
        """

rule merge_pdfs: # optional
    input:
        pdfs = expand(rules.plot_cell_reads.output.pdf, run=runs)
    output:
        pdf = outdir + "/splitted_cell_reads.all_runs.pdf"
    shell:
        """
        merge_pdf.py `dirname {input.pdfs[0]}`/*.pdf {output.pdf}
        """
    
rule trim_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        fq = outdir + "/trimmed/{run}/{cell}.fastq.gz"
    log:
        outdir + "/trimmed/{run}/{cell}.log"
    params:
        fq = rules.split_reads.output.out + "/succeed/{cell}.fastq.gz"
    threads:
        4
    shell:
        """
        ./scripts/demux/trim_reads.py {params.fq} {output.fq} &> {log}
        """