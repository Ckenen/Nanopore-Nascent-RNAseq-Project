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
        expand(outdir + "/trimmed/{run_cell}", run_cell=run_cells),
        # outdir + "/summary_of_trimming.tsv",

rule get_barcodes:
    input:
        fa = BARCODE_FASTA
    output:
        fa = outdir + "/barcodes/{run}.fa",
        tsv = outdir + "/barcodes/{run}.tsv"
    run:
        import subprocess
        d = dat[dat["Run"] == wildcards.run]
        assert len(d) > 0
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
        8
    shell:
        """
        ./scripts/demux/split_reads.py {input} {output} &> {log}
        pigz -p {threads} {output}/*/*.fastq
        """
    
rule trim_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        out = directory(outdir + "/trimmed/{run}/{cell}")
    log:
        outdir + "/trimmed/{run}/{cell}.log"
    params:
        fq = rules.split_reads.output.out + "/succeed/{cell}.fastq.gz"
    threads:
        4
    shell:
        """
        ./scripts/demux/trim_reads.py {params.fq} {output.out} &> {log}
        pigz -p {threads} {output.out}/trimmed.fastq
        """

rule summary_of_trimming:
    input:
        expand(outdir + "/trimmed/{run_cell}", run_cell=run_cells)
    output:
        tsv = outdir + "/summary_of_trimming.tsv"
    run:
        with open(output.tsv, "w+") as fw:
            for i, f in enumerate(sorted(input)):
                lines = open(f + "/stats.tsv").readlines()
                if i == 0:
                    fw.write(lines[0])
                fw.write(lines[1])
                    
