#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "data/datasets"
OUTDIR = "results/demux"

rule all:
    input:
        expand(OUTDIR + "/barcodes/{run}.fa", run=RUNS),
        expand(OUTDIR + "/fbilr/{run}.tsv.gz", run=RUNS),
        expand(OUTDIR + "/fbilr/{run}.stats.tsv.gz", run=RUNS),
        expand(OUTDIR + "/splitted/{run}", run=RUNS),
        expand(OUTDIR + "/trimmed/{run_cell}", run_cell=RUN_CELLS),
        OUTDIR + "/summary_of_trimming.tsv",

rule get_barcodes:
    input:
        fa = config["BARCODES"]
    output:
        fa = OUTDIR + "/barcodes/{run}.fa",
        tsv = OUTDIR + "/barcodes/{run}.tsv"
    run:
        import subprocess
        d = DAT[DAT["Run"] == wildcards.run]
        assert len(d) > 0
        bcs = ["Bar%d" % bc for bc in sorted(set(d["Barcode"]))]
        cmd = "samtools faidx %s %s > %s" % (input.fa, " ".join(bcs), output.fa)
        subprocess.check_call(cmd, shell=True)
        with open(output.tsv, "w+") as fw:
            for cell, bc in d[["Cell", "Barcode"]].values:
                fw.write("%s\tBar%s\n" % (cell, bc))

rule fbilr:
    input:
        fq = INDIR + "/{run}.fastq.gz",
        fa = rules.get_barcodes.output.fa
    output:
        txt = OUTDIR + "/fbilr/{run}.tsv.gz"
    log:
        OUTDIR + "/fbilr/{run}.log"
    threads:
        THREADS
    shell:
        """(
        fbilr -t {threads} -m PE -b {input.fa} {input.fq} \
            | gzip -c > {output.txt} ) &> {log}
        """

rule stat_matrix2:
    input:
        txt = rules.fbilr.output.txt
    output:
        txt = OUTDIR + "/fbilr/{run}.stats.tsv.gz"
    shell:
        """
        ./scripts/demux/stat_matrix2.sh {input.txt} | gzip -c > {output.txt}
        """

rule split_reads:
    input:
        fq = INDIR + "/{run}.fastq.gz",
        mtx = rules.fbilr.output.txt,
        txt = rules.get_barcodes.output.tsv
    output:
        out = directory(OUTDIR + "/splitted/{run}")
    log:
        OUTDIR + "/splitted/{run}.log"
    threads:
        THREADS
    shell:
        """
        ./scripts/demux/split_reads.py {input} {output} &> {log}
        pigz -p {threads} {output}/*/*.fastq
        """
    
rule trim_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        out = directory(OUTDIR + "/trimmed/{run}/{cell}")
    log:
        OUTDIR + "/trimmed/{run}/{cell}.log"
    params:
        fq = rules.split_reads.output.out + "/succeed/{cell}.fastq.gz"
    shell:
        """
        ./scripts/demux/trim_reads.py {params.fq} {output.out} &> {log}
        gzip {output.out}/trimmed.fastq
        """

rule summary_of_trimming:
    input:
        expand(OUTDIR + "/trimmed/{run_cell}", run_cell=RUN_CELLS)
    output:
        tsv = OUTDIR + "/summary_of_trimming.tsv"
    run:
        with open(output.tsv, "w+") as fw:
            for i, f in enumerate(sorted(input)):
                lines = open(f + "/stats.tsv").readlines()
                if i == 0:
                    fw.write(lines[0])
                fw.write(lines[1])
                    
