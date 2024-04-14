#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
CELLS = list(DAT[DAT["Run"] == "GSE128273_NASCseq_K562"]["Cell"].values)
RS = ["1", "2"]
OUTDIR = "results/prepare"

### GSE128273

rule all:
    input:
        # expand(OUTDIR + "/download/GSE128273_NASCseq_K562/sra/{cell}.sra", cell=CELLS),
        expand(OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_{r}.fastq.gz", cell=CELLS, r=RS),

rule prefetch:
    output:
        sra = OUTDIR + "/download/GSE128273_NASCseq_K562/sra/{cell}.sra"
    log:
        OUTDIR + "/download/GSE128273_NASCseq_K562/sra/{cell}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch -o {output.sra} {wildcards.cell} &> {log}
        """

rule sra2fq:
    input:
        sra = rules.prefetch.output.sra
    output:
        tmp1 = OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq",
        tmp2 = OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq",
        fq1 = OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq.gz",
        fq2 = OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq.gz"
    log:
        OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}.log"
    conda:
        "sratools"
    threads:
        6
    shell:
        """(
        fasterq-dump --threads {threads} --split-3 --OUTDIR `dirname {output.fq1}` {input.sra}
        pigz -p {threads} -c {output.tmp1} > {output.fq1}
        pigz -p {threads} -c {output.tmp2} > {output.fq2} ) &> {log}
        """
