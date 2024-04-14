#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
RUN_CELLS = RUN_CELLS_PE
RS = ["1", "2"]
OUTDIR = "results/prepare"

rule all:
    input:
        expand(OUTDIR + "/cutadapt/{run_cell}_{r}.fastq.gz", run_cell=RUN_CELLS, r=RS),
        expand(OUTDIR + "/bowtie2/{run_cell}.bam", run_cell=RUN_CELLS),

def get_input_fastqs(run):
    if run == "GSE128273_NASCseq_K562":
        paths = [
            OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq.gz",
            OUTDIR + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq.gz"]
    else:
        paths = [
            "data/datasets/{cell}_1.fastq.gz",
            "data/datasets/{cell}_2.fastq.gz"]
    return paths

def get_cutadapt_linker_params(run):
    if run == "GSE128273_NASCseq_K562":
        return "-g GTGTATAAGAGACAG -g GCAGAGTACGGG -a CTGTCTCTTATACAC -a CCCGTACTCTGC -G GTGTATAAGAGACAG -G GCAGAGTACGGG -A CTGTCTCTTATACAC -A CCCGTACTCTGC"
    else:
        return "-g GTGTATAAGAGACAG -g ATCAACGCAGAGTAC -a CTGTCTCTTATACAC -a GTACTCTGCGTTGAT -G GTGTATAAGAGACAG -G ATCAACGCAGAGTAC -A CTGTCTCTTATACAC -A GTACTCTGCGTTGAT"

rule cutadapt:
    input:
        fqs = lambda wildcards: get_input_fastqs(wildcards.run)
    output:
        fq1 = OUTDIR + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = OUTDIR + "/cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        OUTDIR + "/cutadapt/{run}/{cell}.log"
    conda:
        "cutadapt"
    params:
        linker = lambda wildcards: get_cutadapt_linker_params(wildcards.run)
    threads:
        THREADS
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 {params.linker} \
            -o {output.fq1} -p {output.fq2} {input.fqs} &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = config["RIBO_BOWTIE2_INDEX"]
    output:
        bam = OUTDIR + "/bowtie2/{run}/{cell}.bam"
    log:
        OUTDIR + "/bowtie2/{run}/{cell}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/bowtie2/{run}/{cell}"
    threads:
        THREADS
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal \
            --un-conc-gz {params.prefix}.fastq.gz \
            -x {input.idx}/ref -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam}
        mv {params.prefix}.fastq.1.gz {params.prefix}.1.fastq.gz
        mv {params.prefix}.fastq.2.gz {params.prefix}.2.fastq.gz ) &> {log}
        """