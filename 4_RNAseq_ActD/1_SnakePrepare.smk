#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
RS = ["1", "2"]
INDIR = "data/datasets"
OUTDIR = "results/prepare"

BOWTIE2_RRNA_INDEX = "/home/chenzonggui/species/homo_sapiens/ncbi/bt2_rrna_index"

rule all:
    input:
        expand(OUTDIR + "/fastqc/{sample}_R{r}_fastqc.html", sample=SAMPLES, r=RS),
        expand(OUTDIR + "/cutadapt/{sample}_R{r}.fastq.gz", sample=SAMPLES, r=RS),
        OUTDIR + "/ribo_bt2_index",
        expand(OUTDIR + "/bowtie2/{sample}.bam", sample=SAMPLES, r=RS)


rule fastqc:
    input:
        fq = INDIR + "/{name}.fastq.gz"
    output:
        html = OUTDIR + "/fastqc/{name}_fastqc.html"
    log:
        OUTDIR + "/fastqc/{name}_fastqc.log"
    conda:
        "fastqc"
    shell:
        """
        fastqc -o `dirname {output.html}` {input.fq} &> {log}
        """
    
rule cutadapt:
    input:
        fq1 = INDIR + "/{sample}_R1.fastq.gz",
        fq2 = INDIR + "/{sample}_R2.fastq.gz"
    output:
        fq1 = OUTDIR + "/cutadapt/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/cutadapt/{sample}_R2.fastq.gz"
    conda:
        "cutadapt"
    log:
        OUTDIR + "/cutadapt/{sample}.log"
    threads:
        THREADS
    shell:
        """
        cutadapt --max-n 2 -m 20 -q 30 \
            -j {threads} \
            -a AGATCGGAAGAGCACACGTC \
            -a GATCGGAAGAGCACACGTCT \
            -A AGATCGGAAGAGCGTCGTGT \
            -A GATCGGAAGAGCGTCGTGTA \
            -o {output.fq1} \
            -p {output.fq2} \
            {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2_index:
    input:
        fa = config["RIBO_FASTA"]
    output:
        out = directory(OUTDIR + "/ribo_bt2_index")
    log:
        OUTDIR + "/ribo_bt2_index.log"
    conda:
        "bowtie2"
    threads:
        THREADS
    shell:
        """(
        mkdir -p {output}
        bowtie2-build --threads {threads} \
            {input.fa} {output}/ref ) &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = rules.bowtie2_index.output.out
    output:
        bam = OUTDIR + "/bowtie2/{sample}.bam"
    conda:
        "bowtie2"
    log:
        OUTDIR + "/bowtie2/{sample}.log"
    params:
        prefix = OUTDIR + "/bowtie2/{sample}"
    threads:
        THREADS
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal \
            --un-conc {params.prefix}.fq \
            -x {input.idx}/ref \
            -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -b -u - \
            | samtools sort -@ {threads} - > {output.bam}
        pigz -p {threads} {params.prefix}.1.fq {params.prefix}.2.fq
        samtools index {output.bam} ) &> {log}
        """