#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
rs = ["1", "2"]
indir = "data/datasets"
outdir = "results/prepare"

BOWTIE2_RRNA_INDEX = "/home/chenzonggui/species/homo_sapiens/ncbi/bt2_rrna_index"

rule all:
    input:
        expand(outdir + "/fastqc/{sample}_R{r}_fastqc.html", sample=samples, r=rs),
        # expand(outdir + "/cutadapt/{sample}_R{r}.fastq.gz", sample=samples, r=rs),
        expand(outdir + "/bowtie2/{sample}.bam", sample=samples, r=rs)


rule fastqc:
    input:
        fq = indir + "/{name}.fastq.gz"
    output:
        html = outdir + "/fastqc/{name}_fastqc.html"
    log:
        outdir + "/fastqc/{name}_fastqc.log"
    shell:
        """
        fastqc -o `dirname {output.html}` {input.fq} &> {log}
        """
    
rule cutadapt:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/{sample}_R1.fastq.gz",
        fq2 = outdir + "/cutadapt/{sample}_R2.fastq.gz"
    log:
        outdir + "/cutadapt/{sample}.log"
    threads:
        6
    shell:
        """
        cutadapt --max-n 2 -m 20 -q 30 -j {threads} \
            -a AGATCGGAAGAGCACACGTC -a GATCGGAAGAGCACACGTCT \
            -A AGATCGGAAGAGCGTCGTGT -A GATCGGAAGAGCGTCGTGTA \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = BOWTIE2_RRNA_INDEX
    output:
        bam = outdir + "/bowtie2/{sample}.bam"
    log:
        outdir + "/bowtie2/{sample}.log"
    params:
        prefix = outdir + "/bowtie2/{sample}"
    threads:
        8
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal --un-conc {params.prefix}.fq \
            -x {input.idx}/ref -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -b -u - \
            | samtools sort -@ {threads} - > {output.bam}
        pigz -p {threads} {params.prefix}.1.fq {params.prefix}.2.fq
        samtools index {output.bam} ) &> {log}
        """