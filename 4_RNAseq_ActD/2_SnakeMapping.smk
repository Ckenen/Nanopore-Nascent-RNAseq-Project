#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/prepare/bowtie2"
OUTDIR = "results/mapping"

rule all:
    input:
        OUTDIR + "/merged/human_fly.fa",
        OUTDIR + "/merged/human_fly.gtf",
        OUTDIR + "/star/index",
        expand(OUTDIR + "/star/mapped/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.{species}.bam", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/filtered/{sample}.{species}.flagstat", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/infer_experiment/{sample}.{species}.txt", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/rmdup/{sample}.{species}.bam", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/rmdup/{sample}.{species}.flagstat", sample=SAMPLES, species=SPECIES),

rule merge_fasta:
    input:
        fa1 = get_fasta("human"),
        fa2 = get_fasta("fly")
    output:
        fa = OUTDIR + "/merged/human_fly.fa"
    conda:
        "star"
    log:
       OUTDIR + "/merged/human_fly.fa.log" 
    shell:
        """(
        cat {input.fa1} {input.fa2} > {output.fa}
        samtools faidx {output.fa} ) &> {log}
        """

rule merge_gtf:
    input:
        gtf1 = get_gtf("human"),
        gtf2 = get_gtf("fly")
    output:
        gtf = OUTDIR + "/merged/human_fly.gtf",
        gtf2 = OUTDIR + "/merged/human_fly.gtf.gz"
    conda:
        "star"
    log:
        OUTDIR + "/merged/human_fly.gtf.log" 
    shell:
        """(
        cat {input.gtf1} {input.gtf2} | grep -v '#' \
            | sort -k1,1 -k4,4n -k5,5n > {output.gtf} 
        bgzip -c {output.gtf} > {output.gtf2}
        tabix -p gff {output.gtf2} ) &> {log}
        """

rule star_index:
    input:
        fa = rules.merge_fasta.output.fa,
        gtf = rules.merge_gtf.output.gtf
    output:
        idx = directory(OUTDIR + "/star/index")
    log:
        OUTDIR + "/star/index.log"
    conda:
        "star"
    threads:
        THREADS
    shell:
        """(
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.idx} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} ) &> {log}
        """

rule star_mapping:
    input:
        fq1 = INDIR + "/{sample}.1.fq.gz",
        fq2 = INDIR + "/{sample}.2.fq.gz",
        index = rules.star_index.output.idx
    output:
        out = directory(OUTDIR + "/star/mapped/{sample}")
    log:
        OUTDIR + "/star/mapped/{sample}.log"
    conda:
        "star"
    threads:
        THREADS
    shell:
        """(
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --genomeDir {input.index} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} ) &> {log}
        """

rule filter_and_split: # filter and split
    input:
        rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/filtered/{sample}.{species}.bam"
    log:
        OUTDIR + "/filtered/{sample}.{species}.log"
    conda:
        "star"
    params:
        pattern = lambda wildcards: get_seqname_pattern(wildcards.species)
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} \
            -q 30 \
            -d "NH:1" \
            -f 2 \
            -F 2308 \
            --expr 'rname =~ "{params.pattern}"' \
            -o {output.bam} \
            {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule infer_experiment:
    input:
        bam = rules.filter_and_split.output.bam,
        bed = lambda wildcards: get_gene_bed(wildcards.species)
    output:
        txt = OUTDIR + "/infer_experiment/{sample}.{species}.txt"
    log:
        OUTDIR + "/infer_experiment/{sample}.{species}.log"
    conda:
        "rseqc"
    shell:
        """
        infer_experiment.py \
            -s 2000000 \
            -i {input.bam} \
            -r {input.bed} > {output.txt} 2> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_and_split.output.bam
    output:
        bam = OUTDIR + "/rmdup/{sample}.{species}.bam",
        txt = OUTDIR + "/rmdup/{sample}.{species}_metrics.txt"
    log:
        OUTDIR + "/rmdup/{sample}.{species}.log"
    conda:
        "picard"
    shell:
        """(
        picard MarkDuplicates \
            --REMOVE_DUPLICATES true \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.txt} 
        samtools index {output.bam} ) &> {log}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    conda:
        "star"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """