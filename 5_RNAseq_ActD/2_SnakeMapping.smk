#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/prepare/bowtie2"
outdir = "results/mapping"

rule all:
    input:
        outdir + "/star/index",
        expand(outdir + "/star/mapped/{sample}", sample=samples),
        expand(outdir + "/filtered/{sample}.{species}.bam", sample=samples, species=species_list),
        expand(outdir + "/filtered/{sample}.{species}.flagstat", sample=samples, species=species_list),
        expand(outdir + "/infer_experiment/{sample}.{species}.txt", sample=samples, species=species_list),
        expand(outdir + "/rmdup/{sample}.{species}.bam", sample=samples, species=species_list),
        expand(outdir + "/rmdup/{sample}.{species}.flagstat", sample=samples, species=species_list),

rule star_index:
    input:
        fa1 = FILES["human"]["GENOME_FASTA"],
        fa2 = FILES["fly"]["GENOME_FASTA"],
        gtf1 = FILES["human"]["ANNOTATION_GTF"],
        gtf2 = FILES["fly"]["ANNOTATION_GTF"]
    output:
        fa = outdir + "/star/index.integrated.fa",
        gtf = outdir + "/star/index.integrated.gtf",
        idx = directory(outdir + "/star/index")
    log:
        outdir + "/star/index.log"
    threads:
        24
    shell:
        """(
        cat {input.fa1} {input.fa2} > {output.fa}
        samtools faidx {output.fa}
        cat {input.gtf1} {input.gtf2} | grep -v '#' | sort -k1,1 -k4,4n -k5,5n > {output.gtf}
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.idx} \
            --genomeFastaFiles {output.fa} --sjdbGTFfile {output.gtf} ) &> {log}
        """

# STAR --genomeLoad Remove --genomeDir results/mapping/star/index

rule star_mapping:
    input:
        fq1 = indir + "/{sample}.1.fq.gz",
        fq2 = indir + "/{sample}.2.fq.gz",
        index = rules.star_index.output.idx
    output:
        out = directory(outdir + "/star/mapped/{sample}")
    log:
        outdir + "/star/mapped/{sample}.log"
    threads:
        12
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --genomeDir {input.index} \
            --genomeLoad LoadAndKeep \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule filter_and_split: # filter and split
    input:
        rules.star_mapping.output.out
    output:
        bam1 = outdir + "/filtered/{sample}.human.bam",
        bam2 = outdir + "/filtered/{sample}.fly.bam"
    log:
        outdir + "/filtered/{sample}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" -f 2 -F 2308 --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -o {output.bam1} {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools view -@ {threads} -q 30 -d "NH:1" -f 2 -F 2308 --expr 'rname =~ "^(2L|2R|3L|3R|[4XY])$"' \
            -o {output.bam2} {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index {output.bam1}
        samtools index {output.bam2} ) &> {log}
        """

rule infer_experiment:
    input:
        bam = outdir + "/filtered/{sample}.{species}.bam",
        bed = lambda wildcards: FILES[wildcards.species]["GENE_BED"]
    output:
        txt = outdir + "/infer_experiment/{sample}.{species}.txt"
    shell:
        """
        infer_experiment.py -s 2000000 -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        """

rule rmdup:
    input:
        bam = outdir + "/filtered/{sample}.{species}.bam"
    output:
        bam = outdir + "/rmdup/{sample}.{species}.bam",
        txt = outdir + "/rmdup/{sample}.{species}_metrics.txt"
    log:
        outdir + "/rmdup/{sample}.{species}.log"
    shell:
        """(
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.bam} -O {output.bam} -M {output.txt} 
        samtools index {output.bam} ) &> {log}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """