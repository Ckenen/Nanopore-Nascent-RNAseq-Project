#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "results/prepare/bowtie2"
outdir = "results/mapping"

rule all:
    input:
        expand(outdir + "/star/{sample}", sample=samples),
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/infer_experiment/{sample}.txt", sample=samples),
        expand(outdir + "/marked_duplicates/{sample}.bam", sample=samples),
rule star_mapping:
    input:
        fq1 = indir + "/{sample}.1.fastq.gz",
        fq2 = indir + "/{sample}.2.fastq.gz",
        idx = config["star_idx"]
    output:
        out = directory(outdir + "/star/{sample}")
    log:
        outdir + "/star/{sample}.log"
    threads:
        12
    shell:
        """
        ./scripts/mapping/star.pe.sh {input} {threads} {output} &> {log}
        """


rule filter_bam:
    input:
        bamdir = rules.star_mapping.output
    output:
        bam = outdir + "/filtered/{sample}.bam"
    shell:
        """
        samtools view -@ {threads} -q 30 -d "NH:1" -f 2 -F 2308 --expr 'rname =~ "^chr([0-9]+|[XY])$"' -o {output.bam} {input.bamdir}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """


rule infer_experiment:
    input:
        bam = outdir + "/filtered/{sample}.bam",
        bed = config["gene_bed"]
    output:
        txt = outdir + "/infer_experiment/{sample}.txt"
    shell:
        """
        infer_experiment.py -s 2000000 -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        """

rule mark_duplicates:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/marked_duplicates/{sample}.bam",
        txt = outdir + "/marked_duplicates/{sample}_metrics.txt"
    log:
        outdir + "/marked_duplicates/{sample}.log"
    threads:
        4
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.bam} -M {output.txt} -O {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """