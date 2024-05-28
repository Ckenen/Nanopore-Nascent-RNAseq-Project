#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
RUN_CELLS = RUN_CELLS
INDIR = "results/prepare/bowtie2"
OUTDIR = "results/mapping"

rule all:
    input:
        expand(OUTDIR + "/star/{run_cell}", run_cell=RUN_CELLS[:1]),
        # expand(OUTDIR + "/filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/filtered/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/rmdup/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/rmdup/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/marked_strand/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/marked_strand/{run_cell}.flagstat", run_cell=RUN_CELLS),

def get_input_fastqs(wildcards):
    run, cell = wildcards.run, wildcards.cell
    layout = get_layout(cell)
    if layout == "PE":
        paths = [
            INDIR + "/%s/%s.1.fastq.gz" % (run, cell),
            INDIR + "/%s/%s.2.fastq.gz" % (run, cell)]
    else:
        paths = [INDIR + "/%s/%s.fastq.gz" % (run, cell)]
    return paths

rule star_mapping:
    input:
        fqs = lambda wildcards: get_input_fastqs(wildcards),
        idx = config["GENOME_STAR_INDEX"]
    output:
        out = directory(OUTDIR + "/star/{run}/{cell}")
    conda:
        "star"
    log:
        OUTDIR + "/star/{run}/{cell}.log"
    params:
        prefix = OUTDIR + "/star/{run}/{cell}/{cell}"
    threads:
        THREADS
    shell:
        """(
        mkdir -p {output.out}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fqs} ) &> {log}
        """

def get_filter_flags(wildcards):
    if get_layout(wildcards.cell) == "PE":
        return "-f 2 -F 2308"
    else:
        return "-F 2308"

rule filter_bam:
    input:
        rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/filtered/{run}/{cell}.bam"
    log:
        OUTDIR + "/filtered/{run}/{cell}.log"
    params:
        flag = lambda wildcards: get_filter_flags(wildcards)
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            {params.flag} -o {output.bam} {input}/{wildcards.cell}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/rmdup/{run}/{cell}.bam"
    log:
        OUTDIR + "/rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -t {threads} -r {input.bam} {output.bam} &> {log}
        """

rule mark_strand:
    input:
        bam = rules.rmdup.output.bam,
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bam = OUTDIR + "/marked_strand/{run}/{cell}.bam",
        txt = OUTDIR + "/marked_strand/{run}/{cell}.tsv"
    params:
        layout = lambda wildcards: get_layout(wildcards.cell)
    log:
        OUTDIR + "/marked_strand/{run}/{cell}.log"
    shell:
        """
        nasctools MarkStrand --gene {input.bed} --tag ST \
            --layout {params.layout} --strand U \
            --summary {output.txt} {input.bam} {output.bam} &> {log}
        samtools index {output.bam}
        """

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
