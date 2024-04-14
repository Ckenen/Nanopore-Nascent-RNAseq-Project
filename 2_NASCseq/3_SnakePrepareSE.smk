#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
RUN_CELLS = RUN_CELLS_SE
INDIR = "data/media"
OUTDIR = "results/prepare"

# The FASTQ files of 50uM s4U treatment 3h ofK562 generated in NASC-seq study were obtained from the authors.
# The FASTQ files are single-end.

rule all:
    input:
        expand(OUTDIR + "/cutadapt/{run_cell}.fastq.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/bowtie2/{run_cell}.bam", run_cell=RUN_CELLS),

rule cutadapt:
    input:
        fq = INDIR + "/{cell}.fastq.gz"
    output:
        fq = OUTDIR + "/cutadapt/{run}/{cell}.fastq.gz"
    log:
        OUTDIR + "/cutadapt/{run}/{cell}.log"
    conda:
        "cutadapt"
    threads:
        THREADS
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g GCAGAGTACGGG \
            -a CTGTCTCTTATACAC -a CCCGTACTCTGC \
            -o {output.fq} {input.fq} &> {log}
        """

rule bowtie2:
    input:
        fq = rules.cutadapt.output.fq,
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
        bowtie2 -p {threads} --local --no-unal --un-gz {params.prefix}.fastq.gz \
            -x {input.idx}/ref -U {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """