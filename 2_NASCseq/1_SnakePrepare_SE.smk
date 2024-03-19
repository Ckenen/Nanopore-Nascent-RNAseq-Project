#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
run_cells = run_cells_se
indir = "data/media"
outdir = "results/prepare"
# cells = cells[:1]

rule all:
    input:
        expand(outdir + "/cutadapt/{run_cell}.fastq.gz", run_cell=run_cells),
        expand(outdir + "/bowtie2/{run_cell}.bam", run_cell=run_cells),

rule cutadapt:
    input:
        fq = indir + "/{cell}.fastq.gz"
    output:
        fq = outdir + "/cutadapt/{run}/{cell}.fastq.gz"
    log:
        outdir + "/cutadapt/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g GCAGAGTACGGG \
            -a CTGTCTCTTATACAC -a CCCGTACTCTGC \
            -o {output.fq} {input.fq} &> {log}
        """

rule bowtie2: # remove rRNA
    input:
        fq = rules.cutadapt.output.fq,
        idx = BOWTIE2_RRNA_INDEX
    output:
        bam = outdir + "/bowtie2/{run}/{cell}.bam"
    log:
        outdir + "/bowtie2/{run}/{cell}.log"
    params:
        prefix = outdir + "/bowtie2/{run}/{cell}"
    threads:
        8
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal --un-gz {params.prefix}.fastq.gz -x {input.idx}/ref -U {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """