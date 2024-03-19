#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
run_cells = run_cells_pe
outdir = "results/prepare"

rule all:
    input:
        # expand(outdir + "/download/GSE128273_NASCseq_K562/sra/{cell}.sra", cell=cells_gse127273),
        expand(outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_{r}.fastq.gz", cell=cells_gse127273, r=["1", "2"]),
        expand(outdir + "/cutadapt/{run_cell}_{r}.fastq.gz", run_cell=run_cells, r=["1", "2"]),
        expand(outdir + "/bowtie2/{run_cell}.bam", run_cell=run_cells),

### GSE128273

rule prefetch:
    output:
        sra = outdir + "/download/GSE128273_NASCseq_K562/sra/{cell}.sra"
    log:
        outdir + "/download/GSE128273_NASCseq_K562/sra/{cell}.log"
    shell:
        """
        prefetch -o {output.sra} {wildcards.cell} &> {log}
        """

rule sra2fq:
    input:
        sra = rules.prefetch.output.sra
    output:
        tmp1 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq",
        tmp2 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq",
        fq1 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq.gz",
        fq2 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq.gz"
    log:
        outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}.log"
    shell:
        """(
        fasterq-dump --threads {threads} --split-3 --outdir `dirname {output.fq1}` {input.sra}
        pigz -p {threads} -c {output.tmp1} > {output.fq1}
        pigz -p {threads} -c {output.tmp2} > {output.fq2} ) &> {log}
        """

rule cutadapt_GSE128273:
    input:
        fq1 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_1.fastq.gz",
        fq2 = outdir + "/download/GSE128273_NASCseq_K562/fastq/{cell}_2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}_2.fastq.gz"
    log:
        outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g GCAGAGTACGGG -a CTGTCTCTTATACAC -a CCCGTACTCTGC \
            -G GTGTATAAGAGACAG -G GCAGAGTACGGG -A CTGTCTCTTATACAC -A CCCGTACTCTGC \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

# Tanglab

rule cutadapt_tanglab:
    input:
        fq1 = "data/datasets/{cell}_1.fastq.gz",
        fq2 = "data/datasets/{cell}_2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        outdir + "/cutadapt/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g ATCAACGCAGAGTAC -a CTGTCTCTTATACAC -a GTACTCTGCGTTGAT \
            -G GTGTATAAGAGACAG -G ATCAACGCAGAGTAC -A CTGTCTCTTATACAC -A GTACTCTGCGTTGAT \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

# Common

rule bowtie2: # remove rRNA
    input:
        fq1 = outdir + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/{run}/{cell}_2.fastq.gz",
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
        bowtie2 -p {threads} --local --no-unal --un-conc-gz {params.prefix}.fastq.gz \
            -x {input.idx}/ref -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam}
        mv {params.prefix}.fastq.1.gz {params.prefix}.1.fastq.gz
        mv {params.prefix}.fastq.2.gz {params.prefix}.2.fastq.gz ) &> {log}
        """