#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/demux/trimmed"
OUTDIR = "results/mapping"

rule all:
    input:
        expand(OUTDIR + "/minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/minimap2/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/extract_umi/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/stat_clip/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/remove_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/chrom_reads/{run_cell}.tsv", run_cell=RUN_CELLS),

rule minimap2:
    input:
        fq = INDIR + "/{run}/{cell}/trimmed.fastq.gz",
        mmi = lambda wildcards: get_genome_splice_mmi(wildcards.cell),
        bed = lambda wildcards: get_transcript_bed(wildcards.cell)
    output:
        bam = OUTDIR + "/minimap2/{run}/{cell}.bam"
    log:
        OUTDIR + "/minimap2/{run}/{cell}.log"
    conda:
        "minimap2"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        THREADS
    shell:
        """(
        minimap2 -ax splice -u f \
            -Y --MD \
            -R '{params.rg}' \
            --junc-bed {input.bed} \
            -t {threads} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/filtered/{run}/{cell}.bam"
    log:
        OUTDIR + "/filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule extract_umi:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/extract_umi/{run}/{cell}.bam"
    log:
        OUTDIR + "/extract_umi/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        ./scripts/mapping/extract_umi.py {input.bam} {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule stat_clip:
    input:
        bam = rules.extract_umi.output.bam,
    output:
        bam = OUTDIR + "/stat_clip/{run}/{cell}.bam",
        tsv = OUTDIR + "/stat_clip/{run}/{cell}.tsv"
    log:
        OUTDIR + "/stat_clip/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        nasctools StatClip -c 20 -s {output.tsv} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

# Time-cost: 20221218_BlastocystC69.C33

rule mark_duplicate: 
    input:
        bam = rules.stat_clip.output.bam
    output:
        bam = OUTDIR + "/mark_duplicate/{run}/{cell}.bam",
        tsv = OUTDIR + "/mark_duplicate/{run}/{cell}.tsv"
    log:
        OUTDIR + "/mark_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        ./scripts/mapping/mark_duplicate.py {input.bam} {output.bam} {output.tsv}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/remove_duplicate/{run}/{cell}.bam",
    log:
        OUTDIR + "/remove_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -F 1024 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule stat_chrom_reads:
    input:
        bam = rules.minimap2.output.bam
    output:
        tsv = OUTDIR + "/chrom_reads/{run}/{cell}.tsv"
    threads:
        4
    shell:
        """
        ./scripts/mapping/stat_chrom_reads.sh {input.bam} {output.tsv}
        """

rule flagstat:
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
