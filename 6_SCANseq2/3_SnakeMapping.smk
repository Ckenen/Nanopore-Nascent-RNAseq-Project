#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/mapping"

rule all:
    input:
        expand(outdir + "/minimap2/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/minimap2/{run_cell}.flagstat", run_cell=run_cells),
        # expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/filtered/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/stat_clip/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/stat_clip/{run_cell}.flagstat", run_cell=run_cells),

rule minimap2:
    input:
        fq = "results/demux/trim_polya/{run}/{cell}.fastq.gz",
        mmi = lambda wildcards: FILES[get_species(wildcards.cell)]["GENOME_SPLICE_MMI"],
        bed = lambda wildcards: FILES[get_species(wildcards.cell)]["TRANSCRIPT_BED"]
    output:
        bam = outdir + "/minimap2/{run}/{cell}.bam"
    log:
        outdir + "/minimap2/{run}/{cell}.log"
    params:
        rg = "@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}"
    threads:
        12
    shell:
        """
        minimap2 -a -x splice -u f -Y --MD -R "{params.rg}" --junc-bed {input.bed} -t {threads} {input.mmi} {input.fq} \
            | samtools view --no-PG -@ {threads} -u - \
            | samtools sort --no-PG -@ {threads} -T {output.bam} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = outdir + "/filtered/{run}/{cell}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' -q 30 -m 200 -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/stat_clip/{run}/{cell}.bam",
        txt = outdir + "/stat_clip/{run}/{cell}.tsv"
    log:
        outdir + "/stat_clip/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools StatClip --max-clip 5 -s {output.txt} -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """

        
