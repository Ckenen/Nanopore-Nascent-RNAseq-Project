#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

new_run_cells = []
for run_cell in run_cells:
    run = run_cell.split("/")[0]
    if run != "20220719_K562R3":
        continue
    cell = run_cell.split("/")[1]
    if dat.loc[cell]["Strain"] != "K562":
        continue
    if not np.isnan(dat.loc[cell]["ActD"]):
        continue
    new_run_cells.append(run_cell)
run_cells = new_run_cells

outdir = "results/bcr_abl_fusion"

rule all:
    input:
        expand(outdir + "/minimap2/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),

rule minimap2:
    input:
        fq = "results/demux/trimmed/{run}/{cell}.fastq.gz",
        mmi = outdir + "/grch38_bcr_abl_fusion.mm2.splice.mmi",
        bed = outdir + "/transcripts.bed"
    output:
        bam = outdir + "/minimap2/{run}/{cell}.bam"
    log:
        outdir + "/minimap2/{run}/{cell}.log"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        12
    shell: 
        """(
        minimap2 -ax splice -u f -Y --MD -R '{params.rg}' --junc-bed {input.bed} -t {threads} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
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
        samtools view -@ {threads} -F 2308 -q 30 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """