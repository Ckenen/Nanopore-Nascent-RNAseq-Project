#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
GTFDIR = "results/assembly/stringtie"
OUTDIR = "results/assembly_custom"

rule all:
    input:
        expand(OUTDIR + "/merged/{group}.config.tsv", group=GROUPS),
        expand(OUTDIR + "/merged/{group}.outputs", group=GROUPS),
        expand(OUTDIR + "/sqanti3/{group}", group=GROUPS),
        expand(OUTDIR + "/gtf/{group}.all.gtf", group=GROUPS),
        expand(OUTDIR + "/gtf_full/{group}.gtf", group=GROUPS),
        # expand(OUTDIR + "/sqanti3_new/{run_cell}", run_cell=run_cells_cellline),

rule make_config:
    output:
        tsv = OUTDIR + "/merged/{group}.config.tsv"
    run:
        cells = get_group_cells(wildcards.group)
        with open(output.tsv, "w+") as fw:
            fw.write("Cell\tGtf\n")
            for cell in cells:
                path = GTFDIR + "/%s/%s.gtf" % (cell.split(".")[0], cell)
                fw.write("%s\t%s\n" % (cell, path))

rule merge_isoforms:
    input:
        tsv = rules.make_config.output.tsv
    output:
        out = directory(OUTDIR + "/merged/{group}.outputs")
    shell:
        """
        ./scripts/assembly/merge_isoforms.py {input.tsv} {output.out}
        """

rule sqanti3:
    input:
        gtf1 = rules.merge_isoforms.output.out,
        gtf2 = lambda wildcards: get_annotation_gtf(get_group_cells(wildcards.group)[0]),
        fa = lambda wildcards: get_genome_fasta(get_group_cells(wildcards.group)[0])
    output:
        out = directory(OUTDIR + "/sqanti3/{group}")
    log:
        OUTDIR + "/sqanti3/{group}.log"
    params:
        gtf1 = rules.merge_isoforms.output.out + "/merged.gtf"
    threads:
        8
    shell:
        """
        ./scripts/assembly/run_sqanti3_orf.sh {params.gtf1} {input.gtf2} {input.fa} {threads} {output} &> {log}
        """

rule merge_gtf:
    input:
        sqdir = rules.sqanti3.output.out,
        ref = lambda wildcards: get_annotation_gtf(get_group_cells(wildcards.group)[0]),
    output:
        gtf1 = OUTDIR + "/gtf/{group}.novel.gtf",
        gtf2 = OUTDIR + "/gtf/{group}.all.gtf"
    params:
        tsv = rules.sqanti3.output.out + "/merged_classification.txt",
        gtf = rules.sqanti3.output.out + "/merged_corrected.gtf.cds.gff",
    shell:
        """
        ./scripts/assembly/report_novel_isoforms.py {params.tsv} {params.gtf} | sort -k1,1 -k4,4n > {output.gtf1}
        cat {input.ref} {output.gtf1} | sort -k1,1 -k4,4n > {output.gtf2}
        bgzip -c {output.gtf1} > {output.gtf1}.gz
        tabix -p gff {output.gtf1}.gz
        bgzip -c {output.gtf2} > {output.gtf2}.gz
        tabix -p gff {output.gtf2}.gz
        """

rule add_attributes:
    input:
        gtf = rules.merge_gtf.output.gtf2
    output:
        gtf = OUTDIR + "/gtf_full/{group}.gtf"
    shell:
        """
        ./scripts/assembly/add_attributes.py {input.gtf} {output.gtf}
        bgzip -c {output.gtf} > {output.gtf}.gz
        tabix -p gff {output.gtf}.gz
        """

def get_group_novel_gtf(cell):
    group = DAT[DAT["Cell"] == cell]["group"].values[0]
    if group == "K562":
        return OUTDIR + "/gtf/K562.all.gtf"
    elif group == "mESC":
        return OUTDIR + "/gtf/mESC.all.gtf"
    else:
        return OUTDIR + "/gtf/MouseBlastocyst.all.gtf"
    
rule sqanti3_new:
    input:
        gtf1 = GTFDIR + "/{run}/{cell}.gtf",
        gtf2 = lambda wildcards: get_group_novel_gtf(wildcards.cell),
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(OUTDIR + "/sqanti3_new/{run}/{cell}")
    log:
        OUTDIR + "/sqanti3_new/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """