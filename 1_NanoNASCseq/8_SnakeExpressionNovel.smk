#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
gtfdir = "results/expression/collapsed"
beddir = "results/expression/collapsed"
outdir = "results/expression_novel"
run_cells = run_cells_cellline + run_cells_blastocyst

rule all:
    input:
        # expand(outdir + "/sqanti3/{run_cell}", run_cell=run_cells),
        # expand(outdir + "/quant_genes/min_read_1_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/quant_genes/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/quant_genes/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/isoform_category/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_isoforms/min_read_1_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_isoforms/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_isoforms/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=run_cells),

def get_novel_gtf(cell):
    species = get_species(cell)
    strain = get_strain(cell)
    if strain == "K562":
        return "results/assembly_custom/gtf/K562.all.gtf"
    elif strain == "mESC":
        return "results/assembly_custom/gtf/mESC.all.gtf"
    elif "blast" in strain.lower() and species == "Mouse":
        return "results/assembly_custom/gtf/MouseBlastocyst.all.gtf"
    else:
        print("Error strain:", strain)
        assert False

rule sqanti3:
    input:
        gtf1 = gtfdir + "/{run}/{cell}.gtf",
        gtf2 = lambda wildcards: get_novel_gtf(wildcards.cell),
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(outdir + "/sqanti3/{run}/{cell}")
    log:
        outdir + "/sqanti3/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """

rule quant_genes:
    input:
        sqanti3_out = outdir + "/sqanti3/{run}/{cell}",
        event_tsv = "results/mismatch/ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = "results/expression/expressed_alleles/{run}/{cell}.tsv"
    output:
        tsv = outdir + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        outdir + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    params:
        sqanti3_tsv = outdir + "/sqanti3/{run}/{cell}/{cell}_classification.txt"
    shell:
        """
        ./scripts/expression/quant_genes.py {params.sqanti3_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """

rule stat_isoform_category:
    input:
        gtf = lambda wildcards: get_novel_gtf(wildcards.cell),
        bed = beddir + "/{run}/{cell}.bed.gz"
    output:
        tsv = outdir + "/isoform_category/{run}/{cell}.tsv"
    log:
        outdir + "/isoform_category/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/expression/stat_isoform_category.py {input} {output} &> {log}
        """

rule quant_isoforms:
    input:
        stat_tsv = outdir + "/isoform_category/{run}/{cell}.tsv",
        event_tsv = "results/mismatch/ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = "results/expression/expressed_alleles/{run}/{cell}.tsv"
    output:
        tsv = outdir + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        outdir + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """
