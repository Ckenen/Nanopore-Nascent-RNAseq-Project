#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_duplicate"
outdir = "results/expression"
#run_cells = run_cells_cellline
#run_cells.remove("20220729_K562R3/20220729_K562R3.C07")

rule all:
    input:
        # expand(outdir + "/expressed_alleles/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/collapsed/{run_cell}.gtf", run_cell=run_cells),
        # expand(outdir + "/sqanti3/{run_cell}", run_cell=run_cells),
        expand(outdir + "/quant_genes/min_read_1_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_genes/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_genes/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/isoform_category/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_isoforms/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/quant_isoforms/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/novel/sqanti3/{run_cell}", run_cell=run_cells),

rule stat_expressed_alleles:
    input:
        bam = indir + "/{run}/{cell}.bam",
        vcf = lambda wildcards: get_snp_vcf(wildcards.cell)
    output:
        tsv = outdir + "/expressed_alleles/{run}/{cell}.tsv"
    log:
        outdir + "/expressed_alleles/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/stat_expressed_alleles.py {input} > {output.tsv} 2> {log}
        """

rule collapse:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        bed = temp(outdir + "/collapsed/{run}/{cell}.bed"),
        bed_gz = outdir + "/collapsed/{run}/{cell}.bed.gz",
        gtf = outdir + "/collapsed/{run}/{cell}.gtf",
        gtf_gz = outdir + "/collapsed/{run}/{cell}.gtf.gz"
    log:
        outdir + "/collapsed/{run}/{cell}.log"
    shell:
        """(
        ./scripts/expression/collapse_umi.py {input.bam} | sort -k1,1 -k2,2n -k3,3n > {output.bed}
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed -f {output.bed_gz}
        ./scripts/expression/bed2gtf.py {output.bed} | sort -k1,1 -k4,4n > {output.gtf}
        bgzip -c {output.gtf} > {output.gtf_gz}
        tabix -p gff -f {output.gtf_gz} ) &> {log}
        """

rule sqanti3:
    input:
        gtf1 = rules.collapse.output.gtf,
        gtf2 = lambda wildcards: get_annotation_gtf(wildcards.cell),
        fa = lambda wildcards: get_genome_fasta(wildcards.cell)
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

# Quantify genes

rule quant_genes:
    input:
        sqanti3_out = rules.sqanti3.output.out,
        event_tsv = "results/mismatch/ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = rules.stat_expressed_alleles.output.tsv
    output:
        tsv = outdir + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        outdir + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    params:
        sqanti3_tsv = rules.sqanti3.output.out + "/{cell}_classification.txt"
    shell:
        """
        ./scripts/expression/quant_genes.py {params.sqanti3_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """

rule stat_isoform_category:
    input:
        gtf = lambda wildcards: get_annotation_gtf(wildcards.cell),
        bed = rules.collapse.output.bed_gz
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
        stat_tsv = rules.stat_isoform_category.output.tsv,
        event_tsv = "results/mismatch/ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = rules.stat_expressed_alleles.output.tsv
    output:
        tsv = outdir + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        outdir + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """

# Novel isoforms

# def get_novel_gtf(cell):
#     species = get_species(cell)
#     strain = get_strain(cell)
#     if strain == "K562":
#         return "results/assembly/custom_merge/gtf/K562.all.gtf"
#     elif strain == "mESC":
#         return "results/assembly/custom_merge/gtf/mESC.all.gtf"
#     elif "blast" in strain.lower() and species == "Mouse":
#         return "results/assembly/custom_merge/gtf/MouseBlastocyst.all.gtf"
#     else:
#         print("Error strain:", strain)
#         assert False

# rule sqanti3_novel:
#     input:
#         gtf1 = rules.collapse.output.gtf,
#         gtf2 = lambda wildcards: get_novel_gtf(wildcards.cell),
#         fasta = lambda wildcards: FILES[get_species(wildcards.cell)]["GENOME_FASTA"]
#     output:
#         out = directory(outdir + "/novel/sqanti3/{run}/{cell}")
#     log:
#         outdir + "/novel/sqanti3/{run}/{cell}.log"
#     params:
#         sqanti3 = "/home/chenzonggui/software/SQANTI3-4.2"
#     threads:
#         4
#     shell:
#         """
#         ./scripts/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
#         """

