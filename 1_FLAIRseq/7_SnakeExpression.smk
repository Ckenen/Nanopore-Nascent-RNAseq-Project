#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/mark_duplicate"
OUTDIR = "results/expression"

rule all:
    input:
        expand(OUTDIR + "/expressed_alleles/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/collapsed/{run_cell}.gtf", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/sqanti3/{run_cell}", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/quant_genes/min_read_1_min_tc_1/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/quant_genes/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/quant_genes/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/isoform_category/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/quant_isoforms/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/quant_isoforms/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=RUN_CELLS),

rule stat_expressed_alleles:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        vcf = lambda wildcards: get_snp_vcf(wildcards.cell)
    output:
        tsv = OUTDIR + "/expressed_alleles/{run}/{cell}.tsv"
    log:
        OUTDIR + "/expressed_alleles/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/stat_expressed_alleles.py {input} > {output.tsv} 2> {log}
        """

rule collapse:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        bed = OUTDIR + "/collapsed/{run}/{cell}.bed",
        bed_gz = OUTDIR + "/collapsed/{run}/{cell}.bed.gz",
        gtf = OUTDIR + "/collapsed/{run}/{cell}.gtf",
        gtf_gz = OUTDIR + "/collapsed/{run}/{cell}.gtf.gz"
    log:
        OUTDIR + "/collapsed/{run}/{cell}.log"
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
        out = directory(OUTDIR + "/sqanti3/{run}/{cell}")
    log:
        OUTDIR + "/sqanti3/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """

rule quant_genes:
    input:
        sqanti3_out = rules.sqanti3.output.out,
        event_tsv = "results/mismatch/ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = rules.stat_expressed_alleles.output.tsv
    output:
        tsv = OUTDIR + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
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
        tsv = OUTDIR + "/isoform_category/{run}/{cell}.tsv"
    log:
        OUTDIR + "/isoform_category/{run}/{cell}.log"
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
        tsv = OUTDIR + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """