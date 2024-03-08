#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/filtered"
outdir = "results/assembly"

rule all:
    input:
        expand(outdir + "/stringtie/{run_cell}.gtf.gz", run_cell=run_cells),
        expand(outdir + "/sqanti3/{run_cell}", run_cell=run_cells),

# StringTie

rule stringtie:
    input:
        bam = indir + "/{run}/{cell}.bam",
        gtf = lambda wildcards: FILES[get_species(wildcards.cell)]["ANNOTATION_GTF"]
    output:
        gtf = outdir + "/stringtie/{run}/{cell}.gtf.gz"
    log:
        outdir + "/stringtie/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        stringtie {input.bam} -G {input.gtf} --fr -L | bedtools sort | awk '$7!="."' | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf} ) &> {log}
        """

# SQANTI3

rule sqanti3:
    input:
        gtf_que = outdir + "/stringtie/{run}/{cell}.gtf.gz",
        gtf_ref = lambda wildcards: FILES[get_species(wildcards.cell)]["ANNOTATION_GTF"],
        fasta = lambda wildcards: FILES[get_species(wildcards.cell)]["GENOME_FASTA"]
    output:
        out = directory(outdir + "/sqanti3/{run}/{cell}")
    log:
        outdir + "/sqanti3/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        set +u; source activate SQANTI3.env
        gzip -d -c {input.gtf_que} > {output.out}.gtf
        ~/software/SQANTI3-4.2/sqanti3_qc.py --cpus {threads} --skipORF --report pdf -d {output.out} \
            {output.out}.gtf {input.gtf_ref} {input.fasta}
        ~/software/SQANTI3-4.2/sqanti3_RulesFilter.py --report skip \
            {output.out}/{wildcards.cell}_classification.txt \
            {output.out}/{wildcards.cell}_corrected.fasta \
            {output.out}/{wildcards.cell}_corrected.gtf 
        rm {output.out}/*.genePred
        rm {output.out}/*.fasta
        rm {output.out}/*.gff
        rm {output.out}/*.gtf
        rm {output.out}/*_junctions.txt
        rm {output.out}/*params.txt
        rm -r {output.out}/GMST
        rm -r {output.out}/RTS
        rm {output.out}.gtf
        gzip {output.out}/*.txt ) &> {log}
        """