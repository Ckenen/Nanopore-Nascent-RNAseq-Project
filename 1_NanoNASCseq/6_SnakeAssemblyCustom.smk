#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
gtfdir = "results/assembly/stringtie"
outdir = "results/assembly_custom"
# strains = strains[:2]

rule all:
    input:
        expand(outdir + "/merged/{strain}.config.tsv", strain=strains),
        expand(outdir + "/merged/{strain}.outputs", strain=strains),
        expand(outdir + "/sqanti3/{strain}", strain=strains),
        expand(outdir + "/gtf/{strain}.all.gtf", strain=strains),
        # expand(outdir + "/sqanti3_new/{run_cell}", run_cell=run_cells_cellline),

rule make_config:
    output:
        tsv = outdir + "/merged/{strain}.config.tsv"
    run:
        cells = get_strain_cells(wildcards.strain)
        with open(output.tsv, "w+") as fw:
            fw.write("Cell\tGtf\n")
            for cell in cells:
                path = gtfdir + "/%s/%s.gtf" % (cell.split(".")[0], cell)
                fw.write("%s\t%s\n" % (cell, path))

rule merge_isoforms:
    input:
        tsv = rules.make_config.output.tsv
    output:
        out = directory(outdir + "/merged/{strain}.outputs")
    shell:
        """
        ./scripts/assembly/merge_isoforms.py {input.tsv} {output.out}
        """

rule sqanti3:
    input:
        gtf1 = rules.merge_isoforms.output.out,
        gtf2 = lambda wildcards: get_annotation_gtf(get_strain_cells(wildcards.strain)[0]),
        fa = lambda wildcards: get_genome_fasta(get_strain_cells(wildcards.strain)[0])
    output:
        out = directory(outdir + "/sqanti3/{strain}")
    log:
        outdir + "/sqanti3/{strain}.log"
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
        ref = lambda wildcards: get_annotation_gtf(get_strain_cells(wildcards.strain)[0]),
    output:
        gtf1 = outdir + "/gtf/{strain}.novel.gtf",
        gtf2 = outdir + "/gtf/{strain}.all.gtf"
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

def get_strain_novel_gtf(cell):
    strain = dat[dat["Cell"] == cell]["Strain"].values[0]
    if strain == "K562":
        return outdir + "/gtf/K562.all.gtf"
    elif strain == "mESC":
        return outdir + "/gtf/mESC.all.gtf"
    else:
        return outdir + "/gtf/MouseBlastocyst.all.gtf"
    

rule sqanti3_new:
    input:
        gtf1 = gtfdir + "/{run}/{cell}.gtf",
        gtf2 = lambda wildcards: get_strain_novel_gtf(wildcards.cell),
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(outdir + "/sqanti3_new/{run}/{cell}")
    log:
        outdir + "/sqanti3_new/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """