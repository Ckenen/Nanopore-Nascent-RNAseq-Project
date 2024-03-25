#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
bamdir = "results/mapping/stat_clip"
outdir = "results/assembly"

rule all:
    input:
        expand(outdir + "/stringtie/{run_cell}.gtf", run_cell=run_cells),
        expand(outdir + "/sqanti3/{run_cell}", run_cell=run_cells_cellline),
        #expand(outdir + "/tama/bed12/{run_cell}.bed", run_cell=run_cells),
        #expand(outdir + "/tama/merged/{cell_line}.filelist.txt", cell_line=cell_lines),
        #expand(outdir + "/tama/merged/{cell_line}.outputs", cell_line=cell_lines),
        #expand(outdir + "/tama/merged/{cell_line}.gtf", cell_line=cell_lines),
        #expand(outdir + "/tama/sqanti3/{cell_line}", cell_line=cell_lines),

rule stringtie:
    input:
        bam = bamdir + "/{run}/{cell}.bam",
        gtf = lambda wildcards: get_annotation_gtf(wildcards.cell)
    output:
        gtf = outdir + "/stringtie/{run}/{cell}.gtf"
    log:
        outdir + "/stringtie/{run}/{cell}.log"
    shell:
        """(
        stringtie {input.bam} -G {input.gtf} --fr -L \
            | sort -k1,1 -k4,4n | awk '$7!="."' > {output.gtf} ) &> {log}
        """

rule sqanti3:
    input:
        gtf1 = rules.stringtie.output.gtf,
        gtf2 = lambda wildcards: get_annotation_gtf(wildcards.cell),
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(outdir + "/sqanti3/{run}/{cell}")
    log:
        outdir + "/sqanti3/{run}/{cell}.log"
    conda:
        "SQANTI3.env"
    threads:
        4
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """

rule tama_format_gtf_to_bed12_stringtie:
    input:
        gtf = rules.stringtie.output.gtf
    output:
        bed = outdir + "/tama/bed12/{run}/{cell}.bed"
    conda:
        "py27"
    shell:
        """
        python {tamadir}/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
            {input} {output} &> /dev/null
        """

# def get_tama_bed(cell_line):
#     paths = []
#     for cell in get_cell_line_cells(cell_line):
#         paths.append(outdir + "/tama/bed12/%s/%s.bed" % (cell.split(".")[0], cell))
#     return paths

# rule make_filelist:
#     input:
#         beds = lambda wildcards: get_tama_bed(wildcards.cell_line)
#     output:
#         txt = outdir + "/tama/merged/{cell_line}.filelist.txt"
#     run:
#         with open(output.txt, "w+") as fw:
#             for bed in input.beds:
#                 cell = bed.split("/")[-1][:-4]
#                 if len(open(bed).readlines()) > 1000:
#                     fw.write("%s\tno_cap\t1,1,1\t%s\n" % (bed, cell))

# rule tama_merge:
#     input:
#         txt = rules.make_filelist.output.txt
#     output:
#         out = directory(outdir + "/tama/merged/{cell_line}.outputs")
#     log:
#         outdir + "/tama/merged/{cell_line}.log"
#     shell:
#         """
#         set +u; source activate py27
#         mkdir {output.out}
#         python {tamadir}/tama_merge.py -d merge_dup -f {input.txt} \
#             -p {output.out}/{wildcards.cell_line} > /dev/null 2> {log}
#         """

# rule tama_bed2gtf:
#     input:
#         outdir + "/tama/merged/{cell_line}.outputs"
#     output:
#         gtf = outdir + "/tama/merged/{cell_line}.gtf"
#     params:
#         bed = outdir + "/tama/merged/{cell_line}.outputs/{cell_line}.bed"
#     shell:
#         """
#         set +u; source activate py27
#         python {tamadir}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
#             {params.bed} {output.gtf} > /dev/null 2>&1
#         """

# rule tama_sqanti3:
#     input:
#         gtf1 = rules.tama_bed2gtf.output.gtf,
#         gtf2 = lambda wildcards: get_annotation_gtf(get_cell_line_cells(wildcards.cell_line)[0]),
#         fasta = lambda wildcards: get_genome_fasta(get_cell_line_cells(wildcards.cell_line)[0])
#     output:
#         out = directory(outdir + "/tama/sqanti3/{cell_line}")
#     log:
#         outdir + "/tama/sqanti3/{cell_line}.log"
#     threads:
#         8
#     shell:
#         """
#         ./scripts/assembly/run_sqanti3.sh {input} {threads} {output} &> {log}
#         """