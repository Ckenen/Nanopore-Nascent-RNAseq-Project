#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/demux"


rule all:
    input:
        expand(outdir + "/read_count/{run}.txt", run=runs),
        expand(outdir + "/barcode/{run}_bc5.fasta", run=runs),
        expand(outdir + "/barcode/{run}_bc3.fasta", run=runs),
        # expand(outdir + "/nanoplexer_5/{run}", run=runs),
        # expand(outdir + "/nanoplexer_3/{run}", run=runs),
        # expand(outdir + "/pychopper/{run_cell}.fastq", run_cell=run_cells),
        # expand(outdir + "/cutadapt/{run_cell}.fastq.gz", run_cell=run_cells),
        expand(outdir + "/trim_polya/{run_cell}.fastq.gz", run_cell=run_cells),


rule get_read_count:
    input:
        fq = "results/prepare/rename/{run}.fastq.gz"
    output:
        txt = outdir + "/read_count/{run}.txt"
    shell:
        """
        zcat {input.fq} | wc -l | awk '{{print $1/4}}' > {output.txt}
        """

def get_barcode_list_5(run):
    tmp = table[table["Library"] == run]
    return ["Bar%d" % bc for bc in sorted(set(tmp["Barcode.5"]))]

def get_barcode_list_3(run):
    tmp = table[table["Library"] == run]
    return ["Bar%d" % bc for bc in sorted(set(tmp["Barcode.3"]))]

rule get_barcode:
    input:
        fa = "../1_NanoNASCseq/data/20210831_barcode_corrected_bar80.fasta"
    output:
        fa5 = outdir + "/barcode/{run}_bc5.fasta",
        fa3 = outdir + "/barcode/{run}_bc3.fasta"
    params:
        barcode_list_5 = lambda wildcards: get_barcode_list_5(wildcards.run),
        barcode_list_3 = lambda wildcards: get_barcode_list_3(wildcards.run)
    shell:
        """
        samtools faidx {input.fa} {params.barcode_list_5} > {output.fa5}
        samtools faidx {input.fa} {params.barcode_list_3} > {output.fa3}
        """

rule nanoplexer_5:
    input:
        fq = "results/prepare/rename/{run}.fastq.gz",
        fa = rules.get_barcode.output.fa5
    output:
        out = directory(outdir + "/nanoplexer_5/{run}")
    log:
        outdir + "/nanoplexer_5/{run}.log"
    threads:
        24
    shell:
        """
        nanoplexer -t {threads} -b {input.fa} -p {output} {input.fq} &> {log}
        rm {output}/unclassified.fastq
        """

rule nanoplexer_3:
    input:
        fqs = rules.nanoplexer_5.output,
        fa = rules.get_barcode.output.fa3
    output:
        out = directory(outdir + "/nanoplexer_3/{run}")
    log:
        outdir + "/nanoplexer_3/{run}.log"
    threads:
        24
    shell:
        """
        (
        mkdir {output}
        for fq in {input.fqs}/Bar*.fastq; do
            out={output}/`basename $fq .fastq`
            nanoplexer -t {threads} -b {input.fa} -p $out $fq &> ${{out}}.log
            rm ${{out}}/unclassified.fastq
        done ) &> {log}
        """

# pychopper

def get_cell_barcode_5(wildcards):
    return "Bar%s" % table[table["Rename"] == wildcards.cell]["Barcode.5"].values[0]

def get_cell_barcode_3(wildcards):
    return "Bar%s" % table[table["Rename"] == wildcards.cell]["Barcode.3"].values[0]

def get_cell_barcode_sequence_5(wildcards):
    return barcodes["Bar%s" % table[table["Rename"] == wildcards.cell]["Barcode.5"].values[0]][0]

def get_cell_barcode_sequence_3(wildcards):
    return barcodes["Bar%s" % table[table["Rename"] == wildcards.cell]["Barcode.3"].values[0]][0]

rule pychopper:
    input:
        fqs = rules.nanoplexer_3.output,
        txt = "data/primer_config.txt"
    output:
        txt = outdir + "/pychopper/{run}/{cell}.primers.txt",
        fq = outdir + "/pychopper/{run}/{cell}.fastq",
        pdf = outdir + "/pychopper/{run}/{cell}.pdf"
    log:
        outdir + "/pychopper/{run}/{cell}.log"
    params:
        bc_5 = lambda wildcards: get_cell_barcode_5(wildcards),
        bc_3 = lambda wildcards: get_cell_barcode_3(wildcards),
        bc_seq_5 = lambda wildcards: get_cell_barcode_sequence_5(wildcards),
        bc_seq_3 = lambda wildcards: get_cell_barcode_sequence_3(wildcards)
    threads:
        8
    shell:
        """
        echo ">MySSP" > {output.txt}
        echo "{params.bc_seq_5}AAGCAGTGGTATCAACGCAGAGTACATGGG" >> {output.txt}
        echo ">MyVNP" >> {output.txt}
        echo "TCAGACGTGTGCTCTTCCGATC{params.bc_seq_3}" >> {output.txt}
        pychopper -m edlib -r {output.pdf} -b {output.txt} -c {input.txt} -t {threads} \
            {input.fqs}/{params.bc_5}/{params.bc_3}.fastq {output.fq} &> {log}
        """

# cutadapt (deprecated)

# rule cutadapt:
#     input:
#         fq = rules.pychopper.output.fq
#     output:
#         fq = outdir + "/cutadapt/{run}/{cell}.fastq.gz",
#         fq2 = outdir + "/cutadapt/{run}/{cell}_untrimmed.fastq.gz"
#     log:
#         outdir + "/cutadapt/{run}/{cell}.log"
#     threads:
#         12
#     shell:
#         """
#         cutadapt -j {threads} -a "A{{15}}" -m 200 -e 0.2 \
#             --untrimmed-output {output.fq2} -o {output.fq} {input.fq} &> {log}
#         """

# remove polyA by custom script

rule trim_polya_and_extract_umi:
    input:
        fq = rules.pychopper.output.fq
    output:
        fq = temp(outdir + "/trim_polya/{run}/{cell}.fastq"),
        fq2 = outdir + "/trim_polya/{run}/{cell}.fastq.gz"
    log:
        outdir + "/trim_polya/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/trim_polya_and_extract_umi.py {input.fq} {output.fq} &> {log}
        pigz -p {threads} -c {output.fq} > {output.fq2}
        """
