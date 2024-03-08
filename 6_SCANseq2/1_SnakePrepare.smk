#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/prepare"
srr_list = [
    "SRR19279042", # UMI_200
    "SRR19279044", # UMI_100
    "SRR19279045", # 4CL
    # "SRR19279046", # 9CL_Mix
    "SRR19279047", # 9CL
    # "SRR19279040", # IGG24
    # "SRR19279041", # IGG6
    # "SRR19279043" # IGG48
]

rule all:
    input:
        # expand(outdir + "/sra/{srr}.sra", srr=srr_list),
        expand(outdir + "/fastq/{srr}.fastq.gz", srr=srr_list),
        expand(outdir + "/rename/{run}.fastq.gz", run=runs),

rule prefetch:
    output:
        sra = outdir + "/sra/{srr}.sra"
    log:
        outdir + "/sra/{srr}.log"
    shell:
        """
        prefetch --max-size 200000000 -o {output} {wildcards.srr} &> {log}
        """

rule fasterq_dump:
    input:
        sra = rules.prefetch.output.sra
    output:
        fq = outdir + "/fastq/{srr}.fastq.gz"
    log:
        outdir + "/fastq/{srr}.log"
    threads:
        12
    shell:
        """(
        fasterq-dump -e {threads} -O `dirname {output.fq}` {input.sra}
        pigz -p {threads} `dirname {output.fq}`/`basename {output.fq} .gz` ) &> {log}
        """

def get_srr_fastq(wildcards):
    srr = infos[infos["Library Name"] == wildcards.run]["Run"].values[0]
    return outdir + "/fastq/%s.fastq.gz" % srr

rule make_linker:
    input:
        fq = lambda wildcards: get_srr_fastq(wildcards)
    output:
        fq = outdir + "/rename/{run}.fastq.gz"
    shell:
        """
        ln -s `readlink -f {input.fq}` {output.fq}
        """