#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
samples = [
    "20221128_K562_Actd_0h_rep1", "20221128_K562_Actd_0h_rep2",
    "20221128_K562_Actd_3h_rep1", "20221128_K562_Actd_3h_rep2",
    "20221128_K562_Actd_6h_rep1", "20221128_K562_Actd_6h_rep2"]
indir = "results/mapping/rmdup"
outdir = "results/snps"

rule all:
    input:
        expand(outdir + "/vcfs/{sample}.vcf.gz", sample=samples[:1]),
        expand(outdir + "/phased_bam/{sample}.bam", sample=samples),
        expand(outdir + "/counts/{sample}.tsv", sample=samples),

rule call_het_snps:
    input:
        bam = indir + "/{sample}.human.bam"
    output:
        vcf = temp(outdir + "/vcfs/{sample}.vcf"),
        vcf_gz = outdir + "/vcfs/{sample}.vcf.gz"
    log:
        outdir + "/vcfs/{sample}.log"
    threads:
        24
    shell:
        """(
        ./scripts/call_het_snps_from_rnaseq.py {input.bam} {wildcards.sample} {threads} {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz} ) &> {log}
        """

rule haplotag:
    input:
        vcf = outdir + "/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz",
        bam = indir + "/{sample}.human.bam"
    output:
        bam = outdir + "/phased_bam/{sample}.bam"
    log:
        outdir + "/phased_bam/{sample}.log"
    shell:
        """
        whatshap haplotag --ignore-read-groups -o {output.bam} {input.vcf} {input.bam} &> {log}
        samtools index {output.bam}
        """

rule stat_hp_reads:
    input:
        vcf = outdir + "/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz",
        bam = outdir + "/phased_bam/{sample}.bam"
    output:
        tsv = outdir + "/counts/{sample}.tsv"
    log:
        outdir + "/counts/{sample}.log"
    shell:
        """
        ./scripts/stat_hp_reads.py {input.vcf} {input.bam} {output.tsv} &> {log}
        """