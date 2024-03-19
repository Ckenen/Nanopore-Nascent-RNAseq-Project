#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/marked_strand"
outdir = "results/mismatch"

rule all:
    input:
        #expand(outdir + "/events/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/marked_nascent/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/ratio/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/pc/{run_cell}.tsv", run_cell=run_cells),
        #expand(outdir + "/signal2noise/{run_cell}.tsv", run_cell=run_cells),
        #expand(outdir + "/estimate_pc/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/expression/fpkm/{run_cell}.tsv", run_cell=run_cells[:1]),

# Events

rule get_events:
    input:
        bam = indir + "/{run}/{cell}.bam",
        bed = SNP_BED_GZ
    output:
        bam = outdir + "/events/{run}/{cell}.bam"
    log:
        outdir + "/events/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools GetEvent --threads {threads} --snp {input.bed} \
            {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_nascent:
    input:
        bam = rules.get_events.output.bam
    output:
        bam = outdir + "/marked_nascent/{run}/{cell}.bam"
    log:
        outdir + "/marked_nascent/{run}/{cell}.log"
    params:
        layout = lambda wildcards: get_layout(wildcards.cell)
    shell:
        """
        nasctools MarkNascent --platform NGS --layout {params.layout} \
            {input.bam} {output.bam} &> {log}
        samtools index {output.bam}
        """

rule report_mismatch:
    input:
        bam = rules.mark_nascent.output.bam
    output:
        tsv = outdir + "/ratio/{run}/{cell}.tsv"
    log:
        outdir + "/ratio/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools ReportMismatch --threads {threads} --strand TAG \
            --strand-tag ST {input.bam} {output.tsv} &> {log}
        """

rule estimate_pc:
    input:
        bam = rules.mark_nascent.output.bam,
        tsv1 = rules.report_mismatch.output.tsv,
        tsv2 = "results/mismatch/Estimate.Pe.Model.K562.PE.tsv"
    output:
        txt = outdir + "/pc/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/estimate_pc.ngs.py {input} > {output}
        """

# rule stat_signal_to_noise:
#     input:
#         bam = rules.mark_nascent.output.bam,
#         txt = rules.report_mismatch.output.txt
#     output:
#         txt = outdir + "/signal2noise/{run}/{cell}.tsv"
#     log:
#         outdir + "/signal2noise/{run}/{cell}.log"
#     shell:
#         """
#         ../public/scripts/mismatch/estimate_pe_pc.py \
#             {input} {output} &> {log}
#         """

# rule estimate_pc:
#     input:
#         bam = outdir + "/marked_nascent/{run}/{cell}.bam",
#         txt = outdir + "/ratio/{run}/{cell}.tsv"
#     output:
#         tsv = outdir + "/estimate_pc/{run}/{cell}.tsv"
#     log:
#         outdir + "/estimate_pc/{run}/{cell}.log"
#     run:
#         import pandas as pd
#         m = pd.read_csv(input.txt, sep="\t", header=0, index_col=0)
#         import json
#         model = json.load(open("pe_model.json"))
#         r = 0
#         for t in model:
#             r += m.loc[t]["Ratio[NoSNP]"] * model[t]["K"] * model[t]["W"]
#         p_e = r
#         # print(p_e)
#         import subprocess
#         subprocess.check_call("./scripts/estimate_pc.ngs.py %s %s %s &> %s" % (input.bam, p_e, output.tsv, log), shell=True)

# FPKM

rule calculate_fpkm:
    input:
        bam = rules.mark_nascent.output.bam,
        bed = TRANSCRIPT_BED_GZ,
        txt = ANNOTATION_TSV
    output:
        txt = outdir + "/expression/fpkm/{run}/{cell}.tsv"
    log:
        outdir + "/expression/fpkm/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM --threads {threads} --layout PE --strand TAG --strand-tag ST --nascent \
            --annotation {input.txt} {input.bam} {input.bed} {output.txt} &> {log}
        """