#!/usr/bin/env runsnakemake
import glob
cells = []
for path in glob.glob("results/linkage_consensus/*.ratio.tsv"):
    cells.append(path.split("/")[-1][:-10])
# print(len(cells))

rule all:
    input:
        expand("results/linkage_consensus/{cell}.pc.tsv", cell=cells)

rule estimate_pc:
    input:
        tsv1 = "results/Estimate.Pe.Model.linkage_consensus.K562.tsv",
        tsv2 = "results/linkage_consensus/{cell}.ratio.tsv",
        tsv3 = "results/linkage_consensus/{cell}.events.tsv",
    output:
        tsv = "results/linkage_consensus/{cell}.pc.tsv"
    shell:
        """
        ../../1_FLAIRseq/scripts/signal2noise/estimate_pc.py -m long -e {input.tsv3} {input.tsv1} {input.tsv2} > {output}
        """
