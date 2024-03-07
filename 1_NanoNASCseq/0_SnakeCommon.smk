#!/usr/bin/env runsnakemake
import os
import numpy as np
import pandas as pd

configfile: "config.yaml"
runs = config["runs"]
threads = 8

## NanoNASCseq.xlsx

dat = pd.read_excel("data/NanoNASCseq.xls")
dat.index = dat["Cell"]
dat["RunCell"] = ["%s/%s" % (run, cell) for run, cell in dat[["Run", "Cell"]].values]

path = "results/mapping/all_cells_final_reads.tsv"
if os.path.exists(path):
    d = pd.read_csv(path, sep="\t", index_col=0)
    d["UniqReads"] = d["Reads"] - d["DuplicateReads"]
    d["UMIs"] = d["UniqReads"]
    dat = dat.merge(d, left_index=True, right_index=True, how="left")

infos = dat
def get_species(cell):
    return infos.loc[cell]["Species"]
def get_strain(cell):
    return infos.loc[cell]["Strain"]
def get_cell_barcode(cell):
    return "Bar%d" % infos.loc[cell]["Barcode"].values[0]
def get_run_barcodes(run):
    return ["Bar%d" % x for x in infos[infos["Run"] == run]["Barcode"]]

# Filtered by run
dat_selected = dat[[r in runs for r in dat["Run"]]]
dat_k562 = dat_selected[dat_selected["Strain"] == "K562"]
dat_mesc = dat_selected[dat_selected["Strain"] == "mESC"]
dat_blastocyst = dat_selected[["blastocyst" not in x.lower() for x in dat_selected["Strain"]]]
run_cells_selected = list(dat_selected["RunCell"])
run_cells_k562 = list(dat_k562["RunCell"])
run_cells_mesc = list(dat_mesc["RunCell"])
run_cells_blastocyst = list(dat_blastocyst["RunCell"])
run_cells = run_cells_selected
#run_cells.remove("20221218_BlastocystC69/20221218_BlastocystC69.C33")
# run_cells = run_cells_k562 + run_cells_mesc
run_cells_cellline = run_cells_k562 + run_cells_mesc
#print("Selected runs: %d" % len(runs))
#print("Selected cells: %d" % len(run_cells_selected))
#print("Final cells: %s" % len(run_cells))


def get_estimate_pe_model(cell):
    strain = dat[dat["Cell"] == cell]["Strain"].values[0]
    if strain == "K562":
        return "results/mismatch/Estimate.Pe.Model.consensus.K562.tsv"
    elif strain == "mESC":
        return "results/mismatch/Estimate.Pe.Model.consensus.mESC.tsv"
    else:
        assert False

## Files

FILES = {
    "Human": {
        "GENOME_FASTA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        "GENOME_SPLICE_MMI": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.mm2.splice.mmi",
        "TRANSCRIPTOME_FASTA": "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.fa",
        "TRANSCRIPTOME_MMI": "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.mm2.map-ont.mmi",
        "TID2GID_TSV": "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/transcript_id_to_gene_id.tsv",
        "TRANSCRIPT_BED": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed",
        "ANNOTATION_GTF": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf",
        "SNP_BED_GZ": "/home/chenzonggui/species/homo_sapiens/hg38/snp151.3.lite.bed.gz",
        "SNP_VCF_GZ": "../4_RNAseq_ActD/results/snps/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz",
    },
    "Mouse": {
        "GENOME_FASTA": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa",
        "GENOME_SPLICE_MMI": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.mm2.splice.mmi",
        "TRANSCRIPTOME_FASTA": "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.fa",
        "TRANSCRIPTOME_MMI": "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.mm2.map-ont.mmi",
        "TID2GID_TSV": "/date/chenzonggui/species/mus_musculus/GRCm38.p6/transcript_id_to_gene_id.tsv",
        "TRANSCRIPT_BED": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.bed",
        "ANNOTATION_GTF": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.sorted.gtf",
        "SNP_BED_GZ": "/home/chenzonggui/species/mus_musculus/mm10/snp142.3.lite.bed.gz",
        "SNP_VCF_GZ": "/home/chenzonggui/species/mus_musculus/mgp/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz"
    }
}
