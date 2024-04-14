#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import os
import numpy as np
import pandas as pd

RUNS = config["RUNS_K562"] + config["RUNS_MESC"] + config["RUNS_HUMAN_BLASTOCYST"] + config["RUNS_MOUSE_BLASTOCYST"]
THREADS = 8
DAT = pd.read_excel(config["TABLE"])
DAT.index = DAT["Cell"]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]

# path = "results/mapping/all_cells_final_reads.tsv"
# if os.path.exists(path):
#     d = pd.read_csv(path, sep="\t", index_col=0)
#     d["UniqReads"] = d["Reads"] - d["DuplicateReads"]
#     d["UMIs"] = d["UniqReads"]
#     dat = dat.merge(d, left_index=True, right_index=True, how="left")

def get_species(cell):
    return DAT.loc[cell]["Species"]
def get_group(cell):
    return DAT.loc[cell]["Group"]
def get_cell_line(cell):
    return DAT.loc[cell]["CellLine"]
def get_label(cell):
    return DAT.loc[cell]["Label"]

DAT_SELECTED = DAT[DAT["Run"].isin(RUNS)]
DAT_K562 = DAT_SELECTED[DAT_SELECTED["Group"] == "K562"]
DAT_MESC = DAT_SELECTED[DAT_SELECTED["Group"] == "mESC"]
DAT_BLASTOCYST = DAT_SELECTED[DAT_SELECTED["Group"] == "MouseBlastocyst"]
RUN_CELLS_SELECTED = list(DAT_SELECTED["RunCell"])
RUN_CELLS_K562 = list(DAT_K562["RunCell"])
RUN_CELLS_MESC = list(DAT_MESC["RunCell"])
RUN_CELLS_BLASTOCYST = list(DAT_BLASTOCYST["RunCell"])
RUN_CELLS = RUN_CELLS_SELECTED
RUN_CELLS_CELLLINE = RUN_CELLS_K562 + RUN_CELLS_MESC

GROUPS = ["K562", "mESC", "MouseBlastocyst"]
TAMA_ROOT = "/lustre/grp/tfclab/chenzg/software/tama" # Root directory of TAMA

def get_group_cells(group):
    d = DAT
    d = d[(d["Group"] == group) & (d["Time"] == 3) & (d["ActD"].isna())]
    if group == "K562":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 50)]
    elif group == "mESC":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    elif group == "MouseBlastocyst":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    else:
        assert False
    return list(sorted(d["Cell"]))

def get_estimate_pe_model(cell):
    return "reports/Estimate.Pe.Model.consensus.%s.tsv" % get_cell_line(cell)

def get_chroms(cell):
    species = get_species(cell)
    if species == "Human":
        return ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
    elif species == "Mouse":
        return ["chr%d" % i for i in range(1, 20)] + ["chrX", "chrY"]

def get_genome_fasta(cell):
    return config["%s_GENOME_FASTA" % get_species(cell).upper()]

def get_genome_splice_mmi(cell):
    return config["%s_SPLICE_MMI" % get_species(cell).upper()]

# def get_transcriptome_fasta(cell):
#     species = get_species(cell)
#     if species == "Human":
#         return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.fa"
#     elif species == "Mouse":
#         return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.fa"
#     else:
#         assert False

# def get_transcriptome_mmi(cell):
#     species = get_species(cell)
#     if species == "Human":
#         return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.mm2.map-ont.mmi"
#     elif species == "Mouse":
#         return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.mm2.map-ont.mmi"
#     else:
#         assert False

# def get_tid_to_gid_tsv(cell):
#     species = get_species(cell)
#     if species == "Human":
#         return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/transcript_id_to_gene_id.tsv"
#     elif species == "Mouse":
#         return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/transcript_id_to_gene_id.tsv"
#     else:
#         assert False

def get_transcript_bed(cell):
    return config["%s_TRANSCRIPT_BED" % get_species(cell).upper()]

def get_annotation_gtf(cell):
    return config["%s_ANNOTATION_GTF_GZ" % get_species(cell).upper()]

def get_snp_bed(cell):
    return config["%s_SNP_BED_GZ" % get_species(cell).upper()]

def get_snp_vcf(cell):
    return config["%s_SNP_VCF_GZ" % get_species(cell).upper()]
