#!/usr/bin/env runsnakemake
import os
import numpy as np
import pandas as pd
configfile: "config.yaml"
runs = config["runs"]
threads = 8
dat = pd.read_excel("data/NanoNASCseq.xls")
dat.index = dat["Cell"]
dat["RunCell"] = ["%s/%s" % (run, cell) for run, cell in dat[["Run", "Cell"]].values]

# path = "results/mapping/all_cells_final_reads.tsv"
# if os.path.exists(path):
#     d = pd.read_csv(path, sep="\t", index_col=0)
#     d["UniqReads"] = d["Reads"] - d["DuplicateReads"]
#     d["UMIs"] = d["UniqReads"]
#     dat = dat.merge(d, left_index=True, right_index=True, how="left")

def get_species(cell):
    return dat.loc[cell]["Species"]
def get_strain(cell):
    return dat.loc[cell]["Strain"]
def get_cell_line(cell):
    return dat.loc[cell]["CellLine"]
def get_label(cell):
    return dat.loc[cell]["Label"]

dat_selected = dat[[r in runs for r in dat["Run"]]]
dat_k562 = dat_selected[dat_selected["Strain"] == "K562"]
dat_mesc = dat_selected[dat_selected["Strain"] == "mESC"]
dat_blastocyst = dat_selected[["blastocyst" in x.lower() for x in dat_selected["Strain"]]]
run_cells_selected = list(dat_selected["RunCell"])
run_cells_k562 = list(dat_k562["RunCell"])
run_cells_mesc = list(dat_mesc["RunCell"])
run_cells_blastocyst = list(dat_blastocyst["RunCell"])
run_cells = run_cells_selected
run_cells_cellline = run_cells_k562 + run_cells_mesc

strains = ["K562", "mESC", "MouseBlastocyst"]
cell_lines = ["K562", "mESC", "MouseBlastocyst"]
tamadir = "/home/chenzonggui/software/tama" # Root directory of TAMA


def get_strain_cells(strain):
    import pandas as pd
    d = pd.read_excel("data/NanoNASCseq_summary_selected.xls")
    if strain == "K562":
        d = d[(d["Strain"] == "K562") & ((d["s4U"] == 0) | (d["s4U"] == 50)) & (d["Time"] == 3) & (d["ActD"].isna())]
    elif strain == "mESC":
        d = d[(d["Strain"] == "mESC") & ((d["s4U"] == 0) | (d["s4U"] == 400)) & (d["Time"] == 3) & (d["ActD"].isna())]
    elif strain == "MouseBlastocyst":
        d = d[["blast" in x.lower() for x in d["Strain"]]]
        d = d[d["Species"] == "Mouse"]
        d = d[((d["s4U"] == 0) | (d["s4U"] == 400)) & (d["Time"] == 3) & (d["ActD"].isna())]
    else:
        assert False
    return list(d["Cell"])

def get_cell_line_cells(cell_line):
    d = dat
    d = d[(d["CellLine"] == cell_line) & (d["Time"] == 3) & (d["ActD"].isna())]
    if cell_line == "K562":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 50)]
    elif cell_line == "mESC":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    elif cell_line == "MouseBlastocyst":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    else:
        assert False
    return list(sorted(d["Cell"]))

def get_estimate_pe_model(cell):
    strain = dat[dat["Cell"] == cell]["Strain"].values[0]
    if strain == "K562":
        return "data/Estimate.Pe.Model.consensus.K562.tsv"
    elif strain == "mESC":
        return "data/Estimate.Pe.Model.consensus.mESC.tsv"
    else:
        assert False

def get_chroms(cell):
    species = get_species(cell)
    if species == "Human":
        return ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
    elif species == "Mouse":
        return ["chr%d" % i for i in range(1, 20)] + ["chrX", "chrY"]
    else:
        assert False

def get_genome_fasta(cell):
    species = get_species(cell)
    if species == "Human":
        return "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa"
    else:
        assert False

def get_genome_splice_mmi(cell):
    species = get_species(cell)
    if species == "Human":
        return "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.mm2.splice.mmi"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.mm2.splice.mmi"
    else:
        assert False

def get_transcriptome_fasta(cell):
    species = get_species(cell)
    if species == "Human":
        return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.fa"
    elif species == "Mouse":
        return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.fa"
    else:
        assert False

def get_transcriptome_mmi(cell):
    species = get_species(cell)
    if species == "Human":
        return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.mm2.map-ont.mmi"
    elif species == "Mouse":
        return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.mm2.map-ont.mmi"
    else:
        assert False

def get_tid_to_gid_tsv(cell):
    species = get_species(cell)
    if species == "Human":
        return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/transcript_id_to_gene_id.tsv"
    elif species == "Mouse":
        return "/date/chenzonggui/species/mus_musculus/GRCm38.p6/transcript_id_to_gene_id.tsv"
    else:
        assert False

def get_transcript_bed(cell):
    species = get_species(cell)
    if species == "Human":
        return "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.bed"
    else:
        assert False

def get_annotation_gtf(cell):
    species = get_species(cell)
    if species == "Human":
        return "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.sorted.gtf"
    else:
        assert False

def get_snp_bed(cell):
    species = get_species(cell)
    if species == "Human":
        return "/home/chenzonggui/species/homo_sapiens/hg38/snp151.3.lite.bed.gz"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/mm10/snp142.3.lite.bed.gz"
    else:
        assert False

def get_snp_vcf(cell):
    species = get_species(cell)
    if species == "Human":
        return "../5_RNAseq_ActD/results/snps/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz"
    elif species == "Mouse":
        return "/home/chenzonggui/species/mus_musculus/mgp/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz"
    else:
        assert False

