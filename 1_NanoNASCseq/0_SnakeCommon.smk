#!/usr/bin/env runsnakemake
import os
import numpy as np
import pandas as pd
runs = [
    ## Test
    # "20220130_Cell8",
    # "20220316_Cell40",
    # "20220413_Cell26",
    # "20220506_Cell25",
    # "20220520_Cell16",
    # "20220525_Run1",
    # "20220525_Run2",
    # "20220601_Cell22",
    # "20220607_Cell40",
    # "20220615_Cell68",
    # "20220621_R10Test",

    ## K562
    "20220719_K562ActD3h00uM", # s4U: 0uM3h, ActD: 3h
    "20220719_K562ActD3h50uM", # s4U: 50uM3h, ActD: 3h
    "20220719_K562ActD6h00uM", # s4U: 0uM3h, ActD: 6h
    "20220719_K562ActD6h50uM", # s4U: 50uM3h, ActD: 6h
    "20220719_K562R1",
    "20220719_K562R2",
    "20220719_K562R3",
    "20220719_K562R4",
    "20220729_K562R1",
    "20220729_K562R2",
    "20220729_K562R3", # Contain some mouse blastocyst cells
    "20221014_K562R1", # R10.4.1
    "20221014_K562R2", # R10.4.1

    ## mESC
    "20220818_mESCR1",
    "20220818_mESCR2",
    "20220818_mESCR3", # Merged to 20220818_mESCR3M
    "20220818_mESCR3M",
    "20220818_mESCR3P2", # Merged to 20220818_mESCR3M
    
    ## MouseBlastocyst (All R10.4.1)
    "20220719_Embryo",
    "20220729_Embryo",
    "20220729_EmbryoR10",
    "20220808_Blastocyst",
    "20220817_Embryo",
    "20220820_Embryo", # contain some mESC cells
    "20220823_MidBlast",
    "20220824_MidBlast",
    "20220831_MidBlast",
    "20220902_Blastocyst", # Merged to 20220902_BlastocystM
    "20220902_BlastocystM",
    "20220902_BlastocystS1", # Merged to 20220902_BlastocystM
    "20220903_Blastocyst", # Merged to 20220903_BlastocystM
    "20220903_BlastocystM",
    "20220903_BlastocystS1", # Merged to 20220903_BlastocystM
    "20220924_Blastocyst",
    "20221008_Blastocyst",
    "20221018_Blastocyst",
    "20221019_BlastocystS1",
    "20221031_Blastocyst",
    "20221031_Blastocyst2",
    "20221120_Blastocyst",
    "20221122_Blastocyst",
    "20221126_Blastocyst",
    "20221127_Blastocyst",
    "20221128_Blastocyst",
    "20221129_Blastocyst",
    "20221207_Blastocyst",
    "20221207_Blastocyst2",
    "20221207_Blastocyst3",
    "20221217_BlastocystC70",
    "20221217_BlastocystC73",
    "20221217_BlastocystC74",
    "20221217_BlastocystC75",
    "20221217_BlastocystC76",
    "20221218_BlastocystC64",
    "20221218_BlastocystC67",
    "20221218_BlastocystC69",
    "20221222_BlastocystC68",
    "20221222_BlastocystC71",
    "20221222_BlastocystC72",
    "20230101_BlastocystC77",
    "20230101_BlastocystC78",
    "20230103_BlastocystC79",
    "20230103_BlastocystC80",
    "20230103_BlastocystC81",
    "20230105_BlastocystC83",
    "20230108_BlastocystC82",
    "20230220_BlastocystC85",
    "20230220_BlastocystC86",
    "20230220_BlastocystC87",
    "20230220_BlastocystC88",

    ## HumanBlastocyst
    "20220912_HumanBlastR1",
    "20220912_HumanBlastR2",
]

# runs = ["20220818_mESCR1"]

BARCODE_FASTA = "data/nanopore_96_barcodes.fa"
THREADS = 8
# threads = 8
dat = pd.read_excel("data/NanoNASCseq_All.xlsx")
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

cell_lines = ["K562", "mESC", "MouseBlastocyst"]
tamadir = "/lustre/grp/tfclab/chenzg/software/tama" # Root directory of TAMA

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
    cell_line = get_cell_line(cell)
    if cell_line == "K562":
        return "data/Estimate.Pe.Model.consensus.K562.tsv"
    elif cell_line == "mESC":
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
        return "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa"
    else:
        assert False

def get_genome_splice_mmi(cell):
    species = get_species(cell)
    if species == "Human":
        return "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.mm2.splice.mmi"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/GRCm38.canonical.mm2.splice.mmi"
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
        return "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.bed"
    else:
        assert False

def get_annotation_gtf(cell):
    species = get_species(cell)
    if species == "Human":
        return "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.sorted.gtf"
    else:
        assert False

def get_snp_bed(cell):
    species = get_species(cell)
    if species == "Human":
        return "/lustre/grp/tfclab/chenzg/species/homo_sapiens/hg38/snp151.3.lite.bed.gz"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/mm10/snp142.3.lite.bed.gz"
    else:
        assert False

def get_snp_vcf(cell):
    species = get_species(cell)
    if species == "Human":
        return "../5_RNAseq_ActD/results/snps/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz"
    elif species == "Mouse":
        return "/lustre/grp/tfclab/chenzg/species/mus_musculus/mgp/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz"
    else:
        assert False

