#!/usr/bin/env runsnakemake
import pandas as pd
infos = pd.read_csv("data/SraRunTable.SCANseq2.csv")
srr_list = list(infos["Run"])

table = pd.read_excel("data/Supplementary_Table_S1_Summary_of_sequenced_cells.xlsx")
runs = ["UMI_100", "UMI_200", "9CL", "4CL"]
table_selected = table[[lib in runs for lib in table["Library"]]]
# print(table)
# print(table.columns)

run_cells = []
for lib, name in table_selected[["Library", "Rename"]].values:
    run_cells.append("%s/%s" % (lib, name))
print("Cells:", len(run_cells))

barcodes = dict()
from Bio import SeqIO
for read in SeqIO.parse("../nanopore_nascent_rna_seq/data/20210831_barcode_corrected_bar80.fasta", "fasta"):
    seq1 = str(read.seq)
    seq2 = str(read.seq.reverse_complement())
    barcodes[read.name] = [seq1, seq2]
# print(barcodes)

def get_species(cell):
    return table[table["Rename"] == cell]["Organism"].values[0]

FILES = {
    "Human": {
        "GENOME_FASTA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        "GENOME_SPLICE_MMI": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.mm2.splice.mmi",
        "TRANSCRIPT_BED": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed",
        "ANNOTATION_GTF": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf",
    },
    "Mouse": {
        "GENOME_FASTA": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa",
        "GENOME_SPLICE_MMI": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.mm2.splice.mmi",
        "TRANSCRIPT_BED": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.transcripts.bed",
        "ANNOTATION_GTF": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/gencode.vM25.annotation.sorted.gtf",
    }
}


