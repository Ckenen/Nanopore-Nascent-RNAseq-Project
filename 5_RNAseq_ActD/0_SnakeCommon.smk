#!/usr/bin/env runsnakemake

samples = [
    "20221128_K562_Actd_0h_rep1", # identical total RNA
    "20221128_K562_Actd_0h_rep2",
    "20221128_K562_Actd_3h_rep1",
    "20221128_K562_Actd_3h_rep2",
    "20221128_K562_Actd_6h_rep1",
    "20221128_K562_Actd_6h_rep2",
    "20221205_K562_Actd_0h_rep1", # identical cell number (deprecated)
    "20221205_K562_Actd_0h_rep2",
    "20221205_K562_Actd_3h_rep1",
    "20221205_K562_Actd_3h_rep2",
    "20221205_K562_Actd_6h_rep1",
    "20221205_K562_Actd_6h_rep2",
]

species_list = ["human", "fly"]

FILES = {
    "human": {
        "GENOME_FASTA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.primary_assembly.genome.fa",
        "ANNOTATION_GTF": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf",
        "GENE_BED": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.genes.bed",
        "TRANSCRIPT_BED_GZ": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz",
        "ANNOTATION_TSV": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.tsv"
    },
    "fly": {
        "GENOME_FASTA": "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.dna_sm.toplevel.fa",
        "ANNOTATION_GTF": "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.gtf",
        "GENE_BED": "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.transcripts.bed",
        "TRANSCRIPT_BED_GZ": "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.transcripts.bed.gz",
        "ANNOTATION_TSV": "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.annotation.tsv"
    }
}
