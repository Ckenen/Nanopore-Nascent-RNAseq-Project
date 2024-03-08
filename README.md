# Workstation of nanopore-based nascent RNA sequencing

## Sumamry

| Directory | Description |
| :-- | :-- |
| 0_Analysis | Source code of advanced analysis and plotting figures. |
| 1_NanoNASCseq | Workflow of the basic analysis of NanoNASC-seq. |
| 2_NASCseq | Workflow of the basic analysis of NASC-seq. |
| 3_NASCseq_SE | Workflow of the basic analysis of NASC-seq (single-end). |
| 4_RNAseq | Workflow of the analysis of RNA-seq. |
| 5_RNAseq_ActD | Workflow of the analysis of RNA-seq under ActD treatment. |
| 6_SCANseq2 | Workflow of the analysis of SCAN-seq2. |
| 7_scCOLORseq | Workflow of the analysis of scCOLOR-seq. |


### 2. Ribosomal RNA

    mkdir -p common/ncbi
    cd common/ncbi

    # Download ribosomal RNA and generated bowtie2 index
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna_from_genomic.fna.gz
    zcat GCF_000001405.39_GRCh38.p13_rna_from_genomic.fna.gz | grep '>' | grep "gbkey=rRNA" > all.rrna.txt
    ./make_ribosomal_rna.py
    mkdir -p bt2_rrna_index/ref
    bowtie2-build ribosomal_rna.fasta bt2_rrna_index/ref
    

### 3. Single nucleotide polymorphism (SNP)

    # cd to the directory.
    mkdir -p common/ucsc
    cd common/uscs

    # Download the SNPs annotation from UCSC.
    # snp151.txt.gz contains 683,635,300 lines, and snp151Common.txt.gz contains 15,175,044 lines.
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz
    
    # Generated a BED6 format file that represents the details information of SNPs. This will be helpful on the IGV.
    gzip -dc snp151.txt.gz | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$4,$5":"$10":"$11":"$12,".",$7}' | bgzip -c > snp151.6.bed.gz
    tabix -p bed snp151.6.bed.gz

    # Generated a BED3 format file that only contains chrom, start, and end of single base conversion. This will be helpful to annotated SNP-incuded mismatch in bam.
    gzip -dc snp151.txt.gz | awk -v FS='\t' -v OFS='\t' '{if($12=="single"){print $2,$3,$4}}' | awk '$3-$2==1' | bgzip -c > snp151.3.bed.gz
    tabix -p bed snp151.3.bed.gz


## Sequencing Data

