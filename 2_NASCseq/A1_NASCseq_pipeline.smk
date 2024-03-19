#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

if False:
    # paired-end, 50 uM, 1 h
    dat = dat[(dat["Run"] == "GSE128273_NASCseq_K562") & (dat["s4U"] == 50) & (dat["Time"] == 1)]
    run_cells = []
    for run, cell in dat[["Run", "Cell"]].values:
        run_cells.append("%s/%s" % (run, cell))
    run_cells = run_cells[:2]
    run_cells = ["GSE128273_NASCseq_K562/SRR8724040"]
outdir = "results/reproduced"

rule all:
    input:
        #expand(outdir + "/trimmedFastqFiles/{run_cell}", run_cell=run_cells),
        #outdir + "/reference/hg38_ercc.fa",
        #outdir + "/reference/hg38_ercc.gtf",
        #outdir + "/reference/strandedness.csv",
        #outdir + "/reference/hg38_ercc.star.index",
        #expand(outdir + "/bamFiles/aligned_bam/{run_cell}", run_cell=run_cells),
        #expand(outdir + "/bamFiles/duplRemoved_bam/{run_cell}", run_cell=run_cells),
        expand(outdir + "/bamFiles/annotated_bam/{run_cell}", run_cell=run_cells),
        #expand(outdir + "/bamFiles/annotated_sorted_bam/{run_cell}", run_cell=run_cells),
        #expand(outdir + "/bamFiles/tagged_bam/{run_cell}", run_cell=run_cells),
        # outdir + "/QC/vcfFilter",
        #expand(outdir + "/bamFiles/filteredTagged_bam/{run_cell}", run_cell=run_cells),
        #expand(outdir + "/QC/errorRates/{run_cell}_ErrorRates.csv", run_cell=run_cells),
        expand(outdir + "/outfiles/pkl_files/{run_cell}_prepared.pkl", run_cell=run_cells),


def get_fastqs(wildcards):
    fq1 = "data/datasets/%s_1.fastq.gz" % wildcards.cell
    fq2 = "data/datasets/%s_2.fastq.gz" % wildcards.cell
    if not os.path.exists(fq1):
        fq1 = "results/prepare/download/GSE128273_NASCseq_K562/fastq/%s_1.fastq.gz" % wildcards.cell
        fq2 = "results/prepare/download/GSE128273_NASCseq_K562/fastq/%s_2.fastq.gz" % wildcards.cell
    return [fq1, fq2]

rule trim_galore: # Encounter error in tanglab3 and tanglab4 nodes
    input:
        fqs = lambda wildcards: get_fastqs(wildcards)
    output:
        out = directory(outdir + "/trimmedFastqFiles/{run}/{cell}")
    log:
        outdir + "/trimmedFastqFiles/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        trim_galore --nextera --length 20 -j {threads} -o {output.out} --paired --dont_gzip {input.fqs} &> {log}
        """

rule make_reference:
    input:
        hg_fa = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        hg_gtf = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf",
        ercc_fa = "data/ERCC92.fa",
        ercc_gtf = "data/ERCC92.gtf"
    output:
        fa = outdir + "/reference/hg38_ercc.fa",
        gtf = outdir + "/reference/hg38_ercc.gtf"
    shell:
        """
        cat {input.hg_fa} {input.ercc_fa} > {output.fa}
        samtools faidx {output.fa}
        cat {input.hg_gtf} {input.ercc_gtf} | sort -k1,1 -k4,4n > {output.gtf}
        """

rule make_strandness:
    input:
        gtf = outdir + "/reference/hg38_ercc.gtf"
    output:
        txt = outdir + "/reference/strandedness.csv"
    shell:
        """
        python ./scripts/make_strandedness.py {input.gtf} {output.txt}
        """

rule build_star_index:
    input:
        fa = outdir + "/reference/hg38_ercc.fa",
        gtf = outdir + "/reference/hg38_ercc.gtf"
    output:
        out = directory(outdir + "/reference/hg38_ercc.star.index")
    log:
        outdir + "/reference/hg38_ercc.star.log"
    threads:
        48
    shell:
        """
        mkdir -p {output.out}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.out} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} &> {log}
        """

rule star_mapping:
    input:
        fqs = outdir + "/trimmedFastqFiles/{run}/{cell}",
        idx = outdir + "/reference/hg38_ercc.star.index"
    output:
        out = directory(outdir + "/bamFiles/aligned_bam/{run}/{cell}")
    log:
        outdir + "/bamFiles/aligned_bam/{run}/{cell}.log"
    params:
        fq1 = outdir + "/trimmedFastqFiles/{run}/{cell}/{cell}_1_val_1.fq",
        fq2 = outdir + "/trimmedFastqFiles/{run}/{cell}/{cell}_2_val_2.fq",
        prefix = outdir + "/bamFiles/aligned_bam/{run}/{cell}/{cell}_"
    threads:
        12
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq1} {params.fq2} \
            --alignSJoverhangMin 1000 \
            --alignSJDBoverhangMin 1 \
            --bamRemoveDuplicatesType UniqueIdentical \
            --outFilterMismatchNoverReadLmax 1 \
            --outFilterMismatchNmax 10 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterMatchNminOverLread 0.66 \
            --outSAMattributes MD \
            --outSAMtype BAM SortedByCoordinate \
            --scoreDelOpen -10000 \
            --scoreInsOpen -10000 \
            --outSAMmultNmax 1 \
            --genomeLoad LoadAndKeep \
            --limitBAMsortRAM 150000000000 &> {log}
        """

# STAR --genomeLoad Remove --genomeDir results/reproduced/reference/hg38_ercc.star.index

rule mark_duplicates:
    input:
        bamdir = outdir + "/bamFiles/aligned_bam/{run}/{cell}"
    output:
        out = directory(outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}")
    log:
        outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}.log"
    params:
        inbam = outdir + "/bamFiles/aligned_bam/{run}/{cell}/{cell}_Aligned.sortedByCoord.out.bam",
        outbam = outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}/{cell}_removeDupl.bam",
        txt = outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}/{cell}_removeDuplMetrics.txt"
    threads:
        8
    shell:
        """
        mkdir -p {output}
        picard MarkDuplicates -I {params.inbam} -O {params.outbam} -M {params.txt} --REMOVE_DUPLICATES true &> {log}
        samtools index -@ {threads} {params.outbam}
        """

rule annotate:
    input:
        bamdir = outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}",
        gtf = outdir + "/reference/hg38_ercc.gtf"
    output:
        out = directory(outdir + "/bamFiles/annotated_bam/{run}/{cell}")
    log:
        outdir + "/bamFiles/annotated_bam/{run}/{cell}.log"
    params:
        inbam = outdir + "/bamFiles/duplRemoved_bam/{run}/{cell}/{cell}_removeDupl.bam",
        outbam = outdir + "/bamFiles/annotated_bam/{run}/{cell}/{cell}_removeDupl.bam.featureCounts.bam"
    shell:
        """
        mkdir -p {output.out}
        Rscript ../public/NASC-seq/scripts/RsubReadFeatureCounts.czg_revised.R {params.inbam} {input.gtf} {output.out} &> {log}
        """

rule annotate_sort:
    input:
        bamdir = outdir + "/bamFiles/annotated_bam/{run}/{cell}"
    output:
        directory(outdir + "/bamFiles/annotated_sorted_bam/{run}/{cell}")
    params:
        inbam = outdir + "/bamFiles/annotated_bam/{run}/{cell}/{cell}_removeDupl.bam.featureCounts.bam",
        outbam = outdir + "/bamFiles/annotated_sorted_bam/{run}/{cell}/{cell}_sorted.bam"
    shell:
        """
        mkdir -p {output}
        samtools sort -o {params.outbam} {params.inbam}
        samtools index {params.outbam}
        """

rule conversiontag:
    input:
        bamdir = outdir + "/bamFiles/annotated_sorted_bam/{run}/{cell}",
        txt = outdir + "/reference/strandedness.csv"
    output:
        out = directory(outdir + "/bamFiles/tagged_bam/{run}/{cell}")
    log:
        outdir + "/bamFiles/tagged_bam/{run}/{cell}.log"
    params:
        inbam = outdir + "/bamFiles/annotated_sorted_bam/{run}/{cell}/{cell}_sorted.bam",
        outbam = outdir + "/bamFiles/tagged_bam/{run}/{cell}/{cell}_PositionTagged.bam",
        txt = outdir + "/bamFiles/tagged_bam/{run}/{cell}/{cell}_PosTag.csv"
    shell:
        """
        mkdir -p {output}
        set +u; source activate py27
		python ../public/NASC-seq/scripts/ConvperPos.czg_revised.py {params.inbam} {params.outbam} {params.txt} {wildcards.cell} {input.txt} &> {log}
        """

rule vcfFilter:
    input:
        expand(outdir + "/bamFiles/tagged_bam/{run_cell}", run_cell=run_cells)
    output:
        out = directory(outdir + "/QC/vcfFilter")
    shell:
        """
        mkdir -p {output}
        Rscript ../public/NASC-seq/scripts/vcfFilter.R results/reproduced/bamFiles/tagged_bam/GSE128273_NASCseq_K562 {output} 0.4 5
        """

rule tagFilter:
    input:
        bamdir = outdir + "/bamFiles/tagged_bam/{run}/{cell}",
        #snpdir = outdir + "/QC/vcfFilter"
    output:
        out = directory(outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}")
    log:
        outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}.log"
    params:
        # txt = outdir + "/QC/vcfFilter/posfile.csv",
        txt = "../public/NASC-seq/data/posfile.chr.csv",
        inbam = outdir + "/bamFiles/tagged_bam/{run}/{cell}/{cell}_PositionTagged.bam",
        outbam = outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}/{cell}_taggedFiltered.bam"
    shell:
        """
        mkdir -p {output.out}
        set +u; source activate py27
        python2.7 ../public/NASC-seq/scripts/filter_reads_paired.czg_revised.py {params.inbam} {params.outbam} {params.txt} &> {log}
        samtools index {params.outbam}
        """

# rule cellQC:
#     input:
#         expand(outdir + "/bamFiles/annotated_bam/{run_cell}", run_cell=run_cells)
#     output:
#         out = directory(outdir + "/QC/featureCount_QC")
#     shell:
#         """
#         Rscript ../public/NASC-seq/scripts/cellQC.R results/reproduced/bamFiles/annotated_bam/GSE128273_NASCseq_K562 {output}
#         """

# rule cellFilter:
#     input:
#     output:
#     shell:

rule calculatePE:
    input:
        bamdir = outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}"
    output:
        txt1 = outdir + "/QC/errorRates/{run}/{cell}_ErrorRates.csv",
        txt2 = outdir + "/QC/p_e/{run}/{cell}_p_e.txt"
    params:
        inbam = outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}/{cell}_taggedFiltered.bam"
    shell:
        """
        set +u; source activate py27
        python2.7 ../public/NASC-seq/scripts/PECalculation_GJH.py {params.inbam} {output.txt1} {output.txt2}
        """

rule prepareData:
    input:
        bamdir = outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}",
        txt = outdir + "/QC/p_e/{run}/{cell}_p_e.txt"
    output:
        pkl = outdir + "/outfiles/pkl_files/{run}/{cell}_prepared.pkl"
    log:
        outdir + "/outfiles/pkl_files/{run}/{cell}_prepared.log"
    params:
        inbam = outdir + "/bamFiles/filteredTagged_bam/{run}/{cell}/{cell}_taggedFiltered.bam"
    shell:
        """
        python3 ../public/NASC-seq/scripts/prepare_pickles.py {params.inbam} {output.pkl} `cat {input.txt}` {wildcards.cell} &> {log}
        """

