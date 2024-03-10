#!/bin/bash
gtf1=$1 # query.gtf
gtf2=$2 # reference.gtf
fasta=$3 # reference.fasta
threads=$4 # threads
outdir=$5 # output directory
name=`basename ${gtf1} .gtf`
set +u; source activate SQANTI3.env
SQANTI3_ROOT="/home/chenzonggui/software/SQANTI3-4.2"
# set +e
${SQANTI3_ROOT}/sqanti3_qc.py -d ${outdir} --cpus ${threads} --skipORF --report skip ${gtf1} ${gtf2} ${fasta} && \
${SQANTI3_ROOT}/sqanti3_RulesFilter.py --report skip ${outdir}/${name}_classification.txt ${outdir}/${name}_corrected.fasta ${outdir}/${name}_corrected.gtf
ls -1 -d ${outdir}/* | grep -v 'classification.txt$' | xargs rm -rf
# exit 0