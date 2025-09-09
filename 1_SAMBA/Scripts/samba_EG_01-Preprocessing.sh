#!/bin/bash

## Set env
DIR=01_SAMBA
DIR_IN=${DIR}/Input
DIR_OUT=${DIR}/Output
DIR_LIB=reference


## Setup env
mkdir ${DIR}
cd ${DIR}
mkdir Input counts trimmedFastq fastQC


## Get Fastq QC info
module load FastQC
fastqc ${DIR_PR}/rawFastq/*.fastq.gz -o ${DIR}/fastQC


## Trim reads and filter
module load cutadapt
TRIM_FWD=tcttgtggaaaggacgaaacaccg
TRIM_REV=gttttagagctagaaatagcaagt
for FQ_IN in $( ls ${DIR_PR}/rawFastq/*_R1_001.fastq.gz )
do 
FQ_OUT=${DIR}/trimmedFastq/$( basename ${FQ_IN} | sed 's/_L001_R1_001.fastq.gz//g' | cut -d'_' -f2 - )
cutadapt -j 0 -e 0.10 -m 15 --discard-untrimmed -g ${TRIM_FWD} ${FQ_IN} -o ${FQ_OUT}_partial-trim.fq.gz &> ${FQ_OUT}_partial-trim.log
cutadapt -j 0 -e 0.10 -m 15 --discard-untrimmed -a ${TRIM_REV} ${FQ_OUT}_partial-trim.fq.gz -o ${FQ_OUT}_trim.fq.gz &> ${FQ_OUT}_trim.log
done



## Align reads
module load Bowtie
for FQ_IN in $( ls ${DIR}/trimmedFastq/*_trim.fq.gz )
do 
OUT=${DIR}/counts/$( basename ${FQ_IN} | sed 's/_trim.fq.gz//g' )
bowtie -v 0 -m 1 --threads 20 --suppress 4,5,6,7 --chunkmbs 2000 --best ${LIBRARY} -q ${FQ_IN} 2> ${OUT}.log > ${OUT}_align.bwt
cut -f3 ${OUT}_align.bwt | sort | uniq -c | gawk '{$1=$1;print $2,$1}' > ${OUT}_count.txt
done



## Create sgRNA count matrix
RESULTS=${DIR}/counts/CountTable_RT.txt
TMP=${DIR}/counts/tmp.txt
for CT in ${DIR}/counts/*_count.txt
do
if [ -s ${RESULTS} ]
then
join -1 1 -2 1 ${RESULTS} ${CT} > ${TMP}
else
cat ${CT} > ${TMP}
fi
cat ${TMP} > ${RESULTS}
done

echo sgRNA $( basename -a ${DIR}/counts/*_count.txt | sed 's/_count.txt//g' | paste -sd' ' ) > ${TMP}
cat ${RESULTS} >> ${TMP}
mv ${TMP} ${RESULTS}






