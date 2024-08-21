#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=4:mem=128gb

SRR_id= ##

wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xvf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.1.0-centos_linux64/bin
which fastq-dump
fastq-dump SRR1067$SRR_id

module load anaconda3/personal

source activate star

mkdir SRR1067${SRR_id}_quant
cd SRR1067${SRR_id}_quant

STAR --outSAMtype None \
        --sjdbGTFfile $HOME/MSc/project/annotations/Homo_sapiens.GRCh38.111.gtf  \
        --quantMode GeneCounts \
        --runThreadN 1 \
        --genomeDir $HOME/MSc/project/genome_index \
        --readFilesIn $TMPDIR/SRR1067${SRR_id}.fastq \
        --outSAMunmapped Within \
        --outSAMattributes Standard

cp -r $TMPDIR/SRR1067${SRR_id}_quant $HOME/MSc/project/GSE141910/
