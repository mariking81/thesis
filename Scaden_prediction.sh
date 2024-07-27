#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=32gb

module load anaconda3/personal
source activate attempt2

scaden predict $HOME/MSc/project/deconvolution/deconvolution_model/GSE141910_genematrix.txt --model_dir $HOME/MSc/project/deconvolution/training/

cp -r $TMPDIR $HOME/
