#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=32gb

module load anaconda3/personal
source activate attempt2

prediction=$HOME/MSc/project/deconvolution/deconvolution_model/GSE141910_genematrix.txt
training=$HOME/MSc/project/deconvolution/simulated/LV_py.h5ad

scaden process $training $prediction

cp -r $TMPDIR $HOME/MSc/project/deconvolution/processed
