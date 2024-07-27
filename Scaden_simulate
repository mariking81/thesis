#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=4:mem=256gb

module load anaconda3/personal
source activate attempt2

scaden simulate --cells 500 --n_samples 2000 --data $HOME/MSc/project/deconvolution/scaden_sim --data-format 'h5ad' --pattern "*.h5ad"  --out $HOME/MSc/project/deconvolution/simulated/
