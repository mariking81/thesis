#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=128:mem=256gb:ngpus=1:gpu_type=RTX6000

module load anaconda3/personal
source activate attempt2

processed_data=$HOME/MSc/project/processed/processed.h5ad

scaden train $processed_data --model_dir $HOME/MSc/project/deconvolution/training
