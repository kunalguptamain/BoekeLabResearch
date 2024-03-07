#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --gres=gpu:3
#SBATCH --partition=gpu4_short
#SBATCH -c4
#SBATCH --mem=300G
#SBATCH -o test_%j.out
#SBATCH -e test_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kg3163@nyu.edu

source /gpfs/data/proteomics/projects/jack/.bashrc
conda activate kunal_diffusion
python3 diffusion_chip_seq.py
