#!/bin/bash
#SBATCH --job-name=cavity-model
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=120:00:00
#SBATCH --partition=sbinlab_gpu
#SBATCH --mem=20G

python3 -u run_pipeline.py
