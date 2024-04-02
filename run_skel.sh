#!/bin/bash
#SBATCH -J job
#SBATCH -o out
#SBATCH -e err
#SBATCH -n 1
#SBATCH --mem 10000
#SBATCH -t 24:00:00

# submit from ~/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/

salloc --account=def-sjsmith -n 1 -t 05:00:00 --mem 6000
