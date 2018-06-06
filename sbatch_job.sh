#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -p 1day
#SBATCH -N 4
#SBATCH -n 50
#SBATCH -c 2
#SBATCH -C avx2
#SBATCH --mem 60G
#SBATCH -o ziqi.out
#SBATCH -e ziqi.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ziqi.deng@embl.de
#SBATCH --job-name=Algo
#SBATCH -C net10G
#SBATCH --switches=1


module load pymol Biopython OpenBabel Rosetta
time -p python main_final.py chk1 2ydi 2phk peptide_candidate_test