#!/bin/bash

#SBATCH -J FTAN4ALL
#SBATCH -o FTAN4ALL_%j.out
#SBATCH -e FTAN4ALL_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=5:00:00
#SBATCH --mem=MaxMemPerNode
#SBATCH --mail-user=hongda.wang@colorado.edu
#SBATCH --mail-type=ALL

cd /projects/howa1663/Code/pyaftan
python run_FTAN_4_JdF_4All.py /lustre/janus_scratch/howa1663/CC_JdF/Test_NewStack
