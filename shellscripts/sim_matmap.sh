#!/bin/bash
#SBATCH --mem=16GB
#SBATCH --partition=ultrasound
#SBATCH -o sim_matmap.out
#SBATCH --exclude=dcc-ultrasound-01
#SBATCH --cpus-per-task=4
#SBATCH --time=40:00:00
                        
date
hostname
echo "Setting up working directory for task..."
echo "Load matlab"
module load Matlab/R2022b
echo "Running scripts"
matlab -nodisplay -nosplash -singleCompThread -r "runfullwave;exit;"
echo "Exiting..."
