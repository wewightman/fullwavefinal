#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --partition=ultrasound
#SBATCH -o sim_beamform.out
#SBATCH --exclude=dcc-ultrasound-01
#SBATCH --cpus-per-task=8
#SBATCH --time=40:00:00
                        
date
hostname

echo "Running script..."
conda run -n data3d python -c "import prepostroutines as ppr; ppr.postbeamform()"
echo "Exiting..."

date
