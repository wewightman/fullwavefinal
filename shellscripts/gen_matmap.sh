#!/bin/bash
#SBATCH --mem=8GB
#SBATCH --partition=ultrasound
#SBATCH -o gen_matmap.out
#SBATCH --exclude=dcc-ultrasound-01
#SBATCH --cpus-per-task=8
#SBATCH --time=40:00:00
                        
date
hostname

echo "Activating python..."
source /hpc/group/ultrasound/wew12/modules/miniconda3/envs/pyus-310/bin/activate
echo "Running script..."
conda run -n data3d python -c "import prepostroutines as ppr; ppr.pre_genmatmap()"
echo "Exiting..."
