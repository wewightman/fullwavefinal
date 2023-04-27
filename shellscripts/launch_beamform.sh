#!/bin/bash
#SBATCH --mem=8GB
#SBATCH --partition=ultrasound
#SBATCH -o /work/wew12/fullwave/launch_beamform.out
#SBATCH --exclude=dcc-ultrasound-01
#SBATCH --cpus-per-task=8
#SBATCH --time=40:00:00
                        
date
hostname

echo "Running script..."
conda run -n data3d python -c "import launch; launch.launch_beamform()"
echo "Exiting..."

date
