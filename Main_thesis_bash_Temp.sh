#!/bin/bash -l
#SBATCH --job-name="Sim"
#SBATCH --account="u12"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lodovico.mazzei@usi.ch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=multithread

module load cray-R

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

for j in {1..100}; do
    echo "running ${j}"
    j=${j} Rscript Main_Thesis_Temp_Cluster.R &
done
wait
