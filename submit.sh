#!/bin/bash
#SBATCH --account=rrg-sorensen
#SBATCH --time=00:01:00
#SBATCH --mem=100M
#SBATCH --output=paramsFIRSTRUN.out
#SBATCH --array=1-19

./Schrodinger.exe $SLURM_ARRAY_TASK_ID


#    # SBATCH --cpus-per-task=1
#    # SBATCH --cpus-per-task=1
#    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#    ../KGLadder in0.42

