#!/bin/bash
#SBATCH --job-name=collective_decision
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=cpu_intel
#SABTCH --time=120:30:00
#SABTCH --output=%j.log

# load the environment
export SLURM_CPUS_PER_NODE
export PATH=/share/apps/MATLAB/R2025a/bin:${PATH}

matlab -nodisplay -nosplash -nodesktop -logfile my.log -r "run('./flu_velocity_force_dependent.m');exit;"
