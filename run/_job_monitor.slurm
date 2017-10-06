#!/bin/bash
#SBATCH -D /BIGDATA/sjtu_hsun_1/
#SBATCH -J msd-monitor
#SBATCH -o _job_slurm.out
#SBATCH -e _job_slurm.err
#SBATCH -p docker_128
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclusive


cd /HOME/sjtu_hsun_1/Project/msd-server/run

hostname
./monitor.py

