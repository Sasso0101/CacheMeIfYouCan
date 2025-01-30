#!/bin/bash

#SBATCH --job-name=bfs_runtime_tests
#SBATCH --output=sout/leonardo_bfs_runtime_tests_%j.out
#SBATCH --error=sout/leonardo_bfs_runtime_tests__%j.err

# https://wiki.u-gov.it/confluence/display/SCAIUS/Booster+Section#BoosterSection-JobManagingandSLURMScheduler
#SBATCH --partition=boost_usr_prod
#SBATCH --account=IscrC_OMG-24
#SBATCH --time=00:30:00
#SBATCH --qos=boost_qos_dbg

#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB

# Import modules
module load gcc/12.2.0 cmake/3.27.7

cd /leonardo/home/userexternal/lpichett/thomas_bfs
scripts/run_runtime_tests.sh