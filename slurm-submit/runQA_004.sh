#!/bin/bash
#SBATCH  --job-name=runQA_004
#SBATCH  --array=0-100
#SBATCH  -o out-slurm/runQA_004_%a
#SBATCH  -e log/runQA_004_%a.err
#SBATCH  --time=00:30:00
#SBATCH  --mem=2G
#SBATCH   -q primary

# This is a SLURM submission script for `bbc_check_v2`
# For MB data there are 26 (--array=0-25) input files
# For HT data there are 120 HT (--array=0-119) input files


trigger=500004
name=runQA
n_events=-1

root_in="./production_pAu200_2015_newTriggers/pAu_2015_${SLURM_ARRAY_TASK_ID}.root"
root_out="out/runQA/500004/runQA_004_${SLURM_ARRAY_TASK_ID}"
# echo $name_out

source set_envpaths.sh
echo ./bin/runQA ${n_events} ${trigger} ${root_in} ${root_out} 
     ./bin/${name} ${n_events} ${trigger} ${root_in} ${root_out} 
