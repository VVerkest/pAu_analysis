#!/bin/bash
#SBATCH  --job-name=bbc_check_s001
#SBATCH  --array=0-100
#SBATCH  -o out-slurm/bbc_check_v2_s001_%a
#SBATCH  -e log/bbc_check_v2_%a.err
#SBATCH  --time=00:30:00
#SBATCH  --mem=2G
#SBATCH   -q primary

# This is a SLURM submission script for `bbc_check_v2`
# For MB data there are 26 (--array=0-25) input files
# For HT data there are 120 HT (--array=0-119) input files


trigger=500001
name="bbc_check_v2"
n_events=-1

in_dir="./production_pAu200_2015/HT/"
# pick one input file per array-member in SLURM
# in_files=(`ls -1 $in_dir`)
root_in="./production_pAu200_2015_newTriggers/pAu_2015_${SLURM_ARRAY_TASK_ID}.root"
root_out=out/${name}/trig_${trigger}_${SLURM_ARRAY_TASK_ID}
# echo $name_out

source set_envpaths.sh
echo ./bin/${name} ${n_events} ${trigger} ${root_in} ${root_out} 
     ./bin/${name} ${n_events} ${trigger} ${root_in} ${root_out} 
