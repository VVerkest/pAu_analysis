#!/bin/bash
#SBATCH  --job-name=adc_count
#SBATCH  --array=0-119
#SBATCH  -o log/adc_count/%j__%a.log
#SBATCH  -e log/adc_count/%j__%a.err
#SBATCH  --time=00:30:00
#SBATCH  --mem=2G
#SBATCH   -q primary

# This is a SLURM submission script for `adc_count`
# For MB data there are 26 (--array=0-25) input files
# For HT data there are 120 HT (--array=0-119) input files


name="adc_count"
n_events=-1

in_dir="./production_pAu200_2015/HT/"
# pick one input file per array-member in SLURM
in_files=(`ls -1 $in_dir`)
root_in="${in_dir}${in_files[${SLURM_ARRAY_TASK_ID}]}"
root_out=out/${name}/${in_files[${SLURM_ARRAY_TASK_ID}]%".root"}__${SLURM_ARRAY_TASK_ID}.root
# echo $name_out

source set_envpaths.sh
echo ./bin/${name} ${root_in} ${root_out} ${n_events}
./bin/${name} ${root_in} ${root_out} ${n_events}
