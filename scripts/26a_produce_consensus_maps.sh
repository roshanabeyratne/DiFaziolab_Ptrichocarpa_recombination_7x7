#!/bin/bash
#PBS -N consensus_maps
#PBS -l nodes=1:ppn=1,pvmem=40gb
#PBS -t 18
#PBS -m n
#PBS -M noemail@hpc.wvu.edu
#PBS -o /group/difazio/populus/gatk-7x7-Stet14/error_files
#PBS -e /group/difazio/populus/gatk-7x7-Stet14/error_files
#PBS -q standby


# STARTTING JOB 
date
cd ${PBS_O_WORKDIR}

source /shared/software/miniconda3/etc/profile.d/conda.sh
conda activate r_3.5.1_onemap_2


# VARIABLE CREATION
chromNames=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 )

### RUNNING script
# echo ${chromNames[${PBS_ARRAYID}]}

Rscript 26b_produce_consensus_maps.R ${chromNames[${PBS_ARRAYID}]}

# FINISH JOB
conda deactivate
date
