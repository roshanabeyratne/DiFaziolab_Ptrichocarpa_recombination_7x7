#!/bin/bash
#PBS -N VQSR_combine_allel_contrib
#PBS -l nodes=1:ppn=1,pvmem=40gb
#PBS -t 1-266
#PBS -m n
#PBS -M noemail@hpc.wvu.edu
#PBS -o /group/difazio/populus/gatk-7x7-Stet14/error_files
#PBS -e /group/difazio/populus/gatk-7x7-Stet14/error_files
#PBS -q standby


# STARTTING JOB 
date
cd ${PBS_O_WORKDIR}
module load lang/R/3.5.2

# VARIABLE CREATION
Nchrom=19
Token=$((${PBS_ARRAYID}-1))
parIdx=$((${Token}/${Nchrom}))
chromIdx=$((${Token}%${Nchrom}))

chromNames=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 )
parentName=( 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073 )

### RUNNING script
# echo ${chromNames[${chromIdx}]} ${parentName[${parIdx}]}

Rscript 13b_VQSR_combine_focal_parent_allele_contribution_across_fullsib_famillies.R ${chromNames[${chromIdx}]} ${parentName[${parIdx}]}

# FINISH JOB
date
