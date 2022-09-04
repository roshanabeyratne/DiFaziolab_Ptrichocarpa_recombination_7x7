#!/bin/bash
#PBS -N 10kb_neg_CO_MEME_iter
#PBS -l nodes=1:ppn=1,pvmem=40gb
#PBS -t 3-10 
#PBS -m n
#PBS -M noemail@hpc.wvu.edu
#PBS -o /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -e /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -q genomics_core


# STARTTING JOB 
date
cd ${PBS_O_WORKDIR}
first_num=${PBS_ARRAYID}
second_num=`echo ${PBS_ARRAYID}-1| bc`

# RUNNING MEME
module load genomics/meme/5.3.0
streme --oc le_10kb_neg_neg${first_num} --objfun de --dna --p neg_con_CO_le_10k${first_num}.fasta --n neg_con_CO_le_10k${second_num}.fasta

# FINISH JOB
date
