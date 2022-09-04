#!/bin/bash
#PBS -N 10kb_neg_CO_MEME_iter
#PBS -l nodes=1:ppn=1,pvmem=40gb
#PBS -t 5-10 
#PBS -m n
#PBS -M noemail@hpc.wvu.edu
#PBS -o /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -e /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -q genomics_core


# STARTTING JOB 
date
cd ${PBS_O_WORKDIR}
module load genomics/bioperl

# INPUT FILES CREATION
#cat filtered_CO_le_10k_df.csv | cut -d ',' -f 1,4,5,6 | tr ',' '\t' | perl -e 'while(<>){chomp; if($_!~/^CHROM/){@var=split; print "$var[0]\t$var[2]\t$var[3]\t1\n";}}' > 10kb_subseqfile_CO.txt
#perl 24b_bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference 10kb_subseqfile_CO.txt > filtered_CO_le_10k.fasta

cat neg_con_CO_le_10k_df${PBS_ARRAYID}.csv | cut -d ',' -f 1,2,3 | tr ',' '\t' | perl -e 'while(<>){chomp; if($_!~/^CHROM/){@var=split; print "$var[0]\t$var[1]\t$var[2]\t1\n";}}' > 10kb_neg_con_subseqfile_CO${PBS_ARRAYID}.txt
perl 24b_bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference 10kb_neg_con_subseqfile_CO${PBS_ARRAYID}.txt > neg_con_CO_le_10k${PBS_ARRAYID}.fasta

# RUNNING MEME
#module load singularity/2.5.2
#singularity shell /shared/software/containers/Meme-5.3.0.simg  
module load genomics/meme/5.3.0
streme --oc le_10kb${PBS_ARRAYID} --objfun de --dna --p filtered_CO_le_10k.fasta --n neg_con_CO_le_10k${PBS_ARRAYID}.fasta

# FINISH JOB
date
