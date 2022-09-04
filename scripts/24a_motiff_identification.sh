#!/bin/bash
#PBS -N 20kb_CO_MEME
#PBS -l nodes=1:ppn=1,pvmem=40gb
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
cat co.detailed.df.csv| cut -d ',' -f 1,4,5,6 | tr ',' '\t' | perl -e 'while(<>){chomp; if($_!~/^CHROM/){@var=split; if($var[1] <=  20000 && $var[1] !=0){print "$var[0]\t$var[2]\t$var[3]\t1\n";}}}' > 20kb_subseqfile_detailed_CO.txt
perl 24b_bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference 20kb_subseqfile_detailed_CO.txt > 20kb_CO_regions_7x7.fasta

# RUNNING MEME
#module load singularity/2.5.2
#singularity shell /shared/software/containers/Meme-5.3.0.simg  
module load genomics/meme/5.3.0
streme --oc 20kb --objfun de --dna --p 20kb_CO_regions_7x7.fasta

# FINISH JOB
date
