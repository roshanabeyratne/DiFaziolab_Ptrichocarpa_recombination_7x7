#!/bin/bash
#PBS -N 30kb_hotspots_MEME
#PBS -l nodes=1:ppn=1,pvmem=40gb
#PBS -t 1-10 
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
cat 30kb_hotspots_bonferronni_cutoff.csv| cut -d ',' -f 1,2,3| tr ',' '\t'| perl -e 'while(<>){chomp; print"$_\t1\n";}' > ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_hotspot_subseqfile.txt
cat ./30kb_hotspots_iterations/negative_ctrl_list/30kb_negative_control_${PBS_ARRAYID}.csv| cut -d ',' -f 1,2,3| tr ',' '\t'| perl -e 'while(<>){chomp; print"$_\t1\n";}' > ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_negative_subseqfile.txt
perl 24b_bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_hotspot_subseqfile.txt > ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_hotspot_regions.fasta
perl 24b_bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_negative_subseqfile.txt > ./30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_negative_regions.fasta

# RUNNING MEME
#module load singularity/2.5.2
#singularity shell /shared/software/containers/Meme-5.3.0.simg  
module load genomics/meme/5.3.0
# streme --oc 30kb_hotspots_iterations/30kb_${PBS_ARRAYID} --objfun de --dna --p 30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_hotspot_regions.fasta --n 30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_negative_regions.fasta
streme --oc 30kb_hotspots_iterations/30kb_${PBS_ARRAYID} --objfun de --dna --p 30kb_hotspots_iterations/30kb_${PBS_ARRAYID}/30kb_negative_regions.fasta

# FINISH JOB
date
