#!/bin/bash
#PBS -N 30kb_hotspots_MEME_w_nc
#PBS -l nodes=1:ppn=1,pvmem=40gb 
#PBS -m n
#PBS -M noemail@hpc.wvu.edu
#PBS -o /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -e /group/difazio/populus/gatk-7x7-Stet14/MEME/error_files
#PBS -q standby

# STARTTING JOB 
date
cd ${PBS_O_WORKDIR}
module load genomics/bioperl

# CREATE FILE WITH HOTSPOT START-END COORDINATES
cat input_files/30kb_hotspots_bonferronni_cutoff_post_filtering.csv | cut -d ',' -f 1,2,3| tr ',' '\t'| perl -e 'while(<>){chomp; print"$_\t1\n";}' > input_files/30kb_hotspot_subseqfile.txt
cat input_files/30kb_negative_control.csv | cut -d ',' -f 1,2,3| tr ',' '\t'| perl -e 'while(<>){chomp; print"$_\t1\n";}' > input_files/30kb_negative_control_subseqfile.txt


# EXTRACT THE FASTA FORMAT OF THE SEQUENCES FOR THE ABOVE COORDINATES
perl bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference input_files/30kb_hotspot_subseqfile.txt > input_files/30kb_hotspot_regions.fasta
perl bpsubseq_file.pl soft_link_to_Stet14_V2_Male_Reference input_files/30kb_negative_control_subseqfile.txt > input_files/30kb_nc_regions.fasta


# RUNNING MEME

#module load singularity/2.5.2
#singularity shell /shared/software/containers/Meme-5.3.0.simg  
module load genomics/meme/5.3.0

# MEME W/O NEGATIVE CONTROL SET
# streme --oc output_files/30kb_hotspots --objfun de --dna --p input_files/30kb_hotspot_regions.fasta

# MEME WITH NEGATIVE CONTROL SET
streme --oc output_files/30kb_hotspots_w_nc --objfun de --dna --p input_files/30kb_hotspot_regions.fasta --n input_files/30kb_nc_regions.fasta


# FINISH JOB
date
