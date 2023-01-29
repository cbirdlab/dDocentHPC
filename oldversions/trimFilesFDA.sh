#!/bin/bash

#rename files for ddocent
#rename R1_001 F *R1_001.fastq.gz   
#rename R2_001 R *R2_001.fastq.gz

#SBATCH --job-name=trimFDA
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=cbird@tamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent/2.24
module load parallel/20160722
module load cutadapt/20151012
module load fastqc/20151106
module load trim_galore/0.4.0


#save file names into text file
echo `date` " Getting filenames"
date
ls *F.fq.gz | sed 's/\.F\.fq\.gz//' > filenames.txt
mkdir assembly
mkdir mapping


#Trim for assembly

echo `date` " Trimming F and R reads for assembly. Read Length 146bp "
cat filenames.txt | parallel "trim_galore --paired --length 146 -q 10 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 8 -e .2 {}.F.fq.gz {}.R.fq.gz --output_dir assembly"  
rename F_val_1 r1 assembly/*_1.fq.gz
rename R_val_2 r2 assembly/*_2.fq.gz

#echo `date` " Trimming F reads for assembly, only required for ezRAD PE data"
#cat filenames.txt | parallel "trim_galore --length 141 -q 0 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG --stringency 8 -e 0.2 --clip_R1 5 assembly/{}F_val_1.fq.gz --output_dir assembly"

#Clean up files after trimming for assembly
mkdir assembly/trimreports
mv assembly/*report.txt assembly/trimreports
#rm assembly/*val_1.fq.gz
#rename _F_val_1_trimmed .r1 assembly/*trimmed.fq.gz

#Trim for mapping
echo `date` " Trimming F and R reads for mapping. Read length 146bp"
cat filenames.txt | parallel "trim_galore --paired --retain_unpaired --length 50 -q 15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 5 -e .2 {}.F.fq.gz {}.R.fq.gz --output_dir mapping"

mkdir mapping/trimreports
mv mapping/*report.txt mapping/trimreports
mkdir mapping/unpairedreads
mv mapping/*unpaired* mapping/unpairedreads
rename F_val_1 R1 mapping/*_1.fq.gz
rename R_val_2 R2 mapping/*_2.fq.gz


