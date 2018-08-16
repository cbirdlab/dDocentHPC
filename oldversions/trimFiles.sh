#!/bin/bash

#rename files for ddocent
#rename R1_001 F *R1_001.fastq.gz   
#rename R2_001 R *R2_001.fastq.gz

#SBATCH --job-name=trimAckiss
#SBATCH --time=96:00:00
#SBATCH -p cbirdq
#SBATCH --nodes=1
#SBATCH --mail-user=cbird@tamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent/2.24
module load parallel/20160722
module load cutadapt/20151012
module load fastqc/20151106
module load trim_galore/0.4.0


#srun -n 1 bash /work/GenomicSamples/cbird/Riginos/dDocent224tref /work/GenomicSamples/cbird/Riginos config.trim

#bash ddocent/2.24 /work/GenomicSamples/cbird/Riginos/config.trim

#save file names into text file
echo -n "Getting filenames"
date
ls *F.fastq.gz | sed 's/\F.fastq\.gz//' > filenames.txt
#mkdir assembly
#mkdir mapping


#Trim for assembly
#seq 1 1 > loop1.txt

#echo -n "Trimming F and R reads for assembly"
#date
#parallel "trim_galore --paired --length 101 -q 10 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 8 -e .2 {2}F.fastq.gz {2}R.fastq.gz --output_dir assembly" :::: loop1.txt filenames.txt  

#echo -n "Trimming F reads for assembly, only required for ezRAD PE data"
#date
#parallel "trim_galore --length 96 -q 0 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG --stringency 8 -e 0.2 --clip_R1 5 assembly/{2}F_val_1.fq.gz --output_dir assembly" :::: loop1.txt filenames.txt

#Clean up files after trimming for assembly
mkdir assembly/trimreports
mv assembly/*report.txt assembly/trimreports
rm assembly/*val_1.fq.gz
rename _F_val_1_trimmed .r1 assembly/*trimmed.fq.gz
rename _R_val_2 .r2 assembly/*_2.fq.gz


#Trim for mapping
#echo -n "Trimming F and R reads for mapping"
#date
#parallel "trim_galore --paired --retain_unpaired --length 50 -q 15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 5 -e .2 {2}F.fastq.gz {2}R.fastq.gz --output_dir mapping" :::: loop1.txt filenames.txt

mkdir mapping/trimreports
mv mapping/*report.txt mapping/trimreports
mkdir mapping/unpairedreads
mv mapping/*unpaired* mapping/unpairedreads
rename _F_val_1 .R1 mapping/*_1.fq.gz
rename _R_val_2 .R2 mapping/*_2.fq.gz


