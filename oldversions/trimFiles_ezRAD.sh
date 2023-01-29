#!/bin/bash

#SBATCH --job-name=trim_ezRAD
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=USER@tamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent/2.24
module load parallel/20160722
module load trim_galore/0.4.0

#save file names into text file
echo -n "Getting filenames"
date
ls *F.fq.gz | sed 's/\F.fq\.gz//' > filenames.txt

#mkdir assembly
#mkdir mapping


#Trim for assembly
seq 1 1 > loop1.txt

echo -n "Trimming F and R reads for assembly"
date
parallel "trim_galore --paired --length 100 -q 10 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 8 -e .2 {2}F.fq.gz {2}R.fq.gz --output_dir assembly" :::: loop1.txt filenames.txt

echo -n "Trimming F reads for assembly, only required for ezRAD PE data"
date
parallel "trim_galore --length 95 -q 0 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG --stringency 8 -e 0.2 --clip_R1 5 assembly/{2}F_val_1.fq.gz --output_dir assembly" :::: loop1.txt filenames.txt

#Clean up files after trimming for assembly
mkdir assembly/trimreports
mv assembly/*report.txt assembly/trimreports
rm assembly/*val_1.fq.gz
rename F_val_1_trimmed r1 assembly/*trimmed.fq.gz
rename R_val_2 r2 assembly/*_2.fq.gz


#Trim for mapping
echo -n "Trimming F and R reads for mapping"
date
parallel "trim_galore --paired --retain_unpaired --length 50 -q 15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 5 -e .2 {2}F.fq.gz {2}R.fq.gz --output_dir mapping" :::: loop1.txt filenames.txt

mkdir mapping/trimreports
mv mapping/*report.txt mapping/trimreports
mkdir mapping/unpairedreads
mv mapping/*_unpaired* mapping/unpairedreads
rename F_val_1 R1 mapping/* 1.fq.gz
rename R_val_2 R2 mapping/* 2.fq.gz


