#!/bin/bash -l
#script to count the number of contigs in each BAM file

#to run: 
# sbatch cntDepthByContig-BAM.bash 5.5

#SBATCH --job-name=${1}
#SBATCH --time=96:00:00
#SBATCH --nodes=1

enable_lmod


module load samtools
module load parallel

CUTOFFS=$1

ls *$CUTOFFS*-R[AG][W.]*bam | parallel --no-notice -k "echo -n {}' ' && samtools view {} | cut -f3 | sort | uniq | wc -l" > DepthByContig-BAM.stats
NumRows=$(seq 1 $(cat DepthByContig-BAM.dat | wc -l))
NumRows2=$((NumRows/2))
paste $(seq 1 $NumRows2) $(grep 'RAW.bam' DepthByContig-BAM.dat) $(grep 'RG.bam' DepthByContig-BAM.dat) > DepthByContig-BAM.dat