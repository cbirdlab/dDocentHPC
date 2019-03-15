#!/bin/bash -l

#SBATCH --job-name=BAMdat
#SBATCH --time=96:00:00
#SBATCH --nodes=1

#script to count the number of contigs in each BAM file

#to run: 
# sbatch cntDepthByContig-BAM.bash 32 5.5

enable_lmod
module load samtools
module load parallel

THREADS=$1
CUTOFFS=$2

#ls *$CUTOFFS*-R[AG][W.]*bam | parallel --no-notice -k "echo -n {}' ' && samtools view {} | cut -f3 | sort | uniq | wc -l" > numContigs-BAM.stats
#NumRows=$(cat numContigs-BAM.stats | wc -l)
#NumRows2=$((NumRows/2))
#paste <(seq 1 $NumRows2) <(grep 'RAW.bam' numContigs-BAM.stats) <(grep 'RG.bam' numContigs-BAM.stats) > numContigs-BAM.dat

#ls *$CUTOFFS*-R[AG][W.]*bam | parallel --no-notice -k "samtools view {} | cut -f3 | uniq -c | sed 's/^  *//g' | tr ' ' '\t' "

ls *$CUTOFFS*-R[AG][W.]*bam | parallel --no-notice -j $THREADS "samtools view {} | cut -f3 | uniq -c | sed 's/^  *//g' | tr ' ' '\t' | sed 's/$/\t{}/g' > {}.dat"

ls *$CUTOFFS*-R[AG][W.]*bam.dat | parallel --no-notice -k -j $THREADS "wc -l {}" > numContigsBAM.$CUTOFFS.dat

ls *$CUTOFFS-R[AG][W.]*bam.dat | parallel --no-notice -k -j $THREADS "cat {}" > DepthByContigBAM.$CUTOFFS.dat

