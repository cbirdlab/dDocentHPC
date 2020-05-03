#!/bin/bash -l

# this script will remove the first four bp from fwd and reverse reads, to remove ezRAD overhangs

THREADS=16

fqCLIPPER(){
fileNAME=$1
paste <(zcat $fileNAME.gz | paste - - - - | cut -f1) \
  <(zcat $fileNAME.gz | paste - - - - | cut -f2 | cut -c5-) \
  <(zcat $fileNAME.gz | paste - - - - | cut -f3) \
  <(zcat $fileNAME.gz | paste - - - - | cut -f4 | cut -c5-) | tr "\t" "\n" | \
  gzip > $fileNAME.gz2
}
export -f fqCLIPPER
parallel --record-env
ls *fq.gz | sed 's/\.gz//' | parallel --env _ --no-notice -j $THREADS "fqCLIPPER {}"
