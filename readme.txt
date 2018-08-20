dDocentHPC.bash [function] [config file]   -   a wrapper script to process GBS and RAD data

dDocentHPC.bash is a hard fork of Dr. Jon Puritz's dDocent wrapper bash script (ddocent.com).  dDocentHPC is designed to be run
without interaction and functions more like a typical unix/linux command line program.  Settings are defined in a config file and then
dDocentHPC is run from the commandline.  dDocentHPC also includes updated algorithms to take advantage of parallel processing.
The resulting vcf files can be filtered with fltrVCF, a separate script.

dDocentHPC Functions: trimFQ, mkREF, mkBAM, fltrBAM, mkVCF

  trimFQ uses trimmomatic to trim fastq files for de novo reference creation (mkREF) and mapping reads to the reference (mkBAM)

  mkREF follows description for de novo reference assembly on ddocent.com

  mkBAM uses bwa mem to map reads to reference genome and outputs raw, unfiltered bam files.
  
  fltrBAM uses samtools view to filter the BAM files.  This is only enabled in PE mode, presently.
  
  mkVCF uses freebayes to genotype individuals or allelotype pools.

Example:
  bash dDocentHPC.bash trimFQ config.4.all
  


