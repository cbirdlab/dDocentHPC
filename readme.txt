dDocentHPC.bash [function] [config file]   -   a wrapper script to process GBS and RAD data

dDocentHPC.bash is a hard fork of Dr. Jon Puritz's dDocent wrapper bash script (ddocent.com).  dDocentHPC is designed to be run
without interaction and functions more like a typical unix/linux command line program.  Settings are defined in a config file 
and then dDocentHPC is run from the commandline.  dDocentHPC also includes updated algorithms to take advantage of parallel 
processing. The resulting vcf files can be filtered with fltrVCF, a separate script.


dDocentHPC Functions: trimFQ, mkREF, mkBAM, fltrBAM, mkVCF

  trimFQ uses trimmomatic to trim fastq files for de novo reference creation (mkREF) and mapping reads to the reference (mkBAM). several fold speedup over 2.2.12

  mkREF follows description for de novo reference assembly on ddocent.com .  several fold speedup over 2.2.12

  mkBAM uses bwa mem to map reads to reference genome and outputs raw, unfiltered bam files.
  
  fltrBAM uses samtools view to filter the BAM files.  This is only enabled in PE mode, presently.
  
  mkVCF uses freebayes to genotype individuals or allelotype pools.

 
Quick Start
	1. Create a project directory of any name that has zipped FASTQ files with following naming convention: 
		PopSamp_IndivID.F.fq.gz
		PopSamp_IndivID.R.fq.gz
	2. Copy the scripts and config file to the project directory
	3.	Run the functions in order, as they are listed above
		a. trimFQ is run from the project directory and creates two additional directories names: mkREF and mkBAM
		b. mkREF should be run from inside the mkREF directory
			i. you will want to run mkREF several times with different cutoff values to identify the best combo
		c. mkBAM, fltrBAM, and mkVCF should be run from inside the mkBAM directory


Example, running dDocentHPC.bash on a workstation:

	bash dDocentHPC.bash trimFQ config.4.all

 
Example SLURM script, running dDocentHPC.bash on a remote HPC:

	#!/bin/bash

	#SBATCH --job-name=trim_ref
	#SBATCH --time=48:00:00
	#SBATCH -p normal
	#SBATCH --nodes=1

	module load ddocent

	#this is an example sbatch script to run dDocentHPC on a slurm supercomputer

	#this will trim the fq.gz files using the settings in config.4.all
	#it is assumed that the directory you run this script from has the 
	#fq.gz files

	bash dDocentHPC.bash trimFQ config.4.all

	#this will assemble the fq.gz files in the mkREF directory that was
	#created by the trim function.  We have to cd to the correct directory
	#then execute dDocentHPC

	cd mkREF
	bash ../dDocentHPC.bash mkREF ../config.4.all

