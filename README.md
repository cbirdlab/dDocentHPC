## dDocentHPC: A pipeline to trim, assemble, map, and genotype reduced representation genomic data

---
### !!! _*new configuration file format `config.6.*` must be used with present version of `dDocentHPC.bash`*_ !!!
	* all fltrBAM options work now
		* remove mapped reads orphaned by fltrBAM after mapping
		* remove mapped reads with excessive soft clipping
		* filter mapped reads by an alignment score threshold that varies with read length for each and every read
	* fltrBAM speed should be increased
	* mkBAM updated to use bwa-meme and automatically adjust to the amount of ram and number of threads on your system
		* this should increase speed
		* if you have <16GB of ram, you might get a memory error

---

dDocentHPC is a hard fork of [Dr. Jon Puritz's dDocent wrapper bash script] (ddocent.com).  dDocentHPC is designed to be run without interaction and functions more like a typical unix/linux commandline program.  Settings are defined in a config file 
and then dDocentHPC is run from the commandline.  dDocentHPC also includes updated algorithms to take advantage of parallel processing. The resulting vcf files can be filtered with [fltrVCF](https://github.com/cbirdlab/fltrVCF), a separate script.

---

## dDocentHPC Functions: trimFQ, mkREF, mkBAM, fltrBAM, mkVCF

  `trimFQ` uses trimmomatic to trim fastq files for de novo reference creation (mkREF) and mapping reads to the reference (mkBAM).
		several fold speedup
		can also use `trimFQref` or `trimFQmap` to only trim fastq files for reference creation or mapping, respectively

  `mkREF` follows description for de novo reference assembly on ddocent.com .  several fold speedup in dDocentHPC

  `mkBAM` uses bwa mem to map reads to reference genome and outputs raw, unfiltered bam files.
  
  `fltrBAM` uses samtools view to filter the BAM files.  This is only enabled in PE mode, presently.
  
  `mkVCF` uses freebayes to genotype individuals or allelotype pools. By default, only SNPs and INDELS are called, not MNPs
		when freebayes calls MNPs, it causes problems downstream with filtering the vcf files with vcftools and vcflib

---
 
## Quick Start

0. [Install all dependencies here using version 2.7.8](https://anaconda.org/bioconda/ddocent)

   ```bash
   conda install -c bioconda ddocent=2.7.8
   # there are unresolved issues if you use 2.9.4
   ```
   
	* Alternatively, [Install all dependencies here](https://github.com/jpuritz/dDocent#installing) or try [here](https://www.ddocent.com/downloads/)
	* you also need [bwa-meme](https://github.com/kaist-ina/BWA-MEME#install-option-1-bioconda)
	* I suggest downloading [anaconda](https://www.anaconda.com/products/distribution) and loading all dependencies in an environment called ddocent.  See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for details on creating and managing environments.
	
1. Create a project directory of any name that has zipped FASTQ files with following naming convention: 
	`PopSamp_IndivID.F.fq.gz`
	`PopSamp_IndivID.R.fq.gz`
		
2. Clone the dDocentHPC repository to your local directory, at the same hierarchical level as the project dir
	```bash
 	git clone https://github.com/cbirdlab/dDocentHPC
	```
  
3. Either add the dDocentHPC repo to your PATH or copy the scripts and config file to the project directory
	
4. Run the functions (trimFQ, mkREF, mkBAM, fltrBAM, and mkVCF) in order, as they are listed above

	a. An example SLURM file, dDocentHPC.sbatch, is provided as a guide for how to run on an HPC
		
	b. dDocentHPC.bash trimFQ is run from the project directory and creates two additional directories names: mkREF and mkBAM
	
	c. dDocentHPC.bash mkREF should be run from inside the mkREF directory
	
		i. you will want to run mkREF several times with different "cutoff" values to identify the best combo
		
	d. mkBAM, fltrBAM, and mkVCF should be run from inside the mkBAM directory
	
5. I strongly recommend that you look at the resulting files and output to determine if things worked as you intended. I have spent countless hours doing this on my projects and have adjusted the default settings in the config file accordingly. However, new projects can throw you a curve ball and the settings may need to be changed from the default values.
		
	a. After every run, read through the output of dDocentHPC to check for errors.  On an hpc, this will be the slurm*.out file. On a linux workstation, output will be printed to screen unless you add a redirect to a file when you run the dDocentHPC.bash script
	
	b. After filtering, run FASTQC and the cntREADS.sbatch script to visualize the results of the filtering
	
	c. After making the reference genome, check the PEAR output for the proportion of overlapping reads. View the scatter plots that help you to decide the cutoff values and adjust your cutoffs accordingly. View the fasta file in an alignment editor, such as seaview. If there are a lot of blocks of very similar sequences, increase the cutoff values.
	
	d. After mapping the reads to the reference genome and filtering them, visualize a sampling of BAM files with IGV before and after the filtering.  If you don't like what you see, adjust the settings in the config file.  I also like to look directly at the BAM files using samtools view.  If you don't understand what you're looking at in a BAM file, then download the [SAM format specification](https://github.com/samtools/hts-specs).  Make sure the reads that you want to filter are being filtered.  Adjust the settings as neccessary.
	
	e. After genotyping, visualize the VCF file.  [The VCF file format specification can be downloaded from here](https://github.com/samtools/hts-specs).  I like to select a sampling of loci and cross reference the VCF against the BAM files visualized in IGV.  Questions to ask: 
	
		i. Are the positions that you think should be called when viewing the BAM file actually called in the VCF?
		
		ii. Is indvididual n genotyped correctly at position k?
		
		iii. Are there poorly mapped reads in the BAM that are causing erroneous variant calls in the VCF.
			
6. Check out [fltrVCF](https://github.com/cbirdlab/fltrVCF) to continue processing the VCF file.

---

## Example, running dDocentHPC.bash on a workstation:
	```bash
	bash dDocentHPC.bash trimFQ config.4.all > trimFQ.out
	```
---
 
## Example SLURM script, running dDocentHPC.bash on a remote HPC:
	```bash
	#!/bin/bash

	#SBATCH --job-name=trimFQ
	#SBATCH --time=96:00:00
	#SBATCH -p normal
	#SBATCH --nodes=1
	#SBATCH -o trimFQ-mkREF-%j.out

	module load ddocent

	#this is an example sbatch script to run dDocentHPC on a slurm supercomputer
	#it is assumed that the files to be trimmed, dDocentHPC.bash, and config.4.all are in the working directory

	#this will trim the fq.gz files using the settings in config.4.all
	#it is assumed that the directory you run this script from has the 
	#fq.gz files

	bash dDocentHPC.bash trimFQ config.4.all

	#this will assemble the fq.gz files in the mkREF directory that was
	#created by the trim function.  We have to cd to the correct directory
	#then execute dDocentHPC

	cd mkREF
	bash ../dDocentHPC.bash mkREF ../config.4.all
	```

---

## Citation

If you use `dDocentHPC` in a publication, please cite the following sources:

Biesack, E. E., Dang, B. T., Ackiss, A. S., Bird, C. E., Chheng, P., Phounvisouk, L., Truong, O.T. & Carpenter, K. E. (2020). Evidence for population genetic structure in two exploited Mekong River fishes across a natural riverine barrier. Journal of Fish Biology, 1-12. 

Bird, C.E. (InsertYearCloned) dDocentHPC. A pipeline to trim, assemble, map, and genotype reduced representation genomic data.  https://github.com/cbirdlab/dDocentHPC

Puritz, J., Hollenbeck, C. M. & Gold, J. R. (2014). dDocent: a RADseq, variant-calling pipeline designed for population genomics of non-model organisms. PeerJ, 2, e314v1.

