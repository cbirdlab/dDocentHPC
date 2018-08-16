#!/usr/bin/env bash

##########dDocent 2.24##########
# v7.0 Jason updated v6.6individuals
# v7.1 RMH updated v7.0; passed all new variables in config file to the function; reordered them to match config file to preserve my sanity

#This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD (SE or PE), ezRAD (SE or PE) data, or SE RAD data.
#It requires that your raw data are split up by tagged individual and follow the naming convention of:

#Pop_Sample1.F.fq and Pop_Sample1.R.fq

#Prints out title and contact info
echo
echo -e "\n dDocent 2.2.16 Hacked for HPC by cbird@tamucc.edu \n"
#echo -e "Contact jpuritz@gmail.com with any problems \n\n "


#dDocent can now accept a configuration file instead of running interactively
#Checks if a configuration file is being used, if not asks for user input
#ceb: I moved the loading of the config file here because i was able to change the config file
#and run the program faster than it took to load the config settings
if [ -n "$1" ]; then
	echo " "
	echo -n "Files output to: "
	pwd
	echo " "
	echo "****************************************************************************************************"
	echo -e "Reading config file \n"
	cat $1
	echo "****************************************************************************************************"
	
	CONFIG=$1
	NUMProc=$(grep -A1 Processor $CONFIG | tail -1)
	MAXMemory=$(grep -A1 Memory $CONFIG | tail -1)
	TRIM=$(grep -A1 '^Trim' $CONFIG | tail -1)
	TRIM_LENGTH_ASSEMBLY=$(grep -A1 TRIM_LENGTH_ASSEMBLY $CONFIG | tail -1)
	SEED_ASSEMBLY=$(grep -A1 SEED_ASSEMBLY $CONFIG | tail -1)
	PALIMDROME_ASSEMBLY=$(grep -A1 PALIMDROME_ASSEMBLY $CONFIG | tail -1)
	SIMPLE_ASSEMBLY=$(grep -A1 SIMPLE_ASSEMBLY $CONFIG | tail -1)
	windowSize_ASSEMBLY=$(grep -A1 windowSize_ASSEMBLY $CONFIG | tail -1)
	windowQuality_ASSEMBLY=$(grep -A1 windowQuality_ASSEMBLY $CONFIG | tail -1)
	TRAILING_ASSEMBLY=$(grep -A1 TRAILING_ASSEMBLY $CONFIG | tail -1)
	TRIM_LENGTH_MAPPING=$(grep -A1 TRIM_LENGTH_MAPPING $CONFIG | tail -1)
	LEADING_MAPPING=$(grep -A1 LEADING_MAPPING $CONFIG | tail -1)
	TRAILING_MAPPING=$(grep -A1 TRAILING_MAPPING $CONFIG | tail -1)
	FixStacks=$(grep -A1 FixStacks $CONFIG | tail -1)
	ASSEMBLY=$(grep -A1 '^Assembly' $CONFIG | tail -1)
	ATYPE=$(grep -A1 Type $CONFIG | tail -1)
	simC=$(grep -A1 Clustering_Similarity $CONFIG | tail -1)
	HPC=$(grep -A1 '^HPC' $CONFIG | tail -1)
	MANCUTOFF=$(grep -A1 Manually $CONFIG | tail -1)
	if [ "$MANCUTOFF" = "no" ]; then
		CUTOFF=$(grep -A1 ^Cutoff1 $CONFIG | tail -1)
		CUTOFF2=$(grep -A1 ^Cutoff2 $CONFIG | tail -1)
		echo ""; echo `date` " CUTOFF=" $CUTOFF
		echo ""; echo `date` " CUTOFF2=" $CUTOFF2
	else
		echo ""; echo `date` " This is the HPC version of dDocent.  Manual cutoffs are not possible.  Please change config file settings."
		echo ""; echo `date` " Aborting dDocent 4 HPC"
		exit
	fi
	MAP=$(grep -A1 Mapping_Reads $CONFIG | tail -1)
	optA=$(grep -A1 _Match $CONFIG | tail -1)
	optB=$(grep -A1 MisMatch $CONFIG | tail -1)
	optO=$(grep -A1 Gap $CONFIG | tail -1)
	MAPPING_MIN_ALIGNMENT_SCORE=$(grep -A1 '^Mapping_Minimum_Alignment_Score' $CONFIG | tail -1)
	MAPPING_CLIPPING_PENALTY=$(grep -A1 '^Mapping_Clipping_Penalty' $CONFIG | tail -1)
	MAPPING_MIN_QUALITY=$(grep -A1 '^Mapping_Min_Quality' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F=$(grep -A1 '^Filter_BAM_Files' $CONFIG | tail -1)
	SNP=$(grep -A1 SNP $CONFIG | tail -1)
	POOLS=$(grep -A1 'Is the data pooled' $CONFIG | tail -1)
	POOL_PLOIDY_FILE=$(grep -A1 '^If the data is pooled, will you provide a copy number variation' $CONFIG | tail -1)
	PLOIDY=$(grep -A1 '^If no cnv file is provided, then what is the ploidy of the samples' $CONFIG | tail -1)
	BEST_N_ALLELES=$(grep -A1 '^Use_Best_N_Alleles' $CONFIG | tail -1)
	MIN_MAPPING_QUAL=$(grep -A1 '^Minimum_Mapping_Quality' $CONFIG | tail -1)
	MIN_BASE_QUAL=$(grep -A1 '^Minimum_Base_Quality' $CONFIG | tail -1)
	HAPLOTYPE_LENGTH=$(grep -A1 '^Haplotype_Length' $CONFIG | tail -1)
	MIN_REPEAT_ENTROPY=$(grep -A1 '^Min_Repeat_Entropy' $CONFIG | tail -1)
	MIN_COVERAGE=$(grep -A1 '^Min_Coverage' $CONFIG | tail -1)
	MIN_ALT_FRACTION=$(grep -A1 '^Min_Alternate_Fraction' $CONFIG | tail -1)
	MIN_ALT_COUNT=$(grep -A1 '^Min_Alternate_Count' $CONFIG | tail -1)
	MIN_ALT_TOTAL=$(grep -A1 '^Min_Alternate_Total' $CONFIG | tail -1)
	READ_MAX_MISMATCH_FRACTION=$(grep -A1 '^Read_Max_Mismatch_Fraction' $CONFIG | tail -1)
	MAIL=$(grep -A1 Email $CONFIG | tail -1)


else
	#GetInfo 
	echo ""; echo `date` " This is the HPC version of dDocent.  A config file must be specified"
	echo ""; echo `date` " Aborting dDocent 4 HPC"
	exit
fi





###Code to check for the required software for dDocent

echo "Checking for required software"
DEP=(freebayes mawk bwa samtools vcftools rainbow gnuplot gawk seqtk cd-hit-est bamToBed bedtools coverageBed parallel vcfcombine bamtools pearRM)
NUMDEP=0
for i in "${DEP[@]}"
do
	if which $i &> /dev/null; then
		foo=0
	else
    		echo `date` " The dependency" $i "is not installed or is not in your" '$PATH'"."
    		NUMDEP=$((NUMDEP + 1))
	fi
done

if find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null| grep -q 'trim' ; then
	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
else
    echo `date` " The dependency trimmomatic is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
fi
	
if find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | grep -q 'Tru' ; then
	ADAPTERS=$(find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | head -1)
else
    echo "The file listing adapters (included with trimmomatic) is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
fi

SAMV1=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c1)
SAMV2=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c3)
	if [ "$SAMV1"  -ge "1" ]; then
		if [ "$SAMV2"  -lt "3" ]; then
        	echo "The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 1.3.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
		fi
	
	else
		    echo "The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 1.3.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
	fi

RAINV=(`rainbow | head -1 | cut -f2 -d' ' `)	
	if [[ "$RAINV" != "2.0.2" && "$RAINV" != "2.0.3" && "$RAINV" != "2.0.4" ]]; then
		echo "The version of Rainbow installed in your" '$PATH' "is not optimized for dDocent."
		echo -en "\007"
		echo -en "\007"
		echo -en "\007"
		echo "Is the version of rainbow installed newer than 2.0.2?  Enter yes or no."
		read TEST
		if [ "$TEST" != "yes" ]; then 
			echo "Please install a version newer than 2.0.2"
			exit 1
		fi
    fi
FREEB=(`freebayes | grep -oh 'v[0-9].*' | cut -f1 -d "." | sed 's/v//' `)	
	if [ "$FREEB" != "1" ]; then
		echo "The version of FreeBayes installed in your" '$PATH' "is not optimized for dDocent."
		echo "Please install at least version 1.0.0"
		exit 1
    fi         	
VCFTV=$( vcftools | grep VCF | sed -e 's/VCFtools (v//' | sed -e 's/)//' )   
#	if [ "$VCFTV" != "0.1.11" ]; then
	# if [ "$VCFTV" != "0.1.11" || "$VCFTV" != "(0.1.14" ]; then
	 if [ "$VCFTV" != "0.1.11" && "$VCFTV" != "0.1.14" ]; then #RMH...still says missing `]' for some reason
	echo "The version of VCFtools installed in your" '$PATH' "is not optimized for dDocent."
		echo "Please install only version 0.1.11"
		exit 1 
	fi
BWAV=$( bwa 2>&1 | mawk '/Versi/' | sed 's/Version: //g' | sed 's/0.7.//g' | sed 's/-.*//g' | cut -c 1-2 )
	if [ "$BWAV" -lt "13" ]; then
        	echo "The version of bwa installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 0.7.13"
        	exit 1
	fi
# BTC=$( bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." )
	# if [ "$BTC" != "2.23" || "$BTC" != "2.26" ]; then
		# echo "The version of bedtools installed in your" '$PATH' "is not optimized for dDocent."
		# echo "Please install only version 2.23.0"
		# exit 1
	# fi
BTC=$( bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." | sed 's/2\.//g' )
	if [ "$BTC" -ge "26" ]; then
		BEDTOOLSFLAG="NEW"
		elif [ "$BTC" == "23" ]; then
		BEDTOOLSFLAG="OLD"
		elif [ "$BTC" != "23" ]; then
		echo "The version of bedtools installed in your" '$PATH' "is not optimized for dDocent."
		echo "Please install version 2.23.0 or version 2.26.0 and above"
		exit 1	
	fi
	
if ! awk --version | fgrep -v GNU &>/dev/null; then
	 awk=gawk
else
	 awk=awk
fi


if [ $NUMDEP -gt 0 ]; then
	echo -e "\nPlease install all required software before running dDocent again."
	exit 1
else
	echo -e "\nAll required software is installed!"
fi

#This code checks for individual fastq files follow the correct naming convention and are gziped
TEST=$(ls *.r1.fq 2> /dev/null | wc -l )
echo TEST=$TEST
#TEST2=$(ls ref.*.fq 2> /dev/null | wc -l )
#echo TEST2=$TEST2
#TEST=$(($TEST - $TEST2))
#echo TEST-TEST2=$TEST
if [ "$TEST" -gt 0 ]; then
	echo -e "\ndDocent is now configured to work on compressed sequence files.  Please run gzip to compress your files."
	echo "This is as simple as 'gzip *.fq'"
	echo "Please rerun dDocent after compressing files."
	echo $TEST
	exit 1
fi



###############################################################################################
#Wrapper for main program functions.  This allows the entire file to be read first before execution
main(){
	echo "";echo `date` " Begin main ddocent function"
	##########User Input Section##########
	#This code gets input from the user and assigns variables
	######################################
	#RMH added ALL variables
	echo "NUMProc=$1"
	echo "MAXMemory=$2"
	echo "TRIM=$3"
	echo "TRIM_LENGTH_ASSEMBLY=$4"
	echo "SEED_ASSEMBLY=$5"
	echo "PALIMDROME_ASSEMBLY=$6"
	echo "SIMPLE_ASSEMBLY=$7"
	echo "windowSize_ASSEMBLY=$8"
	echo "windowQuality_ASSEMBLY=$9"
	echo "TRAILING_ASSEMBLY=${10}"
	echo "TRIM_LENGTH_MAPPING=${11}"
	echo "LEADING_MAPPING=${12}"
	echo "TRAILING_MAPPING=${13}"
	echo "FixStacks=${14}"
	echo "ASSEMBLY=${15}"
	echo "ATYPE=${16}"
	echo "simC=${17}"
	echo "HPC=${18}"
	echo "MANCUTOFF=${19}"
	echo "CUTOFF=${20}"
	echo "CUTOFF2=${21}"
	echo "MAP=${22}"
	echo "optA=${23}"
	echo "optB=${24}"
	echo "optO=${25}"
	echo "MAPPING_MIN_ALIGNMENT_SCORE=${26}"
	echo "MAPPING_CLIPPING_PENALTY=${27}"
	echo "MAPPING_MIN_QUALITY=${28}"
	echo "SAMTOOLS_VIEW_F=${29}"
	echo "SNP=${30}"
	echo "POOLS=${31}"
	echo "POOL_PLOIDY_FILE=${32}"
	echo "PLOIDY=${33}"
	echo "BEST_N_ALLELES=${34}"
	echo "MIN_MAPPING_QUAL=${35}"
	echo "MIN_BASE_QUAL=${36}"
	echo "HAPLOTYPE_LENGTH=${37}"
	echo "MIN_REPEAT_ENTROPY=${38}"
	echo "MIN_COVERAGE=${39}"
	echo "MIN_ALT_FRACTION=${40}"
	echo "MIN_ALT_COUNT=${41}"
	echo "MIN_ALT_TOTAL=${42}"
	echo "READ_MAX_MISMATCH_FRACTION=${43}"
	echo "MAIL=${44}"
	
	#Load config arguments into variables
	# RMH added ALL variables
	NUMProc=$1
	MAXMemory=$2
	TRIM=$3
	TRIM_LENGTH_ASSEMBLY=$4
	SEED_ASSEMBLY=$5
	PALIMDROME_ASSEMBLY=$6
	SIMPLE_ASSEMBLY=$7
	windowSize_ASSEMBLY=$8
	windowQuality_ASSEMBLY=$9
	TRAILING_ASSEMBLY=${10}
	TRIM_LENGTH_MAPPING=${11}
	LEADING_MAPPING=${12}
	TRAILING_MAPPING=${13}
	FixStacks=${14}
	ASSEMBLY=${15}
	ATYPE=${16}
	simC=${17}
	HPC=${18}
	MANCUTOFF=${19}
	CUTOFF=${20}
	CUTOFF2=${21}
	MAP=${22}
	optA=${23}
	optB=${24}
	optO=${25}
	MAPPING_MIN_ALIGNMENT_SCORE=${26}
	MAPPING_CLIPPING_PENALTY=${27}
	MAPPING_MIN_QUALITY=${28}
	SAMTOOLS_VIEW_F=${29}
	SNP=${30}
	POOLS=${31}
	POOL_PLOIDY_FILE=${32}
	PLOIDY=${33}
	BEST_N_ALLELES=${34}
	MIN_MAPPING_QUAL=${35}
	MIN_BASE_QUAL=${36}
	HAPLOTYPE_LENGTH=${37}
	MIN_REPEAT_ENTROPY=${38}
	MIN_COVERAGE=${39}
	MIN_ALT_FRACTION=${40}
	MIN_ALT_COUNT=${41}
	MIN_ALT_TOTAL=${42}
	READ_MAX_MISMATCH_FRACTION=${43}
	MAIL=${44}


	#$$$$$$$$$$$$$$$$$$$$$CEB$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#Count number of individuals in current directory

	#this accounts for separate naming conventions
	#untouched files are F and R
	#files trimmed for assembly are r1 r2
	#files trimmed for mapping are R1 R2
		echo "";echo `date` " cbird has hacked ddocent, it will only eat particular files for particular tasks"
		echo "                              untouched files for trimming are *.F.fq.gz and *.R.fq.gz"
		echo "                              files trimmed for assembly are *r1.fq.gz *r2.fq.gz"
		echo "                              files trimmed for mapping are *R1.fq.gz *R2.fq.gz"
		
	if [ "$TRIM" == "yes" ]; then
		echo "";echo `date` " Trimming file extensions selected"
		F="F"
		Fwild="*.F.fq.gz"
		Fsed=".F.fq.gz"
		R="R"
		Rwild="*.R.fq.gz"
		Rsed=".R.fq.gz"
	elif [ "$ASSEMBLY" == "yes" ]; then
		echo "";echo `date` " Assembly file extensions selected"
		F="r1"
		Fwild="*.r1.fq.gz"
		Fsed=".r1.fq.gz"
		R="r2"
		Rwild="*.r2.fq.gz"
		Rsed=".r2.fq.gz"
		echo "";echo `date` " Assembly file extensions selected"
	elif [ "$MAP" == "yes" ]; then
		echo "";echo `date` " Mapping file extensions selected"
		F="R1"
		Fwild="*.R1.fq.gz"
		Fsed=".R1.fq.gz"
		R="R2"
		Rwild="*.R2.fq.gz"
		Rsed=".R2.fq.gz"
	elif [ "$SNP" == "yes" ]; then
		echo "";echo `date` " F1 R1 Genotyping file extensions selected"
		F="R1"
		Fwild="*.R1.fq.gz"
		Fsed=".R1.fq.gz"
		R="R2"
		Rwild="*.R2.fq.gz"
		Rsed=".R2.fq.gz"
	else
		echo "";echo `date` " Couldn't decide which file extensions were selected"
		F="r1"
		Fwild="*.r1.fq.gz"
		Fsed=".r1.fq.gz"
		R="r2"
		Rwild="*.r2.fq.gz"
		Rsed=".r2.fq.gz"
	fi
	echo "";echo `date` " End extensions selected" $R
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$CEB$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	NumInd=$(ls $Fwild | wc -l)
	NumInd=$(($NumInd - 0))

	#echo "";echo " $NumInd libraries counted"

	#Create list of sample names
	if [ ! -s "namelist.$CUTOFF.$CUTOFF2" ];then
		ls $Fwild > namelist.$CUTOFF.$CUTOFF2
		sed -i'' -e "s/$Fsed//g" namelist.$CUTOFF.$CUTOFF2
	else
		echo "";echo `date` " The namelist file already exists and was not recreated. If you added more *.fq.gz files to this directory, you should delete namelist so that it will be recreated."
	fi
	if [ ! -s "namelistfr.$CUTOFF.$CUTOFF2" ];then
		#echo $Fwild
		ls $Fwild > namelistfr.$CUTOFF.$CUTOFF2
		ls $Rwild >> namelistfr.$CUTOFF.$CUTOFF2
		sed -i'' -e "s/\.fq\.gz//g" namelistfr.$CUTOFF.$CUTOFF2
	else
		echo "";echo `date` " The namelistfr file already exists and was not recreated. If you added more *.fq.gz files to this directory, you should delete namelist so that it will be recreated."
	fi
	#Create an array of sample names
	#NUMNAMES=$(mawk '/_/' namelist.$CUTOFF.$CUTOFF2 | wc -l)
	NUMNAMES=$(grep -c '^' namelist.$CUTOFF.$CUTOFF2)
	echo ""; echo `date` " NUMNAMES=" $NUMNAMES
	echo ""; echo `date` " NUMIND=" $NumInd

	if [ "$NUMNAMES" -eq "$NumInd" ]; then
		NAMES=( `cat "namelist.$CUTOFF.$CUTOFF2" `)
		echo " ";echo " The samples being processed are:"
		cat namelist.$CUTOFF.$CUTOFF2 | parallel "echo {}"
		echo ""
	else
		echo " "
		echo "Individuals do not follow the dDocent naming convention."
		echo "Please rename individuals to: Locality_Individual.F.fq.gz"
		echo "For example: LocA_001.F.fq.gz"
		echo "cbird has hacked ddocent, it will only eat particular files for particular tasks"
		echo "untouched files for trimming are *.F.fq.gz and *.R.fq.gz"
		echo "files trimmed for assembly are *r1.fq.gz *r2.fq.gz"
		echo "files trimmed for mapping are *R1.fq.gz *R2.fq.gz"
		exit 1
	fi

	#Sets a start time variable
	STARTTIME=$(date)

	#ceb_I'm adding an if statement so that this only runs if the user wants it to
	if [ "$FixStacks" == "yes" ]; then
		echo "";echo `date` " Fixing Stacks anomalies"

		#STACKS adds a strange _1 or _2 character to the end of processed reads, this looks for checks for errant characters and replaces them.
		#This functionality is now parallelized and will run if only SE sequences are used.

		STACKS=$(cat namelist.$CUTOFF.$CUTOFF2| parallel -j $NUMProc --no-notice "zcat {}$Fsed | head -1" | mawk '!/\/1/' | wc -l)
		FB1=$(( $NUMProc / 2 ))
		if [ $STACKS -gt 0 ]; then
			echo " ";echo `date` " Removing the _1 character and replacing with /1 in the name of every sequence"
			cat namelist.$CUTOFF.$CUTOFF2 | parallel -j $FB1 --no-notice "zcat {}.$F.fq.gz | sed -e 's:_1$:/1:g' > {}.$F.fq"
			rm $Fwild
			cat namelist.$CUTOFF.$CUTOFF2 | parallel -j $FB1 --no-notice "gzip {}.$F.fq"
		fi

		if [ -f "${NAMES[@]:(-1)}"$Rsed ]; then
			
			STACKS=$(cat namelist.$CUTOFF.$CUTOFF2| parallel -j $NUMProc --no-notice "zcat {}$Rsed | head -1" | mawk '!/\/2/' | wc -l)

			if [ $STACKS -gt 0 ]; then
				echo"";echo `date` " Removing the _2 character and replacing with /2 in the name of every sequence"
				cat namelist.$CUTOFF.$CUTOFF2 | parallel -j $FB1 --no-notice "zcat {}.$R.fq.gz | sed -e 's:_2$:/2:g' > {}.$R.fq"
				rm $Rwild
				cat namelist.$CUTOFF.$CUTOFF2 | parallel -j $FB1 --no-notice "gzip {}.$R.fq"
			fi
		fi
	fi

	##Section of logic statements that dictates the order and function of processing the pipeline

	if [[ "$TRIM" == "yes" && "$ASSEMBLY" == "yes" ]]; then
		echo " ";echo `date` "Trimming reads and simultaneously assembling reference sequences" 
					
			#$$$$$$$$$$$$ceb$$$$$$$$$$$$
		#I made two functions for trimming, TrimReads trims for Mapping, TrimReadsRef trims for reference assembly.  I have modified below
		TrimReadsRef & 2> trimref.log
		TrimReads & 2> trim.log
			Assemble
			#setupRainbow 2> rainbow.log
			wait
	fi

	if [[ "$TRIM" == "yes" && "$ASSEMBLY" != "yes" ]]; then

			TrimReadsRef 2> trimref.log
		TrimReads 2> trim.log
	fi                
					
	if [[ "$TRIM" != "yes" && "$ASSEMBLY" == "yes" ]]; then                
			Assemble
			#setupRainbow 2> rainbow.log
	fi

	#Checks to see if reads will be mapped.
	if [ "$MAP" != "no" ]; then
		echo " ";echo `date` " Using BWA to map reads."
		if [ reference.$CUTOFF.$CUTOFF2.fasta -nt reference.$CUTOFF.$CUTOFF2.fasta.fai ]; then
			samtools faidx reference.$CUTOFF.$CUTOFF2.fasta
			bwa index reference.$CUTOFF.$CUTOFF2.fasta &> index.$CUTOFF.$CUTOFF2.log
		fi
	#dDocent now checks for trimmed read files before attempting mapping
			if [[ "$MAP" != "no" && ! -f "${NAMES[@]:(-1)}"$Rsed ]]; then
				echo "dDocent cannot locate trimmed reads files"
				echo "Please rerun dDocent with quality trimming"
				exit 1
			fi
	#This next section of code checks to see if the reference was assembled by dDocent 
	#and if so, modifies the expected insert length distribution for BWA's metric for proper pairing
			echo ""
			echo -n `date` " was ref assembled by dDocent?"
			if head -1 reference.$CUTOFF.$CUTOFF2.fasta | grep -e 'dDocent' reference.$CUTOFF.$CUTOFF2.fasta 1>/dev/null; then
				echo " Yes. Modifying expected insert length for bwa"
				rm lengths.$CUTOFF.$CUTOFF2.txt &> /dev/null
				for i in "${NAMES[@]}";
					do
					if [ -f "$i.$R.fq.gz" ]; then
						#echo ""
						#echo troubleshooting
						#echo $i
						#echo $R
						#echo $i.$R.fq.gz
						zcat $i.$R.fq.gz | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
					fi
				done	
				if [ -f "lengths.$CUTOFF.$CUTOFF2.txt" ]; then
					MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.$CUTOFF.$CUTOFF2.txt| head -1)
					INSERT=$(($MaxLen * 2 ))
					INSERTH=$(($INSERT + 100 ))
					INSERTL=$(($INSERT - 100 ))
					SD=$(($INSERT / 5))
				fi
				#BWA for mapping for all samples.  As of version 2.0 can handle SE or PE reads by checking for PE read files
				echo ""
				echo `date` " Run bwa mem on dDocent files"
				for i in "${NAMES[@]}"
				do
					if [ -f "$i.R2.fq.gz" ]; then
						bwa mem reference.$CUTOFF.$CUTOFF2.fasta $i.R1.fq.gz $i.R2.fq.gz -L $MAPPING_CLIPPING_PENALTY -I $INSERT,$SD,$INSERTH,$INSERTL -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFF.$CUTOFF2.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q $MAPPING_MIN_QUALITY -f 3 -F $SAMTOOLS_VIEW_F -SbT reference.$CUTOFF.$CUTOFF2.fasta - > $i.$CUTOFF.$CUTOFF2.bam 2>$i.$CUTOFF.$CUTOFF2.bam.log
					else
						bwa mem reference.$CUTOFF.$CUTOFF2.fasta $i.R1.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFF.$CUTOFF2.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q $MAPPING_MIN_QUALITY -f 3 -F $SAMTOOLS_VIEW_F -SbT reference.$CUTOFF.$CUTOFF2.fasta - > $i.$CUTOFF.$CUTOFF2.bam 2>$i.$CUTOFF.$CUTOFF2.bam.log
					fi
					echo ""
					echo `date` " run samtools sort" $i
					samtools sort -@$NUMProc $i.$CUTOFF.$CUTOFF2.bam -o $i.$CUTOFF.$CUTOFF2.bam 
					mv $i.$CUTOFF.$CUTOFF2.bam $i.$CUTOFF.$CUTOFF2-RG.bam
					echo ""
					echo `date` " run samtools insert" $i
					samtools index $i.$CUTOFF.$CUTOFF2-RG.bam
				done
			else
				echo " No. Expected insert length not modified for BWA"
				echo ""
				echo `date` " Run bwa mem on non-dDocent files"
				for i in "${NAMES[@]}"
				do
					if [ -f "$i.R2.fq.gz" ]; then
						bwa mem reference.$CUTOFF.$CUTOFF2.fasta $i.R1.fq.gz $i.R2.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFF.$CUTOFF2.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q $MAPPING_MIN_QUALITY -f 3 -F $SAMTOOLS_VIEW_F -SbT reference.$CUTOFF.$CUTOFF2.fasta - > $i.$CUTOFF.$CUTOFF2.bam 2>$i.$CUTOFF.$CUTOFF2.bam.log
					else
						bwa mem reference.$CUTOFF.$CUTOFF2.fasta $i.R1.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFF.$CUTOFF2.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q $MAPPING_MIN_QUALITY -f 3 -F $SAMTOOLS_VIEW_F -SbT reference.$CUTOFF.$CUTOFF2.fasta - > $i.$CUTOFF.$CUTOFF2.bam 2>$i.$CUTOFF.$CUTOFF2.bam.log
					fi
					samtools sort -@$NUMProc $i.$CUTOFF.$CUTOFF2.bam -o $i.$CUTOFF.$CUTOFF2.bam 
					mv $i.$CUTOFF.$CUTOFF2.bam $i.$CUTOFF.$CUTOFF2-RG.bam
					samtools index $i.$CUTOFF.$CUTOFF2-RG.bam
				done
			fi
	fi

	##Creating mapping intervals if needed, CreateIntervals function is defined later in script
	#If mapping is being performed, intervals are created automatically

	if [ "$MAP" != "no" ]; then
		echo "";echo `date` " Create Intervals"
		ls *.$CUTOFF.$CUTOFF2-RG.bam >bamlist.$CUTOFF.$CUTOFF2.list
		CreateIntervals 
	fi

	##SNP Calling Section of code

	if [ "$SNP" != "no" ]; then
		#Create list of BAM files
		ls *.$CUTOFF.$CUTOFF2-RG.bam >bamlist.$CUTOFF.$CUTOFF2.list
		#If mapping is not being performed, but intervals do not exist they are created
		if [[ "$MAP" == "no" && ! -f "cat.$CUTOFF.$CUTOFF2-RRG.bam" ]]; then
			echo "";echo `date` " recreating intervals 1"
			CreateIntervals 
		fi
		#Check for runs from older versions to ensure the recreation of cat-RRG.bam
		if [[ "$MAP" == "no" && -f "map.$CUTOFF.$CUTOFF2.bed" ]]; then
			echo ""; echo `date` " recreating intervals 2"
			CreateIntervals 
		fi
		#Check to make sure interval files have been created
		if [[ "$MAP" == "no" && ! -f "mapped.$CUTOFF.$CUTOFF2.bed" ]]; then
			echo ""; echo `date` " recreating intervals again"
			bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam > map.$CUTOFF.$CUTOFF2.bed
			bedtools merge -i map.$CUTOFF.$CUTOFF2.bed > mapped.$CUTOFF.$CUTOFF2.bed
			#rm map.$CUTOFF.$CUTOFF2.bed
		fi
		#This code estimates the coverage of reference intervals and removes intervals in 0.01% of depth
		#This allows genotyping to be more effecient and eliminates extreme copy number loci from the data
		echo "";echo `date` " Estimate coverage of ref intervals & remove extreme copy number loci"
		if [ "cat.$CUTOFF.$CUTOFF2-RRG.bam" -nt "cov.$CUTOFF.$CUTOFF2.stats" ]; then
			#coverageBed -abam cat.$CUTOFF.$CUTOFF2-RRG.bam -b mapped.$CUTOFF.$CUTOFF2.bed -counts > cov.$CUTOFF.$CUTOFF2.stats
			if [ "$BEDTOOLSFLAG" == "OLD" ]; then
				coverageBed -abam cat.$CUTOFF.$CUTOFF2-RRG.bam -b mapped.$CUTOFF.$CUTOFF2.bed -counts > cov.$CUTOFF.$CUTOFF2.stats
			else
				bedtools coverage -b cat.$CUTOFF.$CUTOFF2-RRG.bam -a mapped.$CUTOFF.$CUTOFF2.bed -counts -sorted > cov.$CUTOFF.$CUTOFF2.stats
			fi			
		fi
		
		if head -1 reference.$CUTOFF.$CUTOFF2.fasta | grep -e 'dDocent' reference.$CUTOFF.$CUTOFF2.fasta 1>/dev/null; then
		
			DP=$(mawk '{print $4}' cov.$CUTOFF.$CUTOFF2.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFF.$CUTOFF2.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		else
			DP=$(mawk '{print $4}' cov.$CUTOFF.$CUTOFF2.stats | sort -rn | perl -e '$d=.00005;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFF.$CUTOFF2.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		fi
		echo ""
		echo `date` " Making the bed files"
		mawk -v x=$DP '$4 < x' cov.$CUTOFF.$CUTOFF2.stats | sort -V -k1,1 -k2,2 | mawk -v x1="$CUTOFF.$CUTOFF2" -v cutoff=$CC 'BEGIN{i=1} 
		{
			len=$3-$2;lc=len*$4;cov = cov + lc
			if ( cov < cutoff) {x="mapped."i"."x1".bed";print $1"\t"$2"\t"$3 > x}
			else {i=i+1; x="mapped."i"."x1".bed"; print $1"\t"$2"\t"$3 > x; cov=0}
		}' 
		#orig code
		#FB2=$(( $NUMProc / 4 ))
		
		#ceb modified
		FB2=1
		

		#Creates a population file to use for more accurate genotype calling
		if [ ! -s popmap.$CUTOFF.$CUTOFF2 ]; then
			echo ""
			echo `date` " Create popmap"
		
			cut -f1 -d "_" namelist.$CUTOFF.$CUTOFF2 > p.$CUTOFF.$CUTOFF2
			paste namelist.$CUTOFF.$CUTOFF2 p.$CUTOFF.$CUTOFF2 > popmap.$CUTOFF.$CUTOFF2
			rm p.$CUTOFF.$CUTOFF2
			cat popmap.$CUTOFF.$CUTOFF2
		fi
	###############################################################################################	

	###############################################################################################
	###New implementation of SNP calling here to save on memory	

			
		#ceb inserting this to separate splitting of bam files from freebays routine
		call_variants(){
			if [ ! -s split.$1.bam ]; then
				samtools view -@$3 -b -1 -L mapped.$1.bed -o split.$1.bam cat.$2-RRG.bam
			fi
			if [ ! -s split.$1.bam.bai ]; then
				samtools index split.$1.bam
			fi
			# #determine if ploidy or a cnv file is given
			#[[ $4 =~ ^-?[0-9]+$ ]] && ploidy="integer"
			
			# #if $POOLS == no then
			# if [ "${15}" == "no" ]; then
				# #echo Pools are a no ${15}
				# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -V -p ${$16} -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage ${10} -F ${11} -C ${12} -G ${13} -z ${14} --populations popmap.$2 
			# #if $POOLS == true but no pool_ploidy file is provided then
			# elif [ "$4" == "no" ]; then
				# #echo Pools, but no ploidy file $4
				# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -J -w -a -V -p ${$16} -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage ${10} -F ${11} -C ${12} -G ${13} -z ${14} --populations popmap.$2 
			# #if $POOLS == true and a pool_ploidy file is provided then
			# else
				# #echo Pools and a ploidy file $4 ${15}
				# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -J -w -a -V --cnv-map $4  #-n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage ${10} -F ${11} -C ${12} -G ${13} -z ${14} --populations popmap.$2 
			# fi
		}
		export -f call_variants	
		
		echo
		echo `date` Splitting BAM File
		ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --no-notice call_variants {} $CUTOFF.$CUTOFF2 $FB2 $POOL_PLOIDY_FILE $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION $POOLS $PLOIDY 
		
		
		#
		# [[ $POOL_PLOIDY_FILE =~ ^-?[0-9]+$ ]] && ploidy="integer"
		# echo "ploidy is $ploidy"
		# echo $POOL_PLOIDY_FILE 
		# if [ "$POOLS" == "no" ]; then
			# #echo Pools are a no ${15}
		echo
		echo `date` Calling variants using freebayes			
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -V -p $PLOIDY -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2 
		
		# #if $POOLS == true but no pool_ploidy file is provided then
		# elif [ "$ploidy" == "integer" ]; then
			# #echo Pools, but no ploidy file $4
			
			# ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J -w -a -V -p $PLOIDY -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2
		
		# #if $POOLS == true and a pool_ploidy file is provided then
		# elif [ "$ploidy" != "integer" ]; then
			# # # #echo Pools and a ploidy file $4 ${15}
			
			echo
			echo $POOL_PLOIDY_FILE   #for freebayes to call allelotypes, it requires a copy number variation file. this is that file.
			ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J -w -a -V -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2 --cnv-map $POOL_PLOIDY_FILE 
			
			#ceb test  ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J --cnv-map $POOL_PLOIDY_FILE 
			# # ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf > index.txt
			# # # cat index.txt | while read h
			# # # do
				# # # freebayes -b split.$h.bam -t mapped.$h.bed -v raw.$h.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J -w -a -V --cnv-map $POOL_PLOIDY_FILE #-n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2 
			# # # done <index.txt
			# # bamindex=( `cat "index.txt" `)
			# # for i in "${bamindex[@]}"
			# # do
				# # echo $i
				# freebayes -b split.$i.bam -t mapped.$i.bed -v raw.$i.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J -w -a -V --cnv-map $POOL_PLOIDY_FILE #-n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2 
			# # done
		# fi

			# ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf > index.txt
			# bamindex=( `cat "index.txt" `)
			# for i in "${bamindex[@]}"
			# do
				# echo $i
				# freebayes -b split.$i.bam -t mapped.$i.bed -v raw.$i.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -J -w -a -V --cnv-map $POOL_PLOIDY_FILE #-n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2 
			# done

		
		
		# #ceb: this is the original routine for calling variants in individuals
		# call_genos(){
			
			
			# #samtools view -@$3 -b -1 -L mapped.$1.bed -o split.$1.bam cat.$2-RRG.bam
			# #samtools index split.$1.bam
			
			# #purtiz's freebayes for individuals
			# #freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -m 5 -q 5 --haplotype-length 3 --min-repeat-entropy 1 -V --populations popmap.$2 -n 10
			# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -V -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage $10 -F $11 -C $12 -G $13 -z $14 --populations popmap.$2 
			# #rm split.$1.bam*
		# }
		# export -f call_genos
		
		# #ceb: I made this to handle variant calling in pools
		# call_allelotypes(){
			# #samtools view -@$3 -b -1 -L mapped.$1.bed -o split.$1.bam cat.$2-RRG.bam
			# #samtools index split.$1.bam
			
			# #bird's freebayes for pools, attempt 1
			# #if [ "$POOL_PLOIDY_FILE" == "no" ]; then
				# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -K -w -a -V -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage $10 -F $11 -C $12 -G $13 -z $14 --populations popmap.$2 
																												# #-m 20 -q 20 --min-repeat-entropy 1 --min-coverage 5 -C 2 -G 3 -F 0.025 -z 0.20 -w -a -K -p 20 --haplotype-length 125 -n 0 --populations popmap.$2 
			# #else
				# #freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -J -w -a -V --cnv-map $4 -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage $10 -F $11 -C $12 -G $13 -z $14 --populations popmap.$2 
																													# #-m 20 -q 20 --min-repeat-entropy 1 --min-coverage 5 -C 2 -G 3 -F 0.025 -z 0.20 -w -a -K -p 20 --haplotype-length 125 -n 0 --populations popmap.$2 	
			# #fi
			# #rm split.$1.bam*
		# }
		# export -f call_allelotypes
		# call_allelotypes2(){
			# #samtools view -@$3 -b -1 -L mapped.$1.bed -o split.$1.bam cat.$2-RRG.bam
			# #samtools index split.$1.bam
			
			# #bird's freebayes for pools, attempt 1
			# #if [ "$POOL_PLOIDY_FILE" == "no" ]; then
				# #freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -K -w -a -V -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage $10 -F $11 -C $12 -G $13 -z $14 --populations popmap.$2 
																												# #-m 20 -q 20 --min-repeat-entropy 1 --min-coverage 5 -C 2 -G 3 -F 0.025 -z 0.20 -w -a -K -p 20 --haplotype-length 125 -n 0 --populations popmap.$2 
			# #else
				# freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.$2.fasta -J -w -a -V --cnv-map $4 -n $5 -m $6 -q $7 -E $8 --min-repeat-entropy $9 --min-coverage $10 -F $11 -C $12 -G $13 -z $14 --populations popmap.$2 
																													# #-m 20 -q 20 --min-repeat-entropy 1 --min-coverage 5 -C 2 -G 3 -F 0.025 -z 0.20 -w -a -K -p 20 --haplotype-length 125 -n 0 --populations popmap.$2 	
			# #fi
			# #rm split.$1.bam*
		# }
		# export -f call_allelotypes2


		#CEB decide wether to call genotypes or allelotypes based upon config file
		#CEB if config files says genotypes then...
		#if [ "$POOLS" == "no" ]; then
			#echo "";echo `date` " Using GNU Parallel & FreeBayes to call genotypes in diploid individuals"
			
			#this line executes the call_genos function and passes 2 variables to it
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice call_genos {} $CUTOFF.$CUTOFF2 $FB2 $POOL_PLOIDY_FILE $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -m 5 -q 5 --haplotype-length 3 --min-repeat-entropy 1 -V --populations popmap.$CUTOFF.$CUTOFF2 -n 10
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL --haplotype-length $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY -V --populations popmap.$CUTOFF.$CUTOFF2 -n $BEST_N_ALLELES
			
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -V -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2
			

		#CEB elseif config file specifies pools then...
		#elif [ "$POOL_PLOIDY_FILE" == "no" ]; then
			
			#echo "";echo `date` " Preparing CNV files for freebayes"
			
			#rm cnv.*.$CUTOFF.$CUTOFF2.bed
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' > bednames.txt
			#awk '{print $1}' popmap.$CUTOFF.$CUTOFF2 > poolnames.txt
			#awk '{print $2}' pool_ploidy.txt > ploidy.txt
			#ceb make sure there is a file named pool_ploidy.txt that has the following format: pool ploidy
			#ceb make cnv file to feed pool level ploidy to freebayes, this can be used and adapted if different contigs have different ploidy
			#parallel 'sed "s/$/\t{2}/" mapped.{1}.bed >> cnv.{1}.bed ' :::: bednames.txt $POOL_PLOIDY_FILE
			#ceb sort the cnv file
			#ls cnv.*.$CUTOFF.$CUTOFF2.bed | parallel 'sort {} -o {}' 
			
			#echo "";echo `date` " Using GNU Parallel & FreeBayes to call allelotypes in pools or polyploid individuals"
			#this line executes the call_genos function and passes 2 variables to it
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice call_allelotypes {} $CUTOFF.$CUTOFF2 $FB2 $POOL_PLOIDY_FILE $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFF.$CUTOFF2.fasta -K -w -a -V -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $MIN_ALT_COUNT -G $MIN_ALT_TOTAL -z $READ_MAX_MISMATCH_FRACTION --populations popmap.$CUTOFF.$CUTOFF2
		
		#else
			#echo "";echo `date` " Using GNU Parallel & FreeBayes to call allelotypes in pools or polyploid individuals"
			#ls mapped.*.$CUTOFF.$CUTOFF2.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --progress --memfree $MAXMemory -j $NUMProc --no-notice call_allelotypes2 {} $CUTOFF.$CUTOFF2 $FB2 $POOL_PLOIDY_FILE $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION
		#fi
	
		#CEBendif
	####	
		#ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --memfree $MAXMemory -j $FB1 --no-notice --delay 1 freebayes -L bamlist.list -t mapped.{}.bed -v raw.{}.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10
		#ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --memfree $MAXMemory -j $FB1 --no-notice "samtools view -b -L mapped.{}.bed | freebayes -c -t mapped.{}.bed -v raw.{}.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10"


		#rm mapped.*.$CUTOFF.$CUTOFF2.bed 

		mv raw.1.$CUTOFF.$CUTOFF2.vcf raw.01.$CUTOFF.$CUTOFF2.vcf
		mv raw.2.$CUTOFF.$CUTOFF2.vcf raw.02.$CUTOFF.$CUTOFF2.vcf
		mv raw.3.$CUTOFF.$CUTOFF2.vcf raw.03.$CUTOFF.$CUTOFF2.vcf
		mv raw.4.$CUTOFF.$CUTOFF2.vcf raw.04.$CUTOFF.$CUTOFF2.vcf
		mv raw.5.$CUTOFF.$CUTOFF2.vcf raw.05.$CUTOFF.$CUTOFF2.vcf
		mv raw.6.$CUTOFF.$CUTOFF2.vcf raw.06.$CUTOFF.$CUTOFF2.vcf
		mv raw.7.$CUTOFF.$CUTOFF2.vcf raw.07.$CUTOFF.$CUTOFF2.vcf
		mv raw.8.$CUTOFF.$CUTOFF2.vcf raw.08.$CUTOFF.$CUTOFF2.vcf
		mv raw.9.$CUTOFF.$CUTOFF2.vcf raw.09.$CUTOFF.$CUTOFF2.vcf

		vcfcombine raw.*.$CUTOFF.$CUTOFF2.vcf | sed -e 's/	\.\:/	\.\/\.\:/g' > TotalRawSNPs.$CUTOFF.$CUTOFF2.vcf

		if [ ! -d "raw.$CUTOFF.$CUTOFF2.vcf" ]; then
			mkdir raw.$CUTOFF.$CUTOFF2.vcf
		fi

		mv raw.*.$CUTOFF.$CUTOFF2.vcf ./raw.$CUTOFF.$CUTOFF2.vcf

		echo "";echo `date` " Using VCFtools to parse SNPS.vcf for SNPs that are called in at least 90% of individuals"
		vcftools --vcf TotalRawSNPs.$CUTOFF.$CUTOFF2.vcf --geno 0.9 --out Final.$CUTOFF.$CUTOFF2 --counts --recode --non-ref-af 0.001 --max-non-ref-af 0.9999 --mac 1 --minQ 30 --recode-INFO-all &>VCFtools.$CUTOFF.$CUTOFF2.log
		
		echo "";echo `date` " Renaming vcf output"
		mv Final.recode.vcf Final.recode.$CUTOFF.$CUTOFF2.vcf
		mv Final.frq.count Final.frq.$CUTOFF.$CUTOFF2.count
	fi

	##Checking for possible errors

	if [ "$MAP" != "no" ]; then
	ERROR1=$(mawk '/developer/' bwa* | wc -l 2>/dev/null) 
	fi
	ERROR2=$(mawk '/error/' *.bam.log | wc -l 2>/dev/null)
	ERRORS=$(($ERROR1 + $ERROR2))

	#Move various log files to own directory
	if [ ! -d "logfiles" ]; then
	mkdir logfiles
	fi
	mv *.log log ./logfiles 2> /dev/null

	#Sending a completion email

	if [ $ERRORS -gt 0 ]; then
			echo -e "dDocent has finished with errors in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\nPlease check log files\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/Final.log` "\n\ndDocent 2.24 \nThe 'd' is silent, hillbilly." | mailx -s "dDocent has finished with ERRORS!" $MAIL
	else
			echo -e "dDocent has finished with an analysis in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/Final.log` "\n\ndDocent 2.24 \nIt is pronounced Dee-Docent, professor." | mailx -s "dDocent has finished" $MAIL
	fi


	#Creates (or appends to) a dDcoent run file recording variables
	#RMH added ALL variables
	echo "Variables used in dDocent Run at" $STARTTIME >> dDocent.runs
	echo "Number of Processors" >> dDocent.runs
	echo $NUMProc >> dDocent.runs
	echo "Maximum Memory" >> dDocent.runs
	echo $MAXMemory >> dDocent.runs
	echo "Trimming" >> dDocent.runs
	echo $TRIM >> dDocent.runs
	echo "TRIM_LENGTH_ASSEMBLY" >> dDocent.runs
	echo $TRIM_LENGTH_ASSEMBLY >> dDocent.runs
	echo "SEED_ASSEMBLY" >> dDocent.runs
	echo $SEED_ASSEMBLY >> dDocent.runs
	echo "PALIMDROME_ASSEMBLY" >> dDocent.runs
	echo $PALIMDROME_ASSEMBLY >> dDocent.runs
	echo "SIMPLE_ASSEMBLY" >> dDocent.runs
	echo $SIMPLE_ASSEMBLY >> dDocent.runs
	echo "windowSize_ASSEMBLY" >> dDocent.runs
	echo $windowSize_ASSEMBLY >> dDocent.runs
	echo "windowQuality_ASSEMBLY" >> dDocent.runs
	echo $windowQuality_ASSEMBLY >> dDocent.runs
	echo "TRAILING_ASSEMBLY" >> dDocent.runs
	echo $TRAILING_ASSEMBLY >> dDocent.runs
	echo "TRIM_LENGTH_MAPPING" >> dDocent.runs
	echo $TRIM_LENGTH_MAPPING >> dDocent.runs
	echo "LEADING_MAPPING" >> dDocent.runs
	echo $LEADING_MAPPING >> dDocent.runs
	echo "TRAILING_MAPPING" >> dDocent.runs
	echo $TRAILING_MAPPING >> dDocent.runs
	echo "FixStacks" >> dDocent.runs
	echo $FixStacks >> dDocent.runs
	echo "Assembly?" >> dDocent.runs
	echo $ASSEMBLY >> dDocent.runs
	echo "Type_of_Assembly" >> dDocent.runs
	echo $ATYPE >> dDocent.runs
	echo "Clustering_Similarity%" >> dDocent.runs
	echo $simC >> dDocent.runs
	echo "HPC" >> dDocent.runs
	echo $HPC >> dDocent.runs
	echo "MANCUTOFF" >> dDocent.runs
	echo $MANCUTOFF >> dDocent.runs	
	echo "CUTOFF" >> dDocent.runs
	echo $CUTOFF >> dDocent.runs
	echo "CUTOFF2" >> dDocent.runs 
	echo $CUTOFF2 >> dDocent.runs
	echo "Mapping_Reads?" >> dDocent.runs
	echo $MAP >> dDocent.runs
	echo "Mapping_Match_Value" >> dDocent.runs
	echo $optA >> dDocent.runs
	echo "Mapping_MisMatch_Value" >> dDocent.runs
	echo $optB >> dDocent.runs
	echo "Mapping_GapOpen_Penalty" >> dDocent.runs
	echo $optO >> dDocent.runs
	echo "MAPPING_MIN_ALIGNMENT_SCORE" >> dDocent.runs
	echo $MAPPING_MIN_ALIGNMENT_SCORE >> dDocent.runs
	echo "MAPPING_CLIPPING_PENALTY" >> dDocent.runs
	echo $MAPPING_CLIPPING_PENALTY >> dDocent.runs
	echo "MAPPING_MIN_QUALITY" >> dDocent.runs
	echo $MAPPING_MIN_QUALITY >> dDocent.runs
	echo "SAMTOOLS_VIEW_F" >> dDocent.runs
	echo $SAMTOOLS_VIEW_F >> dDocent.runs
	echo "Calling_SNPs?" >> dDocent.runs
	echo $SNP >> dDocent.runs
	echo "POOLS" >> dDocent.runs
	echo $POOLS >> dDocent.runs
	echo "POOL_PLOIDY_FILE" >> dDocent.runs
	echo $POOL_PLOIDY_FILE >> dDocent.runs
	echo "PLOIDY" >> dDocent.runs
	echo $PLOIDY >> dDocent.runs
	echo "BEST_N_ALLELES" >> dDocent.runs
	echo $BEST_N_ALLELES >> dDocent.runs
	echo "MIN_MAPPING_QUAL" >> dDocent.runs
	echo $MIN_MAPPING_QUAL >> dDocent.runs
	echo "MIN_BASE_QUAL" >> dDocent.runs
	echo $MIN_BASE_QUAL >> dDocent.runs
	echo "HAPLOTYPE_LENGTH" >> dDocent.runs
	echo $HAPLOTYPE_LENGTH >> dDocent.runs
	echo "MIN_REPEAT_ENTROPY" >> dDocent.runs
	echo $MIN_REPEAT_ENTROPY >> dDocent.runs
	echo "MIN_COVERAGE" >> dDocent.runs
	echo $MIN_COVERAGE >> dDocent.runs
	echo "MIN_ALT_FRACTION" >> dDocent.runs
	echo $MIN_ALT_FRACTION >> dDocent.runs
	echo "MIN_ALT_COUNT" >> dDocent.runs
	echo $MIN_ALT_COUNT >> dDocent.runs
	echo "MIN_ALT_TOTAL" >> dDocent.runs
	echo $MIN_ALT_TOTAL >> dDocent.runs
	echo "READ_MAX_MISMATCH_FRACTION" >> dDocent.runs
	echo $READ_MAX_MISMATCH_FRACTION >> dDocent.runs
	echo "Email" >> dDocent.runs
	echo $MAIL >> dDocent.runs

}
###############################################################################################


##Function definitions


###############################################################################################
#Function for trimming reads using trimmomatic
TrimReadsRef () { 
	echo " "
	echo "Trimming reads for reference genome"
		date
for i in "${NAMES[@]}";
	do
        zcat $i.$F.fq.gz | head -2 | tail -1 >> lengths_ref.txt
        done	
        MLen=$(mawk '{ print length() | "sort -rn" }' lengths_ref.txt| head -1)
        MLen=$(($MLen / 2))
	TW="MINLEN:$MLen"
	
for i in "${NAMES[@]}"
do
#echo "Trimming Sample $i"
if [ -f $i.$R.fq.gz ]; then

	#$$$$$$$ceb$$$$$$$$$$$$$
	#adjusted settings to retain as many reads as possible
	java -jar $TRIMMOMATIC PE -threads $NUMProc -phred33 $i.F.fq.gz $i.R.fq.gz $i._r1.fq.gz $i.unpairedF.fq.gz $i.r2.fq.gz $i.unpairedR.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY TRAILING:$TRAILING_ASSEMBLY SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY CROP:$TRIM_LENGTH_ASSEMBLY MINLEN:$TRIM_LENGTH_ASSEMBLY  &> $i.pe.trim.log

	#$$$$$$ceb$$$$$
	#make read1 shorter than read 2 by clipping off first few bp
	java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 $i._r1.fq.gz $i.r1.fq.gz HEADCROP:4 &> $i.se.trim.log

else java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 $i.F.fq.gz $i.r1.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY TRAILING:$TRAILING_ASSEMBLY SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY CROP:$TRIM_LENGTH_ASSEMBLY MINLEN:$TRIM_LENGTH_ASSEMBLY &> $i.trim.log
fi 
done
mkdir unpaired_ref &>/dev/null
mv *unpaired*.gz ./unpaired_ref &>/dev/null	
#RMH housekeeping
rm *_r1.fq.gz
mv *r1.fq.gz assembly
mv *r2.fq.gz assembly
}
###############################################################################################

##############################################################################################
#Function for trimming reads using trimmomatic
TrimReads () { 
	echo " "
	echo "Trimming reads for mapping"
		date
for i in "${NAMES[@]}";
	do
        zcat $i.$F.fq.gz | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
        done	
        MLen=$(mawk '{ print length() | "sort -rn" }' lengths.$CUTOFF.$CUTOFF2.txt| head -1)
        MLen=$(($MLen / 2))
	TW="MINLEN:$MLen"
	
for i in "${NAMES[@]}"
do
#echo "Trimming Sample $i"
if [ -f $i.R.fq.gz ]; then
java -jar $TRIMMOMATIC PE -threads $NUMProc -phred33 $i.F.fq.gz $i.R.fq.gz $i.R1.fq.gz $i.unpairedF.fq.gz $i.R2.fq.gz $i.unpairedR.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY LEADING:$LEADING_MAPPING TRAILING:$TRAILING_MAPPING SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY MINLEN:$TRIM_LENGTH_MAPPING &> $i.trim.log
else java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 $i.F.fq.gz $i.R1.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY LEADING:$LEADING_MAPPING TRAILING:$TRAILING_MAPPING SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY MINLEN:$TRIM_LENGTH_MAPPING &> $i.trim.log
fi 
done
mkdir unpaired &>/dev/null
mv *unpaired*.gz ./unpaired &>/dev/null	
#RMH housekeeping
mv *R1.fq.gz mapping
mv *R2.fq.gz mapping
}
###############################################################################################



###############################################################################################
#Main function for assembly
Assemble()
{
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
SED1='s/^[ 	]*//'
SED2='s/ /	/g'
FRL=$(zcat ${NAMES[0]}.r1.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)

	

#$$$$$$$$$$$CEB$$$$$$$$$$$$$$$$
#this block of code will either align or concatenate r1 & r2 seqs and create a file of unique sequences 
	echo " "; echo `date` " Reference Genome Assembly Step 1, Select Cutoffs"


	if [ ${NAMES[@]:(-1)}.r1.fq.gz -nt ${NAMES[@]:(-1)}.uniq.seqs ];then
		echo " ";echo `date` " the *.fq.gz files are newer than the *uniq.seqs files or the *.uniq.seq files do not exist"
		
		if [ ! -s ${NAMES[@]:(-1)}.uniq.seqs ];then
			echo "";echo `date` " make *.uniq.seqs files"		
			if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
			#If PE assembly, creates a concatenated file of every unique for each individual in parallel
				#CEB this would be faster if fwd & rev could be created in 1 line
				#cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "zcat {}.r1.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
				#cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "zcat {}.r2.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
				#here is my try at making this go faster
				cat namelistfr.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "zcat {}.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.XXX"
				rename r1.XXX forward *r1.XXX
				rename r2.XXX reverse *r2.XXX
				if [ "$ATYPE" = "RPE" ]; then
					echo "";echo `date` " RPE assembly"
					cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNNNNNNNNNNNN/' | sort | uniq -c -w $FRL| sed -e '$SED1' | sed -e '$SED2' > {}.uniq.seqs"
				else
					echo "";echo `date` " PE assembly"
					cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNNNNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
				fi
				rm *.forward
				rm *.reverse
			fi
			if [ "$ATYPE" == "SE" ]; then
				echo "";echo `date` " SE assembly"
			#if SE assembly, creates files of every unique read for each individual in parallel
				cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "zcat {}.r1.fq.gz | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
			fi
			if [ "$ATYPE" == "OL" ]; then
			#If OL assembly, dDocent assumes that the majority of PE reads will overlap, so the software PEAR is used to merge paired reads into single reads
				echo "";echo `date` " OL assembly"
				for i in "${NAMES[@]}";
					do
					#zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
					zcat $i$Rsed | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
				done	
				MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.$CUTOFF.$CUTOFF2.txt| head -1)
				LENGTH=$(( $MaxLen / 3))
				echo "";echo `date` " OL assembly: PEAR "
				for i in "${NAMES[@]}"
					do
					#pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH 
					pearRM -f $i$Fsed -r $i$Rsed -o $i -j $NUMProc -n $LENGTH -p 0.0001
				done
				echo "";echo `date` " OL assembly: create *.uniq.seqs "
				cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
			fi
		else
			echo "";echo "***************************************************************************"
			echo `date` " if you want to recreate *uniq.seqs files then you need to delete the present *.uniq.seq files "
			echo "***************************************************************************"
			echo " "
		fi
	fi
#end block of code that will either align or concatenate R1 & R2 seqs and create a file of unique sequences 
#$$$$$$$$$$$$CEB$$$$$$$$$$$$$$$$



#Create a data file with the number of unique sequences and the number of occurrences
#if the creation of uniqseq.data were done in parallel, it would be much faster.  Evan?
#what happens now is the uniqseq file from each individual is combined into 1, then the 1 file is repeatedly queried to make uniqseqs.data which
#which is a simple summary 

if [ -s "uniq.seqs.gz" ]; then
	echo "";echo `date` " uniq.seqs.gz exists"
	if [ uniq.seqs.gz -nt uniq.seqs ]; then
		echo "";echo `date` " unzip uniq.seqs.gz"
		gunzip uniq.seqs.gz 2>/dev/null
	fi
fi

if [ ! -e "uniqseq.data" ]; then
	#echo "";echo `date` " initiating make uniq.seqs"
	#cat *.uniq.seqs > uniq.seqs
	
	# if uniqseq.data does not exist or does not have a file size greater than zero then do this
	#this removes singletons,should speed up making uniqseq.data since singletons aren't used in uniq.seqs
	#first step: remove singletons from *.uniq.seqs files in parallel
	echo "";echo `date` " removing singletons from *.uniq.seqs"	
	ls *.uniq.seqs | parallel "grep -v '^1[[:space:]]' {} > {}.temp "
	ls *.uniq.seqs | parallel "mv {}.temp {} "
	
	#second step 
	#echo "";echo `date` " make uniq.seqs"
	#cat *.uniq.seqs | grep -v "^1[[:space:]]" > uniq.seqs
	
fi

if [ ! -e "uniqseq.data" ]; then
	echo "";echo `date` " make uniqseq.data"
	seq 2 20 > pfile
	#cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
	ls *.uniq.seqs > uniqseqfilenames.txt

	echo "";echo `date` " counting up uniqseqs"
	#CEB i changed the following line to overwrite the previous results rather than append
	parallel --no-notice "(echo -n -e {2}'\t' && grep -c {2} {1}) > {1}.{2}.test " :::: uniqseqfilenames.txt pfile
	
	#Evan, here, we should count up the total number of lines in the uniq.seqs files, then subtract the total from the 2-20 counts
	#that will give us  the number greater than 20, which can be appended to the *.data files
	#ls *.uniq.seqs | parallel --no-notice "sed -n '$=' {} > {}.tot.test"

	ls *.uniq.seqs | parallel --no-notice "cat {}.*.test | sort -g > {}.data"
	
	echo "";echo `date` " combining uniqseq counts"
	cat *uniq.seqs.data | awk '{sums[$1] += $2;} END { for (i in sums) print i " " sums[i]; }' | sort -g > uniqseq.predata

	echo "";echo `date` " calculating >=X for uniqseq.data"
	#need to make function of awk statement to put inside parallel
	cntUniq() {
	 awk -v var=$1 'FNR >= var {f+=$2}END{print f}' uniqseq.predata >> grtrthans
	}
	rm -rf grtrthans
	#need to do this
	export -f cntUniq
	#subtract 1 from each number in pfile and feed to parallel to run function
	awk '{$1 = $1 - 1; print}' pfile | parallel --no-notice cntUniq
	#make uniqseq.data
	sort -nr grtrthans > greaterthans
	paste pfile greaterthans > uniqseq.data
	#cleanup junk
	rm -rf grtrthans
	rm -rf greaterthans
	rm -rf pfile
	rm -rf *.test
	rm -rf *uniq.seqs.data
	rm -rf *uniqseq.predata

	
	
	#cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
	
fi

#Plot graph of above data
	echo "";echo `date` " plot graph for cutoff 1"; echo ""

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [1:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
echo "";echo `date` " CEB has modified dDocent here such that there is no user input"
echo "";echo `date` " Use the graph above to set CUTOFF1 in the config file"
echo "";echo `date` " CUTOFF1 is currently set to $CUTOFF"
echo "";echo `date` " The graph below is dependent upon the value of CUTOFF1"
if [ "$MANCUTOFF" = "yes" ];then
	echo -e "Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly"
	read CUTOFF
fi

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if [ "$ATYPE" == "RPE" ]; then
	if [ ! -s "uniqCperindv.$CUTOFF" ];then
		echo "";echo `date` " make uniqCperindv.$CUTOFF  tradRAD, not ezRAD & ddRAD"
		parallel --no-notice -j $NUMProc mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | sort | uniq -c -w $FRL | sed -e 's/^[ 	]*//' | sed -e 's/ /	/g' > uniqCperindv.$CUTOFF
	fi
else
	if [ ! -s "uniqCperindv.$CUTOFF" ];then
		echo "";echo `date` " make uniqCperindv.$CUTOFF  ezRAD & ddRAD, not tradRAD"
		#this one takes a long time. Can it be sped up? Evan?
		parallel --no-notice -j $NUMProc mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv.$CUTOFF
		# ceb,5/16/2018 beginning attempt to speed up creation of uniqCperindv.$CUTOFF
		#the grep removes singletons 
		#ls *uniq.seqs | parallel 'cut -f2' | sort | uniq -c | grep -v '^1[[:space:]]' > uniqCperindv.$CUTOFF
	fi
fi
if [ "$NumInd" -gt 10 ]; then
	NUM=$(($NumInd / 2))
else
	NUM=$NumInd
fi

if [ ! -s "uniqseq.$CUTOFF.peri.data" ];then
	echo "";echo `date` " make uniqseq.$CUTOFF.peri.data"
	#for ((i = 2; i <= $NUM; i++)); do
	#	echo $i >> ufile.$CUTOFF
	#done
	#l

	#cat ufile.$CUTOFF
	# this one is slow.  Can it be sped up?
	seq 2 $NUM | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv.$CUTOFF | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.$CUTOFF.peri.data
	#cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '$1 >= x' uniqCperindv.$CUTOFF | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.$CUTOFF.peri.data
	#rm ufile.$CUTOFF
fi

#Plot graph of above data
	echo "";echo `date` " plot graph for cutoff 2"
	echo ""
	export PlotFile=uniqseq.$CUTOFF.peri.data
gnuplot << \EOF 
filename=system("echo $PlotFile")
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot filename with lines notitle
pause -1
EOF

echo -en "\007"
echo -en "\007"
echo -en "\007"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
echo "";echo `date` " CEB has modified dDocent here such that there is no user input"
if [ "$MANCUTOFF" = "yes" ]; then
	echo -e "Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals"
	read CUTOFF2
fi
echo "";echo `date` " Once you have used the graphs to decide upon cutoffs, adjust the config file so that dDocent will complete reference genome assembly"
if [ "$HPC" = "yes" ]; then
	echo "";echo `date` " Stopping because config file is set to HPC, Get Graphs for cutoffs, then stop? yes"
	exit
fi
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Prints instructions on how to move analysis to background and disown process
echo "";echo `date` " Assembly Phase 2"
echo "";echo `date` " At this point, all configuration information has been entered and dDocent may take several hours to run." 
#echo "It is recommended that you move this script to a background operation and disable terminal input and output."
#echo "All data and logfiles will still be recorded."
#echo "To do this:"
#echo "Press control and Z simultaneously"
#echo "Type 'bg' without the quotes and press enter"
#echo "Type 'disown -h' again without the quotes and press enter"
#echo ""
#echo "Now sit back, relax, and wait for your analysis to finish."

#Now that data cutoffs have been chosen, reduce data set to specified set of unique reads, convert to FASTA format,
#and remove reads with substantial amounts of adapters

echo "";echo `date` " mawking"
mawk -v x=$CUTOFF2 '$1 >= x' uniqCperindv.$CUTOFF > uniq.k.$CUTOFF.c.$CUTOFF2.seqs
cut -f2 uniq.k.$CUTOFF.c.$CUTOFF2.seqs > totaluniqseq.$CUTOFF.$CUTOFF2
mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq.$CUTOFF.$CUTOFF2 > uniq.$CUTOFF.$CUTOFF2.full.fasta
LENGTH=$(mawk '!/>/' uniq.$CUTOFF.$CUTOFF2.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(($LENGTH * 3 / 4))
$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.$CUTOFF.$CUTOFF2.full.fasta > uniq.$CUTOFF.$CUTOFF2.fq
#$$$$$ceb$$$$$$$$$$$
#added cp line and turned off the trimmomatic
cp uniq.$CUTOFF.$CUTOFF2.fq uniq.$CUTOFF.$CUTOFF2.fq1
#java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.$CUTOFF.$CUTOFF2.fq uniq.$CUTOFF.$CUTOFF2.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.$CUTOFF.$CUTOFF2.fq1 > uniq.$CUTOFF.$CUTOFF2.fasta
mawk '!/>/' uniq.$CUTOFF.$CUTOFF2.fasta > totaluniqseq.$CUTOFF.$CUTOFF2
rm uniq.$CUTOFF.$CUTOFF2.fq*

#If this is a PE assembly
if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	echo ""
	echo "begin PE Assembly"
	#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
	echo ""
	echo -n "Reads are first clustered using only the Forward reads using CD-hit instead of rainbow "
	date
	sed -e 's/NNNNNNNNNNNNNNNNNNNN/	/g' uniq.$CUTOFF.$CUTOFF2.fasta | cut -f1 > uniq.$CUTOFF.$CUTOFF2.F.fasta
	CDHIT=$(python -c "print max("$simC" - 0.1,0.8)")
	cd-hit-est -i uniq.$CUTOFF.$CUTOFF2.F.fasta -o xxx.$CUTOFF.$CUTOFF2 -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.$CUTOFF.$CUTOFF2.log
	mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.$CUTOFF.$CUTOFF2.clstr | sed 's/[>dDococent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids.$CUTOFF.$CUTOFF2
	paste sort.contig.cluster.ids.$CUTOFF.$CUTOFF2 totaluniqseq.$CUTOFF.$CUTOFF2 > contig.cluster.totaluniqseq.$CUTOFF.$CUTOFF2
	sort -k2,2 -g contig.cluster.totaluniqseq.$CUTOFF.$CUTOFF2 | sed -e 's/NNNNNNNNNNNNNNNNNNNN/	/g' > rcluster.$CUTOFF.$CUTOFF2
	#CD-hit output is converted to rainbow format
	echo ""
	echo -n "Running rainbow div "
	date
	rainbow div -i rcluster.$CUTOFF.$CUTOFF2 -o rbdiv.$CUTOFF.$CUTOFF2.out -f 0.5 -K 10
	echo ""
	echo -n "Running rainbow merge "
	date
	rainbow merge -i rbdiv.$CUTOFF.$CUTOFF2.out -a -o rbasm.$CUTOFF.$CUTOFF2.out -N10000 -l 20 -f 0.75 -r 2 -R10000 
	echo ""
	echo `date` " Selecting contigs"

	#This AWK code replaces rainbow's contig selection perl script
	cat rbasm.$CUTOFF.$CUTOFF2.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
		if (NR == 1) e=$2;
		else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNNNNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/C/) clus=$2;
		else if ($1 ~/L/) len=$2;
		else if ($1 ~/S/) seq=$2;
		else if ($1 ~/N/) freq=$2;
		else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
		else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
		else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
        }' > rainbow.$CUTOFF.$CUTOFF2.fasta

	seqtk seq -r rainbow.$CUTOFF.$CUTOFF2.fasta > rainbow.$CUTOFF.$CUTOFF2.RC.fasta
	mv rainbow.$CUTOFF.$CUTOFF2.RC.fasta rainbow.$CUTOFF.$CUTOFF2.fasta
	
	echo ""
	echo -n "Check for overlap in paired end reads with Pear"
	date
	#The rainbow assembly is checked for overlap between newly assembled Forward and Reverse reads using the software PEAR
	sed -e 's/NNNNNNNNNNNNNNNNNNNN/	/g' rainbow.$CUTOFF.$CUTOFF2.fasta | cut -f1 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.$CUTOFF.$CUTOFF2.F.fq
	sed -e 's/NNNNNNNNNNNNNNNNNNNN/	/g' rainbow.$CUTOFF.$CUTOFF2.fasta | cut -f2 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.$CUTOFF.$CUTOFF2.R.fq

	seqtk seq -r ref.$CUTOFF.$CUTOFF2.R.fq > ref.$CUTOFF.$CUTOFF2.RC.fq
	mv ref.$CUTOFF.$CUTOFF2.RC.fq ref.$CUTOFF.$CUTOFF2.R.fq
	LENGTH=$(mawk '!/>/' rainbow.$CUTOFF.$CUTOFF2.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
	echo length is $LENGTH
	LENGTH=$(( $LENGTH * 5 / 4))
	

	pearRM -f ref.$CUTOFF.$CUTOFF2.F.fq -r ref.$CUTOFF.$CUTOFF2.R.fq -o overlap.$CUTOFF.$CUTOFF2 -p 0.0001 -j $NUMProc -n $LENGTH -v 20

	rm ref.$CUTOFF.$CUTOFF2.F.fq ref.$CUTOFF.$CUTOFF2.R.fq
	
	
	echo ""
	echo -n "More mawking"
	date
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.$CUTOFF.$CUTOFF2.assembled.fastq > overlap.$CUTOFF.$CUTOFF2.fasta
	mawk '/>/' overlap.$CUTOFF.$CUTOFF2.fasta > overlap.$CUTOFF.$CUTOFF2.loci.names
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.$CUTOFF.$CUTOFF2.unassembled.forward.fastq > other.$CUTOFF.$CUTOFF2.F
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.$CUTOFF.$CUTOFF2.unassembled.reverse.fastq > other.$CUTOFF.$CUTOFF2.R
	paste other.$CUTOFF.$CUTOFF2.F other.$CUTOFF.$CUTOFF2.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNNNNNNNNNNNN/g' > other.$CUTOFF.$CUTOFF2.FR

	cat other.$CUTOFF.$CUTOFF2.FR overlap.$CUTOFF.$CUTOFF2.fasta > totalover.$CUTOFF.$CUTOFF2.fasta

	rm *.F *.R
fi

if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" ]]; then
	cp uniq.$CUTOFF.$CUTOFF2.fasta totalover.$CUTOFF.$CUTOFF2.fasta
fi
	echo ""
	echo -n "More CD-HITting"
	date
cd-hit-est -i totalover.$CUTOFF.$CUTOFF2.fasta -o reference.$CUTOFF.$CUTOFF2.fasta.original -M 0 -T 0 -c $simC

echo ""
echo -n "sed command "
date
sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.$CUTOFF.$CUTOFF2.fasta.original > reference.$CUTOFF.$CUTOFF2.fasta

echo ""
echo -n "samtools faidx "
date
samtools faidx reference.$CUTOFF.$CUTOFF2.fasta
echo ""
echo -n "bwa index reference "
date
bwa index reference.$CUTOFF.$CUTOFF2.fasta
	echo ""
	echo -n "End Assembly of Reference Genome"
	date
	echo ""
}
###############################################################################################


###############################################################################################
##Create alignment intervals
##This takes advantage of the fact that RAD loci are very discrete.  Instead of calculating intervals for every BAM file,
##this function merges all BAM files together and removes duplicates.  This overall BAM file 
##is used to create a single list of intervals, saving a large amount of computational time.

CreateIntervals()
{

echo ""
echo `date` " samtools merge"
samtools merge -@$NUMProc -b bamlist.$CUTOFF.$CUTOFF2.list -f cat.$CUTOFF.$CUTOFF2-RRG.bam &>/dev/null
echo ""
echo `date` " samtools index"
samtools index cat.$CUTOFF.$CUTOFF2-RRG.bam 
wait
echo ""
echo `date` " bamToBed"
# bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam > map.$CUTOFF.$CUTOFF2.bed
# echo ""
# echo `date` " bedtools merge"
# bedtools merge -i map.$CUTOFF.$CUTOFF2.bed > mapped.$CUTOFF.$CUTOFF2.bed
# rm map.$CUTOFF.$CUTOFF2.bed

bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam | bedtools merge > mapped.$CUTOFF.$CUTOFF2.bed

}
###############################################################################################


###############################################################################################
#This checks that dDocent has detected the proper number of individuals and exits if incorrect
GetInfo(){
echo "$NumInd individuals are detected. Is this correct? Enter yes or no and press [ENTER]"

read Indcorrect

if [ "$Indcorrect" == "no" ]; then
        echo "Please double check that all fastq files are named Ind01.F.fq.gz and Ind01.R.fq.gz"
        exit 1
elif [ "$Indcorrect" == "yes" ]; then
            echo "Proceeding with $NumInd individuals"
else
        echo "Incorrect Input"
        exit 1
fi

#Tries to get number of processors, if not asks user
NUMProc=( `grep -c ^processor /proc/cpuinfo 2> /dev/null` ) 
NUMProc=$(($NUMProc + 0)) 

echo "dDocent detects $NUMProc processors available on this system."
echo "Please enter the maximum number of processors to use for this analysis."
        read NUMProc
        
if [ $NUMProc -lt 1 ]; then
        echo "Incorrect. Please enter the number of processing cores on this computer"
        read NUMProc
fi                
if [ $NUMProc -lt 1 ]; then
        echo "Incorrect input, exiting"
        exit 1
fi

#Tries to get maximum system memory, if not asks user
MAXMemory=$(($(grep -Po '(?<=^MemTotal:)\s*[0-9]+' /proc/meminfo | tr -d " ") / 1048576))G

echo "dDocent detects $MAXMemory maximum memory available on this system."
echo "Please enter the maximum memory to use for this analysis. The size can be postfixed with 
K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 
1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
echo "For example, to limit dDocent to ten gigabytes, enter 10G or 10g"
        read MAXMemory

while [[ -z $MAXMemory ]];
	do
	echo "Incorrect input"
	echo -e "Please enter the maximum memory to use for this analysis. The size can be postfixed with K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
	echo -e "This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0"
	echo -e "Then press [ENTER]"
	read MAXMemory
	done

#Asks if user wants to trim reads.  This allows this part of the pipeline to be skipped during subsequent analyses
echo -e "\nDo you want to quality trim your reads?" 
echo "Type yes or no and press [ENTER]?"

read TRIM

#Asks if user wants to perform an assembly.  This allows this part of the pipeline to be skipped during subsequent analyses

echo -e "\nDo you want to perform an assembly?"
echo "Type yes or no and press [ENTER]?"

read ASSEMBLY

if [ "$ASSEMBLY" == "no" ]; then
        echo -e "\nReference contigs need to be in a file named reference.$CUTOFF.$CUTOFF2.fasta\n"
        sleep 1
else
	echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
	echo -e "Then press [ENTER]"
	read ATYPE

	while [[ $ATYPE != "SE" && $ATYPE != "PE" && $ATYPE != "OL" && $ATYPE != "RPE" ]];
	do
	echo "Incorrect input"
	echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
	echo -e "Then press [ENTER]"
	read ATYPE
	done
fi
#If performing de novo assembly, asks if the user wants to enter a different -c value
if [ "$ASSEMBLY" == "yes" ]; then
    echo "Reads will be assembled with Rainbow"
    echo "CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa."
    echo "Would you like to enter a new c parameter now? Type yes or no and press [ENTER]"
    read optC
    if [ "$optC" == "no" ]; then
            echo "Proceeding with default 0.9 value."
            simC=0.9
        elif [ "$optC" == "yes" ]; then
            echo "Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)"
            read newC
            simC=$newC
        else
            echo "Incorrect input. Proceeding with the default value."
            simC=0.9
        fi
fi

#Asks if user wants to map reads and change default mapping variables for BWA
echo "Do you want to map reads?  Type yes or no and press [ENTER]"
read MAP
if [ "$MAP" == "no" ]; then
        echo "Mapping will not be performed"
        optA=1
    	optB=4
    	optO=6
        else
                echo "BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa."
                echo "Would you like to enter a new parameters now? Type yes or no and press [ENTER]"
                read optq

        if [ "$optq" == "yes" ]; then
        echo "Please enter new value for A (match score).  It should be an integer.  Default is 1."
        read newA
        optA=$newA
                echo "Please enter new value for B (mismatch score).  It should be an integer.  Default is 4."
        read newB
        optB=$newB
                echo "Please enter new value for O (gap penalty).  It should be an integer.  Default is 6."
        read newO
        optO=$newO
        else
                echo "Proceeding with default values for BWA read mapping."
                optA=1
                optB=4
                optO=6
        fi
fi

#Does user wish to call SNPs?
echo "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
read SNP

while [[ $SNP != "yes" && $SNP != "no" ]];
	do
	echo "Incorrect input"
	echo -e "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
	read SNP
	done

#Asks user for email address to notify when analysis is complete
echo ""
echo "Please enter your email address.  dDocent will email you when it is finished running."
echo "Don't worry; dDocent has no financial need to sell your email address to spammers."
read MAIL
echo ""
echo ""

if [ "$ASSEMBLY" == "no" ]; then
#Prints instructions on how to move analysis to background and disown process
echo "At this point, all configuration information has been entered and dDocent may take several hours to run." 
echo "It is recommended that you move this script to a background operation and disable terminal input and output."
echo "All data and logfiles will still be recorded."
echo "To do this:"
echo "Press control and Z simultaneously"
echo "Type 'bg' without the quotes and press enter"
echo "Type 'disown -h' again without the quotes and press enter"
echo ""
echo "Now sit back, relax, and wait for your analysis to finish."
fi

if [ "$ASSEMBLY" == "yes" ]; then
echo "dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background."
fi
}


#Actually starts program
if [ -n "$1" ]; then
	#main $1
	# RMH added ALL variables below; main $NUMProc $MAXMemory $TRIM $FixStacks $ASSEMBLY $ATYPE $simC $MAP $optA $optB $optO $SNP $MAIL $HPC $MANCUTOFF $CUTOFF $CUTOFF2 $TRIM_LENGTH_ASSEMBLY $TRIM_LENGTH_MAPPING $SEED_ASSEMBLY $PALIMDROME_ASSEMBLY $SIMPLE_ASSEMBLY $windowSize_ASSEMBLY $windowQuality_ASSEMBLY
	main $NUMProc $MAXMemory $TRIM $TRIM_LENGTH_ASSEMBLY $SEED_ASSEMBLY $PALIMDROME_ASSEMBLY $SIMPLE_ASSEMBLY $windowSize_ASSEMBLY $windowQuality_ASSEMBLY $TRAILING_ASSEMBLY $TRIM_LENGTH_MAPPING $LEADING_MAPPING $TRAILING_MAPPING $FixStacks $ASSEMBLY $ATYPE $simC $HPC $MANCUTOFF $CUTOFF $CUTOFF2 $MAP $optA $optB $optO $MAPPING_MIN_ALIGNMENT_SCORE $MAPPING_CLIPPING_PENALTY $MAPPING_MIN_QUALITY $SAMTOOLS_VIEW_F $SNP $POOLS $POOL_PLOIDY_FILE $PLOIDY $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION $MAIL
else
	#main
	echo ""; echo `date` " This is the HPC version of dDocent.  A config file must be specified"
	echo ""; echo `date` " Aborting dDocent 4 HPC"
	exit
fi

#Compress Large Leftover files
echo ""
echo -n "Compress large files"
date
gzip -f concat.fasta concat.seq rcluster.$CUTOFF.$CUTOFF2 rbdiv.$CUTOFF.$CUTOFF2.out rbasm.$CUTOFF.$CUTOFF2.out rainbow.$CUTOFF.$CUTOFF2.fasta reference.$CUTOFF.$CUTOFF2.fasta.original uniq.seqs uniq.$CUTOFF.$CUTOFF2.fasta totaluniqseq.$CUTOFF.$CUTOFF2 uniq..$CUTOFF.$CUTOFF2F.fasta uniq.$CUTOFF.$CUTOFF2.RC.fasta 2> /dev/null &
###############################################################################################

