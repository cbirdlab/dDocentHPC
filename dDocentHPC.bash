#!/usr/bin/env bash
VERSION=4.1
#This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD (SE or PE), ezRAD (SE or PE) data, or SE RAD data.
#It requires that your raw data are split up by tagged individual and follow the naming convention of:

#Pop_Sample1.F.fq and Pop_Sample1.R.fq

#Prints out title and contact info
echo; echo -e "\n* dDocentHPC v$VERSION Forked by cbird@tamucc.edu * \n"
#echo -e "Contact jpuritz@gmail.com with any problems \n\n "


#Determine which functions to run
#demultiplexFQ trimFQ	mkREF	mkBAM		fltrBAM		mkVCF	fltrVCF
if [ -n "$1" ]; then
	#getting functions from command line
	FUNKTION=$(echo $1)
	echo; echo "Running dDocentHPC $FUNKTION..."
else
	echo ""; echo `date` "ERROR:		dDocentHPC must be run with 2 arguments, "
	echo "			dDocentHPC.bash [function] [config file]"
	echo "			functions: trimFQ, mkREF, mkBAM, fltrBAM, mkVCF, fltrVCF"
	exit
fi


#dDocent can now accept a configuration file instead of running interactively
if [ -n "$2" ]; then
	echo " "
	echo `date` "Files output to: " $(pwd)
	
	echo ""; echo `date` "Reading config file... "
	echo " "
	#echo "BEGIN*CONFIG*FILE*****************************************************************************************"
#	cat $1
	cat $2
	#echo "*********END*CONFIG*FILE*******************************************************************************************"
	
	echo ""; echo `date` "Reading in variables from config file..."
	CONFIG=$2
	NUMProc=$(grep -A1 Processor $CONFIG | tail -1)
	MAXMemory=$(grep -A1 Memory $CONFIG | tail -1)
	#TRIM=$(grep -A1 Trim $CONFIG | tail -1)
	TRIM="RemoveThisVar"
	TRIM_LENGTH_ASSEMBLY=$(grep 'trimmomatic MINLEN (integer, mkREF only)' $CONFIG | awk '{print $1;}')
	SEED_ASSEMBLY=$(grep 'trimmomatic ILLUMINACLIP:<seed mismatches> (integer)' $CONFIG | awk '{print $1;}')
	PALIMDROME_ASSEMBLY=$(grep 'trimmomatic ILLUMINACLIP:<palindrome clip thresh> (integer)' $CONFIG | awk '{print $1;}')
	SIMPLE_ASSEMBLY=$(grep 'trimmomatic ILLUMINACLIP:<simple clip thresh> (integer)' $CONFIG | awk '{print $1;}')
	windowSize_ASSEMBLY=$(grep 'trimmomatic SLIDINGWINDOW:<windowSize> (integer)' $CONFIG | awk '{print $1;}')
	windowQuality_ASSEMBLY=$(grep 'trimmomatic SLIDINGWINDOW:<windowQuality> (integer)' $CONFIG | awk '{print $1;}')
	TRAILING_ASSEMBLY=$(grep 'trimmomatic SLIDINGWINDOW:<windowQuality> (integer)' $CONFIG | awk '{print $1;}')
	TRIM_LENGTH_MAPPING=$(grep 'trimmomatic MINLEN (integer, mkBAM only)' $CONFIG | awk '{print $1;}')
	LEADING_MAPPING=$(grep 'trimmomatic LEADING:<quality> (integer, mkBAM only)' $CONFIG | awk '{print $1;}')
	TRAILING_MAPPING=$(grep 'trimmomatic TRAILING:<quality> (integer, mkBAM only)' $CONFIG | awk '{print $1;}')
	HEADCROP=$(grep 'HEADCROP:<length> (integer, only Read1 for ezRAD)' $CONFIG | awk '{print $1;}')
	
	FixStacks=$(grep -A1 FixStacks $CONFIG | tail -1)
	ASSEMBLY=$(grep -A1 '^Assembly' $CONFIG | tail -1)
	ATYPE=$(grep -A1 Type $CONFIG | tail -1)
	simC=$(grep -A1 Simi $CONFIG | tail -1)
	HPC=$(grep -A1 HPC $CONFIG | tail -1)
	MANCUTOFF=$(grep -A1 Manually $CONFIG | tail -1)
	if [ "$MANCUTOFF" = "no" ]; then
		CUTOFF=$(grep -A1 ^Cutoff1 $CONFIG | tail -1)
		CUTOFF2=$(grep -A1 ^Cutoff2 $CONFIG | tail -1)
		#echo ""; echo `date` " CUTOFF=" $CUTOFF
		#echo ""; echo `date` " CUTOFF2=" $CUTOFF2
	else
		echo ""; echo `date` " ERROR: This is the HPC version of dDocent.  Manual cutoffs are not possible.  Please change config file settings."
		echo "  Aborting dDocentHPC"
		exit
	fi
	MAP=$(grep -A1 Mapping_Reads $CONFIG | tail -1)
	optA=$(grep -A1 _Match $CONFIG | tail -1)
	optB=$(grep -A1 MisMatch $CONFIG | tail -1)
	optO=$(grep -A1 Gap $CONFIG | tail -1)
	MAPPING_MIN_ALIGNMENT_SCORE=$(grep -A1 '^Mapping_Minimum_Alignment_Score' $CONFIG | tail -1)
	MAPPING_CLIPPING_PENALTY=$(grep -A1 '^Mapping_Clipping_Penalty' $CONFIG | tail -1)
	
	
	FILTERMAP=$(grep -A1 'Filter_Mapped_Read_Alignments' $CONFIG | tail -1)
	MAPPING_MIN_QUALITY=$(grep -A1 '^Mapping_Min_Quality' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F4=$(grep -A1 '^Remove_unmapped_reads' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F8=$(grep -A1 '^Remove_read_pair_if_one_is_unmapped' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F256=$(grep -A1 '^Remove_secondary_alignments' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F512=$(grep -A1 '^Remove_reads_not_passing_platform_vendor_filters' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F1024=$(grep -A1 '^Remove_PCR_or_optical_duplicates' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F2048=$(grep -A1 '^Remove_supplementary_alignments' $CONFIG | tail -1)
	SAMTOOLS_VIEW_f2=$(grep -A1 '^Keep_only_properly_aligned_read_pairs' $CONFIG | tail -1)

	if [ "$SAMTOOLS_VIEW_F4" == "yes" ]; then
		F4=4
	else
		F4=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F8" == "yes" ]; then
		F8=8
	else
		F8=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F256" == "yes" ]; then
		F256=256
	else
		F256=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F512" == "yes" ]; then
		F512=512
	else
		F512=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F1024" == "yes" ]; then
		F1024=1024
	else
		F1024=0
	fi

	if [ "$SAMTOOLS_VIEW_F2048" == "yes" ]; then
		F2048=2048
	else
		F2048=0
	fi
	
	SAMTOOLS_VIEW_F=$(($F4+$F8+$F256+$F512+$F1024+$F2048)) 
	
	if [ "$SAMTOOLS_VIEW_f2" == "yes" ]; then
		SAMTOOLS_VIEW_f=2
	else
		SAMTOOLS_VIEW_f=0
	fi
	
	SAMTOOLS_VIEW_Fcustom=$(grep -A1 '^Custom_samtools_view_F_bit_value' $CONFIG | tail -1)
	SAMTOOLS_VIEW_fcustom=$(grep -A1 '^Custom_samtools_view_f_bit_value' $CONFIG | tail -1)
	SOFT_CLIP_CUT=$(grep -A1 '^Remove_reads_with_excessive_soft_clipping' $CONFIG | tail -1)
	SOFT_CLIP_CUTOFF=$((($SOFT_CLIP_CUT+9)/10))
	FILTER_MIN_AS=$(grep -A1 '^Remove_reads_with_alignment_score_below' $CONFIG | tail -1)
	FILTER_ORPHANS=$(grep -A1 '^Remove_reads_orphaned_by_filters' $CONFIG | tail -1)

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
	FREEBAYES_Q=$(grep 'freebayes -Q --mismatch-base-quality-threshold' $CONFIG | awk '{print $1;}')
	FREEBAYES_U=$(grep 'freebayes -U --read-mismatch-limit' $CONFIG | awk '{print $1;}')
	MAIL=$(grep -A1 Email $CONFIG | tail -1)

else
	echo ""; echo `date` "ERROR:		dDocentHPC must be run with 2 arguments, "
	echo "			dDocentHPC.bash [function] [config file]"
	echo "			functions: trimFQ, mkREF, mkBAM, fltrBAM, mkVCF, fltrVCF"
	echo "			the default config file is called config.*"
	exit
fi


####################################################################################################
###Function to check for the required software for dDocent
####################################################################################################

CheckDependencies(){
	echo ""; echo " Running CheckDependencies Function..."
	DEP=$1
	EXIT_IF_DEP_MISSING=$2
	dDocentFUNCTION=$3
	NUMDEP=0
	for i in "${DEP[@]}"
	do
		
		if which $i &> /dev/null; then
			echo ""; echo "  The dependency $i is installed!"
		else
			echo ""; echo "  The dependency" $i "is not installed or is not in your" '$PATH'"."
			NUMDEP=$(($NUMDEP + 1))
		fi
		
		if [ "$i" == "freebayes" ]; then
			FREEB=(`freebayes | grep -oh 'v[0-9].*' | cut -f1 -d "." | sed 's/v//' `)	
			if [ "$FREEB" != "1" ]; then
				echo "  The version of FreeBayes installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install at least version 1.0.0"
				NUMDEP=$(($NUMDEP + 1))
			fi  
			
		#elif [ "$i" == "mawk" ]; then
			
		elif [ "$i" == "bwa" ]; then
			BWAV=$( bwa 2>&1 | mawk '/Versi/' | sed 's/Version: //g' | sed 's/0.7.//g' | sed 's/-.*//g' | cut -c 1-2 )
			if [ "$BWAV" -lt "13" ]; then
					echo "  The version of bwa installed in your" '$PATH' "is not optimized for dDocent."
					echo "  Please install at least version 0.7.13"
					NUMDEP=$(($NUMDEP + 1))
			fi
		
		elif [ "$i" == "samtools" ]; then
			SAMV1=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c1)
			SAMV2=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c3)
			if [ "$SAMV1"  -ge "1" ]; then
				if [ "$SAMV2"  -lt "3" ]; then
					echo "  The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
					echo "  Please install at least version 1.3.0"
				NUMDEP=$(($NUMDEP + 1))
				fi
			else
				echo "  The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install at least version 1.3.0"
				NUMDEP=$(($NUMDEP + 1))
			fi
		
		elif [ "$i" == "vcftools" ]; then
			VCFTV=$(vcftools | grep 'VCF' | sed -e 's/VCFtools (//' | sed -e 's/)//')   
			if [[ "$VCFTV" != 0.1.[1-9][0-9] ]]; then 
				echo "  The version of VCFtools installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install only version 0.1.11"
				NUMDEP=$(($NUMDEP + 1))
			fi

		elif [ "$i" == "rainbow" ]; then
			RAINV=(`rainbow | head -1 | cut -f2 -d' ' `)	
			if [[ "$RAINV" != "2.0.2" && "$RAINV" != "2.0.3" && "$RAINV" != "2.0.4" ]]; then
				echo "  The version of Rainbow installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install a version newer than 2.0.2"
				NUMDEP=$(($NUMDEP + 1))
			fi
		#elif [ "$i" == "gnuplot" ]; then
		
		#elif [ "$i" == "gawk" ]; then
		
		#elif [ "$i" == "seqtk" ]; then
		
		#elif [ "$i" == "cd-hit-est" ]; then
		
		#elif [ "$i" == "bamToBed" ]; then
		
		elif [ "$i" == "bedtools" ]; then
			BTC=$(bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." )
			if [[ "$BTC" != 2.[2-9][0-9] ]]; then
				echo "  The version of bedtools installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install only version 2.xx.0"
				NUMDEP=$(($NUMDEP + 1))
			fi
			BTC=$( bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." | sed 's/2\.//g' )
			if [ "$BTC" -ge "26" ]; then
				BEDTOOLSFLAG="NEW"
				elif [ "$BTC" == "23" ]; then
				BEDTOOLSFLAG="OLD"
				elif [ "$BTC" != "23" ]; then
				echo "  The version of bedtools installed in your" '$PATH' "is not optimized for dDocent."
				echo "  Please install version 2.23.0 or version 2.26.0 and above"
				NUMDEP=$(($NUMDEP + 1))	
			fi
		
		#elif [ "$i" == "coverageBed" ]; then
		
		#elif [ "$i" == "parallel" ]; then
		
		#elif [ "$i" == "vcfcombine" ]; then
		
		#elif [ "$i" == "bamtools" ]; then
		
		#elif [ "$i" == "pearRM" ]; then
		
		elif [ "$i" == "trimmomatic" ]; then
			if find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null| grep -q 'trim' ; then
				TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
			else
				echo " The dependency trimmomatic is not installed or is not in your" '$PATH'"."
				NUMDEP=$(($NUMDEP + 1))
			fi
				
			if find ${PATH//:/ } -maxdepth 2 -name TruSeq3-PE-2.fa 2> /dev/null | grep -q 'Tru' ; then
				ADAPTERS=$(find ${PATH//:/ } -maxdepth 2 -name TruSeq3-PE-2.fa 2> /dev/null | head -1)
			else
				echo " The file listing adapters (included with trimmomatic) is not installed or is not in your" '$PATH'"."
				NUMDEP=$(($NUMDEP + 1))
			fi
		fi
	done
	
	if [[ $NUMDEP -gt 0 && $EXIT_IF_DEP_MISSING == "TRUE" ]]; then
		echo -e "\n ERROR: Please install all required software then try again."
		exit 1
	elif [[ $NUMDEP -gt 0 && $EXIT_IF_DEP_MISSING == "FALSE" ]]; then
		echo -e "\n ERROR: Some software is not installed but $dDocentFUNCTION .  Please install all required software for full functionality"
	else
		echo -e "\n All dependencies are installed and up to date!"
	fi
}
echo ""; echo `date` "Checking for all required dDocent software..."
DEP=(freebayes mawk bwa samtools vcftools rainbow gnuplot gawk seqtk cd-hit-est bamToBed bedtools coverageBed parallel vcfcombine bamtools pearRM)
CheckDependencies $DEP "FALSE" "dDocentHPC will continue to run with limited functionality"

if ! awk --version | fgrep -v GNU &>/dev/null; then
	 awk=gawk
else
	 awk=awk
fi


#This code checks for individual fastq files follow the correct naming convention and are gziped
TEST=$(ls *.fq 2> /dev/null | wc -l )

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
	echo "";echo ""; echo `date` " Begin main ddocent function"
	##########User Input Section##########
	#This code gets input from the user and assigns variables
	######################################
	#RMH added ALL variables
	echo "";echo " Loading variable values..."
	echo "  FUNKTION=	${55}"
	echo "  NUMProc=	$1"
	echo "  MAXMemory=	$2"
	echo "  TRIM=	$3"
	echo "  TRIM_LENGTH_ASSEMBLY=	$4"
	echo "  SEED_ASSEMBLY=	$5"
	echo "  PALIMDROME_ASSEMBLY=	$6"
	echo "  SIMPLE_ASSEMBLY=	$7"
	echo "  windowSize_ASSEMBLY=	$8"
	echo "  windowQuality_ASSEMBLY=	$9"
	echo "  TRAILING_ASSEMBLY=	${10}"
	echo "  TRIM_LENGTH_MAPPING=	${11}"
	echo "  LEADING_MAPPING=	${12}"
	echo "  TRAILING_MAPPING=	${13}"
	echo "  HEADCROP=	${56}"
	echo "  FixStacks=	${14}"
	echo "  ASSEMBLY=	${15}"
	echo "  ATYPE=	${16}"
	echo "  simC=	${17}"
	echo "  HPC=	${18}"
	echo "  MANCUTOFF=	${19}"
	echo "  CUTOFF=	${20}"
	echo "  CUTOFF2=	${21}"
	echo "  MAP=	${22}"
	echo "  optA=	${23}"
	echo "  optB=	${24}"
	echo "  optO=	${25}"
	echo "  MAPPING_MIN_ALIGNMENT_SCORE=	${26}"
	echo "  MAPPING_CLIPPING_PENALTY=	${27}"
	
	echo "  FILTERMAP=	${45}"
	
	echo "  MAPPING_MIN_QUALITY=	${28}"
	
	echo "  SAMTOOLS_VIEW_f=	${46}"
	
	echo "  SAMTOOLS_VIEW_F=	${29}"
	
	echo "  SAMTOOLS_VIEW_Fcustom=	${47}"
	echo "  SAMTOOLS_VIEW_fcustom=	${48}"
	echo "  SOFT_CLIP_CUTOFF=	${49}"
	echo "  FILTER_MIN_AS=	${52}"
	echo "  FILTER_ORPHANS=	${50}"
	
	echo "  SNP=	${30}"
	echo "  POOLS=	${31}"
	echo "  POOL_PLOIDY_FILE=	${32}"
	echo "  PLOIDY=	${33}"
	echo "  BEST_N_ALLELES=	${34}"
	echo "  MIN_MAPPING_QUAL=	${35}"
	echo "  MIN_BASE_QUAL=	${36}"
	echo "  HAPLOTYPE_LENuuGTH=	${37}"
	echo "  MIN_REPEAT_ENTROPY=	${38}"
	echo "  MIN_COVERAGE=	${39}"
	echo "  MIN_ALT_FRACTION=	${40}"
	echo "  MIN_ALT_COUNT=	${41}"
	echo "  MIN_ALT_TOTAL=	${42}"
	echo "  READ_MAX_MISMATCH_FRACTION=	${43}"
	echo "  FREEBAYES_Q=	${53}"
	echo "  FREEBAYES_U=	${54}"
	echo "  MAIL=	${44}"
	echo "  BEDTOOLSFLAG=	${51}"
	
	#Load config arguments into variables
	# RMH added ALL variables
	FUNKTION=${55}
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
	HEADCROP=${56}
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
	FILTERMAP=${45}
	MAPPING_MIN_QUALITY=${28}
	SAMTOOLS_VIEW_f=${46}
	SAMTOOLS_VIEW_F=${29}
	SAMTOOLS_VIEW_Fcustom=${47}
	SAMTOOLS_VIEW_fcustom=${48}
	SOFT_CLIP_CUTOFF=${49}
	FILTER_MIN_AS=${52}
	FILTER_ORPHANS=${50}
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
	FREEBAYES_Q=${53}
	FREEBAYES_U=${54}
	MAIL=${44}
	BEDTOOLSFLAG=${51}

	CUTOFFS=$CUTOFF.$CUTOFF2
	
	#$$$$$$$$$$$$$$$$$$$$$CEB$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#Count number of individuals in current directory

	#this alerts user to fastq file naming conventions
		echo "";echo " The HPC version of dDocent will only digest files with particular extensions for particular tasks"
		echo "  untouched files for trimming must be *.F.fq.gz and *.R.fq.gz"
		echo "  files trimmed for assembly must be *r1.fq.gz *r2.fq.gz"
		echo "  files trimmed for mapping must be *R1.fq.gz *R2.fq.gz"
		
#	if [ "$TRIM" == "yes" ]; then
	if [ "$FUNKTION" == "trimFQ" ]; then
		F="F"
		Fwild="*.F.fq.gz"
		Fsed=".F.fq.gz"
		R="R"
		Rwild="*.R.fq.gz"
		Rsed=".R.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
#	elif [ "$ASSEMBLY" == "yes" ]; then
	elif [ "$FUNKTION" == "mkREF" ]; then
		F="r1"
		Fwild="*.r1.fq.gz"
		Fsed=".r1.fq.gz"
		R="r2"
		Rwild="*.r2.fq.gz"
		Rsed=".r2.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
#	elif [ "$MAP" == "yes" ]; then
	elif [ "$FUNKTION" == "mkBAM" ]; then
		F="R1"
		Fwild="*.R1.fq.gz"
		Fsed=".R1.fq.gz"
		R="R2"
		Rwild="*.R2.fq.gz"
		Rsed=".R2.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
#	elif [ "$FILTERMAP" == "yes" ]; then
	elif [ "$FUNKTION" == "fltrBAM" ]; then
		F="R1"
		Fwild="*-RAW.bam"
		Fsed=".${CUTOFFS}-RAW.bam"
		R="R2"
		Rwild="*.R2.fq.gz"
		Rsed=".R2.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
#	elif [ "$SNP" == "yes" ]; then
	elif [ "$FUNKTION" == "mkVCF" ]; then
		F="R1"
		Fwild="*-RG.bam"
		Fsed=".${CUTOFFS}-RG.bam"
		R="R2"
		Rwild="*.R2.fq.gz"
		Rsed=".R2.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
	else
		echo "";echo "  Couldn't decide which file extensions were selected"
		F="r1"
		Fwild="*.r1.fq.gz"
		Fsed=".r1.fq.gz"
		R="r2"
		Rwild="*.r2.fq.gz"
		Rsed=".r2.fq.gz"
		echo "";echo "  extensions selected: $Fwild $Rwild"
	fi

	NumInd=$(ls $Fwild | wc -l)
	NumInd=$(($NumInd - 0))

	#Create list of sample names
	if [ ! -s "namelist.$CUTOFF.$CUTOFF2" ];then
		ls $Fwild > namelist.$CUTOFFS
		sed -i'' -e "s/$Fsed//g" namelist.$CUTOFFS
	else
		echo "";echo " The namelist file already exists and was not recreated. "
		echo "  If you experience errors, you should delete the namelist file."
	fi

	#Create an array of sample names
	#NUMNAMES=$(mawk '/_/' namelist.$CUTOFFS | wc -l)
	NUMNAMES=$(grep -c '^' namelist.$CUTOFFS)

	if [ "$NUMNAMES" -eq "$NumInd" ]; then
		NAMES=( `cat "namelist.$CUTOFFS" `)
		echo " ";echo " The samples being processed are:"
		cat namelist.$CUTOFFS | parallel --no-notice "echo '  '{}"
		echo ""
	else
		echo ""; echo " "`date` " ERROR: NUMNAMES=$NUMNAMES NUMIND=$NumInd"
		echo "  Individuals do not follow the dDocent naming convention."
		echo "  Please rename individuals to: Locality_Individual.F.fq.gz"
		echo "  For example: LocA_001.F.fq.gz"
		exit 1
	fi

	#Sets a start time variable
	STARTTIME=$(date)

	##Section of logic statements that dictates the order and function of processing the pipeline
	
	if [ "$FixStacks" == "yes" ]; then
		FIXSTACKS
	fi
	
#	if [ "$TRIM" == "yes" ]; then
	if [ "$FUNKTION" == "trimFQ" ]; then
		echo " ";echo " "`date` "Trimming reads " 
		TrimReadsRef #& 2> ./mkREF/trimref.log
		TrimReads #& 2> ./mkBAM/trim.log
	fi

#	elif [ "$ASSEMBLY" == "yes" ]; then
	if [ "$FUNKTION" == "mkREF" ]; then
		Assemble
	fi

#	elif [ "$MAP" == "yes" ]; then
	if [ "$FUNKTION" == "mkBAM" ]; then
		MAP2REF $CUTOFF.$CUTOFF2 $NUMProc $CONFIG $ATYPE $Rsed NAMES
	fi

#	elif [ "$FILTERMAP" == "yes" ]; then
	if [ "$FUNKTION" == "fltrBAM" ]; then
		FILTERBAM $CUTOFFS $NUMProc $CONFIG $ATYPE
	fi
	
#	elif [ "$SNP" == "yes" ]; then
	if [ "$FUNKTION" == "mkVCF" ]; then
		GENOTYPE $CUTOFFS $NUMProc $CONFIG
	fi

	##Checking for possible errors

#	if [ "$MAP" != "no" ]; then
	if [ "$FUNKTION" == "mkBAM" ]; then
		ERROR1=$(mawk '/developer/' bwa* | wc -l 2>/dev/null) 
	fi
	#CEB this was causing an error, there's no bam.log created
	#	ERROR2=$(mawk '/error/' *.bam.log | wc -l 2>/dev/null)
	#ERRORS=$(($ERROR1 + $ERROR2))

	#Move various log files to own directory
	if [ ! -d "logfiles" ]; then
		mkdir logfiles
	fi
	mv *.log log ./logfiles 2> /dev/null

	#Sending a completion email

	#CEB this was causing error: mawk: cannot open ./logfiles/Final.log (No such file or directory)
	#if [ $ERRORS -gt 0 ]; then
	#			echo -e "dDocent has finished with errors in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\nPlease check log files\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/Final.log` "\n\ndDocent 2.24 \nThe 'd' is silent, hillbilly." | mailx -s "dDocent has finished with ERRORS!" $MAIL
	#else
	#		echo -e "dDocent has finished with an analysis in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/Final.log` "\n\ndDocent 2.24 \nIt is pronounced Dee-Docent, professor." | mailx -s "dDocent has finished" $MAIL
	#fi


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
	echo "HEADCROP" >> dDocent.runs
	echo $HEADCROP >> dDocent.runs
	echo "FixStacks" >> dDocent.runs
	echo $FixStacks >> dDocent.runs
	echo "Assembly?" >> dDocent.runs
#	echo $ASSEMBLY >> dDocent.runs
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
#	echo $MAP >> dDocent.runs
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
#	echo $SNP >> dDocent.runs
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




##############################################################################################
#Function for correcting stacks file errors (antiquated?)
###############################################################################################

FIXSTACKS () {
	echo "";echo " "`date` " Fixing Stacks anomalies"

	#STACKS adds a strange _1 or _2 character to the end of processed reads, this looks for checks for errant characters and replaces them.
	#This functionality is now parallelized and will run if only SE sequences are used.

	STACKS=$(cat namelist.$CUTOFFS| parallel -j $NUMProc --no-notice "zcat {}$Fsed | head -1" | mawk '!/\/1/' | wc -l)
	FB1=$(( $NUMProc / 2 ))
	if [ $STACKS -gt 0 ]; then
		echo "  Removing the _1 character and replacing with /1 in the name of every sequence"
		cat namelist.$CUTOFFS | parallel -j $FB1 --no-notice "zcat {}.$F.fq.gz | sed -e 's:_1$:/1:g' > {}.$F.fq"
		rm $Fwild
		cat namelist.$CUTOFFS | parallel -j $FB1 --no-notice "gzip {}.$F.fq"
	fi

	if [ -f "${NAMES[@]:(-1)}"$Rsed ]; then
		
		STACKS=$(cat namelist.$CUTOFFS| parallel -j $NUMProc --no-notice "zcat {}$Rsed | head -1" | mawk '!/\/2/' | wc -l)

		if [ $STACKS -gt 0 ]; then
			echo "  Removing the _2 character and replacing with /2 in the name of every sequence"
			cat namelist.$CUTOFFS | parallel -j $FB1 --no-notice "zcat {}.$R.fq.gz | sed -e 's:_2$:/2:g' > {}.$R.fq"
			rm $Rwild
			cat namelist.$CUTOFFS | parallel -j $FB1 --no-notice "gzip {}.$R.fq"
		fi
	fi
}



##############################################################################################
#Function for trimming reads using trimmomatic
###############################################################################################

TrimReadsRef () { 
	echo " "
	echo `date` "Trimming reads for reference genome"

	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
	ADAPTERS=$(find ${PATH//:/ } -maxdepth 2 -name TruSeq3-PE-2.fa 2> /dev/null | head -1)
	
	echo "TRIMMOMATIC=	$TRIMMOMATIC"
	echo "ADAPTERS=	$ADAPTERS"
	echo "NUMProc=	$NUMProc"
	echo "SEED_ASSEMBLY=	$SEED_ASSEMBLY"
	echo "PALIMDROME_ASSEMBLY=	$PALIMDROME_ASSEMBLY"
	echo "SIMPLE_ASSEMBLY=	$SIMPLE_ASSEMBLY"
	echo "TRAILING_ASSEMBLY=	$TRAILING_ASSEMBLY"
	echo "windowSize_ASSEMBLY=	$windowSize_ASSEMBLY"
	echo "windowQuality_ASSEMBLY=	$windowQuality_ASSEMBLY"
	echo "TRIM_LENGTH_ASSEMBLY=	$TRIM_LENGTH_ASSEMBLY"
	echo "F=	$F"
	echo "R=	$R"

	if [ ! -d "mkREF" ]; then mkdir mkREF &>/dev/null; fi
	if [ ! -d "./mkREF/unpaired" ]; then mkdir ./mkREF/unpaired &>/dev/null; fi
	if [ ! -d "./mkREF/logs" ]; then mkdir ./mkREF/logs &>/dev/null; fi
	

	if [ -f "${NAMES[1]}".$R.fq.gz ]; then
		Proc=$((NUMProc/4))
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $Proc "java -jar $TRIMMOMATIC PE -threads 4 -phred33 {}.F.fq.gz {}.R.fq.gz ./mkREF/{}.r1.fq.gz ./mkREF/unpaired/{}.unpairedF.fq.gz ./mkREF/{}.r2.fq.gz ./mkREF/unpaired/{}.unpairedR.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY HEADCROP:$HEADCROP TRAILING:$TRAILING_ASSEMBLY SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY CROP:$TRIM_LENGTH_ASSEMBLY MINLEN:$TRIM_LENGTH_ASSEMBLY  &> ./mkREF/logs/{}.trim.log"
	else 
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $NUMProc "java -jar $TRIMMOMATIC SE -threads 1 -phred33 {}.F.fq.gz ./mkREF/{}.r1.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY HEADCROP:$HEADCROP TRAILING:$TRAILING_ASSEMBLY SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY CROP:$TRIM_LENGTH_ASSEMBLY MINLEN:$TRIM_LENGTH_ASSEMBLY &> ./mkREF/logs/{}.trim.log"
	fi
	
	if [[ $HEADCROP != 0 ]]; then
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $NUMProc "java -jar $TRIMMOMATIC SE -threads 1 -phred33 ./mkREF/{}.r1.fq.gz ./mkREF/{}.r1.fq.gz.ezRAD HEADCROP:$HEADCROP &> ./mkREF/logs/{}.trim.log"
		mkdir ./mkREF/unheadcropped
		ls ./mkREF/*r1.fq.gz | parallel --no-notice "mv {} ./mkREF/unheadcropped"
		rename .fq.gz.ezRAD .fq.gz ./mkREF/*.fq.gz.ezRAD
	fi
	
	#RMH housekeeping
#	rm *_r1.fq.gz
}


##############################################################################################
#Function for trimming reads using trimmomatic
###############################################################################################

TrimReads () { 
	echo " "
	echo; echo `date` "Trimming reads for mapping"

	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
	ADAPTERS=$(find ${PATH//:/ } -maxdepth 2 -name TruSeq3-PE-2.fa 2> /dev/null | head -1)
	
	if [ ! -d "mkBAM" ]; then mkdir mkBAM &>/dev/null; fi
	if [ ! -d "./mkBAM/unpaired" ]; then mkdir ./mkBAM/unpaired &>/dev/null; fi
	if [ ! -d "./mkBAM/logs" ]; then mkdir ./mkBAM/logs &>/dev/null; fi

	if [ -f "${NAMES[1]}".$R.fq.gz ]; then
		Proc=$((NUMProc/4))
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $Proc "java -jar $TRIMMOMATIC PE -threads 4 -phred33 {}.F.fq.gz {}.R.fq.gz ./mkBAM/{}.R1.fq.gz ./mkBAM/unpaired/{}.unpairedF.fq.gz ./mkBAM/{}.R2.fq.gz ./mkBAM/unpaired/{}.unpairedR.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY HEADCROP:$HEADCROP LEADING:$LEADING_MAPPING TRAILING:$TRAILING_MAPPING SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY MINLEN:$TRIM_LENGTH_MAPPING &> ./mkBAM/logs/{}.trim.log"
	else 
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $NUMProc "java -jar $TRIMMOMATIC SE -threads 1 -phred33 {}.F.fq.gz ./mkBAM/{}.R1.fq.gz ILLUMINACLIP:$ADAPTERS:$SEED_ASSEMBLY:$PALIMDROME_ASSEMBLY:$SIMPLE_ASSEMBLY HEADCROP:$HEADCROP LEADING:$LEADING_MAPPING TRAILING:$TRAILING_MAPPING SLIDINGWINDOW:$windowSize_ASSEMBLY:$windowQuality_ASSEMBLY MINLEN:$TRIM_LENGTH_MAPPING &> ./mkBAM/logs/{}.trim.log"
	fi 
	
	if [[ $HEADCROP != 0 ]]; then
		echo "${NAMES[@]}" | sed 's/ /\n/g' | parallel --no-notice -j $NUMProc "java -jar $TRIMMOMATIC SE -threads 1 -phred33 ./mkREF/{}.r1.fq.gz ./mkREF/{}.r1.fq.gz.ezRAD HEADCROP:$HEADCROP &> ./mkREF/logs/{}.trim.log"
		mkdir ./mkREF/unheadcropped
		ls ./mkREF/*r1.fq.gz | parallel --no-notice "mv {} ./mkREF/unheadcropped"
		rename .fq.gz.ezRAD .fq.gz ./mkREF/*.fq.gz.ezRAD
	fi
}


###############################################################################################
#Main function for assembly
###############################################################################################

Assemble(){
	AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
	AWK2='!/>/'
	AWK3='!/NNN/'
	PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
	SED1='s/^[ 	]*//'
	SED2='s/ /	/g'
	FRL=$(zcat ${NAMES[0]}.r1.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)

	if [ ! -s "namelistfr.$CUTOFFS" ];then
		ls $Fwild > namelistfr.$CUTOFFS
		ls $Rwild >> namelistfr.$CUTOFFS
		sed -i'' -e "s/\.fq\.gz//g" namelistfr.$CUTOFFS
	else
		echo "";echo " The namelistfr file already exists and was not recreated."
		echo "  If you experience errors, you should delete the namelistfr file."
	fi

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
						cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' | sort | uniq -c -w $FRL| sed -e '$SED1' | sed -e '$SED2' > {}.uniq.seqs"
					else
						echo "";echo `date` " PE assembly"
						cat namelist.$CUTOFF.$CUTOFF2 | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
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
					# for i in "${NAMES[@]}";
						# do
						# #zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
						# zcat $i$Rsed | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt
					# done	
					#parallel code ceb aug 2018
					parallel --no-notice -j $NUMProc "zcat {}$Rsed | head -2 | tail -1 >> lengths.$CUTOFF.$CUTOFF2.txt" ::: "${NAMES[@]}"
					
					MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.$CUTOFF.$CUTOFF2.txt| head -1)
					LENGTH=$(( $MaxLen / 3))
					echo "";echo `date` " OL assembly: PEAR "
					# for i in "${NAMES[@]}"
						# do
						# #pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH 
						# pearRM -f $i$Fsed -r $i$Rsed -o $i -j $NUMProc -n $LENGTH -p 0.0001
					# done
					#parallel code ceb aug 2018; need to check and see how many threads pear can use
					parallel --no-notice -j "$NUMProc pearRM -f {}$Fsed -r {}$Rsed -o {} -j $NUMProc -n $LENGTH -p 0.0001" ::: "${NAMES[@]}"
					
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
		sed -e 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/	/g' uniq.$CUTOFF.$CUTOFF2.fasta | cut -f1 > uniq.$CUTOFF.$CUTOFF2.F.fasta
		CDHIT=$(python -c "print max("$simC" - 0.1,0.8)")
		cd-hit-est -i uniq.$CUTOFF.$CUTOFF2.F.fasta -o xxx.$CUTOFF.$CUTOFF2 -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.$CUTOFF.$CUTOFF2.log
		mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.$CUTOFF.$CUTOFF2.clstr | sed 's/[>dDococent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids.$CUTOFF.$CUTOFF2
		paste sort.contig.cluster.ids.$CUTOFF.$CUTOFF2 totaluniqseq.$CUTOFF.$CUTOFF2 > contig.cluster.totaluniqseq.$CUTOFF.$CUTOFF2
		sort -k2,2 -g contig.cluster.totaluniqseq.$CUTOFF.$CUTOFF2 | sed -e 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/	/g' > rcluster.$CUTOFF.$CUTOFF2
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
			else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
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
		sed -e 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/	/g' rainbow.$CUTOFF.$CUTOFF2.fasta | cut -f1 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.$CUTOFF.$CUTOFF2.F.fq
		sed -e 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/	/g' rainbow.$CUTOFF.$CUTOFF2.fasta | cut -f2 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.$CUTOFF.$CUTOFF2.R.fq

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
		paste other.$CUTOFF.$CUTOFF2.F other.$CUTOFF.$CUTOFF2.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' > other.$CUTOFF.$CUTOFF2.FR

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
#Map Reads 2 Reference
###############################################################################################

MAP2REF(){

	CUTOFFS=$1
	NUMProc=$2
	CONFIG=$3
	ATYPE=$4
	Rsed=$5
	names=$6[@]
	NAMES=("${!names}")

	optA=$(grep -A1 _Match $CONFIG | tail -1)
	optB=$(grep -A1 MisMatch $CONFIG | tail -1)
	optO=$(grep -A1 Gap $CONFIG | tail -1)
	MAPPING_MIN_ALIGNMENT_SCORE=$(grep -A1 '^Mapping_Minimum_Alignment_Score' $CONFIG | tail -1)
	MAPPING_CLIPPING_PENALTY=$(grep -A1 '^Mapping_Clipping_Penalty' $CONFIG | tail -1)

	echo " ";echo `date` " Using BWA to map reads."
	if [ reference.$CUTOFFS.fasta -nt reference.$CUTOFFS.fasta.fai ]; then
		samtools faidx reference.$CUTOFFS.fasta
		bwa index reference.$CUTOFFS.fasta &> index.$CUTOFFS.log
	fi
	#dDocent now checks for trimmed read files before attempting mapping
	if [ ! -f "${NAMES[@]:(-1)}"$Rsed ]; then
		echo "dDocent cannot locate trimmed reads files"
		echo "Please rerun dDocent with quality trimming"
		exit 1
	fi
	#This next section of code checks to see if the reference was assembled by dDocent 
	#and if so, modifies the expected insert length distribution for BWA's metric for proper pairing
	echo ""
	echo -n `date` " was ref assembled by dDocent?"
	if head -1 reference.$CUTOFFS.fasta | grep -e 'dDocent' reference.$CUTOFFS.fasta 1>/dev/null; then
		echo " Yes. Modifying expected insert length for bwa"
		rm lengths.$CUTOFFS.txt &> /dev/null
		for i in "${NAMES[@]}";
			do
			if [ -f "$i.$R.fq.gz" ]; then
				#echo ""
				#echo troubleshooting
				#echo $i
				#echo $R
				#echo $i.$R.fq.gz
				zcat $i.$R.fq.gz | head -2 | tail -1 >> lengths.$CUTOFFS.txt
			fi
		done	
		if [ -f "lengths.$CUTOFFS.txt" ]; then
			MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.$CUTOFFS.txt| head -1)
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
				#bwa mem reference.$CUTOFFS.fasta $i.R1.fq.gz $i.R2.fq.gz -L $MAPPING_CLIPPING_PENALTY -I $INSERT,$SD,$INSERTH,$INSERTL -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFFS.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q $MAPPING_MIN_QUALITY -f 3 -F $SAMTOOLS_VIEW_F -SbT reference.$CUTOFFS.fasta - > $i.$CUTOFFS.bam 2>$i.$CUTOFFS.bam.log
				#CEB: updated to output minimially-filtered bamfiles here
				bwa mem reference.$CUTOFFS.fasta $i.R1.fq.gz $i.R2.fq.gz -L $MAPPING_CLIPPING_PENALTY -I $INSERT,$SD,$INSERTH,$INSERTL -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFFS.log | samtools view -@$NUMProc -SbT reference.$CUTOFFS.fasta - > $i.$CUTOFFS.bam 2>$i.$CUTOFFS.bam.log
				
			else
				bwa mem reference.$CUTOFFS.fasta $i.R1.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFFS.log | samtools view -@$NUMProc -SbT reference.$CUTOFFS.fasta - > $i.$CUTOFFS.bam 2>$i.$CUTOFFS.bam.log
			fi
			echo ""
			echo `date` " run samtools sort" $i
			samtools sort -@$NUMProc $i.$CUTOFFS.bam -o $i.$CUTOFFS.bam 
			mv $i.$CUTOFFS.bam $i.$CUTOFFS-RAW.bam
			echo ""
			echo `date` " run samtools index" $i
			samtools index $i.$CUTOFFS-RAW.bam
		done
	else
		echo " No. Expected insert length not modified for BWA"
		echo ""
		echo `date` " Run bwa mem on non-dDocent files"
		for i in "${NAMES[@]}"
		do
			if [ -f "$i.R2.fq.gz" ]; then
				bwa mem reference.$CUTOFFS.fasta $i.R1.fq.gz $i.R2.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFFS.log | samtools view -@$NUMProc -SbT reference.$CUTOFFS.fasta - > $i.$CUTOFFS.bam 2>$i.$CUTOFFS.bam.log
			else
				bwa mem reference.$CUTOFFS.fasta $i.R1.fq.gz -L $MAPPING_CLIPPING_PENALTY -t $NUMProc -a -M -T $MAPPING_MIN_ALIGNMENT_SCORE -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.$CUTOFFS.log | samtools view -@$NUMProc -SbT reference.$CUTOFFS.fasta - > $i.$CUTOFFS.bam 2>$i.$CUTOFFS.bam.log
			fi
			samtools sort -@$NUMProc $i.$CUTOFFS.bam -o $i.$CUTOFFS.bam 
			mv $i.$CUTOFFS.bam $i.$CUTOFFS-RAW.bam
			samtools index $i.$CUTOFFS-RAW.bam
		done
	fi
}

###############################################################################################
#FILTER BAM FILES
###############################################################################################

function FILTERBAM(){
	CUTOFFS=$1
	NUMProc=$2
	CONFIG=$3
	ATYPE=$4

	MAPPING_MIN_QUALITY=$(grep -A1 '^Mapping_Min_Quality' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F4=$(grep -A1 '^Remove_unmapped_reads' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F8=$(grep -A1 '^Remove_read_pair_if_one_is_unmapped' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F256=$(grep -A1 '^Remove_secondary_alignments' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F512=$(grep -A1 '^Remove_reads_not_passing_platform_vendor_filters' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F1024=$(grep -A1 '^Remove_PCR_or_optical_duplicates' $CONFIG | tail -1)
	SAMTOOLS_VIEW_F2048=$(grep -A1 '^Remove_supplementary_alignments' $CONFIG | tail -1)
	SAMTOOLS_VIEW_f2=$(grep -A1 '^Keep_only_properly_aligned_read_pairs' $CONFIG | tail -1)
	if [ "$SAMTOOLS_VIEW_F4" == "yes" ]; then
		F4=4
	else
		F4=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F8" == "yes" ]; then
		F8=8
	else
		F8=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F256" == "yes" ]; then
		F256=256
	else
		F256=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F512" == "yes" ]; then
		F512=512
	else
		F512=0
	fi
	
	if [ "$SAMTOOLS_VIEW_F1024" == "yes" ]; then
		F1024=1024
	else
		F1024=0
	fi

	if [ "$SAMTOOLS_VIEW_F2048" == "yes" ]; then
		F2048=2048
	else
		F2048=0
	fi
	
	SAMTOOLS_VIEW_F=$(($F4+$F8+$F256+$F512+$F1024+$F2048)) 
	
	if [ "$SAMTOOLS_VIEW_f2" == "yes" ]; then
		SAMTOOLS_VIEW_f=2
	else
		SAMTOOLS_VIEW_f=0
	fi
	
	SAMTOOLS_VIEW_Fcustom=$(grep -A1 '^Custom_samtools_view_F_bit_value' $CONFIG | tail -1)
	SAMTOOLS_VIEW_fcustom=$(grep -A1 '^Custom_samtools_view_f_bit_value' $CONFIG | tail -1)
	SOFT_CLIP_CUT=$(grep -A1 '^Remove_reads_with_excessive_soft_clipping' $CONFIG | tail -1)
	SOFT_CLIP_CUTOFF=$((($SOFT_CLIP_CUT+9)/10))
	FILTER_MIN_AS=$(grep -A1 '^Remove_reads_with_alignment_score_below' $CONFIG | tail -1)
	FILTER_ORPHANS=$(grep -A1 '^Remove_reads_orphaned_by_filters' $CONFIG | tail -1)

	echo "";echo " "`date` "Filtering raw BAM Files"
	if [ "$ATYPE" == "PE" ]; then 	#paired end alignments
		#Filter 1: remove reads based on samtools flags
			echo "";echo "  "`date` " Applying Filter 1: removing paired reads mapping to different contigs, secondary, and supplementary alignments"
			BITS=$(($SAMTOOLS_VIEW_F+$SAMTOOLS_VIEW_f2))
			BITScustom=$(($SAMTOOLS_VIEW_Fcustom+$SAMTOOLS_VIEW_fcustom))
			#if [ $FILTER_MIN_AS < 100 ]; then
				#FILTER_MIN_AS_TENS=$((($FILTER_MIN_AS+9)/10))
				if [[ "$BITS" != "0" && "$BITScustom" != "0" ]]; then
					ls *$CUTOFFS-RAW.bam | sed 's/\-RAW.bam//g' | parallel --no-notice -j $NUMProc "samtools view -h -q $MAPPING_MIN_QUALITY -F $SAMTOOLS_VIEW_F -f $SAMTOOLS_VIEW_f2 {}-RAW.bam | samtools view -Sh1 -F $SAMTOOLS_VIEW_Fcustom -f SAMTOOLS_VIEW_fcustom - -o {}-RG.bam "
				elif [[ "$BITS" != 0 && "$BITScustom" == 0 ]]; then
					ls *$CUTOFFS-RAW.bam | sed 's/\-RAW.bam//g' | parallel --no-notice -j $NUMProc "samtools view -h1 -q $MAPPING_MIN_QUALITY -F $SAMTOOLS_VIEW_F -f SAMTOOLS_VIEW_f2 {}-RAW.bam -o {}-RG.bam "
				elif [[ "$BITS" == 0 && "$BITScustom" != 0 ]]; then
					ls *$CUTOFFS-RAW.bam | sed 's/\-RAW.bam//g' | parallel --no-notice -j $NUMProc "samtools view -h1 -q $MAPPING_MIN_QUALITY -F $SAMTOOLS_VIEW_Fcustom -f SAMTOOLS_VIEW_fcustom {}-RAW.bam -o {}-RG.bam "
				elif [[ "$BITScustom" == 0 && "$BITS" == 0 ]]; then
					ls *$CUTOFFS-RAW.bam | sed 's/\-RAW.bam//g' | parallel --no-notice -j $NUMProc "samtools view -h1 -q $MAPPING_MIN_QUALITY {}-RAW.bam -o {}-RG.bam "
				fi
			#else
			
			#fi
		
		#Filter 2: remove reads with excessive soft clipping and orphans
			if [[ "$SOFT_CLIP_CUTOFF" != "no" || "$FILTER_ORPHANS" != "no" ]]; then
				#Function for filtering BAM files
				SoftClipOrphanFilter(){
					if [[ "$2" != "no" && "$3" == "yes" ]]; then
						samtools view $1 | mawk '!/(\t([$2-9].|[1-9][0-9][0-9])S|([$2-9].|[1-9][0-9][0-9])S\t)/' | awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] > 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' | cat <(samtools view -H $1) - | samtools view -S1T $4 - | samtools sort - -o $1 
					elif [[ "$2" != "no" && "$3" == "no" ]]; then
						samtools view -h $1 | mawk '!/(\t([$2-9].|[1-9][0-9][0-9])S|([$2-9].|[1-9][0-9][0-9])S\t)/' | samtools view -S1 - -o $1
					elif [[ "$2" == "no" && "$3" == "yes" ]]; then
						samtools view $1 | awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] > 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' | cat <(samtools view -H $1) - | samtools view -S1T $4 - | samtools sort - -o $1
					fi
				}
				export -f SoftClipOrphanFilter
				
				echo "";echo "  "`date` " Applying Filter 2: removing excessively soft clipped reads (and their mates)"
				echo "   "`date` " SOFT_CLIP_CUTOFF is $SOFT_CLIP_CUTOFF * 10"
				ls *$CUTOFFS-RG.bam | parallel --no-notice -j $NUMProc "SoftClipOrphanFilter {} $SOFT_CLIP_CUTOFF $FILTER_ORPHANS reference.$CUTOFFS.fasta $FILTER_MIN_AS"
			fi
		
		#Index the filtered bam files 
			echo "";echo "  "`date` " Indexing the filtered BAM files"
			ls *$CUTOFFS-RG.bam | parallel --no-notice -j $NUMProc "samtools index {}" 
		
	#elif [ "$ATYPE" == "OL" ]; then					#single end alignments, -f2 turned off
		
	# elif [ "$ATYPE" == "RPE" ]; then					#single end alignments
		
	# elif [ "$ATYPE" == "SE" ]; then					#single end alignments	
		
	fi
}

###############################################################################################
#Call SNPs
###############################################################################################

function GENOTYPE(){
	echo ""; echo `date` "Genotyping initiated..."	
	CUTOFFS=$1
	NUMProc=$2
	CONFIG=$3
	
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
	#MIN_ALT_COUNT=$(grep -A1 '^Min_Alternate_Count' $CONFIG | tail -1)
	#MIN_ALT_TOTAL=$(grep -A1 '^Min_Alternate_Total' $CONFIG | tail -1)
	#READ_MAX_MISMATCH_FRACTION=$(grep -A1 '^Read_Max_Mismatch_Fraction' $CONFIG | tail -1)

	FREEBAYES_z=$(grep 'freebayes -z --read-max-mismatch-fraction' $CONFIG | awk '{print $1;}')
	FREEBAYES_C=$(grep 'freebayes -C --min-alternate-count' $CONFIG | awk '{print $1;}')
	FREEBAYES_3=$(grep 'freebayes ~3 ~~min-alternate-qsum' $CONFIG | awk '{print $1;}'); FREEBAYES_3=$((FREEBAYES_3 * FREEBAYES_C))
	FREEBAYES_G=$(grep 'freebayes -G --min-alternate-total' $CONFIG | awk '{print $1;}')
	FREEBAYES_Q=$(grep 'freebayes -Q --mismatch-base-quality-threshold' $CONFIG | awk '{print $1;}')
	FREEBAYES_U=$(grep 'freebayes -U --read-mismatch-limit' $CONFIG | awk '{print $1;}')
	FREEBAYES_DOLLAR=$(grep 'freebayes -\$ --read-snp-limit' $CONFIG | awk '{print $1;}')
	FREEBAYES_e=$(grep 'freebayes -e --read-indel-limit' $CONFIG | awk '{print $1;}')

	FREEBAYES_w=$(grep 'freebayes -w --hwe-priors-off' $CONFIG | awk '{print $1;}'); if [ $FREEBAYES_w == "no" ]; then FREEBAYES_w=""; else FREEBAYES_w="-w "; fi
	FREEBAYES_V=$(grep 'freebayes -V --binomial-obs-priors-off' $CONFIG | awk '{print $1;}'); if [ $FREEBAYES_V == "no" ]; then FREEBAYES_V=""; else FREEBAYES_V="-V "; fi
	FREEBAYES_a=$(grep 'freebayes -a --allele-balance-priors-off' $CONFIG | awk '{print $1;}'); if [ $FREEBAYES_a == "no" ]; then FREEBAYES_a=""; else FREEBAYES_a="-a "; fi
	FREEBAYES_no_partial_observations=$(grep 'freebayes --no-partial-observations' $CONFIG | awk '{print $1;}'); if [ ${FREEBAYES_no_partial_observations} == "no" ]; then FREEBAYES_no_partial_observations=""; else FREEBAYES_no_partial_observations="--no-partial-observations "; fi
	
	

	echo ""; echo "	Settings read in from config file:"
	echo "		POOLS=$POOLS"
	echo "		POOL_PLOIDY_FILE=$POOL_PLOIDY_FILE"
	echo "		PLOIDY=$PLOIDY"
	echo "		BEST_N_ALLELES=$BEST_N_ALLELES"
	echo "		MIN_MAPPING_QUAL $MIN_MAPPING_QUAL"
	echo "		MIN_BASE_QUAL $MIN_BASE_QUAL"
	echo "		HAPLOTYPE_LENGTH $HAPLOTYPE_LENGTH"
	echo "		MIN_REPEAT_ENTROPY $MIN_REPEAT_ENTROPY"
	echo "		MIN_COVERAGE $MIN_COVERAGE"
	echo "		MIN_ALT_FRACTION $MIN_ALT_FRACTION"
#	echo "		MIN_ALT_COUNT $MIN_ALT_COUNT"
#	echo "		MIN_ALT_TOTAL $MIN_ALT_TOTAL"
#	echo "		READ_MAX_MISMATCH_FRACTION $READ_MAX_MISMATCH_FRACTION"

	echo "		freebayes -z $FREEBAYES_z"
	echo "		freebayes -C $FREEBAYES_C"
	echo "		freebayes -G $FREEBAYES_G"
	echo "		freebayes -3 $FREEBAYES_3"
	echo "		freebayes -Q $FREEBAYES_Q"
	echo "		freebayes -U $FREEBAYES_U"
	echo "		freebayes -\$ $FREEBAYES_DOLLAR"
	echo "		freebayes -e $FREEBAYES_e"

	echo "		freebayes -w $FREEBAYES_w"
	echo "		freebayes -V $FREEBAYES_V"
	echo "		freebayes -a $FREEBAYES_a"
	if [ $FREEBAYES_w != "" ]; then echo "		freebayes ${FREEBAYES_w}"; fi
	if [ $FREEBAYES_V != "" ]; then echo "		freebayes ${FREEBAYES_V}"; fi
	if [ $FREEBAYES_a != "" ]; then echo "		freebayes ${FREEBAYES_a}"; fi
	if [ $FREEBAYES_no_partial_observations != "" ]; then echo "		freebayes ${FREEBAYES_no_partial_observations}"; fi
	
	echo ""; echo " "`date` " Preparing files for genotyping..."
	
	#evaluate input bam files and process them appropriately
	RGBAM=$(ls *$CUTOFFS-RG.bam | wc -l)
	if [ ! -f split.1.$CUTOFFS.bam ]; then
		if [ $RGBAM == 0 ]; then
			RAWBAM=$(ls *${CUTOFFS}*.bam | wc -l)
			if [ $RAWBAM == 0 ]; then
						echo "";echo "  "`date` " ERROR: Genotyping terminated.  No BAM files found."
						exit
			else
				echo "";echo "  "`date` "BAM files detected without -RG.bam extension.  Making -RG.bam files..."
				ls *${CUTOFFS}*.bam | sed 's/\.bam//g' parallel --no-notice -j $NUMProc "cp {}.bam {}-RG.bam"
			fi
		else
			echo "";echo "  "`date` "Filtered BAM files detected with -RG.bam."
		fi
		
		#determine if individual bam files need to be assembled in to a large bam file
		if [ ! -f "cat.$CUTOFFS-RRG.bam" ]; then
			echo "";echo "  "`date` "-RRG.bam file not detected, the individual bam files will be merged"
			echo "  "`date` " samtools merge"
			ls *.$CUTOFFS-RG.bam > bamlist.$CUTOFFS.list
			samtools merge -@$NUMProc -b bamlist.$CUTOFF.$CUTOFF2.list -f cat.$CUTOFF.$CUTOFF2-RRG.bam &>/dev/null
			echo "  "`date` " samtools index"
			samtools index cat.$CUTOFF.$CUTOFF2-RRG.bam 
			wait
		fi
		if [ "cat.${CUTOFFS}-RRG.bam" -nt mapped.$CUTOFFS.bed ]; then
			if [ ! -f cat.$CUTOFF.$CUTOFF2-RRG.bam.bai ]; then samtools index cat.$CUTOFF.$CUTOFF2-RRG.bam; fi
			bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam | bedtools merge > mapped.$CUTOFF.$CUTOFF2.bed
		fi
		#This code estimates the coverage of reference intervals and removes intervals in 0.01% of depth
		#This allows genotyping to be more effecient and eliminates extreme copy number loci from the data
		echo "  "`date` " Estimating coverage of ref intervals & remove extreme copy number loci..."
		if [ "cat.${CUTOFFS}-RRG.bam" -nt "cov.$CUTOFFS.stats" ]; then
			if [ "$BEDTOOLSFLAG" == "OLD" ]; then
				coverageBed -abam cat.$CUTOFFS-RRG.bam -b mapped.$CUTOFFS.bed -counts > cov.$CUTOFFS.stats
			else
				bedtools coverage -b cat.$CUTOFFS-RRG.bam -a mapped.$CUTOFFS.bed -counts -sorted > cov.$CUTOFFS.stats
			fi			
		fi
		if head -1 reference.$CUTOFFS.fasta | grep -e 'dDocent' reference.$CUTOFFS.fasta 1>/dev/null; then
			DP=$(mawk '{print $4}' cov.$CUTOFFS.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		else
			DP=$(mawk '{print $4}' cov.$CUTOFFS.stats | sort -rn | perl -e '$d=.00005;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		fi

		echo "  "`date` " Making the bed files..."
		mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | sort -V -k1,1 -k2,2 | mawk -v x1="$CUTOFFS" -v cutoff=$CC 'BEGIN{i=1} 
		{
			len=$3-$2;lc=len*$4;cov = cov + lc
			if ( cov < cutoff) {x="mapped."i"."x1".bed";print $1"\t"$2"\t"$3 > x}
			else {i=i+1; x="mapped."i"."x1".bed"; print $1"\t"$2"\t"$3 > x; cov=0}
		}' 


			
		split_bam(){
			if [ ! -s split.$1.bam ]; then samtools view -@ 1 -b -1 -L mapped.$1.bed -o split.$1.bam cat.$2-RRG.bam; fi
			if [ ! -s split.$1.bam.bai ]; then samtools index split.$1.bam; fi
		}
		export -f split_bam	
		echo "  "`date` " Splitting BAM File"
		ls mapped.*.$CUTOFFS.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --no-notice -j $NUMProc "split_bam {} $CUTOFFS"
	else
		echo "";echo "  "`date` "Split BAM files detected and they will be genotyped."
		#This code estimates the coverage of reference intervals and removes intervals in 0.01% of depth
		#This allows genotyping to be more effecient and eliminates extreme copy number loci from the data
		echo "  "`date` " Estimating coverage of ref intervals & remove extreme copy number loci..."
		if [ "cat.${CUTOFFS}-RRG.bam" -nt "cov.$CUTOFFS.stats" ]; then
			if [ "$BEDTOOLSFLAG" == "OLD" ]; then
				coverageBed -abam cat.$CUTOFFS-RRG.bam -b mapped.$CUTOFFS.bed -counts > cov.$CUTOFFS.stats
			else
				bedtools coverage -b cat.$CUTOFFS-RRG.bam -a mapped.$CUTOFFS.bed -counts -sorted > cov.$CUTOFFS.stats
			fi			
		fi
		if head -1 reference.$CUTOFFS.fasta | grep -e 'dDocent' reference.$CUTOFFS.fasta 1>/dev/null; then
			DP=$(mawk '{print $4}' cov.$CUTOFFS.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		else
			DP=$(mawk '{print $4}' cov.$CUTOFFS.stats | sort -rn | perl -e '$d=.00005;@l=<>;print $l[int($d*@l)]')
			CC=$( mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
		fi

		echo "  "`date` " Making the bed files..."
		mawk -v x=$DP '$4 < x' cov.$CUTOFFS.stats | sort -V -k1,1 -k2,2 | mawk -v x1="$CUTOFFS" -v cutoff=$CC 'BEGIN{i=1} 
		{
			len=$3-$2;lc=len*$4;cov = cov + lc
			if ( cov < cutoff) {x="mapped."i"."x1".bed";print $1"\t"$2"\t"$3 > x}
			else {i=i+1; x="mapped."i"."x1".bed"; print $1"\t"$2"\t"$3 > x; cov=0}
		}' 
	fi

	if [ ! -s popmap.$CUTOFFS ]; then
		echo "  "`date` " Creating popmap..."
		cut -f1 -d "_" namelist.$CUTOFFS > p.$CUTOFFS
		paste namelist.$CUTOFFS p.$CUTOFFS > popmap.$CUTOFFS
		rm p.$CUTOFFS
		cat popmap.$CUTOFFS
	fi

	
	if [ "$POOLS" == "no" ]; then
		echo; echo " "`date` " Genotyping individuals of ploidy $PLOIDY using freebayes..."			
		#ls mapped.*.$CUTOFFS.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --no-notice -j $NUMProc "freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFFS.fasta -V -p $PLOIDY -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $FREEBAYES_C -G $FREEBAYES_G -3 $FREEBAYES_3 -e FREEBAYES_e -z $FREEBAYES_z -Q $FREEBAYES_Q -U $FREEBAYES_U --populations popmap.$CUTOFFS "
		ls mapped.*.$CUTOFFS.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --no-notice -j $NUMProc "freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFFS.fasta -p $PLOIDY -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $FREEBAYES_C -G $FREEBAYES_G -3 $FREEBAYES_3 -e $FREEBAYES_e -z $FREEBAYES_z -Q $FREEBAYES_Q -U $FREEBAYES_U -$ $FREEBAYES_DOLLAR --populations popmap.$CUTOFFS ${FREEBAYES_w}${FREEBAYES_V}${FREEBAYES_a}${FREEBAYES_no_partial_observations}"
	elif [ "$POOL_PLOIDY_FILE" == "no" ]; then
		echo; echo " "`date` "Running freebayes on pools of cumulative ploidy ${PLOIDY}..."
		ls mapped.*.$CUTOFFS.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --no-notice "freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFFS.fasta -J -p $PLOIDY -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $FREEBAYES_C -G $FREEBAYES_G -3 $FREEBAYES_3 -e $FREEBAYES_e -z $FREEBAYES_z -Q $FREEBAYES_Q -U $FREEBAYES_U -$ $FREEBAYES_DOLLAR --populations popmap.$CUTOFFS ${FREEBAYES_w}${FREEBAYES_V}${FREEBAYES_a}${FREEBAYES_no_partial_observations}"
	elif [ "$POOL_PLOIDY_FILE" != "no" ]; then
		echo; echo " "`date` "Running freebayes on pools with the following cnv file: ${POOL_PLOIDY_FILE}..."
		ls mapped.*.$CUTOFFS.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --no-notice "freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.$CUTOFFS.fasta -J -n $BEST_N_ALLELES -m $MIN_MAPPING_QUAL -q $MIN_BASE_QUAL -E $HAPLOTYPE_LENGTH --min-repeat-entropy $MIN_REPEAT_ENTROPY --min-coverage $MIN_COVERAGE -F $MIN_ALT_FRACTION -C $FREEBAYES_C -G $FREEBAYES_G -3 $FREEBAYES_3 -e $FREEBAYES_e -z $FREEBAYES_z -Q $FREEBAYES_Q -U $FREEBAYES_U -$ $FREEBAYES_DOLLAR --populations popmap.$CUTOFFS --cnv-map $POOL_PLOIDY_FILE ${FREEBAYES_w}${FREEBAYES_V}${FREEBAYES_a}${FREEBAYES_no_partial_observations}" 
	fi

	echo ""; echo "  "`date` "Cleaning up files..."
	#rm split.*.$CUTOFFS.bam*
	#rm mapped.*.$CUTOFFS.bed 

	mv raw.1.$CUTOFFS.vcf raw.01.$CUTOFFS.vcf
	mv raw.2.$CUTOFFS.vcf raw.02.$CUTOFFS.vcf
	mv raw.3.$CUTOFFS.vcf raw.03.$CUTOFFS.vcf
	mv raw.4.$CUTOFFS.vcf raw.04.$CUTOFFS.vcf
	mv raw.5.$CUTOFFS.vcf raw.05.$CUTOFFS.vcf
	mv raw.6.$CUTOFFS.vcf raw.06.$CUTOFFS.vcf
	mv raw.7.$CUTOFFS.vcf raw.07.$CUTOFFS.vcf
	mv raw.8.$CUTOFFS.vcf raw.08.$CUTOFFS.vcf
	mv raw.9.$CUTOFFS.vcf raw.09.$CUTOFFS.vcf

	echo ""; echo " "`date` "Assembling final VCF file..."
	vcfcombine raw.*.$CUTOFFS.vcf | sed -e 's/	\.\:/	\.\/\.\:/g' > TotalRawSNPs.$CUTOFFS.vcf
	bgzip -@ $NUMProc -c TotalRawSNPs.$CUTOFFS.vcf > TotalRawSNPs.$CUTOFFS.vcf.gz
	tabix -p vcf TotalRawSNPs.$CUTOFFS.vcf.gz

	if [ ! -d "raw.$CUTOFFS.vcf" ]; then
		mkdir raw.$CUTOFFS.vcf
	fi

	mv raw.*.$CUTOFFS.vcf ./raw.$CUTOFFS.vcf

}

###############################################################################################
##Create alignment intervals
##This takes advantage of the fact that RAD loci are very discrete.  Instead of calculating intervals for every BAM file,
##this function merges all BAM files together and removes duplicates.  This overall BAM file 
##is used to create a single list of intervals, saving a large amount of computational time.



CreateIntervals()
{


echo "  "`date` " samtools merge"
samtools merge -@$NUMProc -b bamlist.$CUTOFF.$CUTOFF2.list -f cat.$CUTOFF.$CUTOFF2-RRG.bam &>/dev/null
echo "  "`date` " samtools index"
samtools index cat.$CUTOFF.$CUTOFF2-RRG.bam 
wait

echo "  "`date` " bamToBed"
# bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam > map.$CUTOFF.$CUTOFF2.bed
# echo ""
# echo `date` " bedtools merge"
# bedtools merge -i map.$CUTOFF.$CUTOFF2.bed > mapped.$CUTOFF.$CUTOFF2.bed
# rm map.$CUTOFF.$CUTOFF2.bed

bamToBed -i cat.$CUTOFF.$CUTOFF2-RRG.bam | bedtools merge > mapped.$CUTOFF.$CUTOFF2.bed

}
###############################################################################################


###############################################################################################

# GetInfo(){
	# echo "$NumInd individuals are detected. Is this correct? Enter yes or no and press [ENTER]"

	# read Indcorrect

	# if [ "$Indcorrect" == "no" ]; then
			# echo "Please double check that all fastq files are named Ind01.F.fq.gz and Ind01.R.fq.gz"
			# exit 1
	# elif [ "$Indcorrect" == "yes" ]; then
				# echo "Proceeding with $NumInd individuals"
	# else
			# echo "Incorrect Input"
			# exit 1
	# fi

	# #Tries to get number of processors, if not asks user
	# NUMProc=( `grep -c ^processor /proc/cpuinfo 2> /dev/null` ) 
	# NUMProc=$(($NUMProc + 0)) 

	# echo "dDocent detects $NUMProc processors available on this system."
	# echo "Please enter the maximum number of processors to use for this analysis."
			# read NUMProc
			
	# if [ $NUMProc -lt 1 ]; then
			# echo "Incorrect. Please enter the number of processing cores on this computer"
			# read NUMProc
	# fi                
	# if [ $NUMProc -lt 1 ]; then
			# echo "Incorrect input, exiting"
			# exit 1
	# fi

	# #Tries to get maximum system memory, if not asks user
	# MAXMemory=$(($(grep -Po '(?<=^MemTotal:)\s*[0-9]+' /proc/meminfo | tr -d " ") / 1048576))G

	# echo "dDocent detects $MAXMemory maximum memory available on this system."
	# echo "Please enter the maximum memory to use for this analysis. The size can be postfixed with 
	# K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 
	# 1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
	# echo "For example, to limit dDocent to ten gigabytes, enter 10G or 10g"
			# read MAXMemory

	# while [[ -z $MAXMemory ]];
		# do
		# echo "Incorrect input"
		# echo -e "Please enter the maximum memory to use for this analysis. The size can be postfixed with K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
		# echo -e "This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0"
		# echo -e "Then press [ENTER]"
		# read MAXMemory
		# done

	# #Asks if user wants to trim reads.  This allows this part of the pipeline to be skipped during subsequent analyses
	# echo -e "\nDo you want to quality trim your reads?" 
	# echo "Type yes or no and press [ENTER]?"

	# read TRIM

	# #Asks if user wants to perform an assembly.  This allows this part of the pipeline to be skipped during subsequent analyses

	# echo -e "\nDo you want to perform an assembly?"
	# echo "Type yes or no and press [ENTER]?"

	# read ASSEMBLY

	# if [ "$ASSEMBLY" == "no" ]; then
			# echo -e "\nReference contigs need to be in a file named reference.$CUTOFF.$CUTOFF2.fasta\n"
			# sleep 1
	# else
		# echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
		# echo -e "Then press [ENTER]"
		# read ATYPE

		# while [[ $ATYPE != "SE" && $ATYPE != "PE" && $ATYPE != "OL" && $ATYPE != "RPE" ]];
		# do
		# echo "Incorrect input"
		# echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
		# echo -e "Then press [ENTER]"
		# read ATYPE
		# done
	# fi
	# #If performing de novo assembly, asks if the user wants to enter a different -c value
	# if [ "$ASSEMBLY" == "yes" ]; then
		# echo "Reads will be assembled with Rainbow"
		# echo "CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa."
		# echo "Would you like to enter a new c parameter now? Type yes or no and press [ENTER]"
		# read optC
		# if [ "$optC" == "no" ]; then
				# echo "Proceeding with default 0.9 value."
				# simC=0.9
			# elif [ "$optC" == "yes" ]; then
				# echo "Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)"
				# read newC
				# simC=$newC
			# else
				# echo "Incorrect input. Proceeding with the default value."
				# simC=0.9
			# fi
	# fi

	# #Asks if user wants to map reads and change default mapping variables for BWA
	# echo "Do you want to map reads?  Type yes or no and press [ENTER]"
	# read MAP
	# if [ "$MAP" == "no" ]; then
		# echo "Mapping will not be performed"
		# optA=1
		# optB=4
		# optO=6
	# else
		# echo "BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa."
		# echo "Would you like to enter a new parameters now? Type yes or no and press [ENTER]"
		# read optq

		# if [ "$optq" == "yes" ]; then
			# echo "Please enter new value for A (match score).  It should be an integer.  Default is 1."
			# read newA
			# optA=$newA
					# echo "Please enter new value for B (mismatch score).  It should be an integer.  Default is 4."
			# read newB
			# optB=$newB
					# echo "Please enter new value for O (gap penalty).  It should be an integer.  Default is 6."
			# read newO
			# optO=$newO
		# else
			# echo "Proceeding with default values for BWA read mapping."
			# optA=1
			# optB=4
			# optO=6
		# fi
	# fi

	# #Does user wish to call SNPs?
	# echo "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
	# read SNP

	# while [[ $SNP != "yes" && $SNP != "no" ]];
		# do
		# echo "Incorrect input"
		# echo -e "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
		# read SNP
	# done

	# #Asks user for email address to notify when analysis is complete
	# echo ""
	# echo "Please enter your email address.  dDocent will email you when it is finished running."
	# echo "Don't worry; dDocent has no financial need to sell your email address to spammers."
	# read MAIL
	# echo ""
	# echo ""

	# if [ "$ASSEMBLY" == "no" ]; then
		# #Prints instructions on how to move analysis to background and disown process
		# echo "At this point, all configuration information has been entered and dDocent may take several hours to run." 
		# echo "It is recommended that you move this script to a background operation and disable terminal input and output."
		# echo "All data and logfiles will still be recorded."
		# echo "To do this:"
		# echo "Press control and Z simultaneously"
		# echo "Type 'bg' without the quotes and press enter"
		# echo "Type 'disown -h' again without the quotes and press enter"
		# echo ""
		# echo "Now sit back, relax, and wait for your analysis to finish."
	# fi

	# if [ "$ASSEMBLY" == "yes" ]; then
		# echo "dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background."
	# fi
# }


#Actually starts program
if [ -n "$1" ]; then
	#main $1
	# RMH added ALL variables below; main $NUMProc $MAXMemory $TRIM $FixStacks $ASSEMBLY $ATYPE $simC $MAP $optA $optB $optO $SNP $MAIL $HPC $MANCUTOFF $CUTOFF $CUTOFF2 $TRIM_LENGTH_ASSEMBLY $TRIM_LENGTH_MAPPING $SEED_ASSEMBLY $PALIMDROME_ASSEMBLY $SIMPLE_ASSEMBLY $windowSize_ASSEMBLY $windowQuality_ASSEMBLY
	main $NUMProc $MAXMemory $TRIM $TRIM_LENGTH_ASSEMBLY $SEED_ASSEMBLY $PALIMDROME_ASSEMBLY $SIMPLE_ASSEMBLY $windowSize_ASSEMBLY $windowQuality_ASSEMBLY $TRAILING_ASSEMBLY $TRIM_LENGTH_MAPPING $LEADING_MAPPING $TRAILING_MAPPING $FixStacks $ASSEMBLY $ATYPE $simC $HPC $MANCUTOFF $CUTOFF $CUTOFF2 $MAP $optA $optB $optO $MAPPING_MIN_ALIGNMENT_SCORE $MAPPING_CLIPPING_PENALTY $MAPPING_MIN_QUALITY $SAMTOOLS_VIEW_F $SNP $POOLS $POOL_PLOIDY_FILE $PLOIDY $BEST_N_ALLELES $MIN_MAPPING_QUAL $MIN_BASE_QUAL $HAPLOTYPE_LENGTH $MIN_REPEAT_ENTROPY $MIN_COVERAGE $MIN_ALT_FRACTION $MIN_ALT_COUNT $MIN_ALT_TOTAL $READ_MAX_MISMATCH_FRACTION $MAIL $FILTERMAP $SAMTOOLS_VIEW_f $SAMTOOLS_VIEW_Fcustom $SAMTOOLS_VIEW_fcustom $SOFT_CLIP_CUTOFF $FILTER_ORPHANS $BEDTOOLSFLAG $FILTER_MIN_AS $FREEBAYES_Q $FREEBAYES_U $FUNKTION $HEADCROP
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

