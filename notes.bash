# different lines in ddocent than dDocentHPC

SED1='s/^[ \t]*//'
SED2='s/\s/\t/g'
FRL=$(gunzip -c ${NAMES[0]}.F.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)


special_uniq(){
	mawk -v x=$1 '$1 >= x' $2  |cut -f2 | sed -e 's/NNNNNNNNNN/	/g' | cut -f1 | uniq
}
export -f special_uniq


	