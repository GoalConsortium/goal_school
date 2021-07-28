#!/bin/bash
#check_inputfiles.sh

usage() {
  echo "Required Columns for RNASeq Workflow:"
  echo "SampleID"
  echo "FqR1"
  echo "FqR2: if paired-end sequencing"
  echo "-d  --design file"
  echo "-e  --sequence type (pe or se)"
  echo "Example: bash checkdesignfile.sh -x pe -d designfile.txt"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :d:e: opt
do
    case $opt in
        d) design=$OPTARG;;
        e) rpair=$OPTARG;;
    	h) usage;;
    esac
done

shift $(($OPTIND -1))

baseDir="`dirname \"$0\"`"

if [[ -z $rpair ]]
then
    rpair=$1
fi
if [[ -z $design ]]
then
    design=$2
fi
IFS=$'\t'

perl -pe 's/\r\n*/\n/g' $design |perl -pe ' s/FullPathTo//g' > design.fix.txt

hasfqr1=`grep FqR1 design.fix.txt`
hassid=`grep SampleID design.fix.txt`
ispe=`grep FqR2 design.fix.txt`
p='a b'

if [[ ${rpair} == 'pe' ]] &&  [[ -z $ispe ]]
then
    usage
fi
if [[ -z $hassid ]] &&  [[ -z $hasfqr1 ]]
then
    usage
fi
header="SampleID\tSampleName\tSubjectID\tSampleGroup\tFqR1"
if [[ ${rpair} == 'pe' ]]
then
    header="${header}\tFqR2"
fi
echo -e $header > design.valid.txt

while read i; do
    line=($i)
    colnames=($(head -n 1 design.fix.txt))
    while [ "${#colnames[@]}" -gt 0 ] ; do
	var="$(echo -e "${line[0]}" | tr -d '[:space:]')"
	if [[ ! $colnames[0] =~ 'Fq' ]]
	then
	    var="$(echo -e "${var}" | tr -d '-')"
	fi
	eval ${colnames[0]}="${var}"
	colnames=( "${colnames[@]:1}" )
	line=( "${line[@]:1}" )
    done
    if [[ $SampleID != "SampleID" ]]
    then
	SampleGroup="$(echo -e "$SampleGroup" | tr -d '_')"
	if [[ ! -f $FqR1 ]]
	then
	    continue
	fi
	if [[ $SampleID =~ ^[0-9] ]]
	then
	    SampleID="S${SampleID}"
	fi
	if [[ $FqR1 =~ 'gz' ]]
	then
	    mv ${FqR1} ${SampleID}.R1.fastq.gz
	else
	    mv ${FqR1} ${SampleID}.R1.fastq
	    pigz -f ${SampleID}.R1.fastq
	fi
	FqR1="${SampleID}.R1.fastq.gz"
	fqf1="$FqR1"
	fqf2="$FqR1"
        if [[ -f ${FqR2} ]]
	then 
	    if [[ $FqR2 =~ 'gz' ]]
	    then
	        mv ${FqR2} ${SampleID}.R2.fastq.gz
	    else
	        mv ${FqR2} ${SampleID}.R2.fastq
	        pigz -f ${SampleID}.R2.fastq
	    fi
	    FqR2="${SampleID}.R2.fastq.gz"
            fqf1="${fqf1}\t$FqR2"
	    fqf2="${fqf2},$FqR2"
        fi
	if [[ -z $SampleGroup ]]
	then
	    SampleGroup=`printf "%s\n" ${p[@]} | shuf | head -1`
	fi
	if [[ -z $SampleName ]]
	then
	    SampleName=$SampleID
	fi
	if [[ -z $SubjectID ]]
	then
	    SubjectID=$SampleID
	fi
	echo -e "${SampleID}\t${SampleName}\t${SubjectID}\t${SampleGroup}\t${fqf1}" >> design.valid.txt
	echo -e "${SampleID},${fqf2}"
    fi 
done <design.fix.txt
