#!/bin/bash
#check_inputfiles.sh

fqs=`ls *.f*`

for i in $fqs;
do
   if [[ ${i} == *.fq ]];
   then
	new_name=`echo ${i} | sed -e "s/.fq\$/_good.fastq/"`;
	mv ${i} ${new_name};
	`pigz -f ${new_name}`;
   elif [[ ${i} == *.fastq ]];
   then
	new_name=`echo ${i} | sed -e "s/.fastq\$/_good.fastq/"`;
	mv ${i} ${new_name};
	`pigz -f ${new_name}`;
   elif [[ ${i} == *.fq.gz ]];
   then
	new_name=`echo ${i} | sed -e "s/.fq.gz\$/_good.fastq.gz/"`;
	mv ${i} ${new_name};
   else
	new_name=`echo ${i} | sed -e "s/.fastq.gz\$/_good.fastq.gz/"`;
	mv ${i} ${new_name};
   fi;
done
