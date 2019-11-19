#!/bin/bash

KAAS_dir='/home/alo_mb19/KAAS_result'
mkdir /home/alo_mb19/KAAS_result_cleaned
out_dir='/home/alo_mb19/KAAS_result_cleaned'

array=($KAAS_dir/*)

for f in "${array[@]}";
do 
	path_to_file_no_extension=${f:: -4} #remove the .txt 
	name=${path_to_file_no_extension##*/} #retain part after last /
	#echo "${name}.cleaned.txt"
	#echo $KAAS_dir/"${name}.cleaned.txt"
	
	
	grep '\sK' $f > $out_dir/"${name}.cleaned.txt"


# grep '\sK' SaanichInlet_10m_6_ORFs_ko.txt > SaanichInlet_10m_6_ORFs_ko.cleaned.txt

done
