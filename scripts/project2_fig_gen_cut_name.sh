#!/bin/bash

for f in /home/alo_mb19/High_Qual_MAGs_Prokka_Result_AA/*/*faa
do
	prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
	mag_id=$( echo $f | sed 's/.faa//g' )
	mag_id=${mag_id##*/} # retain file path part after last /
	mag_id=${mag_id/./_}
	echo $prokka_id,$mag_id
done > Prokka_MAG_map.csv
