#!/bin/bash

reads_dir='/projects/micb405/resources/project_2/2019/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/MedQPlus_MAGs'
out_dir='/home/alo_mb19'
gtdbtk_dir='/projects/micb405/resources/project_2/2019/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/gtdbtk_output'

array=($reads_dir/*)
bacteria_count=0
archaea_count=0

#f = path to the file including file name
for f in "${array[@]}"; 
do
	path_to_file_no_extension=${f:: -3} #remove the last 3 characters
	name=${path_to_file_no_extension##*/} #retain part after last /
	echo "$name"

	
	if grep -q "\<$name\>" $gtdbtk_dir/gtdbtk.ar122.*_pplacer.tsv
	then 
		kingdom="Archaea"
		echo "$kingdom"
		archaea_count=$(($archaea_count + 1))
	else
		kingdom="Bacteria"
		echo "$kingdom"
		bacteria_count=$(($bacteria_count + 1))
	fi


prokka \
	--outdir $out_dir/Prokka_output3/$name/ \
	--prefix $name \
	--kingdom $kingdom \
	--cpus 4 \
	$reads_dir/$name.fa
done

echo "$bacteria_count"
echo "$archaea_count"

