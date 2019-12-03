#!/bin/bash

BWA_output_dir='/home/alo_mb19/BWA_output'
home_dir='/home/alo_mb19'
prokka_dir='/home/alo_mb19/High_Qual_MAGs_Prokka_Result_AA'

mkdir /home/alo_mb19/RPKM_outputs

array=($BWA_output_dir/*.sam)

#Note: chnaged fna -> ffn

for f in "${array[@]}"
do
	name=${f##*/} #retain part after last /
	/projects/micb405/resources/project_2/2019/rpkm \
	-c $home_dir/SaanichInlet_MAG_ORFs.ffn \
	-a $f \
	-o $home_dir/RPKM_outputs/"${name}_RPKM.csv" 
done
