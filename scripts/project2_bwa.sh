#/bin/bash


mkdir /home/alo_mb19/BWA_output
out_dir='/home/alo_mb19/BWA_output'
MetaT_reads_dir='/projects/micb405/project2/SaanichInlet_10m/MetaT_reads'


#make index
bwa index -p $out_dir/BWA_index /home/alo_mb19/SaanichInlet_MAG_ORFs.ffn #change fna->ffn
echo "Done making index"

array=($MetaT_reads_dir/*)


#align reads and create SAM files for every cruise at 10m
for f in "${array[@]}";
do
	#name=${f##*/} #retain part after last /
	bwa mem -t 8 -p $out_dir/BWA_index $f > $out_dir/SI072_SaanichInlet_MAG_ORFs.sam 
done

echo "Done creating SAM files"
