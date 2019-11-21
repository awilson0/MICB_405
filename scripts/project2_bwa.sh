#/bin/bash


mkdir /home/alo_mb19/BWA_output
out_dir='/home/alo_mb19/BWA_output'
metatranscriptomes_dir='/projects/micb405/resources/project_2/2019/Metatranscriptomes'


#make index
bwa index -p $out_dir/BWA_index /home/alo_mb19/high_qual_MAGs_all_genes.fna
echo "Done making index"

array=($metatranscriptomes_dir/*_10m.qtrim.artifact.rRNA.clean.fastq*)


#align reads and create SAM files for every cruise at 10m
for f in "${array[@]}";
do
	name=${f##*/} #retain part after last /
	bwa mem -t -2 $out_dir/BWA_index $f > $out_dir/"${name}.sam" 
done

echo "Done creating SAM files"
