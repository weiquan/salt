#!/bin/bash 
sim_read=T #True for sim reads /False for not sim reads
build_index=T
read_len=100

wgsim_folder=../Simulator/wgsim-master
read_folder=../Reads
read_prefix=Read
ref_folder=../Genome
ref=Genome.fa
variants=../Variants/mutations.txt
hapmap=../Variants/hapmap.txt
salt_folder=../../Bin
salt_index_folder=../Index
salt_index=${ref}.salt.index
salt_sam_folder=../SAM


if [ $sim_read == T ];then
  echo ===sim reads===
  $wgsim_folder/wgsim -e 0 -r 0.05 -R 0 -d 500 -s 50 -N 1000000 -1 $read_len -2 $read_len  -h $ref_folder/$ref $read_folder/${read_prefix}1.fq $read_folder/${read_prefix}2.fq >$variants
  echo ===sim reads end===
fi

echo ===convert mutations.txt to hapmap.txt type===
awk '$3 != "-" && $4 != "-" && length($4)==1{if ($3<$4)$4=$3"/"$4; else $4=$4"/"$3; print $1"\t"$2"\t"$4"\t"$3}' $variants >$hapmap
echo ===convert mutations.txt to hapmap.txt end===


echo ===Indexing===
if [ $build_index == T ];then
  $salt_folder/salt-idx -k 19 $ref_folder/$ref $hapmap $salt_index_folder/$salt_index
fi
echo ===Indexing end===

#./run_se_test.sh $salt_folder $salt_index_folder $salt_index $read_folder $read_prefix $salt_sam_folder $wgsim_folder $read_len
./run_pe_test.sh $salt_folder $salt_index_folder $salt_index $read_folder $read_prefix $salt_sam_folder $wgsim_folder $read_len


