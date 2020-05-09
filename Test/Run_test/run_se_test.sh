#!/bin/bash 
echo ======================Single end reads aln test================================================
salt_folder=$1
salt_index_folder=$2
salt_index=$3
read_folder=$4
read_prefix=$5
salt_sam_folder=$6
wgsim_folder=$7
salt_sam=${read_prefix}_SE.sam
read_len=$8
arg="-d -r 1 -l $read_len -n 20 -c -m 500 -t 4" 

echo ===aln reads by salt===
time $salt_folder/salt  $arg  $salt_index_folder/$salt_index $read_folder/${read_prefix}1.fq >$salt_sam_folder/$salt_sam
echo ===aln reads by salt end===


echo ===stat results by wgsim_eval.pl===
$wgsim_folder/wgsim_eval.pl alneval $salt_sam_folder/$salt_sam 2>$salt_sam_folder/${read_prefix}_unmapped.sam
echo ===stat results by wgsim_eval.pl end===
echo ======================Single end reads aln testi finish========================================

