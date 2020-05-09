#!/bin/bash 
echo ======================Paired end reads aln test================================================

salt_folder=$1
salt_index_folder=$2
salt_index=$3
read_folder=$4
read_prefix=$5
salt_sam_folder=$6
wgsim_folder=$7
salt_sam=${read_prefix}_PE.sam
read_len=$8

arg=" -d -p -e -l $read_len -c -a 350 -b 650 -r 5 -t 4"


echo ===aln reads by salt===
time $salt_folder/salt $arg   $salt_index_folder/$salt_index $read_folder/${read_prefix}1.fq $read_folder/${read_prefix}2.fq >$salt_sam_folder/$salt_sam
#time ../../Debug/salt_debug $arg   $salt_index_folder/$salt_index $read_folder/${read_prefix}1.fq $read_folder/${read_prefix}2.fq >$salt_sam_folder/$salt_sam
echo ===aln reads by salt end===


echo ===stat results by wgsim_eval.pl===
$wgsim_folder/wgsim_eval.pl alneval $salt_sam_folder/$salt_sam 2>$salt_sam_folder/${read_prefix}_unmapped.sam
echo ===stat results by wgsim_eval.pl end===
echo ======================Paired end reads aln testi finish========================================
