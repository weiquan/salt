#!/bin/bash 
echo ======================Paired end reads aln test================================================

snpaln_folder=$1
snpaln_index_folder=$2
snpaln_index=$3
read_folder=$4
read_prefix=$5
snpaln_sam_folder=$6
wgsim_folder=$7
snpaln_sam=${read_prefix}_PE.sam
read_len=$8

arg=" -d -p -e -l $read_len -c -a 350 -b 650 -r 5 -t 4"


echo ===aln reads by snpaln===
time $snpaln_folder/snpaln $arg   $snpaln_index_folder/$snpaln_index $read_folder/${read_prefix}1.fq $read_folder/${read_prefix}2.fq >$snpaln_sam_folder/$snpaln_sam
#time ../../Debug/snpaln_debug $arg   $snpaln_index_folder/$snpaln_index $read_folder/${read_prefix}1.fq $read_folder/${read_prefix}2.fq >$snpaln_sam_folder/$snpaln_sam
echo ===aln reads by snpaln end===


echo ===stat results by wgsim_eval.pl===
$wgsim_folder/wgsim_eval.pl alneval $snpaln_sam_folder/$snpaln_sam 2>$snpaln_sam_folder/${read_prefix}_unmapped.sam
echo ===stat results by wgsim_eval.pl end===
echo ======================Paired end reads aln testi finish========================================
