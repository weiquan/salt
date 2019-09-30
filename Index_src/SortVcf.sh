#!/bin/bash 
if [ $# -ne 1 ]; then
	echo "$0 <VCF>"
	exit
fi
fn=$1	
cat=zcat #*.gz use zcat, *.bz2 use bzcat, Uncompressed files use cat
chrome_list=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
sorted_list=''
for chr in "${chrome_list[@]}"; do
	$cat $1|grep -w $chr>$fn.$chr.vcf
       	sort -n -k 2 -t '	' $fn.$chr.vcf >$fn.$chr.sorted
	sorted_list=$sorted_list" "$fn.$chr.sorted 
done
cat $sorted_list >$fn.sorted

for chr in "${chrome_list[@]}"; do
	rm $fn.$chr.vcf $fn.$chr.sorted
done
