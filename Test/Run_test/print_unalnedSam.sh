#!/bin/bash 
usage="$0 <SAM>"

if [ $# != 1 ]; then
    echo $usage
    exit 1;
elif [ ! -f $1 ];then
    echo "[Error]: Can\'t find file $1 !"
    exit 1;
fi


SAM=$1
FLAG_UMAPPED=4
awk -v u=$FLAG_UMAPPED '/^[^@]/{if(and($2, u) != 0) {print $0}}' $1
