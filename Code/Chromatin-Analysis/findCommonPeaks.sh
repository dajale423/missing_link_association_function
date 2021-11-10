#!/bin/bash
'''
author: Sumaiya Nazeen <sumaiya_nazeen@hms.harvard.edu>
Finds common merged peaks between different chromatin marks per tissue type.
'''
dnase=$1
ac=$2
me1=$3
me3=$4
output=$5

bedtools intersect -wao -a ${dnase} -b ${ac} ${me1} ${me3} | awk '{if($6!=".") print $0}' > ${output}
