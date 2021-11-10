#!/bin/bash
'''
author: Sumaiya Nazeen <sumaiya_nazeen@hms.harvard.edu>
Finds common merged peaks between different chromatin marks per tissue type.
'''
ac=$1
me1=$2
me3=$3
output=$4

bedtools intersect -wao -a ${ac} -b ${me1} ${me3} | awk '{if($6!=".") print $0}' > ${output}
