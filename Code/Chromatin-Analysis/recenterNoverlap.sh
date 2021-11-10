input=$1
sizefile=$2
blacklist=$3
outprefix=$4

bedtools slop -i ${input} -b 175 -g ${sizefile} | awk 'OFS="\t"{l=$3-$2; if (l < 500){$2=$2-int((500-l)/2); $3=$3 + int((500-l)/2);} print $0}' | bedtools sort -i stdin -faidx ${sizefile} | bedtools merge -i stdin -c 7 -o mean,count | bedtools intersect -v -wa -a stdin -b ${blacklist} | bedtools sort -i stdin -faidx ${sizefile} | bedtools merge -i stdin -c 4 -o mean > ${outprefix}.recentered.bed;
