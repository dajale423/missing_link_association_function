#!/usr/bin/env python
'''
author: Sumaiya Nazeen <sumaiya_nazeen@hms.harvard.edu>
Finds candidate narrow peaks per chromatin marks per tissue type within +/- 1 Mb of 
the putatively causative genes.
'''

import os
import sys
import glob

def main():
	if len(sys.argv) < 4:
		print("Usage: python find_candidate_peaks.py macs2.narrowPeaks regions.bed output.narrowPeaks")
		exit(1)
	peaksFile = sys.argv[1]
	targetRegionsFile = sys.argv[2]
	outputFile = sys.argv[3]

	tss = {}
	lines = [l.strip() for l in open(targetRegionsFile)]
	for l in lines:
		x = l.split('\t')
		if x[0] not in tss.keys():
			tss[x[0]] = [int(x[1])]
		else:
			tss[x[0]].append(int(x[1]))

	of = open(outputFile, "w")
	with open(peaksFile) as fp:
		l = fp.readline().strip()
		while(l):
			x = l.split('\t')
			chr = x[0]
			st = int(x[1])
			en = int(x[2])
			if chr in tss.keys():
				chk = tss[chr]
				for v in chk:
					if st >= v - 1000000 and st <= v + 1000000:
						of.write(l+'\n')
						break
					elif en >= v - 1000000 and en <= v + 1000000:
						of.write(l+'\n')
						break
			l = fp.readline().strip()			
	of.close()

if __name__=="__main__":
	main()
