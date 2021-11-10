import sys
import os
import numpy as np
import scipy as sp
from scipy.stats import mstats
import glob

def isOverlapped(st,en,tss,ten):
	if (st < tss and en < tss) or st > ten:
		return "No"
	else:
		return "Yes"

def get_neighbors_acts(w, gmap):
	v = w.split(':')
	chr = v[0]
	y = v[1].split('_')
	st = int(y[0])
	en = int(y[1])
	neigh_acts = []
	feats = [x for x,_ in gmap]
	vals = [y for _,y in gmap]
	for i in range(len(feats)):
		f = feats[i].split(':')
		p = f[1].split('_')
		stf = int(p[0])
		enf = int(p[1])
		if chr != f[0] or stf < st-5000000 or stf > en+5000000:
			 continue
		print(w, feats[i])
		neigh_acts.append(vals[i])
	norm_x = sum(neigh_acts)
	return norm_x

def main():
	if len(sys.argv) < 4:
		print("Usage: python abd-compute.py peaks_file gencode_file outprefix")
		exit(1)
	
	peaks_file = sys.argv[1]
	gencode_file = sys.argv[2]
	outprefix = sys.argv[3]
	outfull = outprefix + '_full.tsv'
	outsmall = outprefix + '_top5.tsv'
	
	lines = [l.strip() for l in open(peaks_file)]
	m = {}
	for l in lines:
		x = l.split('\t')
		k = x[0]+':'+x[1]+'_'+x[2]
		if k not in m.keys():
			m[k] = [[float(x[3])],x[4].split(','),[],[],[],[],0]
		if x[5] == '1':
			m[k][2].append(float(x[9]))
			m[k][3] += x[10].split(',')
		elif x[5] == '2':
			m[k][4].append(float(x[9]))
			m[k][5] += x[10].split(',')

	print("Step1")	
	for k in m.keys():
		scores = m[k][0] + m[k][2] + m[k][4]
		m[k][6] = mstats.gmean(scores)
	
	glines = [l.strip() for l in open(gencode_file)]
	gm = {}
	for l in glines:
		x = l.split('\t')
		if x[2] not in gm.keys():
			gm[x[2]] = [(x[0],int(x[4]),int(x[5]))]
		else:
			gm[x[2]].append((x[0],int(x[4]),int(x[5])))
	print("Step2")		
	m2 = {}
	for k in m.keys():
		v = k.split(':')
		chr = v[0]
		y = v[1].split('_')
		st = int(y[0])
		en = int(y[1])
		for gene, tss, ten in gm[chr]:
			dist = abs(tss - ((st+en)/2.0))+1
			if dist <= 100000:
				if k not in m2.keys():
					m2[k] = m[k] + [[],[],[],[],[],[]]
				m2[k][7].append(gene)
				m2[k][8].append(tss)
				m2[k][9].append(isOverlapped(st,en,tss,ten))
				m2[k][10].append(dist)
				m2[k][11].append(m[k][6]*100000./dist)
	
#	for k in m2.keys():
#		genes = m2[k][7]
#		acts = m2[k][11]
#		m3 = {}
#		for i in range(len(genes)):
#			if genes[i] not in m3.keys():
#				m3[genes[i]] = [acts[i]]
#			else:
#				m3[genes[i]].append(acts[i])
#		for i in range(len(acts)):
#			abc = acts[i]/sum(acts)
#			m2[k][12].append(abc)
	print("Step3")
	print(len(m2.keys()))
	m4 = {}
	for k in m2.keys():
		genes = m2[k][7]
		acts = m2[k][11]
		m2[k][12] = [0.0 for _ in acts]
		m3 = {}
		for i in range(len(genes)):
			if genes[i] not in m3.keys():
				m3[genes[i]] = [acts[i]]
			else:
				m3[genes[i]].append(acts[i])
		gs = list(m3.keys())
		acs = [0.0 for _ in gs]
		for i in range(len(gs)):
			acs[i] = max(m3[gs[i]])
			if gs[i] not in m4.keys():
				m4[gs[i]] = [(k, acs[i])]
			else:
				m4[gs[i]].append((k, acs[i]))
			
	
	for g in m4.keys():
		print(g, m4[g])
		for i in range(len(m4[g])):
			w,x = m4[g][i]
			v = w.split(':')
			chr = v[0]
			y = v[1].split('_')
			st = int(y[0])
			en = int(y[1])
			norm_x = get_neighbors_acts(w,m4[g])
			genes = m2[w][7]
			for j in range(len(genes)):
				if genes[j] == g:
					m2[w][12][j] = m2[w][11][j] / norm_x				

	of = open(outfull,"w")
	of2 = open(outsmall,"w")
	of.write("FeatureId\tChr\tStart\tEnd\tH3K27ac_Enrichment\tH3K27ac_Summits\tH3K4me1_Enrichment\tH3K4me1_Summits\tH3K4me3_Enrichment\tH3K4me3_Summits\tActivity_gmean\tGenes\tGeneTSS\tOverlaps\tDistances\tABD_scores\n")
	of2.write("FeatureId\tChr\tStart\tEnd\tH3K27ac_Enrichment\tH3K27ac_Summits\tH3K4me1_Enrichment\tH3K4me1_Summits\tH3K4me3_Enrichment\tH3K4me3_Summits\tActivity_gmean\tGenes\tGeneTSS\tOverlaps\tDistances\tABD_scores\n")
	for k in m2.keys():
		fid = k
		x = k.split(':')
		chr = x[0]
		y = x[1].split('_')
		st = y[0]
		en = y[1]
		s = fid + '\t' + chr + '\t' + st + '\t' + en + '\t' + ';'.join([str(v) for v in m2[k][0]]) + '\t' + ';'.join(m2[k][1]) + '\t' + ';'.join([str(v) for v in m2[k][2]]) + '\t' + ';'.join(m2[k][3]) + '\t' + ';'.join([str(v) for v in m2[k][4]]) + '\t' + ';'.join(m2[k][5]) + '\t' + str(m2[k][6]) + '\t' 
		s_abd = np.sort(m2[k][12])
		rs_abd = s_abd[::-1]
		s_arg = np.argsort(m2[k][12])
		rs_arg = s_arg[::-1]
		all_genes = [m2[k][7][w] for w in rs_arg]
		all_tss = [m2[k][8][w] for w in rs_arg]
		all_ovss = [m2[k][9][w] for w in rs_arg] 
		all_dists = [m2[k][10][w] for w in rs_arg]
		top5_abd = rs_abd[0:5]
		top5_genes = all_genes[0:5]
		top5_tss = all_tss[0:5]
		top5_ovss = all_ovss[0:5]
		top5_dists = all_dists[0:5]
		s2 = s
		s += ';'.join(all_genes) + '\t'
		s += ';'.join([str(v) for v in all_tss]) + '\t'
		s += ';'.join(all_ovss)
		s += '\t' + ';'.join([str(v) for v in all_dists])
		s += '\t' + ';'.join([str(v) for v in rs_abd]) + '\n'
		of.write(s)
		s2 += ';'.join(top5_genes) + '\t'
		s2 += ';'.join([str(v) for v in top5_tss]) + '\t'
		s2 += ';'.join(top5_ovss)
		s2 += '\t' + ';'.join([str(v) for v in top5_dists])
		s2 += '\t' + ';'.join([str(v) for v in top5_abd]) + '\n'
		of2.write(s2)
	
	of.close()
	of2.close()

if __name__=="__main__":
	main()
