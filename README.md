This code is associated with the following manuscript. If you use any part of the source code, please cite us:

"Connally, Noah, Sumaiya Nazeen, Daniel Lee, Huwenbo Shi, John Stamatoyannopoulos, Sung Chun, Chris Cotsapas, Christopher Cassa, and Shamil Sunyaev. 
The missing link between genetic association and regulatory function. medRxiv (2021). https://www.medrxiv.org/content/10.1101/2021.06.08.21258515v2"

### Requirements
	Python >= 3.5.3
	Python packages:
		numpy 1.18.5
		scipy 1.4.1
	R 3.4.3
	R packages: 
		susieR 0.9.0
		matrixcalc 1.0-3
		Matrix 1.2-12
		reticulate 1.16 

	This code has been run with GCC 6.3.0.

### Directory Structure
	Data: contains auxiliary data files needed to run the scripts.
	Code/
		Fine-mapping: contains the R script for running fine-mapping of gwas variants using SuSiE alogorithm.
		Chromatin-Analysis: contains the scripts for computing ABD scores for the gene-feature pairs per tissue type.

### Instructions	
	i)
	ii) Fine-mapping:
		```
		Rscript susie-finemapping.R <z_file> <ld_matrix_npz_file> <N> <out_file>
		```
	iii) Computing ABD scores:
		a. Finding candidate peaks per chromatin mark per tissue type:
		   ```
		   python3 find_candidate_peaks.py input.narrowPeak causativeGenes.bed step1.candidatePeak
		   ```
		b. Recenter and overlap per chromatin mark per tissue type:
		   ```
		   ./recenterNoverlap.sh step1.candidatePeak chr_sizes blacklist-hg19.bed step2.recenteredPeak
		   ```
		c. Find common chromatin peaks between H3K27AC, H3K4me1, and H3K4me3 peaks per tissue type
		   ```
		   ./findCommonPeaks.sh ac.recenteredPeak me1.recenteredPeak me3.recenteredPeak step3.commonPeak
		   ```
		d. Compute activity by distance per tissue type:
		   ```
		   python3 abd-compute.py step3.commonPeak reFlat.gencode.v19 output
		   ```


### Contact: 
