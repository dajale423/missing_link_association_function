#!/usr/bin/env Rscript
#' author: Sumaiya Nazeen
#' Performs fine-mapping of GWAS hits using the SuSiE method described in the following paper. 
#' Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). 
#' A simple new approach to variable selection in regression, with application to genetic fine mapping. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology). https://doi.org/10.1111/rssb.12388
#' 
#' @param z_file Path of a file containing GWAS summary statistics. This file has seven space-separated columns with the following headers: 
#'               rsid chromosome position allele1 allele2 maf beta se
#' @param ld_matrix_npz_file LD matrix stored in the Scipy sparse matrix format as a .npz file 
#' @param N Number of study participants in the original GWAS study
#' @param out_file Path to a file which stores the posterior inclusion probabilities of the GWAS variants.
#' output fine-mapped GWAS variants with posterior inclusion probabilities between 0 and 1.

# import required packages
library(reticulate)
library(susieR)
library(Matrix)
library(matrixcalc)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
	stop("Insufficient arguments.
		Required parameters: 
		z_file
		ld_matrix_npz_file
		N
		out_file"
	, call.=FALSE)
} else {
	z_file <- args[1]
	ld_mat_file <- args[2]
	N <- as.integer(args[3])
	outfile <- args[4]
	
	z <- read.table(z_file, sep=' ', header=T)
	z[is.na(z)] <- 0
	betahat <- z$beta
	sehat <- z$se
	scipy <- import("scipy.sparse")
	sparse_mat<- scipy$load_npz(ld_mat_file)
	ld_mat <- as.matrix(sparse_mat)

	ind <- which(is.na(betahat))
	if (length(ind) != 0){
		betahat <- as.numeric(na.omit(betahat))
		sehat <- as.numeric(na.omit(sehat))
		ld_mat <- ld_mat[-ind,-ind]
		z <- na.omit(z)
	}

	#if(is.positive.semi.definite(ld_mat)){
	#	R <- ld_mat
	#} else{
	#	nc <- nearPD(ld_mat, conv.tol=1e-7, corr=T, keepDiag=T)
	#	R <- matrix(nc$mat, nrow=dim(nc$mat)[1])
	#}
	R <- ld_mat + t(ld_mat)
	diag(R) <- 1
	res <- susie_suff_stat(bhat=betahat,shat=sehat, R=R, n=N, estimate_residual_variance=FALSE, estimate_prior_variance = FALSE)
	pip <- res$pip
	converged <- res$converged
	zres <- cbind(z, pip, converged)
	zresn <- zres[order(-zres$pip, zres$position),]
	write.table(zresn, outfile, sep='\t', quote=FALSE, row.names=FALSE)
}
