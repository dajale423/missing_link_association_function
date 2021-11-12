import pandas as pd
import csv
import sys
import os
import glob

def split_eqtl(eqtl): #add new columns to eqtl_copy
    eqtl_copy = eqtl.copy()
    
    eqtl_copy["zscore"] = eqtl_copy['slope']/eqtl_copy['slope_se']

    variant_id = eqtl_copy["variant_id"].str.split("_", n = 4, expand = True)
    eqtl_copy["chromosome"] = variant_id[0]
    eqtl_copy["position"] = variant_id[1].astype(int)
    eqtl_copy["a1"] = variant_id[2].astype(str)
    eqtl_copy["a2"] = variant_id[3].astype(str)
    eqtl_copy = eqtl_copy.sort_values(by=['position'])
    
    return eqtl_copy

def filter_gwas_rows(gwas, eqtl):
    
    gwas_copy = gwas.copy()
    
    for index, row in gwas_copy.iterrows():
        get_eqtl_row = eqtl.loc[eqtl["position"] == row["pos"]]

        a1_eqtl = list(get_eqtl_row["a1"])
        a2_eqtl = list(get_eqtl_row["a2"])

        if (row["a1"] in a1_eqtl) == False:
            gwas_copy.drop(index, inplace=True)
        elif (row["a2"] in a2_eqtl) == False:
            gwas_copy.drop(index, inplace=True)
    return gwas_copy
    

def filter_eqtl_position(gwas, eqtl):
    overlap_gwas_set = set(gwas["pos"])
    overlap_eqtl_set = set(eqtl["position"])
    
    eqtl_copy = eqtl.copy()
    
    eqtl_copy = eqtl_copy.loc[eqtl_copy["position"].isin(list(overlap_eqtl_set - overlap_gwas_set)) == False]
    
    return eqtl_copy

def filter_eqtl_allele(gwas, eqtl):
    eqtl_copy = eqtl.copy()
    
    for index, row in eqtl_copy.iterrows():
        get_gwas_row = gwas.loc[gwas["pos"] == row["position"]]

        a1_gwas = list(get_gwas_row["a1"])
        a2_gwas = list(get_gwas_row["a2"])
        

        if (row["a1"] in a1_gwas) == False:
            eqtl_copy.drop(index, inplace=True)
        elif (row["a2"] in a2_gwas) == False:
            eqtl_copy.drop(index, inplace=True)

    return eqtl_copy

def main():
    #list of file names    
    new_gwas_directory = "/net/data/GTEx/gwas/ukbb_conditioned_gwas/"
    
    phenotype = sys.argv[1]

    tissue_capitalized = sys.argv[2]
    tissue = tissue_capitalized.lower()
    
    if phenotype == "t2d":
        noah_filename = '/net/data/GTEx/gwas/mahajan_2018_t2d/'        
        regions_file = noah_filename + 'cojo/results_hg19/main_pheno_peaks/mahajan_t2d_hg19.indexSNP.tsv'
        gwas_filename = new_gwas_directory + 'mahajan.t2d.hg19.conditioned.cojo'
    elif phenotype == "bc":
        noah_filename = '/net/data/GTEx/gwas/zhang_2020_bc/'
        regions_file = noah_filename + 'cojo/results_hg19/main_pheno_peaks/zhang_bc_hg19.indexSNP.tsv'
        gwas_filename = new_gwas_directory + 'zhang.bc.hg19.conditioned.cojo'        
    elif phenotype == "cd":
        noah_filename = '/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/'
        regions_file = noah_filename + 'cojo/results_hg19/main_pheno_peaks/liu_cd_hg19.indexSNP.tsv'
        gwas_filename = new_gwas_directory + 'liu.cd.hg19.conditioned.cojo'    
    elif phenotype == "uc":
        noah_filename = '/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/'
        regions_file = noah_filename + 'cojo/results_hg19/main_pheno_peaks/liu_uc_hg19.indexSNP.tsv'
        gwas_filename = new_gwas_directory + 'liu.uc.hg19.conditioned.cojo'    
    
    eqtl_file_directory = '/net/data/GTEx/GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/' + tissue + '/'

    output_filename_gwas_directory = new_gwas_directory + phenotype + '/eCAVIAR_zscore/' + tissue_capitalized + '/'
    output_filename_eqtl_begin = output_filename_gwas_directory + 'eqtl_'
    output_filename_gwas_begin = output_filename_gwas_directory + 'gwas_'

    if not os.path.exists(output_filename_gwas_directory):
        os.mkdir(output_filename_gwas_directory)
    else:
        files = glob.glob(output_filename_gwas_directory + '*')
        for f in files:
            os.remove(f)
#             if os.path.getsize(f) < 1 * 1024:
#                 os.remove(f)
    
    
    #read files
    regions_df = pd.read_csv(regions_file, dtype={'CHR': 'int', 'SNP': 'str', 'BP': 'int', 'STARTBP': 'int', \
                                                  'ENDBP': 'int'} ,sep='\t')
    gwas = pd.read_csv(gwas_filename, sep = '\t')

    #make edits to regions_df 
    regions_df["CHROM_NAME"] = "chr" +  regions_df["CHR"].astype(str)

    #make edits to gwas file
    gwas_chr_header = "Chr"
    gwas_pos_header = 'bp'

    snp = gwas["SNP"].str.split(":", n = 4, expand = True)
    gwas["a2"] = snp[2]
    gwas["a1"] = snp[3]

    gwas = gwas[gwas[gwas_chr_header] != "chr"]
    gwas['pos'] = gwas[gwas_pos_header].astype(int)
    gwas['b'] = gwas['bC'].astype(float)
    gwas['se'] = gwas['bC_se'].astype(float)
    gwas["zscore"] = gwas["b"]/gwas["se"]
    gwas["zscore"] = gwas['zscore'].astype(float)
    
    for chromosome in regions_df["CHR"].unique():    
        eqtl_filename = eqtl_file_directory + tissue_capitalized + '.chr' + str(chromosome) +'.txt'
    
        #read eQTL file
        eqtl = pd.read_csv(eqtl_filename, sep = '\t')

        #make edits to eQTL file
        eqtl = split_eqtl(eqtl)

        chr_gwas = gwas[gwas[gwas_chr_header] == chromosome]

        chr_regions_df = regions_df[regions_df["CHR"] == chromosome]
        startbp_list = [x for x in chr_regions_df['STARTBP']]
        endbp_list = [x for x in chr_regions_df['ENDBP']]
        bp_list = [x for x in chr_regions_df['BP']]

        bp_list = zip(startbp_list, endbp_list, bp_list)

        for bp in bp_list:
            region_eqtl = eqtl[(eqtl["position"] >= bp[0]) & (eqtl["position"] <= bp[1])]
            region_gwas = chr_gwas[(chr_gwas["pos"] >= bp[0]) & (chr_gwas["pos"] <= bp[1])]

            overlap = set(region_eqtl.position).intersection(set(region_gwas.pos))
            overlap_eqtl = region_eqtl[region_eqtl['position'].isin(overlap)]

            for gene in overlap_eqtl["gene_id"].unique():
                #skip if a file exists
                output_filename_eqtl = output_filename_eqtl_begin + str(chromosome) + '_' + str(bp[2]) + '_' + gene + '.Z'
                
                if os.path.isfile(output_filename_eqtl):
                    continue
                
                
                gene_eqtl = overlap_eqtl[overlap_eqtl['gene_id'] == gene]

                overlap_gene = set(gene_eqtl.position).intersection(set(region_gwas.pos))
                overlap_gwas = region_gwas[region_gwas["pos"].isin(overlap_gene)]

                #some filtering for minor/major allele
                overlap_gwas = filter_gwas_rows(overlap_gwas, gene_eqtl)

                if len(gene_eqtl) > len(overlap_gwas):
                    gene_eqtl = filter_eqtl_position(overlap_gwas, gene_eqtl)
                    gene_eqtl = filter_eqtl_allele(overlap_gwas, gene_eqtl)

                zscore_gwas = list(overlap_gwas['zscore'])
                zscore_eqtl = list(gene_eqtl['zscore'])
                name_eqtl = list(gene_eqtl['variant_id'])

                #output
                column_names = ["name", "zscore"]

                output_filename_eqtl = output_filename_eqtl_begin + str(chromosome) + '_' + str(bp[2]) + '_' + gene + '.Z'
                output_filename_gwas = output_filename_gwas_begin  + str(chromosome) + '_' + str(bp[2]) + '_' + gene + '.Z'

                print(output_filename_eqtl)

                data = {'name': name_eqtl, 'zscore': zscore_gwas}

                output_df_gwas = pd.DataFrame(data)
                output_df_gwas.to_csv(output_filename_gwas,index = False, header = None, sep='\t')

                output_df_eqtl = pd.DataFrame(columns = column_names)
                output_df_eqtl["name"] = name_eqtl
                output_df_eqtl["zscore"]= zscore_eqtl
                output_df_eqtl.to_csv(output_filename_eqtl,index = False, header = None, sep='\t')
                
                
if __name__ == "__main__":
    main()

