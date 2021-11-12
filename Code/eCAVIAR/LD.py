import glob
import pandas as pd
import gzip
import numpy as np
from sympy import *
from sympy import roots, solve_poly_system
import sys
import os

#form the location_names, get the LD file names
def round_down(x):
    x = int(x)
    x -= x % 1000000
    return(x)

def main():    
    phenotype = sys.argv[1]
    tissue_capitalized = sys.argv[2]
    
    #new gwas directory for Noah
    new_gwas_directory = "/net/data/GTEx/gwas/ukbb_conditioned_gwas/"
    #make output of LD files for eCAVIAR    
    zfile_directory = new_gwas_directory + phenotype + '/eCAVIAR_zscore/' + tissue_capitalized + '/'
    
    LD_upper_directory = new_gwas_directory + phenotype + '/eCAVIAR_LD/'
    output_directory = new_gwas_directory + phenotype + '/eCAVIAR_LD/' + tissue_capitalized + '/'
        
    LD_directory = '/net/data/GTEx/ukbb_ld/'
    chromosome_select = str(sys.argv[3])

    print(zfile_directory)
    print("chrom: " + chromosome_select)
    
    #make directory if it doesn't exist
    if not os.path.exists(LD_upper_directory):
        os.mkdir(LD_upper_directory)    
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    
    zfiles = glob.glob(zfile_directory + '*.Z')

    location_list = []
    chr_list = []

    for i in zfiles:
        filename = i.split("/")[-1]
        location = filename.split("_")[2]
        chromosome = filename.split("_")[1]
        location_list.append(location)
        chr_list.append(chromosome)

    data = {'chr': chr_list, 'location': location_list, "zscore_file": zfiles}
    zfile_df = pd.DataFrame(data)
    zfile_df = zfile_df.sort_values(['location'])

    result = [round_down(x) for x in zfile_df['location']]

    LD_filenames = [str(i + 1) + '_' + str(i + 3000001) for i in result]
    zfile_df["LD_filename"] = LD_filenames
    zfile_df['LD_filename']=zfile_df[['chr','LD_filename']].values.tolist()
    zfile_df['LD_filename']=zfile_df['LD_filename'].apply('_'.join)
    zfile_df['LD_filename'] = 'chr' + zfile_df['LD_filename'].astype(str)

    zfile_df = zfile_df[zfile_df["zscore_file"].str.contains('eqtl_')]

    #remove files that are already made
    LDfiles = glob.glob(output_directory + '*.LD')

    filename_list = []
    for i in LDfiles:
        filename = i.split("/")[-1]
        filename = filename[:-3]
        filename_list.append(filename)

    for i in filename_list:
        zfile_df["zscore_file"] = zfile_df["zscore_file"].astype(str)
        zfile_df = zfile_df[zfile_df["zscore_file"].str.contains(i) == False]
        
    zfile_df["chr"] = zfile_df["chr"].astype(str)
    zfile_df = zfile_df[zfile_df["chr"] == chromosome_select]
    
    MHC_regions = ['chr6_28000001_31000001', 'chr6_29000001_32000001', 'chr6_30000001_33000001']
#     npz2_regions = ['chr6_31000001_34000001', 'chr6_32000001_35000001', 'chr11_121000001_124000001', 'chr11_48000001_51000001']
    npz2_regions = []

    zfile_df = zfile_df[zfile_df["LD_filename"].isin(MHC_regions) == False]
    zfile_df = zfile_df[zfile_df["LD_filename"].isin(npz2_regions) == False]

    print(len(zfile_df))

    old_LDfilename = ''

    for index, row in zfile_df.iterrows():

        #load zfile
        zscore_file_name = row['zscore_file']
        z_file = pd.read_csv(zscore_file_name , sep= '\t', names = ["name", "zscore"])

        print(zscore_file_name)

        if len(z_file) == 0:
            output_filename = output_directory + zscore_file_name.split("/")[-1][5:-2] + '.LD'
            np.savetxt(output_filename, [], fmt='%6.3f')
            continue

        #add chromosome and position column to zscore dataframe
        variant_id = z_file["name"].str.split("_", n = 4, expand = True)
        z_file["chromosome"] = variant_id[0]
        z_file["position"] = variant_id[1].astype(int)
        z_file["major_allele"] = variant_id[2]
        z_file["minor_allele"] = variant_id[3]

        new_LDfilename = row["LD_filename"]

        #reset LDfilename if min position is less
        lower_bound = int(row["LD_filename"].split("_")[1])
        if min(z_file['position']) < lower_bound:
            new_LDfilename = 'chr'+ row['chr'] + '_' + str(lower_bound - 1000000) + '_' + str(lower_bound + 2000000)
        
        if new_LDfilename in MHC_regions:
            continue
        elif new_LDfilename in npz2_regions:
            continue

        if old_LDfilename != new_LDfilename:
            npz_filename = LD_directory + new_LDfilename + '.npz'
            gz_filename = LD_directory + new_LDfilename + '.gz'

            #load files
            data = np.load(npz_filename)
            names = pd.read_csv(gz_filename, sep= '\t')

        old_LDfilename = new_LDfilename

        #make a smaller list from using maximum and minimum indexes
        minimum_value = min(names.loc[names['position'] >= min(z_file['position'])]['position'])
        minimum = names.loc[names['position'] == minimum_value].index[0]
        maximum_value = max(names.loc[names['position'] <= max(z_file['position'])]['position'])
        maximum = names.loc[names['position'] == maximum_value].index[0]


        small_data_row = data['row'][(data['row'] >= minimum) & (data['row'] <= maximum) & (data['col'] >= minimum) & (data['col'] <= maximum)]
        small_data_col = data['col'][(data['row'] >= minimum) & (data['row'] <= maximum) & (data['col'] >= minimum) & (data['col'] <= maximum)]
        small_data_data = data['data'][(data['row'] >= minimum) & (data['row'] <= maximum) & (data['col'] >= minimum) & (data['col'] <= maximum)]
        small_data_df = pd.DataFrame({'row': small_data_row[:], 'col': small_data_col[:], 'data': small_data_data[:]})

        #do further cutoff from list of actual indexes
        names_zscore = names[names['position'].isin(z_file["position"].unique())]
        names_zscore_new = names[names['position'].isin(z_file["position"].unique())]

        #check if the minor allele actually match
        for index, row in names_zscore.iterrows():
            if z_file[z_file['position'] == row['position']]["minor_allele"].iloc[0] != row["allele2"]:
                names_zscore_new.drop(index, inplace=True)

        #if there is a position in zfile that doesn't exist for the LD file, delete the position in the zscore files
        if len(names_zscore_new) < len(z_file):
            #update for eqtl
            z_file = z_file.loc[z_file['position'].isin(names_zscore_new["position"].unique()) == True]
            z_file[["name", "zscore"]].to_csv(zscore_file_name,index = False, header = None, sep='\t')

            #do the same for gwas
            gwas_zscore_file_name = zfile_directory + 'gwas_' + zscore_file_name.split("/")[-1][5:-2] + '.Z'
            gwas_z_file = pd.read_csv(gwas_zscore_file_name , sep= '\t', names = ["name", "zscore"])

            variant_id = gwas_z_file["name"].str.split("_", n = 4, expand = True)
            gwas_z_file["position"] = variant_id[1].astype(int)

            gwas_z_file = gwas_z_file.loc[gwas_z_file['position'].isin(names_zscore_new["position"].unique()) == True]
            gwas_z_file[["name", "zscore"]].to_csv(gwas_zscore_file_name,index = False, header = None, sep='\t')

        index_list = list(names_zscore_new.index)
        small_data_df = small_data_df[small_data_df['row'].isin(index_list) & small_data_df['col'].isin(index_list)]

        small_data_df_new = small_data_df.copy()
    
        for i in small_data_df["row"].unique():
            length = len(small_data_df[small_data_df["row"] == i]) + len(small_data_df[small_data_df["col"] == i])
            if length != len(small_data_df["row"].unique()) + 1:
                small_data_df_new = small_data_df_new[(small_data_df_new["row"] != i)&(small_data_df_new["col"] != i)]

        names_zscore_new = names_zscore_new[names_zscore_new.index.isin(list(small_data_df_new["row"].unique()))]
        index_list = list(names_zscore_new.index)
        small_data_df = small_data_df_new.copy()

        #if there is a position in zfile that doesn't exist for the LD file, delete the position in the zscore files
        if len(names_zscore_new) < len(z_file):
            #update for eqtl
            z_file = z_file.loc[z_file['position'].isin(names_zscore_new["position"].unique()) == True]
            z_file[["name", "zscore"]].to_csv(zscore_file_name,index = False, header = None, sep='\t')

            #do the same for gwas
            gwas_zscore_file_name = zfile_directory + 'gwas_' + zscore_file_name.split("/")[-1][5:-2] + '.Z'
            gwas_z_file = pd.read_csv(gwas_zscore_file_name , sep= '\t', names = ["name", "zscore"])

            variant_id = gwas_z_file["name"].str.split("_", n = 4, expand = True)
            gwas_z_file["position"] = variant_id[1].astype(int)

            gwas_z_file = gwas_z_file.loc[gwas_z_file['position'].isin(names_zscore_new["position"].unique()) == True]
            gwas_z_file[["name", "zscore"]].to_csv(gwas_zscore_file_name,index = False, header = None, sep='\t')



        #calculate the number of SNPs by usign triangular number solution
        x = symbols('x', positive=True)
        k = len(small_data_df)
        length = solve(x**2 + x - 2*k, x)
        if sympify(length[0]).is_integer:
            length = int(length[0])
        else:
            print(gz_filename + 'does not give an integer valued triangular number solution')
            continue

        arr = np.zeros((length,length))

        row_num = 0
        col_num = 0

        #make the array
        for index, row in small_data_df.iterrows():
            arr[row_num, col_num] = row['data']
            if row_num == col_num:
                col_num = col_num + 1
                row_num = 0
            else:
                row_num = row_num + 1

        LDmatrix_gwas = np.add(arr, arr.transpose())

        output_filename = output_directory + zscore_file_name.split("/")[-1][5:-2] + '.LD'

        np.savetxt(output_filename, LDmatrix_gwas, fmt='%6.3f')

        print(output_filename)
    
if __name__ == "__main__":
    main()
