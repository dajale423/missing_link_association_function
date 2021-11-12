import glob
import pandas as pd
import sys


#parameters
coloc_cutoff = float(sys.argv[1])
output_file = 'BC'+ str(coloc_cutoff) + '.txt'

### mendelian list

filename = "/net/home/dlee/noah/v7_eCAVIAR_GTEx/BC_mendelian_genes.txt"

with open(filename) as f:
    lines = [line.rstrip() for line in f]

BC_mendelian = [i.replace(" ", "") for i in lines]
BC_mendelian = [i.split("|")[-2] for i in BC_mendelian]

#make output of LD files for eCAVIAR
ecaviar_directory = '/net/home/dlee/GTEx/gwas/zhang_2020_bc/eCAVIAR_results/eCAVIAR_output/Breast_Mammary_Tissue/'

efiles = glob.glob(ecaviar_directory + '*col')

efiles_nodir = [i.split("/")[-1] for i in efiles]
position = [int(i.split("_")[1]) for i in efiles_nodir]
genes = [i.split("_")[2] for i in efiles_nodir]

data = {'Filename': efiles, 'Locus': position, 'Gene': genes}
efile_df = pd.DataFrame(data)

efile_df = efile_df.sort_values(['Locus'])


output_list = []

print(len(efile_df["Locus"].unique()))

counter = 0

for i in efile_df["Locus"].unique():
    print(counter)
    
    efile_locus_df = efile_df[efile_df["Locus"] == i]
    
    output_string = str(i) + ': '
    
    output_string = output_string + str(len(efile_locus_df)) + ": "
    
    for index, row in efile_df.iterrows():
        eCAVIAR_output = pd.read_csv(row["Filename"] , sep= '\t')
        maxCLPP = max(eCAVIAR_output["CLPP"])
        maxPCausalset = max(eCAVIAR_output["Prob_in_pCausalSet"])
        
        if maxCLPP > coloc_cutoff:
            output_string = output_string + str(row["Gene"]) + ", "
            
    counter += 1
    
    output_list.append(output_string)

with open(output_file, 'w') as f:
    for item in output_list:
        f.write("%s\n" % item)