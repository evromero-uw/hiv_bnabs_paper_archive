import os
import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
import numpy as np
import pandas as pd

from par_data_class import Pardata


SEL_THRESH = 0.044

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will use all alleles at each site
MULTI_SEG = True

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

# In this file I am going to look at all the sites MPL identified as selected. 

#First I will get the MPL selection coefficients for all sites

mpl_dir = '../../../../results/3BNC/MPL/%s/analysis/'
data_dir = '../../../../data/clyde_westfall_2024_final/3BNC117/%s/835_%s_NT.translated.fasta'
coeff_file = 'total-selection-%s.csv'
out_dir = '../../../../results/pub_figs_07-2025/figure_4/MPL_analysis/'

par_list = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']



# First I need to get the identified nucleotide coordinates for each site
def get_nucleotide_coord(aa_coord):
    """Convert amino acid coordinate to get the nucleotide coordinates of
    the first nucleotide of the codon."""
    if aa_coord % 1 != 0:
        frac_part = aa_coord % 1
        aa_coord = int(aa_coord - frac_part - 1) * 3 + 1
        return aa_coord + frac_part
    return (aa_coord - 1) * 3 + 1

mpl_support_df = []

# Now for each participant, I will get the selection coefficients at these
# sites
for par in par_list:

    # Get the selection coefficients for this participant
    curr_mpl_dir = mpl_dir % par
    curr_coeff_file = coeff_file % par
    mpl_coeff = pd.read_csv(curr_mpl_dir + curr_coeff_file)

    # Now get the sites that have high MPL coefficients
    mpl_coeff['s_MPL'] = [np.round(x, 3) for x in mpl_coeff['s_MPL']]
    mpl_coeff = mpl_coeff[mpl_coeff['s_MPL'] >= SEL_THRESH]
    print(mpl_coeff)

    # Next get the aa coordinates for these sites
    for curr_hxb2 in mpl_coeff['HXB2_index']:
        only_nums = ''.join(filter(str.isdigit, curr_hxb2))
        only_letters = ''.join(filter(str.isalpha, curr_hxb2))
        only_letters = only_letters.lower()
        
        if len(only_letters) == 0:
            aa_coord = int(only_nums)
            if aa_coord % 3 != 0:
                aa_coord = np.ceil(aa_coord/3)
            else:
                aa_coord = aa_coord/3

        else: 
            aa_coord = int(only_nums)
            if aa_coord % 3 != 0:
                aa_coord = np.ceil(aa_coord/3)
            else:
                aa_coord = aa_coord/3
            frac_part = ALPHABET.index(only_letters) / 3
            frac_part = np.ceil(frac_part) / 1000
            aa_coord += frac_part
    

        #Now I will check the site in the participant alignment
        #to see how many sequences it is present in
        inFile = data_dir % (par, par)
        participant_dat = Pardata(inFile, 'clyde2024', par)
        participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

        #Check that the site varies in more than SEQ_NUM_THRESH sequences
        seq_arr = participant_dat.seq_arr
        hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
        aa_coord_arr = hxb2_nuc_coords_hta[aa_coord]
        seq_arr_coord = seq_arr[:, aa_coord_arr]
        seq_num = np.sum(seq_arr_coord != '-')
        arr_values, arr_counts = np.unique(seq_arr_coord[seq_arr_coord != '-'], return_counts=True)
        arr_counts = np.sort(arr_counts)[::-1]
        if len(arr_counts) < 2 or arr_counts[1] < 10:
            print(f"Skipping {curr_hxb2} as it does not vary in enough sequences.")
            continue

        # Save the sites in a DataFrame
        mpl_support_df.append([par, aa_coord, 
                                mpl_coeff[mpl_coeff['HXB2_index'] == curr_hxb2]['s_MPL'].values[0],
                                curr_hxb2])

# Combine all the DataFrames into one
mpl_support_df = pd.DataFrame(mpl_support_df, columns=['participant', 'hxb2_coord_AA', 's_MPL', 'HXB2_index'])
mpl_support_df.drop_duplicates(['participant', 'hxb2_coord_AA'], inplace=True)


mpl_support_df.to_csv(out_dir + 'mpl_support_sites.csv', index=False)