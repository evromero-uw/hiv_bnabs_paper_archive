import os
import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
import numpy as np
import pandas as pd

from par_data_class import Pardata

#This code finds the additional variation sites which passed the
#threshold for the single locus variant of MPL, but not for the multiple locus variant
#It makes a supplemental table of these "additional variation" sites
#and also makes a dataframe for incorporating into the summary plot (panel a)
#of figure 4

SEL_THRESH = 0.044

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will use all alleles at each site
MULTI_SEG = True

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

#In this file I'm going to check out which coefficients MPL adjusted for hitchiking 

#First I will get the MPL selection coefficients for all sites

mpl_dir = '../../../../results/3BNC/MPL/%s/analysis/'
data_dir = '../../../../data/clyde_westfall_2024_final/3BNC117/%s/835_%s_NT.translated.fasta'
coeff_file = 'total-selection-%s.csv'
out_dir = '../../../../results/pub_figs_07-2025/supplemental_tables/mpl_single_locus/'

out_dir_fig4 = '../../../../results/pub_figs_07-2025/figure_4/data/'

par_list = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']

contact_regions_3BNC117 = {'Loop_D': (275, 283),
                    #  'V1': (131, 157),
                    #  'V2': (158, 196),
                    #  'V3': (296, 331),
                    #  'V4': (385, 418),
                     'V5': (460, 470),
                     'CD4_Binding_Loop': (364, 374)}




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
    
    # Get the sites wehre the traditional selection coefficient is above the threshold
    # but the adjusted MPL coefficient is below the threshold
    mpl_coeff = mpl_coeff[mpl_coeff['s_SL'] >= SEL_THRESH]
    mpl_coeff = mpl_coeff[mpl_coeff['s_MPL'] < SEL_THRESH]


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
        if aa_coord not in hxb2_nuc_coords_hta:
            print(f"Skipping {curr_hxb2} as it is not in the alignment.")
            continue
        aa_coord_arr = hxb2_nuc_coords_hta[aa_coord]
        seq_arr_coord = seq_arr[:, aa_coord_arr]
        seq_num = np.sum(seq_arr_coord != '-')
        arr_values, arr_counts = np.unique(seq_arr_coord[seq_arr_coord != '-'], return_counts=True)
        arr_counts = np.sort(arr_counts)[::-1]
        if len(arr_counts) < 2 or arr_counts[1] < 10:
            print(f"Skipping {curr_hxb2} as it does not vary in enough sequences.")
            continue

        #check if the site is in any of the contact regions
        in_contact_region = False
        for region, (start, end) in contact_regions_3BNC117.items():
            if start <= aa_coord <= end:
                in_contact_region = True
                break
        if not in_contact_region:
            print(f"Skipping {curr_hxb2} as it is not in a contact region.")
            continue

        mpl_support_df.append([par, aa_coord, 
                                mpl_coeff[mpl_coeff['HXB2_index'] == curr_hxb2]['s_MPL'].values[0],
                                curr_hxb2, mpl_coeff[mpl_coeff['HXB2_index'] == curr_hxb2]['s_SL'].values[0]])
        
mpl_support_df = pd.DataFrame(mpl_support_df, columns=['participant', 'hxb2_coord_AA', 's_MPL', 'HXB2_index', 's_SL'])
mpl_support_df['hxb2_coord_NT'] = mpl_support_df['hxb2_coord_AA'].apply(get_nucleotide_coord)

#Get only the unique AAs
mpl_support_df = mpl_support_df.drop_duplicates(subset=['participant', 'hxb2_coord_AA'])
mpl_support_df.reset_index(inplace=True, drop=True)
print(mpl_support_df)

#Make a version that will feed into the figure 4 summary plot
mpl_support_fig4 = mpl_support_df.copy()
mpl_support_fig4.drop(columns = ['s_MPL', 's_SL', 'HXB2_index', 'hxb2_coord_NT'], inplace=True)
mpl_support_fig4['method'] = 'linked_variation'
mpl_support_fig4.to_csv(out_dir_fig4 + 'additional_variation.csv', index=False)

#Make a stand alone supplemental table
mpl_support_df.drop(columns = ['hxb2_coord_NT', 'HXB2_index'], inplace=True)
#label any positions aligned to hxb2 gaps
string_labels = []
for curr_str in mpl_support_df['hxb2_coord_AA']:
    if curr_str % 1 == 0:
        string_labels.append(str(int(curr_str)))
    else:
        remainder = curr_str % 1
        remainder = int(np.round(remainder * 1000))  # Convert to integer for display
        alphabet_remainder = ALPHABET[remainder]
        string_labels.append(f"{int(curr_str)}{alphabet_remainder}")
mpl_support_df['hxb2_coord_AA'] = string_labels
mpl_support_df.rename(columns={'participant': 'Participant'}, inplace=True)
mpl_support_df.to_csv(out_dir + 'mpl_single_variant_support.csv', index=False)
print(mpl_support_df)