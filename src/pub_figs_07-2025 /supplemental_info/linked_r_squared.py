import sys
sys.path.append('../../../bin/')
sys.path.append('../../../bin/wrappers/')
import numpy as np
import pandas as pd
from scipy.stats import mode
import r2Analysis as rSquared

from par_data_class import Pardata

#In this file, I am going to calculate the R-squared value for the highly
#linked mutations I found in the 3BNC117 dataset.
supp_table_dif = '../../../results/pub_figs_07-2025/supplemental_tables/'
inDir_3BNC = '../../../data/clyde_westfall_2024_final/3BNC117/'
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24', 'D0', 'W12']


#I think I want to get the dataframe of all sites and separate out
#the additional sites
identified_muts_file = supp_table_dif + '3BNC117_all_data.csv'
identified_muts_df = pd.read_csv(identified_muts_file)
par_list_3BNC = identified_muts_df['participant'].unique().tolist()

linked_sites = identified_muts_df[identified_muts_df['method'] == 'additional_sites']
main_sites = identified_muts_df[identified_muts_df['method'] != 'additional_sites']

possible_genotype_combos = [(False, False), (False, True),
                            (True, False), (True, True)]
possible_genotype_combos = [('W4', 'W4'), ('W4', 'D0'),
                            ('D0', 'W4'), ('D0', 'D0')]
possible_genotype_labels = ['ab', 'aB', 'Ab', 'AB']

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will use all alleles at each site
MULTI_SEG = True
###############################################################################
def get_most_common_char(arr):
    """Get the most common character in an array."""
    unique_chars, counts = np.unique(arr, return_counts=True)
    most_common_char = unique_chars[np.argmax(counts)]
    return most_common_char

def get_most_common_char_non_d0(arr, d0_chars):
    """Get the most common character in an array, excluding day 0 characters."""
    filtered_arr = [char for char in arr if char not in d0_chars]
    if len(filtered_arr) == 0:
        return None  # Return None if no characters are left after filtering
    return get_most_common_char(np.array(filtered_arr))
###############################################################################
#make a place to store the results
all_r2_results = []

#Loop through each participant and calculate the R-squared value for each pair of sites
for curr_par in par_list_3BNC:
    print(curr_par)

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    participant_dat = Pardata(inFile, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

    # Next, I will get a couple of items out of the dataset
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath
    seq_info_df = participant_dat.seq_info_df
    seq_arr = participant_dat.seq_arr

    # Also get the majority allele at day 0 by getting the most common
    # allele at each site
    day_0_df = seq_info_df[seq_info_df['time_label'] == 'D0']
    day_0_seq_arr = seq_arr[day_0_df['seq_index'].values, :]
    day_0_modes = np.apply_along_axis(get_most_common_char, 0, day_0_seq_arr)

    # Get the most common allele at each site at week 4
    week_4_df = seq_info_df[seq_info_df['time_label'] == 'W4']
    week_4_seq_arr = seq_arr[week_4_df['seq_index'].values, :]
    week_4_modes = []
    for i in range(week_4_seq_arr.shape[1]):
        day_0_char = day_0_modes[i]
        week_4_modes.append(get_most_common_char_non_d0(week_4_seq_arr[:, i], day_0_char))
    week_4_modes = np.array(week_4_modes)
        


    # Now I want to get the sequences that circulated from weeks 1-8
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
    seq_arr = seq_arr[seq_info_df['seq_index'].values, :]
    
    # Next I'm going to pull out the alleles at each site
    curr_linked_sites = linked_sites[linked_sites['participant'] == curr_par]
    curr_main_sites = main_sites[main_sites['participant'] == curr_par]

    # Now I will loop through each linked site and calculate the R-squared
    # value for each main site
    
    for linked_site in curr_linked_sites['hxb2_coord_AA'].unique():
        linked_site_index = hxb2_nuc_coords_hta[linked_site]
        linked_site_alleles = seq_arr[:, linked_site_index]
        
        #Label each position as major or minor
        linked_d0 = day_0_modes[linked_site_index]
        linked_w4 = week_4_modes[linked_site_index]
        linked_majority_genotypes = []
        for i in range(len(linked_site_alleles)):
            if linked_site_alleles[i] == linked_d0:
                linked_majority_genotypes.append('D0')
            elif linked_site_alleles[i] == linked_w4:
                linked_majority_genotypes.append('W4')
            else:
                linked_majority_genotypes.append('other')

        #check linkage to the main and additional sites
        possible_links = set(curr_main_sites['hxb2_coord_AA'].unique()).union(set(curr_linked_sites['hxb2_coord_AA'].unique()))
        for main_site in possible_links:
            if main_site == linked_site:
                continue
            main_site_index = hxb2_nuc_coords_hta[main_site]
            main_site_alleles = seq_arr[:, main_site_index]

            main_d0 = day_0_modes[main_site_index]
            main_w4 = week_4_modes[main_site_index]
            main_majority_genotypes = []
            for i in range(len(main_site_alleles)):
                if main_site_alleles[i] == main_d0:
                    main_majority_genotypes.append('D0')
                elif main_site_alleles[i] == main_w4:
                    main_majority_genotypes.append('W4')
                else:
                    main_majority_genotypes.append('other')

            #Join the two arrays together
            combined_genotypes = np.vstack((linked_majority_genotypes, main_majority_genotypes))
            combined_genotypes = combined_genotypes.T
  

            #Count the number of times each combination occurs
            combined_genotypes, combined_counts = np.unique(combined_genotypes, axis=0, return_counts=True)
            count_dict = dict(zip(map(tuple, combined_genotypes), combined_counts))
            
            labeled_count_dict = {}
            for curr_combo in possible_genotype_combos:
                curr_combo_label = possible_genotype_labels[possible_genotype_combos.index(curr_combo)]
                if curr_combo in count_dict:
                    labeled_count_dict[curr_combo_label] = count_dict[curr_combo]
                else:
                    labeled_count_dict[curr_combo_label] = 0

            #Calculate the R-squared value
            d_prime_result = rSquared.calc_D_prime(AB_obs=labeled_count_dict['AB'],
                                                Ab_obs=labeled_count_dict['Ab'],
                                                aB_obs=labeled_count_dict['aB'],
                                                ab_obs=labeled_count_dict['ab'])
            r_squared_result = rSquared.r2(AB_obs=labeled_count_dict['AB'],
                                            Ab_obs=labeled_count_dict['Ab'],
                                            aB_obs=labeled_count_dict['aB'],
                                            ab_obs=labeled_count_dict['ab'])
            all_r2_results.append([curr_par, linked_site, linked_d0, linked_w4, main_site, main_d0, main_w4, d_prime_result, r_squared_result, labeled_count_dict['AB'], labeled_count_dict['Ab'], labeled_count_dict['aB'], labeled_count_dict['ab']])

all_r2_results_df = pd.DataFrame(all_r2_results, columns=['participant', 'linked_site', 'linked_d0', 'linked_w4', 'main_site', 'main_d0', 'main_w4', 'd_prime', 'r_squared', 'AB_count', 'Ab_count', 'aB_count', 'ab_count'])
print(all_r2_results_df)

#Get the max value for each participant and linked site
all_r2_results_df = all_r2_results_df.groupby(['participant', 'linked_site']).agg({
    'd_prime': 'max',
    'r_squared': 'max'
})


print(all_r2_results_df)