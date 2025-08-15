import os
import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd

import data_util as du
import dataset_metadata
from par_data_class import Pardata

#In this plot I am going to try and identify major contributors to escape and
#trace the haplotype lineages

#Of the annotated sites, I am looking at the derived alleles present at high
#frequencies > 3*IQR above the 75th percentile of the distribution of derived
#alleles at each site

#In this file I am going to try putting the haplotype frequencies for the
#explained participants into the dataframs I need to make area plots

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

#Then I am going to use these sites to make the lineage haplotype plots
outFolder = 'inclusive_hap_diagrams/'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/'

inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'
par_list_3BNC = ['2E7', '2E4', '2E5', '2C5', '2E2', '2C1', '2D1', '2E1', '2E3']


time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2']

inDir_fig4_summary = '../../../../results/pub_figs_07-2025/supplemental_tables/3BNC117_all_data.csv'
fig4_summary_df = pd.read_csv(inDir_fig4_summary)

HXB2_RES_POS = dataset_metadata.RESISTANCE_POS_AA_HXB2 

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will only use the two most frequent alleles at each site
MULTI_SEG = True

# Day 0 freq threshold
D0_FREQ_THRESH = 0

HAP_FREQ_THRESH = 0
FREQ_DIFF_THRESH = 0.1

# Require the difference in frequency between to be more than one read
num_seqs_required = 10

################################ Helper Functions #############################
###############################################################################
def sort_escape_alleles(escape_alleles):
    """Takes in a list of escape alleles and sorts them by their HXB2
    coordinates and then by their derived allele. Returns a sorted list.
    """
    number_components = [float(x[:-1]) for x in escape_alleles]
    number_dict = {}
    for escape_allele in escape_alleles:
        curr_number = float(escape_allele[:-1])
        curr_letter = escape_allele[-1]
        if curr_number not in number_dict:
            number_dict[curr_number] = [curr_letter]
        else:
            number_dict[curr_number].append(curr_letter)
    
    number_components = sorted(number_components)
    new_number_components = []
    for curr_number in number_components:
        # If the number is a float, convert it to a string with the letter
        # representation of the fractional part
        if curr_number % 1 != 0:
            frac_part = curr_number % 1
            curr_number = int(curr_number)
            frac_part_str = np.round(frac_part * 1000)
            frac_part_str = ALPHABET[int(frac_part_str)]
            curr_number = str(curr_number) + frac_part_str

        else:
            # If the number is an integer, convert it to a string
            curr_number = int(curr_number)
            curr_number = str(curr_number)
        new_number_components.append(curr_number)
    

    new_escape_alleles = []
    for i in range(len(number_components)):
        curr_number = number_components[i]
        curr_letters = number_dict[curr_number]
        curr_letters = sorted(curr_letters)
        for curr_letter in curr_letters:
            new_escape_alleles.append(new_number_components[i] + curr_letter)
    print('New escape alleles:', new_escape_alleles)
    return new_escape_alleles

###############################################################################
# fig4_summary_df = fig4_summary_df[fig4_summary_df['method'] != 'Glycan variability']
# escape_sites_hxb2 = fig4_summary_df['hxb2_coord_AA'].unique()
# escape_sites_hxb2 = [int(x) for x in escape_sites_hxb2]
# escape_sites_hxb2 = set(escape_sites_hxb2)

all_hap_info_df = []

# Next I will loop through the participants
# Next, I need to get all of the hxb2 positions across all individuals
hxb2_pos_unique = set()
for curr_par in par_list_3BNC:
    escape_sites_hxb2 = fig4_summary_df[fig4_summary_df['participant'] == curr_par]['hxb2_coord_AA'].unique()
    # escape_sites_hxb2 = [int(x) for x in escape_sites_hxb2]
    escape_sites_hxb2 = set(escape_sites_hxb2)
    print('Participant:', curr_par)
    print('Escape sites:', escape_sites_hxb2)
    if len(escape_sites_hxb2) == 0:
        print('No escape sites for participant:', curr_par)
        continue

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    # First, I need to load the data
    participant_dat = Pardata(inFile, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

    # Next, I will get a couple of items out of the dataset
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]

    # Get the unique sites in array coordinates
    escape_sites_arr = [hxb2_nuc_coords_hta[x] for x in escape_sites_hxb2]

    if curr_par == '2D1':
        print(escape_sites_hxb2)


    #Get the segregating sites
    all_seg_freq_dict = du.get_seg_sites(seqArr, set(),
                            allele_freq_thresh = ALLELE_FREQ_THRESH,
                            return_multi = MULTI_SEG,
                            include_drm_seg = True,
                            verbose = False)

    #Get the allele frequencies at each timepoint
    allele_freq_df, all_time_freqs, time_sample_sizes = du.seg_freq_by_timepoint(
                                                        all_seg_freq_dict,
                                                        seqArr,
                                                        seq_info_df)
    time_sample_sizes = seq_info_df.groupby('time_label').size().to_dict()

    #look for alleles which increase in frequency
    increasing_sites_dict = {}
    frequency_diff_df = allele_freq_df.copy()
    #only check timepoints where the frequency increase is more than one sequence 
    #basically ignore sparse sampling times
    frequency_diff_df['sample_size'] = frequency_diff_df['time'].apply(lambda x: time_sample_sizes[x])
    frequency_diff_df = frequency_diff_df[frequency_diff_df['sample_size'] > num_seqs_required]

    #Also ignore week 12 since the bnab is washing out at that point
    frequency_diff_df = frequency_diff_df[frequency_diff_df['time'] != 'W12']
    
    print('Participant:', curr_par)
    print(frequency_diff_df['time'].unique())


    for name, group in frequency_diff_df.groupby(['position', 'allele']):
        if group.shape[0] > 1:
            day_0_freq = group[group['time'] == 'D0']['freqs'].values[0]
            other_freqs = group[group['time'] != 'D0']['freqs'].values
            freq_diff = other_freqs - day_0_freq
            max_freq_diff = np.max(freq_diff)
            
            
            if max_freq_diff > FREQ_DIFF_THRESH:
                if curr_par == '2E3':
                    print(name, max_freq_diff)

                if name[0] not in increasing_sites_dict:
                    increasing_sites_dict[name[0]] = [name[1]]
                else:
                    increasing_sites_dict[name[0]].append(name[1])



    #Identify alleles that are common at day 0
    common_allele_dict = {}
    for curr_site in escape_sites_arr:
        d0_freqs = allele_freq_df.loc[allele_freq_df['time'] == 'D0']
        d0_freqs = d0_freqs[d0_freqs['position'] == curr_site]
        d0_freqs = d0_freqs[d0_freqs['freqs'] > D0_FREQ_THRESH]

        #If the site is not segregating in our dataset, then we will skip it
        if d0_freqs.shape[0] == 0:
            continue
        else:
            #Get the most common allele
            d0_freqs = d0_freqs.sort_values(by = 'freqs', ascending = False)
            common_allele_dict[curr_site] = d0_freqs['allele'].values[0]

    escape_sites_arr = set(common_allele_dict.keys())
    increasing_sites_arr = set(increasing_sites_dict.keys())

    escape_sites_arr = escape_sites_arr.intersection(increasing_sites_arr)
    escape_sites_arr = list(escape_sites_arr)
    
    #Now, I will get the haplotypes at the sites of interest
    escape_hap_arr = seqArr[:, escape_sites_arr]

    
    #Next, I will mark out all of the alleles that are common at day 0
    #or do not increase in frequency
    new_escape_hap_arr = []
    for curr_column in range(escape_hap_arr.shape[1]):
        curr_column_site = escape_sites_arr[curr_column]
        curr_common_allele = common_allele_dict[curr_column_site]


        curr_increasing_allele = increasing_sites_dict[curr_column_site]
        curr_increasing_allele = set(curr_increasing_allele)


        allele_values = escape_hap_arr[:, curr_column]

        new_col = []
        for curr_entry in allele_values:
            if curr_entry in curr_common_allele:
                new_col.append('Common')
            elif curr_entry in curr_increasing_allele:
                new_col.append(curr_entry)
            else:
                new_col.append('-')

        new_escape_hap_arr.append(new_col)
    new_escape_hap_arr = np.array(new_escape_hap_arr).T
    escape_hap_arr = new_escape_hap_arr

    if escape_hap_arr.shape[0] == 0:
        continue
    
    
    #Look at how many rare alleles appear after day 0 for each site
    later_timepoints = seq_info_df[seq_info_df['time_label'] != 'D0']
    later_timepoints = later_timepoints['seq_index'].values
    escape_hap_arr_later = escape_hap_arr[later_timepoints, :]
  
    escape_hap_arr_later = ~np.isin(escape_hap_arr_later, ['Common', '-'])
    num_later_seqs = escape_hap_arr_later.shape[0]
    escape_hap_arr_later = np.sum(escape_hap_arr_later, axis=0)
    escape_hap_arr_later = escape_hap_arr_later / num_later_seqs

   
    #First I will find the array coordinates of the sites to plot
    escape_sites_to_plot_arr = np.argwhere(escape_hap_arr_later > 0)
    escape_sites_to_plot_arr = [escape_sites_arr[x[0]] for x in escape_sites_to_plot_arr]
    #Next I will get their indices in my mini hap array
    escape_sites_to_plot_inds = np.argwhere(np.isin(escape_sites_arr, escape_sites_to_plot_arr))
    escape_sites_to_plot_inds = [x[0] for x in escape_sites_to_plot_inds]

    #Now I will get the HXB2 coordinates of the sites to plot
    escape_to_plot_hxb2 = [hxb2_nuc_coords_ath[x] for x in escape_sites_to_plot_arr]
 
    #Now I will subset the hap array down to the sites to plot
    escape_hap_arr = escape_hap_arr[:, escape_sites_to_plot_inds]


    #Separate the haplotypes by timepoint
    timepoints = seq_info_df['time_label'].unique()
    hap_info_df = []
    for curr_time in timepoints:
        curr_escape_inds = seq_info_df[seq_info_df['time_label'] == curr_time]
        curr_escape_inds = curr_escape_inds['seq_index'].values
        curr_escape_haplotypes = escape_hap_arr[curr_escape_inds, :]


        #Now I will get the unique haplotypes and their counts
        unique_haplotypes, hap_counts = np.unique(curr_escape_haplotypes, 
                                                  axis=0, return_counts=True)

   
        if np.sum(hap_counts) < 10:
            continue
   
        #For each haplotype record the alleles it has
        for i in range(unique_haplotypes.shape[0]):
            #Get the haplotype and its count
            curr_hap = unique_haplotypes[i]
            hap_count = hap_counts[i]
            hap_count = hap_count / np.sum(hap_counts)

            #Now record each of the escape alleles and their identities
            escape_allele_count = 0
            res_region_allele_count = 0
            escape_allele_list = []
            for site in range(len(curr_hap)):
                derived_allele= curr_hap[site]
                if curr_hap[site] == 'Common' or derived_allele == '-':
                    continue
                else:
                    site_hxb2 = hxb2_nuc_coords_ath[escape_sites_to_plot_arr[site]]
                    seqArr_index = escape_sites_to_plot_arr[site]
                    hxb2_index = hxb2_nuc_coords_ath[seqArr_index]
                    mut_str = str(hxb2_index) + derived_allele
                    escape_allele_list.append(mut_str)

            # print(escape_allele_list)
            escape_allele_list = sort_escape_alleles(escape_allele_list)
            hap_info_df.append([curr_time, curr_par, hap_count, 
                                escape_allele_count, escape_allele_list,
                                res_region_allele_count])
    hap_info_df = pd.DataFrame(hap_info_df, columns=['time_label', 
                                                    'participant', 'hap_freq',
                                                    'escape_count', 'escape_alleles',
                                                    'res_region_count'])
    

    #Collapse any duplicate entries which can occur if there were haplotypes with
    #gaps
    hap_info_df['escape_alleles'] = hap_info_df['escape_alleles'].apply(lambda x: tuple(x))
    new_hap_info_df = []
    for name, group in hap_info_df.groupby(['time_label', 'participant', 'escape_alleles']):
        new_freq = np.sum(group['hap_freq'])
        other_vals = group.iloc[0]
        escape_allele_count = other_vals['escape_count']
        res_region_count = other_vals['res_region_count']
        new_hap_info_df.append([name[0], name[1], new_freq, escape_allele_count,
                                name[2], res_region_count])
    hap_info_df = pd.DataFrame(new_hap_info_df,
                                columns=['time_label', 'participant', 'hap_freq',
                                    'escape_count', 'escape_alleles', 'res_region_count'])
    
    #Do a little bookkeeping on the columns to make the saved dataframe nicer
    hap_info_df['time_label_int'] = hap_info_df['time_label'].apply(lambda x: int(x[1:])*7)
    hap_info_df['escape_alleles'] = hap_info_df['escape_alleles'].apply(lambda x: ','.join(x))
    hap_info_df.drop(columns = ['escape_count', 'res_region_count'], inplace = True)        

    #Make sure all combos of haplotype and timepoint are present
    all_times = hap_info_df['time_label'].unique()
    for curr_time in all_times:
        curr_time_df = hap_info_df[hap_info_df['time_label'] == curr_time]
        for curr_hap in hap_info_df['escape_alleles'].unique():
            if curr_hap not in curr_time_df['escape_alleles'].unique():
                recombination = hap_info_df[hap_info_df['escape_alleles'] == curr_hap]
                hap_info_df.loc[len(hap_info_df)] = [curr_time, curr_par, 
                                                     0, curr_hap, 
                                                     int(curr_time[1:])*7]
                

    all_hap_info_df.append(hap_info_df)

all_hap_info_df = pd.concat(all_hap_info_df)
all_hap_info_df.to_csv(outDir + outFolder + 'hap_info_df.csv', index = False)
