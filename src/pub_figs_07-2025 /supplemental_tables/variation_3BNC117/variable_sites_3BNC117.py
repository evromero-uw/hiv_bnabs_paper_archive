import sys
import numpy as np
import pandas as pd
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')

import data_util as du
from par_data_class import Pardata

#This file makes supplemental tables of the ancestral loss analysis in the
#3BNC117 cohort

#In this file I am going to replicate Elena Giorgi's ancestral loss analysis
alignment_file = '../../../../data/clyde_westfall_2024_final/3BNC117/%s/835_%s_NT.translated.fasta'
table_folder = '../../../../results/pub_figs_07-2025/supplemental_tables/variation_3BNC117/'
par_list_3BNC = ['2E7', '2E4', '2E5', '2C5', '2E2', '2C1', '2D1', '2E1', '2E3']
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W12']
ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

VARIABLE_LOOP_COORDS = [(459, 470),(384, 418), (294, 331), (155, 196), (129, 157)]
variable_loop_set = [set(range(start, end + 1)) for start, end in VARIABLE_LOOP_COORDS]
variable_loop_set = variable_loop_set[0].union(*variable_loop_set[1:])


MULTI_SEG = True
ALLELE_FREQ_THRESH = 0

#Alleles must have a frequency of at least this to be considered the day 0 majority
DAY_0_FREQ_MIN = 0.5

MAJORITY_LOSS_THRESH = 0.05

#Make a place to store the losses
majority_loss_df = []
#########################################################################################
for curr_par in par_list_3BNC:
    # First, I need to load the data
    curr_alignment_file = alignment_file % (curr_par, curr_par)
    participant_dat = Pardata(curr_alignment_file, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

    # Next, I will get a couple of items out of the dataset
    seqArr, seq_info_df = du.fasta_to_dataStructs(curr_alignment_file, clyde2024=True)
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath


    # Now I will get the allele frequency trajectories for the focal sites
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
    
    #Loop through each position and get the ancestral loss
    for name, group in allele_freq_df.groupby('position'):
        #Get the allele with max day 0 frequency
        day_0_freqs = group[group['time'] == 'D0']
        
        #Get the amino acid with the highest frequency at day 0
        max_row = day_0_freqs.loc[day_0_freqs['freqs'].idxmax()]
        max_allele = max_row['allele']

        #Check that the majority allele has a frequency above the threshold
        if max_row['freqs'] < DAY_0_FREQ_MIN:
            #Skip the site if there is no majority allele at day 0
            continue

        #Next I need to check how much this allele has changed in frequency
        majority_allele_freqs = group[group['allele'] == max_allele]
        majority_d0 = majority_allele_freqs[majority_allele_freqs['time'] == 'D0']['freqs'].values[0]
        majority_other = majority_allele_freqs[majority_allele_freqs['time'] != 'D0']['freqs'].values
        majority_loss = majority_d0 - majority_other

        #get only the positive losses
        majority_loss = majority_loss[majority_loss > 0]
        
        #calculate the percentage loss
        if len(majority_loss) > 0:
            majority_loss = majority_loss / majority_d0
            majority_loss = np.max(majority_loss)

            #Add the results to the DataFrame
            hxb2_coord = hxb2_nuc_coords_ath[name]
            majority_loss_df.append([curr_par, hxb2_coord, max_allele, majority_loss])


majority_loss_df = pd.DataFrame(majority_loss_df, columns=['participant', 'hxb2_coord_AA', 'majority_allele', 'majority_loss'])
majority_loss_df = majority_loss_df[majority_loss_df['majority_loss'] >= MAJORITY_LOSS_THRESH]
majority_loss_df['Variable Loop'] = majority_loss_df['hxb2_coord_AA'].apply(
    lambda x: int(x) in variable_loop_set)

majority_loss_summary = []
#I want to go through and summarize the losses across participants
for name, group in majority_loss_df.groupby('hxb2_coord_AA'):
    #Get the number of participants that have a loss at this site
    num_participants = group['participant'].nunique()

    par_list = group['participant'].tolist()
    freq_change_list = group['majority_loss'].tolist()
    freq_change_list = [round(x * 100, 1) for x in freq_change_list]
    freq_change_list = [str(x) + '%' for x in freq_change_list]
    variable_flag = group['Variable Loop'].iloc[0]

    
    majority_loss_summary.append([name, num_participants, 
                                  par_list, freq_change_list, variable_flag])
    
majority_loss_summary_df = pd.DataFrame(majority_loss_summary, 
                                        columns=['hxb2_coord_AA', 'num_participants',
                                                 'participant_list',
                                                 'freq_change_list', 'Variable Loop'])
majority_loss_summary_df = majority_loss_summary_df[majority_loss_summary_df['num_participants'] >= 3]
majority_loss_summary_df = majority_loss_summary_df.sort_values(by='num_participants', ascending=False)
majority_loss_summary_df.reset_index(drop=True, inplace=True)

#label any positions aligned to hxb2 gaps
string_labels = []
for curr_str in majority_loss_summary_df['hxb2_coord_AA']:
    if curr_str % 1 == 0:
        string_labels.append(str(int(curr_str)))
    else:
        remainder = curr_str % 1
        remainder = int(np.round(remainder * 1000))  # Convert to integer for display
        alphabet_remainder = ALPHABET[remainder]
        string_labels.append(f"{int(curr_str)}{alphabet_remainder}")
majority_loss_summary_df['hxb2_coord_AA'] = string_labels

majority_loss_summary_df.reset_index(drop=True, inplace=True)
majority_loss_summary_df.rename(columns={'hxb2_coord_AA': 'HXB2 Coordinate',
                                         'num_participants': 'Participant count',
                                         'participant_list': 'Participants with variation',
                                         'freq_change_list': 'Max loss (%) by participant'
                                         }, inplace=True)
majority_loss_summary_df.to_csv(table_folder + '3BNC117_ancestral_loss_summary.csv', index=False)