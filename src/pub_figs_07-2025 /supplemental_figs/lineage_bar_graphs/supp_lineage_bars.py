import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/10-1074/')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import matched_anc
import data_util as du
import dataset_metadata
from matplotlib import rcParams
from par_data_class import Pardata

#The code in this file is for making the bar graphs of lineages and their escape
#mutations in the 10-1074 cohort.

params = {'figure.figsize': (7.08, 8.66), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

inDir = '../../../../data/clyde_westfall_2024_final/10-1074/'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/'
outFolder = 'lineage_bars/'

#Exclude 1HD4k & 1HD10k since they are in the main text
par_list = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD5K', '1HD6K', '1HD7K',
             '1HD11K']


# TIMEPOINT_DICT = {'W1': ('D0'), 'W4' : ('D0', 'W1')}
# REV_TIMEPOINT_DICT = {('D0'): 'W1', ('D0', 'W1'): 'W4'}

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0
# We will only use the two most frequent alleles at each site
MULTI_SEG = True

# Initial frequency threshold
INITIAL_FREQ_THRESH = 0

# Window size for ancestor tracing
WINDOWSIZE = 300

###############################################################################
fig, axs = plt.subplots(4, 2)
axs = axs.flatten()


#Loop through each participant and make the network + highlighter plots
for ax_ind, SAMPLE_PAR in enumerate(par_list):
    curr_ax = axs[ax_ind]
    print(SAMPLE_PAR)
    if SAMPLE_PAR == '1HC2':
        TIMEPOINT_DICT = {'W4': 'D0'}
        REV_TIMEPOINT_DICT = {'D0': 'W4'}
        WEEK_NUM = '4'

    elif SAMPLE_PAR == '1HD6K':
        TIMEPOINT_DICT = {'W8': 'D0'}
        REV_TIMEPOINT_DICT = {'D0': 'W8'}
        WEEK_NUM = '8'

    else:
        TIMEPOINT_DICT = {'W1': ('D0'), 'W4' : ('D0', 'W1')}
        REV_TIMEPOINT_DICT = {('D0'): 'W1', ('D0', 'W1'): 'W4'}
        WEEK_NUM = '4'

    
    timepoint_set = set(TIMEPOINT_DICT.keys()).union(set(TIMEPOINT_DICT.values()))

    #The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir + SAMPLE_PAR + '/885_' + SAMPLE_PAR  + '_NT_filtered.fasta'
    hxb2_res_positions = dataset_metadata.RESISTANCE_POS_NT_HXB2


    #Construct the data object and load the data
    participant_dat = Pardata(inFile, 'clyde2024', SAMPLE_PAR)
    participant_dat.load_data_10_1074(hxb2_res_positions, ALLELE_FREQ_THRESH, MULTI_SEG)
    seq_info_df = participant_dat.seq_info_df
    seq_arr = participant_dat.seq_arr
    seg_freq_dict = participant_dat.seg_freq_dict
    arr_res_set = participant_dat.arr_res_set


    #Slice Window/2 bp on either side of the resistance positions
    arr_res_positions = participant_dat.arr_res_positions
    arr_res_min = min(arr_res_positions[0][0], arr_res_positions[1][0])
    arr_res_max = max(arr_res_positions[0][1], arr_res_positions[1][1])
    arr_res_window_min = arr_res_min - int(WINDOWSIZE/2)
    arr_res_window_max = arr_res_max + int(WINDOWSIZE/2)
    seq_arr = seq_arr[:, arr_res_window_min:arr_res_window_max]


    #Now, for each sequence, we'll get a list of the segregating mutations 
    #Perhaps we should make a dictionary where each key is a sequence 
    #index and each value is a set of minor alleles that sequence contains
    seg_set_dict = du.make_seg_sets(seg_freq_dict, participant_dat.seq_arr, seq_info_df,
                                    arr_res_window_min, arr_res_window_max, arr_res_set)


    #Make a dictionary mapping sequence names to their array indices
    index_to_name = dict(zip(seq_info_df['seq_index'], seq_info_df['orig_name']))

    #Get the hxb2 position of the start of the array now that it's been sliced
    hxb2_res_window_min = participant_dat.hxb2_nuc_coords_ath[arr_res_window_min]
    hxb2_res_window_max = participant_dat.hxb2_nuc_coords_ath[arr_res_window_max]


    #A dictionary where each key is an index at day 0 and each value is all the sequences
    #it is matched with as the closest match at week 1
    closest_seqs = {}
    #The dictionary above with each key and value reversed
    closest_seqs_reversed = {}


    #A dictionary where each key is an index at day 0 and each value is the minimum
    #hamming distance between that sequence and its closest match at week 1
    min_dist_dict = {}

    #Subset the window and mask the resistance positions
    masked_seq_arr, removed_arrs = du.mask_hxb2_coords(seq_arr,
                                        start_pos = hxb2_res_positions[0][0],
                                        end_pos = hxb2_res_positions[0][1],
                                        arr_start_pos = hxb2_res_window_min,
                                        second_site = hxb2_res_positions[1])
    
    #Make a dataframe counting the number of identical sequences at each timepoint
    #This will be used for labeling identical sequences
    labeled_counts = du.label_identical_seqs(masked_seq_arr, seq_info_df)
    labeled_counts['num_identical'] = [len(x) + 1 for x in labeled_counts['identical_seqs']]
    # labeled_counts = labeled_counts[['orig_name', 'seq_index', 'num_identical']]

    for curr_timepoint in TIMEPOINT_DICT.keys():

        #Get the sequences from the current timepoint
        curr_time_info = seq_info_df[seq_info_df['time_label'] == curr_timepoint].copy()


        #Get the sequences from the previous timepoint
        prev_timepoint = TIMEPOINT_DICT[curr_timepoint]
        if isinstance(prev_timepoint, str):
            prev_timepoint = [prev_timepoint]

        prev_time_info = seq_info_df[seq_info_df['time_label'].isin(prev_timepoint)].copy()
        prev_time_info = prev_time_info.reset_index(drop = True)
        prev_time_info = du.label_identical_seqs(masked_seq_arr, prev_time_info, ignore_time = True)

        
        curr_time_info = curr_time_info.reset_index(drop = True)
        curr_time_info = du.label_identical_seqs(masked_seq_arr, curr_time_info, ignore_time = True)
            
        #For each sequence, find it's closest match in the previous timepoint
        for index, row in curr_time_info.iterrows():
            curr_seq_ind = row['seq_index']

            #Find the closest ancestor
            closest_anc_ind, min_hdist, percent_gap = \
                matched_anc.find_closest_anc(masked_seq_arr, curr_seq_ind,
                                            prev_time_info, return_multi = True, ignore_gap = True)
            
            #Get the info for the closest ancestor and save it
            min_dist_dict[curr_seq_ind] = min_hdist

            if isinstance(closest_anc_ind, np.int64):
                closest_anc_ind = [closest_anc_ind]


            closest_seqs_reversed[curr_seq_ind] = closest_anc_ind


            for i in closest_anc_ind:
                if i not in closest_seqs.keys():
                    closest_seqs[i] = [curr_seq_ind]
                else:
                    closest_seqs[i].append(curr_seq_ind)
    
    seq_lineage_dict = matched_anc.count_lineages(closest_seqs)

    #Label all of the sequences with their lineage number
    labeled_counts['lineage_id'] = labeled_counts['seq_index'].map(seq_lineage_dict)
    labeled_counts = labeled_counts[labeled_counts['time_label'].isin(timepoint_set)].copy()
    labeled_counts = labeled_counts.dropna(subset = ['lineage_id'])
    
    #expand the labeled counts dataframe so that each sequence is represented by a row
    print(labeled_counts)
    new_labels = []
    for index, row in labeled_counts.iterrows():
        if row['identical_seqs'] == ():
            new_labels.append([row['seq_index']])
        else:
            new_labels.append([row['seq_index']] + list(row['identical_seqs']))
    labeled_counts['identical_seqs'] = new_labels
    labeled_counts = labeled_counts.explode('identical_seqs')

    #Get the resistance mutations for each sequence
    new_res_mut_list = []
    for index, row in labeled_counts.iterrows():
        if row['res_muts'] == [None]:
            new_res_mut_list.append('Susceptible')
        else:
            new_res_mut_list.append(row['res_muts'][0])

    labeled_counts['res_muts'] = new_res_mut_list
    labeled_counts['lineage_id'] = labeled_counts['lineage_id'].astype('str')
    labeled_counts = labeled_counts[~labeled_counts['time_label'].isin(['D0', 'W1'])]

    labeled_counts.rename(columns = {'res_muts': 'Escape Identity'}, inplace = True)


    sns.histplot(data = labeled_counts, y = 'lineage_id', 
                  hue = 'Escape Identity', multiple = 'stack', shrink = 0.8,
                  palette=dataset_metadata.RESISTANCE_HUE_DICT,
                  ax = curr_ax)
    curr_ax.get_legend().remove()
    if ax_ind in [0, 2, 4, 6]:
        curr_ax.set_ylabel('Day 0 lineage')
    else:
        curr_ax.set_ylabel('')
    curr_ax.set_xlabel('Number of week ' + WEEK_NUM + ' descendants')
    curr_ax.yaxis.set_major_formatter(plt.NullFormatter())
    curr_ax.set_title(SAMPLE_PAR)
    curr_ax.tick_params(left = False)
plt.subplots_adjust(left= 0.08, right = 0.98, wspace = 0.05, hspace = 0.5, top= 0.95, bottom = 0.05)
plt.savefig(outDir + outFolder + 'lineage_escape_supp.png', dpi = 300)
plt.close()