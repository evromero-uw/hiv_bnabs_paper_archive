import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd


import data_util as du

# I am going to look at the number of sequences with each glycosylation pattern
# Elena Giorgi sent me the glycosylation frequency data so I just need to check
# in that

infile = 'PNG.V5.freq.csv'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/glycan_shift_summary/'
glycan_freq_df = pd.read_csv(outDir + infile)
#Remove week 12 too since we think the antibody is decreasing in conversation by then
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24', 'W12']

inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']

sample_times_int = [1, 4, 8]
MAX_DIFF_THRESHOLD = 0.1  # Threshold for maximum difference to consider

###############################################################################
#I need to load the sequence data to get counts of sequences per timepoint
all_time_summary = []
for curr_par in par_list_3BNC:
    print(curr_par)

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    #Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
    time_counts = seq_info_df['time_label'].value_counts()
    par_time_counts = time_counts.to_dict()
    par_time_counts = pd.DataFrame.from_dict(par_time_counts, 
                                                orient='index', columns=['count'])
    par_time_counts['par'] = curr_par

    #give a new index
    par_time_counts.reset_index(inplace=True)
    all_time_summary.append(par_time_counts)

#Concatenate all the time counts
all_time_summary_df = pd.concat(all_time_summary, ignore_index=True)
all_time_summary_df['time_label_int'] = all_time_summary_df['index'].apply(lambda x: int(x[1:]))

#For each participant, get the majority pattern at day 0 
glycan_freq_max_df = []
for name, group in glycan_freq_df.groupby('ID'):
    #Get the maximum frequency
    max_freq = group['0'].idxmax()
    #Get the glycan pattern with that maximum frequency
    max_pattern = group.loc[max_freq]
    glycan_freq_max_df.append(max_pattern)
glycan_freq_max_df = pd.DataFrame(glycan_freq_max_df)

#I want to calculate the largest frequency difference for each glycan over time
for curr_time in sample_times_int:
    curr_time_diff = np.abs(glycan_freq_max_df['0'] - glycan_freq_max_df[str(curr_time)])
    glycan_freq_max_df[str(curr_time) + '_diff'] = curr_time_diff

print(glycan_freq_max_df)

#Mask the time points where the difference is only supported by 1 or zero sequences
for index, row in glycan_freq_max_df.iterrows():
    for curr_time in sample_times_int:
        #Get the time label
        time_label = str(curr_time)
        #Get the count for that time label
        count = all_time_summary_df[(all_time_summary_df['par'] == row['ID']) & 
                                    (all_time_summary_df['time_label_int'] == curr_time)]['count'].values[0]
        count = row[str(curr_time) + '_diff'] * count  # Adjust count by the frequency of the glycan

        #If the count is less than 2, mask the difference
        if count < 2:
            glycan_freq_max_df.at[index, time_label + '_diff'] = np.nan


#Get the maximum frequency difference
glycan_freq_max_df['max_diff'] = glycan_freq_max_df[[str(t) + '_diff' for t in sample_times_int]].max(axis=1)
glycan_freq_max_df['max_diff_time'] = glycan_freq_max_df[[str(t) + '_diff' for t in sample_times_int]].idxmax(axis=1)


#Get only those with a maximum difference greater than 0.1
glycan_freq_max_df = glycan_freq_max_df[glycan_freq_max_df['max_diff'] > MAX_DIFF_THRESHOLD]
print(glycan_freq_max_df)

