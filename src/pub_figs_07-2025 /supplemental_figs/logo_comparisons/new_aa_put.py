import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd

import data_util as du
from par_data_class import Pardata

#In this file I am going to count how many new AAs we identified at putative escape sites
#compared to the single genome sequencing data.

inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'
outDir_3BNC = '../../../../results/pub_figs_07-2025/supplemental_tables/new_putative_AAs/'
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24']
label_list = ['SGS', 'SMRT-UMI']
type_time = ['SGS D0', 'SMRT-UMI D0', 'SGS W4', 'SMRT-UMI W4']

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

# We will only use all alleles at each site
MULTI_SEG = True

unique_escape = []
unique_escape_pars = []

###############################################################################


#Get the putative escape loci
putative_escape_dir = '../../../../results/pub_figs_07-2025/supplemental_tables/3BNC117_all_data.csv'
putative_escape_df = pd.read_csv(putative_escape_dir, index_col=None)

out_data = []

for curr_par in par_list_3BNC:
    print(curr_par)

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    # First, I need to load the data
    participant_dat = Pardata(inFile, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

    # Next, I will get a couple of items out of the dataset
    seq_info_df = participant_dat.seq_info_df
    seqArr = participant_dat.seq_arr
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath
        
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]


    # I also want to get the allele frequencies at each time point
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
    
    # Now I will remove day 0 
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(['D0'])]

    # I need to get the sgs and porpid sequences separately
    seq_info_df['smrt_umi'] = ['KX' not in x for x in seq_info_df['orig_name']]
    seq_info_df['sample_type'] = ['SMRT-UMI' if x else 'SGS' for x in seq_info_df['smrt_umi']]

    # Separate the dataframes based on sample type
    sgs_df = seq_info_df[seq_info_df['sample_type'] == 'SGS']
    sgs_arr = seqArr[sgs_df['seq_index'].values, :]
    if sgs_arr.shape[0] == 0:
        print(f"No SGS sequences for {curr_par}, skipping...")
        continue
    smrt_umi_df = seq_info_df[seq_info_df['sample_type'] == 'SMRT-UMI']
    smrt_umi_arr = seqArr[smrt_umi_df['seq_index'].values, :]

    # Now I will loop through the putative escape sites
    curr_put_escape_df = putative_escape_df[putative_escape_df['participant'] == curr_par]
    curr_put_escape_df = curr_put_escape_df.drop_duplicates(subset=['hxb2_coord_AA'])
    for index, row in curr_put_escape_df.iterrows():
        # Get the putative escape site
        escape_site = row['hxb2_coord_AA']
        escape_site_arr = hxb2_nuc_coords_hta[escape_site]

        #check the sequences at the putative escape site
        sgs_alleles = np.unique(sgs_arr[:, escape_site_arr])
        smrt_umi_alleles = np.unique(smrt_umi_arr[:, escape_site_arr])

        #get the alleles unique to smart umi
        #I checked and there are no alleles that are unique to SGS at these sites
        smrt_umi_unique = np.setdiff1d(smrt_umi_alleles, sgs_alleles)
        if len(smrt_umi_unique) > 0:
            # If there are unique alleles, I will add them to the list
            for curr_aa in smrt_umi_unique:
                unique_escape.append(curr_aa)

                curr_allele_freq = allele_freq_df[allele_freq_df['position'] == hxb2_nuc_coords_hta[escape_site]]
                curr_allele_freq = curr_allele_freq[curr_allele_freq['allele'] == curr_aa].copy()
                curr_allele_freq.reset_index(drop=True, inplace=True)
                #Get the total sequence count from the dictionary for each timepoint
                curr_allele_freq['total_seqs'] = curr_allele_freq['time'].apply(
                    lambda x: time_sample_sizes[x])
                curr_allele_freq['num_seqs'] = curr_allele_freq['freqs'] * curr_allele_freq['total_seqs']

                for index, row in curr_allele_freq.iterrows():
                    if row['freqs'] == 0:
                        continue
                    # If the frequency is greater than 0, I will add it to the output
                    escape_site_hxb2 = hxb2_nuc_coords_ath[row['position']]
                    out_data.append([curr_par, escape_site_hxb2, curr_aa, row['freqs'], row['time'], int(row['num_seqs'])])

            unique_escape_pars.append(curr_par)

#get only unique participants
unique_escape_pars = list(set(unique_escape_pars))
#but keep all additional aas since they are across different sites

out_data = pd.DataFrame(out_data, columns=['Participant', 'Escape Site', 'Unique AA', 'Frequency', 'Timepoint', '# Sequences'])
out_data['# Sequences'] = out_data['# Sequences'].astype(int)
out_data['time_label_int'] = out_data['Timepoint'].apply(lambda x: int(x[1:])*7)
out_data = out_data.sort_values(by=['Participant', 'Escape Site', 'Unique AA', 'time_label_int'], ascending=True)
out_data.drop(columns=['time_label_int'], inplace=True)

#Remove sites only observed once across all timepoints
site_counts = out_data.groupby(['Participant', 'Escape Site', 'Unique AA']).aggregate({'# Sequences': 'sum'}).reset_index()
site_counts.rename(columns={'# Sequences': 'All Time Seqs'}, inplace=True)
out_data = out_data.merge(site_counts, on=['Participant', 'Escape Site', 'Unique AA'])
out_data = out_data[out_data['All Time Seqs'] > 1]
out_data.drop(columns=['All Time Seqs'], inplace=True)

#label any nonint ticks with the 
string_labels = []
for curr_str in out_data['Escape Site']:
    if curr_str % 1 == 0:
        string_labels.append(str(int(curr_str)))
    else:
        remainder = curr_str % 1
        remainder = int(np.round(remainder * 1000))  # Convert to integer for display
        alphabet_remainder = ALPHABET[remainder]
        string_labels.append(f"{int(curr_str)}{alphabet_remainder}")
out_data['Escape Site'] = string_labels



#Get the timepoint with the max frequency for each participant, escape site, and unique AA
freq_times = out_data.groupby(['Participant', 'Escape Site', 'Unique AA'])['Frequency'].idxmax()
freq_times = out_data.loc[freq_times, ['Participant', 'Escape Site', 'Unique AA', 'Timepoint', 'Frequency']]

#Now make the data wide and label the max frequency and timepoint
out_data_wide = out_data.pivot_table(index=['Participant','Escape Site', 'Unique AA'],
                                 columns='Timepoint', values='# Sequences',
                                 aggfunc='sum').reset_index()

out_data_wide = out_data_wide.merge(freq_times, on=['Participant', 'Escape Site', 'Unique AA'], how='left')
out_data_wide.rename(columns={'Frequency': 'Max Frequency',
                              'Timepoint': 'Max Timepoint'}, inplace=True)

out_data_wide.to_csv(outDir_3BNC + 'new_putative_escape_AAs.csv', index=False)

#check the max frequencies of 459 mutations that arent 459g
out_data_459 = out_data[out_data['Escape Site'] == '459a']
out_data_459 = out_data_459[out_data_459['Unique AA'] != 'G']
sum_freqs = out_data_459.groupby(['Participant', 'Timepoint'])['Frequency'].sum().reset_index()
print(sum_freqs)


