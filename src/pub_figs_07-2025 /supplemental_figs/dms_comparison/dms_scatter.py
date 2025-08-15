import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import data_util as du
import dataset_metadata
from matplotlib import rcParams

#This code makes the supplemental scatter plots of DMS escape scores vs. 
# in vivo allele frequencies for the 3BNC117 and 10-1074 cohorts

#slide size parameters
params = {'figure.figsize':(2, 2), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial', 'figure.titlesize': 6, 'axes.titlesize': 6}
rcParams.update(params)

# In this file I am going to be looking at enrichment in the DMS data
# I think it would be good to both use the exact AA mutations and also
# the sites

dms_indir = "../../../../data/Radford2025/"
TIMES_SEEN_FILTER_BF520 = 3
TIMES_SEEN_FILTER_TRO11 = 2

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0
# We will only use the two most frequent alleles at each site
MULTI_SEG = True
#Participant 2D1 only has sequence data out to hxb2 510
seq_end_2D1 = 510

# I need to start with getting a list of all the escape associated sites
inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']
inDir_1074 = '../../../../data/clyde_westfall_2024_final/10-1074/'
par_list_1074 = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD4K', '1HD5K', '1HD6K',
                  '1HD7K', '1HD10K', '1HD11K']
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24']

outFolder = 'dms_comparison/'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/'

###############################################################################
####################### Load and filter the DMS data ##########################
bf520_1074 = pd.read_csv(dms_indir + "BF520/10-1074_mut_effect.csv")
bf520_1074['strain'] = 'BF520'
bf520_1074['antibody'] = '10-1074'
bf520_3BNC = pd.read_csv(dms_indir + "BF520/3BNC117_mut_effect.csv")
bf520_3BNC['strain'] = 'BF520'
bf520_3BNC['antibody'] = '3BNC117'
bf520_3BNC = bf520_3BNC[bf520_3BNC['times_seen'] >= TIMES_SEEN_FILTER_BF520]

tro11_1074 = pd.read_csv(dms_indir + "TRO11/10-1074_mut_effect.csv")
tro11_1074['strain'] = 'TRO11'
tro11_1074['antibody'] = '10-1074'
tro11_3BNC = pd.read_csv(dms_indir + "TRO11/3BNC117_mut_effect.csv")
tro11_3BNC['strain'] = 'TRO11'
tro11_3BNC['antibody'] = '3BNC117'
tro11_3BNC = tro11_3BNC[tro11_3BNC['times_seen'] >= TIMES_SEEN_FILTER_TRO11]

# Concatenate the dataframes to gather all the DMS data together
dms_info_3BNC = pd.concat([bf520_3BNC, tro11_3BNC])
dms_info_3BNC = dms_info_3BNC[['site', 'wildtype', 'mutant', 'mutation',
                               'escape_mean', 'escape_median', 'escape_std',
                               'strain', 'antibody', 'times_seen']]

dms_info_1074 = pd.concat([bf520_1074, tro11_1074])
dms_info_1074 = dms_info_1074[['site', 'wildtype', 'mutant', 'mutation',
                                 'escape_mean', 'escape_median', 'escape_std',
                                    'strain', 'antibody', 'times_seen']]

###############################################################################
########################## Load the in vivo data ##############################
# For this, I need to get every allele and it's maximum change in frequency
allele_info_3BNC = []
for curr_par in par_list_3BNC:
    print(curr_par)

    # The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par  + \
                    '_NT.translated.fasta'
    
    # Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    time_counts = seq_info_df['time_label'].value_counts()
    time_counts = time_counts[time_counts > 10]
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]

    # Get the segregating sites
    all_seg_freq_dict = du.get_seg_sites(seqArr, set(),
                            allele_freq_thresh = ALLELE_FREQ_THRESH,
                            return_multi = MULTI_SEG,
                            include_drm_seg = True,
                            verbose = False)

    # Get the allele frequencies at each timepoint
    allele_freq_df, all_time_freqs, time_sample_sizes = du.seg_freq_by_timepoint(
                                                        all_seg_freq_dict,
                                                        seqArr,
                                                        seq_info_df)
    
    #Make a dataframe of frequency changes from day 0
    for curr_seg_site in all_seg_freq_dict:
        curr_seg_info = allele_freq_df[allele_freq_df['position'] == curr_seg_site]
        curr_seg_hxb2 = hxb2_nuc_coords_ath[curr_seg_site]

        if curr_par == '2D1' and curr_seg_hxb2 > seq_end_2D1:
            continue
        
        #Loop through each allele
        for curr_allele in all_seg_freq_dict[curr_seg_site]:
            #Get the frequency at day 0
            day_0_freq = curr_seg_info[curr_seg_info['time'] == 'D0']
            day_0_freq = day_0_freq[day_0_freq['allele'] == curr_allele[0]]['freqs'].values[0]


            #Get the frequencies at the other timepoints
            curr_allele_info = curr_seg_info[curr_seg_info['allele'] == curr_allele[0]]
            
            #Loop through the timepoints
            for curr_timepoint in curr_allele_info['time'].unique():
                if curr_timepoint == 'D0':
                    continue

                #Get the frequency at the current timepoint
                curr_freq = curr_allele_info[curr_allele_info['time'] == curr_timepoint]['freqs'].values[0]

                #If there is only one read, skip it
                curr_sample_size = time_sample_sizes[curr_timepoint]
                if curr_sample_size < 5:
                    continue

                #Calculate the frequency change
                freq_change = curr_freq - day_0_freq

                #Also check the gap percentage at that site
                gap_freq = seqArr[:, curr_seg_site]
                gap_freq = np.sum(gap_freq == '-')/len(gap_freq)


                #Append the info to the dataframe
                allele_info_3BNC.append([curr_seg_site, curr_seg_hxb2, 
                                        curr_allele[0], curr_par, 
                                        curr_timepoint, freq_change, 
                                        day_0_freq, curr_freq,
                                        gap_freq])

#Make a dataframe from the list
allele_info_3BNC = pd.DataFrame(allele_info_3BNC, columns = ['position', 
                                            'hxb2_coord_AA', 'derived_allele',
                                            'participant', 'timepoint', 
                                            'freq_change', 'day_0_freq',
                                            'curr_freq', 'col_gap_freq'])

#Get the maximum frequency change for each site
max_freq_df = []
for name, group in allele_info_3BNC.groupby(['hxb2_coord_AA', 'derived_allele',
                                             'participant']):
    max_freq = group['freq_change'].max()
    max_freq_row = group[group['freq_change'] == max_freq]
    max_freq_row = max_freq_row.iloc[0]
    max_freq_df.append(max_freq_row)
allele_info_3BNC = pd.DataFrame(max_freq_df)

###############################################################################
######################### Get the sites we can compare ########################
escape_info_df = []

# Get the escape scores from the DMS data
for index, row in allele_info_3BNC.iterrows():
    # Get the escape score for the current site
    coord_str = str(round(row['hxb2_coord_AA']))
    curr_score = dms_info_3BNC[(dms_info_3BNC['site'] == coord_str) &
                            (dms_info_3BNC['mutant'] == row['derived_allele'])]
    
    # If there is no escape score, continue
    if curr_score.empty:
        continue
    
    escape_info = [row['position'], row['hxb2_coord_AA'], row['derived_allele'],
                        row['participant'], row['timepoint'], row['freq_change'],
                        row['day_0_freq'], row['curr_freq'], row['col_gap_freq']]
    for strain_ind, indiv_strain in curr_score.iterrows():
        escape_info_curr = escape_info + [indiv_strain['strain'], indiv_strain['antibody'],
                                        indiv_strain['escape_mean'],
                                        indiv_strain['escape_median'], indiv_strain['escape_std']]
        escape_info_df.append(escape_info_curr)
    
# Make a dataframe from the list
fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
escape_info_df = pd.DataFrame(escape_info_df, columns = ['position', 'hxb2_coord_AA', 'derived_allele',
                                                        'participant', 'timepoint', 'freq_change',
                                                        'day_0_freq', 'curr_freq', 'col_gap_freq',
                                                        'strain', 'antibody', 'escape_mean',
                                                        'escape_median', 'escape_std'])

sns.scatterplot(data=escape_info_df, x='escape_mean', y='freq_change',
                s=15, ax = ax, marker = '.', color = 'black', alpha = 0.5)
plt.subplots_adjust(left=0.2, right=0.95, bottom=0.2, top=0.9)
ax.set_xlim(0, 3.5)
ax.set_ylim(0, 1)
ax.set_xlabel('DMS escape score')
ax.set_ylabel('Maximum frequency\nchange in vivo')
ax.xaxis.labelpad = 1
ax.yaxis.labelpad = 1
ax.set_title('3BNC117')

plt.savefig(outDir + outFolder + 'DMS_escape_vs_freq_change_3BNC117.png', dpi = 300)
plt.close()

###############################################################################
############################# Load the 10-1074 data ###########################
# For this, I need to get every allele and it's maximum change in frequency
allele_info_1074 = []
for curr_par in par_list_1074:
    print(curr_par)

    # The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir_1074 + curr_par + '/885_' + curr_par  + \
                        '_NT_filtered.translated.fasta'
    hxb2_res_positions = dataset_metadata.RESISTANCE_POS_AA_HXB2


    #Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    time_counts = seq_info_df['time_label'].value_counts()
    time_counts = time_counts[time_counts > 50]
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]

    # Get the segregating sites
    all_seg_freq_dict = du.get_seg_sites(seqArr, set(),
                            allele_freq_thresh = ALLELE_FREQ_THRESH,
                            return_multi = MULTI_SEG,
                            include_drm_seg = True,
                            verbose = False)

    # Get the allele frequencies at each timepoint
    allele_freq_df, all_time_freqs, time_sample_sizes = du.seg_freq_by_timepoint(
                                                        all_seg_freq_dict,
                                                        seqArr,
                                                        seq_info_df)
    
    #Make a dataframe of frequency changes from day 0
    for curr_seg_site in all_seg_freq_dict:
        curr_seg_info = allele_freq_df[allele_freq_df['position'] == curr_seg_site]
        curr_seg_hxb2 = hxb2_nuc_coords_ath[curr_seg_site]

        if curr_par == '2D1' and curr_seg_hxb2 > seq_end_2D1:
            continue
        
        #Loop through each allele
        for curr_allele in all_seg_freq_dict[curr_seg_site]:
            #Get the frequency at day 0
            day_0_freq = curr_seg_info[curr_seg_info['time'] == 'D0']
            day_0_freq = day_0_freq[day_0_freq['allele'] == curr_allele[0]]['freqs'].values[0]


            #Get the frequencies at the other timepoints
            curr_allele_info = curr_seg_info[curr_seg_info['allele'] == curr_allele[0]]
            
            #Loop through the timepoints
            for curr_timepoint in curr_allele_info['time'].unique():
                if curr_timepoint == 'D0':
                    continue

                #Get the frequency at the current timepoint
                curr_freq = curr_allele_info[curr_allele_info['time'] == curr_timepoint]['freqs'].values[0]

                #If there is only one read, skip it
                curr_sample_size = time_sample_sizes[curr_timepoint]
                if curr_sample_size < 5:
                    continue

                #Calculate the frequency change
                freq_change = curr_freq - day_0_freq

                #Also check the gap percentage at that site
                gap_freq = seqArr[:, curr_seg_site]
                gap_freq = np.sum(gap_freq == '-')/len(gap_freq)


                #Append the info to the dataframe
                allele_info_1074.append([curr_seg_site, curr_seg_hxb2, 
                                        curr_allele[0], curr_par, 
                                        curr_timepoint, freq_change, 
                                        day_0_freq, curr_freq,
                                        gap_freq])

#Make a dataframe from the list
allele_info_1074 = pd.DataFrame(allele_info_1074, columns = ['position', 
                                            'hxb2_coord_AA', 'derived_allele',
                                            'participant', 'timepoint', 
                                            'freq_change', 'day_0_freq',
                                            'curr_freq', 'col_gap_freq'])

#Get the maximum frequency change for each site
max_freq_df = []
for name, group in allele_info_1074.groupby(['hxb2_coord_AA', 'derived_allele',
                                             'participant']):
    max_freq = group['freq_change'].max()
    max_freq_row = group[group['freq_change'] == max_freq]
    max_freq_row = max_freq_row.iloc[0]
    max_freq_df.append(max_freq_row)
allele_info_1074 = pd.DataFrame(max_freq_df)

###############################################################################
######################### Get the sites we can compare ########################
escape_info_df = []

# Get the escape scores from the DMS data
for index, row in allele_info_1074.iterrows():
    # Get the escape score for the current site
    coord_str = str(round(row['hxb2_coord_AA']))
    curr_score = dms_info_1074[(dms_info_1074['site'] == coord_str) &
                            (dms_info_1074['mutant'] == row['derived_allele'])]
    
    # If there is no escape score, continue
    if curr_score.empty:
        continue
    
    escape_info = [row['position'], row['hxb2_coord_AA'], row['derived_allele'],
                        row['participant'], row['timepoint'], row['freq_change'],
                        row['day_0_freq'], row['curr_freq'], row['col_gap_freq']]
    for strain_ind, indiv_strain in curr_score.iterrows():
        escape_info_curr = escape_info + [indiv_strain['strain'], indiv_strain['antibody'],
                                        indiv_strain['escape_mean'],
                                        indiv_strain['escape_median'], indiv_strain['escape_std']]
        escape_info_df.append(escape_info_curr)
    
# Make a dataframe from the list
fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)        
escape_info_df = pd.DataFrame(escape_info_df, columns = ['position', 'hxb2_coord_AA', 'derived_allele',
                                                        'participant', 'timepoint', 'freq_change',
                                                        'day_0_freq', 'curr_freq', 'col_gap_freq',
                                                        'strain', 'antibody', 'escape_mean',
                                                        'escape_median', 'escape_std'])

sns.scatterplot(data=escape_info_df, x='escape_mean', y='freq_change', 
                s=15, ax = ax, marker = '.', color = 'black', alpha = 0.5)
plt.subplots_adjust(bottom=0.2, wspace=0.2, left = 0.2, right = 0.95)
ax.set_xlim(0, 3.5)
ax.set_ylim(0, 1)
ax.set_xlabel('DMS escape score')
ax.set_ylabel('Maximum frequency\nchange in vivo')
ax.set_title('10-1074')
ax.xaxis.labelpad = 1
ax.yaxis.labelpad = 1

plt.savefig(outDir + outFolder + 'DMS_escape_vs_freq_change_1074.png', dpi = 300)
plt.close()