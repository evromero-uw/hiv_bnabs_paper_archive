import os
import re
import sys
import json
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
import matplotlib.pyplot as plt

import data_util as du
import dataset_metadata
from par_data_class import Pardata


# In this file I am making a diagram which incorporates the 3BNC117
# haplotypes and their frequencies. This version includes all sites
# identified by mpl, multiple encoding analysis, and "additional
# variation"

info_indir = '../../../../results/pub_figs_07-2025/supplemental_figs/inclusive_hap_diagrams/'
hap_info_path = info_indir + "hap_info_df.csv"
inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'

outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/inclusive_hap_diagrams/'

TEXT_FONTSIZE = 6
params = {
          'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

# First I need to get all the 3BNC117 haplotypes and their frequencies
par_list_3BNC = ['2E7', '2E4', '2E5', '2C5', '2E2', '2C1', '2D1', '2E1', '2E3', '2C1']



time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2']

HXB2_RES_POS = dataset_metadata.RESISTANCE_POS_AA_HXB2 

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will only use the two most frequent alleles at each site
MULTI_SEG = True

# Day 0 freq threshold
D0_FREQ_THRESH = 0.05

# Only look at timepoints with > 10 sequences
SAMPLE_SIZE_CUTOFF = 0

HAP_FREQ_THRESH = 0

# a color pallette for the haplotype backbones
color_palette = 'tab20'


###############################################################################
# Make a function that combines rare haplotypes into a single haplotype

def combine_rare_haps(curr_info_df, HAP_FREQ_THRESH):
    """ This function takes a dataframe of haplotype frequencies and combines
    haplotypes that never reach a frequency above HAP_FREQ_THRESH into a single
    haplotype.
    """
    haps_to_combine = []
    for name, group in curr_info_df.groupby(['escape_alleles']):
        max_freq = np.max(group['hap_freq'])
        if max_freq < HAP_FREQ_THRESH:
            haps_to_combine.append(name[0])

    combined_hap_df = []
    for name, group in curr_info_df.groupby(['time_label']):
        curr_combine = group[group['escape_alleles'].isin(haps_to_combine)]

        un_combine = group[~group['escape_alleles'].isin(haps_to_combine)]
        if curr_combine.shape[0] > 0:
            new_freq = np.sum(curr_combine['hap_freq'])
            escape_alleles = 'Rare variants'
            time_label = name[0]
            participant = curr_par
            time_label_int = int(time_label[1:])*7
            combined_hap_df.append([time_label, participant, new_freq, 
                                    escape_alleles, time_label_int, (), ()])
        for index, row in un_combine.iterrows():
            combined_hap_df.append([row['time_label'], row['participant'],
                                    row['hap_freq'], row['escape_alleles'],
                                    row['time_label_int'], row['escape_loci'],
                                    row['escape_AAs']])
    curr_info_df = pd.DataFrame(combined_hap_df, columns=['time_label', 
                                                    'participant', 'hap_freq',
                                                    'escape_alleles',
                                                    'time_label_int',
                                                    'escape_loci', 
                                                    'escape_AAs']) 
    return curr_info_df

###############################################################################
hap_info_df = pd.read_csv(hap_info_path)

#I need to convert the escape mut strings to a list of loci
escape_loci = []
escape_AAs = []
for curr_hap in hap_info_df['escape_alleles']:
    curr_escape_loci = []
    curr_escape_AA = {}
    if type(curr_hap) == str:
        indiv_muts = curr_hap.split(',')
        for curr_mut in indiv_muts:
            only_nums = re.sub(r'[^\d]+', '', curr_mut) 
            only_letters = re.sub(r'[^\D]+', '', curr_mut)
            curr_escape_AA[int(only_nums)] = only_letters
            curr_escape_loci.append(int(only_nums))
        
    curr_escape_loci.sort()
    escape_loci.append(tuple(curr_escape_loci))
    curr_escape_AA = [curr_escape_AA[k] for k in sorted(curr_escape_AA)]
    escape_AAs.append(curr_escape_AA)
hap_info_df['escape_loci'] = escape_loci
hap_info_df['escape_AAs'] = escape_AAs

#Next I will loop through each individual and pair of timepoints
#Loop through each participant and each of their sites
for curr_par in par_list_3BNC:
    currOutDir = outDir + curr_par + '/'

    if not os.path.exists(currOutDir):
        os.mkdir(currOutDir)

    seen_AAs = set()
    par_hap_info_df = hap_info_df[hap_info_df['participant'] == curr_par].copy()
    par_hap_info_df = combine_rare_haps(par_hap_info_df, HAP_FREQ_THRESH)

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    # First, I need to load the data
    participant_dat = Pardata(inFile, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)
    
    # Next, I will get a couple of items out of the dataset
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath

    # Now, I need to get the sample sizes at each timepoint
    time_sample_sizes = {}
    for time in seq_info_df['time_label'].unique():
        sample_size = len(seq_info_df[seq_info_df['time_label'] == time])
        if sample_size > SAMPLE_SIZE_CUTOFF or time == 'D0':
            time_sample_sizes[time] = sample_size

    # I'll use this to filter out timepoints with less than 10 sequences
    par_hap_info_df = par_hap_info_df[par_hap_info_df['time_label'].isin(time_sample_sizes.keys())]

    # Next, I'll assign each haplotype to a color
    unique_haps = par_hap_info_df['escape_alleles'].unique()
    unique_haps = [x for x in unique_haps if str(x) != 'nan']
    color_list = sns.color_palette(color_palette, len(unique_haps) + 1)
    #check if tab blue is in the color list, if so, remove it
    tab_blue = sns.color_palette('tab10')[0]
    if tab_blue in color_list:
        color_list.remove(tab_blue)
    else:
        color_list = color_list[:-1]  # remove the last color if not found
    hap_colors = list(zip(unique_haps, color_list))
    hap_color_dict = dict(hap_colors)
    hap_color_dict[np.nan] = 'gray'
    hap_color_dict['Rare variants'] = 'tab:blue'  # color for rare variants
    par_hap_info_df['color'] = par_hap_info_df['escape_alleles'].map(hap_color_dict)
    

    # Now, I am going to loop through each pair of adjacent timepoints
    time_int_list = par_hap_info_df['time_label_int']
    time_int_list = time_int_list.unique()
    time_int_list.sort()

    num_freq_ax = len(time_int_list)-1
    num_highlighter_ax = 1
    total_plot_num = num_freq_ax + num_highlighter_ax
    


    fig, freq_ax = plt.subplots(1, 1, figsize = (3.8/3, 1.5), dpi=300)

    # ax = ax.flatten()
    # highlighter_ax = []
    # freq_ax = ax[0]
    # highlighter_ax = ax[1]
    

    for i in range(len(time_int_list)-1):
        first_time = time_int_list[i]
        first_time_haps = par_hap_info_df[par_hap_info_df['time_label_int'] == first_time]
        second_time = time_int_list[i+1]
        second_time_haps = par_hap_info_df[par_hap_info_df['time_label_int'] == second_time]

        # Next I need to rank order the haplotypes by their frequencies at each time
        first_time_haps = first_time_haps.sort_values(by='escape_alleles')
        second_time_haps = second_time_haps.sort_values(by='escape_alleles')

        # Move the nan haplotypes and rare haplotypes to the start
        first_time_nan = first_time_haps[first_time_haps['escape_alleles'].isna()]
        first_time_haps_rare = first_time_haps[first_time_haps['escape_alleles'] == 'Rare variants']
        first_time_haps = first_time_haps[first_time_haps['escape_alleles'] != 'Rare variants']
        first_time_haps = first_time_haps[~first_time_haps['escape_alleles'].isna()]
        first_time_haps = pd.concat([first_time_nan, first_time_haps_rare, first_time_haps])


        second_time_nan = second_time_haps[second_time_haps['escape_alleles'].isna()]
        second_time_haps_rare = second_time_haps[second_time_haps['escape_alleles'] == 'Rare variants']
        second_time_haps = second_time_haps[second_time_haps['escape_alleles'] != 'Rare variants']
        second_time_haps = second_time_haps[~second_time_haps['escape_alleles'].isna()]
        second_time_haps = pd.concat([second_time_nan, second_time_haps_rare, second_time_haps])
        
        # I don't need to plot the haplotypes which are not present at either time
        present_time_1 = first_time_haps[first_time_haps['hap_freq'] > 0]
        present_time_1 = present_time_1['escape_alleles'].tolist()

        present_time_2 = second_time_haps[second_time_haps['hap_freq'] > 0]
        present_time_2 = present_time_2['escape_alleles'].tolist()
        all_haps = set(present_time_1 + present_time_2)
        
        first_time_haps = first_time_haps[first_time_haps['escape_alleles'].isin(all_haps)]
        second_time_haps = second_time_haps[second_time_haps['escape_alleles'].isin(all_haps)]

        # Now I want to plot the haplotypes which are present at both timepoints
        curr_y_t1 = 0
        curr_y_t2 = 0
        for index, row in first_time_haps.iterrows():
            curr_hap = row['escape_alleles']
            curr_freq = row['hap_freq']
            next_freq = second_time_haps[second_time_haps['escape_alleles'] == curr_hap]
            if next_freq.shape[0] == 0:
                if np.isnan(curr_hap):
                    next_freq = second_time_haps[second_time_haps['escape_alleles'].isna()]
            next_freq = next_freq['hap_freq'].values[0]

            prev_y_t1 = curr_y_t1
            prev_y_t2 = curr_y_t2
            curr_y_t1 += curr_freq
            curr_y_t2 += next_freq

            #now I need to fill between the previous y and the current y
            freq_ax.fill_between([first_time, second_time], [prev_y_t1, prev_y_t2], [curr_y_t1, curr_y_t2], color=row['color'],
                                alpha=0.5, linewidth=0)
            freq_ax.set_ylim(0, 1)
            freq_ax.set_xlabel('Trial timepoint')
            freq_ax.set_ylabel('Haplotype frequency')

            freq_ax.set_xticks(time_int_list)
            new_ticklabels = []
            for curr_time in time_int_list:
                time_str = int(curr_time/7)
                new_ticklabels.append(str(time_str))
            freq_ax.set_xticklabels(new_ticklabels)
            freq_ax.tick_params(axis='x', pad=1, rotation=45)
            freq_ax.set_title('Participant ' + curr_par, fontsize=TEXT_FONTSIZE)
                

        if i == len(time_int_list)-2:
            final_time_haps = second_time_haps
    
    plt.subplots_adjust(left=0.2, right=0.95, bottom=0.2)
    plt.savefig(currOutDir + curr_par +  '_3BNC117_hap_freqs.png', dpi=300)
    plt.close()


    

    # I need to get the number of escape haplotypes and the number of unique escape loci
    hap_count = len(par_hap_info_df['escape_alleles'].unique())
    if hap_count == 1:
        print('Only one haplotype for participant ' + curr_par + ', skipping highlighter plot')
        continue
    if hap_count == 0:
        print('No haplotypes for participant ' + curr_par + ', skipping highlighter plot')
        continue
    escape_loci = set()
    
    for index, row in par_hap_info_df.iterrows():
        curr_escape_loci = row['escape_loci']
        for curr_loc in curr_escape_loci:
            escape_loci.add(curr_loc)
    sorted_escape_loci = sorted(list(escape_loci))
    

    highlighter_width = len(escape_loci) * 0.2
    if curr_par == '2E1':
        highlighter_width += 0.2  # Add extra width for 2E1 to accommodate the letter labels
    print('Highlighter width: ' + str(highlighter_width))
    if highlighter_width > 3.8/2:
        highlighter_width = 3.8/2
    fig, highlighter_ax = plt.subplots(1, 1, figsize = (highlighter_width, 1.5), dpi=300)


    # This gives me the dimensions of the plot, from here I need to color boxes
    # based on the escape loci and annotate the escape AAs in text

    #Maybe I will make the patches first
    
    par_hap_info_df.drop_duplicates(subset=['escape_alleles'], inplace=True)
    par_hap_info_df.sort_values(by='escape_alleles', inplace=True)
    #Move nan haplotypes to the end
    par_hap_info_df_nan = par_hap_info_df[par_hap_info_df['escape_alleles'].isna()].reset_index(drop=True)
    par_hap_info_df_rare = par_hap_info_df[par_hap_info_df['escape_alleles'] == 'Rare variants'].reset_index(drop=True)

    par_hap_info_df = par_hap_info_df[~par_hap_info_df['escape_alleles'].isna()].reset_index(drop=True)
    par_hap_info_df = par_hap_info_df[par_hap_info_df['escape_alleles'] != 'Rare variants'].reset_index(drop=True)


    par_hap_info_df = pd.concat([par_hap_info_df_nan, par_hap_info_df_rare,  par_hap_info_df], ignore_index=True)
    par_hap_info_df.reset_index(drop=True, inplace=True)
    print(par_hap_info_df)

    # wildtype_row = pd.DataFrame({'escape_alleles': ['nan'], 'color': ['gray']})
    # final_time_haps = pd.concat([final_time_haps, wildtype_row], ignore_index=True)
    # final_time_haps.reset_index(drop=True, inplace=True)

    

    for index, row in par_hap_info_df.iterrows():
        if row['escape_alleles'] == 'Rare variants':
            rect = plt.Rectangle((0, index * 11), len(escape_loci) * 11 -1, 10, linewidth=1, edgecolor='none', facecolor= row['color'], alpha=0.5)
            highlighter_ax.add_patch(rect)
            #highlighter_ax.text(len(escape_loci)*11/2, index * 11 + 5, 'Rare escape haplotypes', fontsize=TEXT_FONTSIZE, color='black', ha='center', va='center')
            continue
        if str(row['escape_alleles']) == 'nan':
            print('in if ' + curr_par)
            rect = plt.Rectangle((0, index * 11), len(escape_loci) * 11 -1, 10, linewidth=1, edgecolor='none', facecolor= row['color'], alpha=0.5)
            highlighter_ax.add_patch(rect)
            #highlighter_ax.text(len(escape_loci)*11/2, index * 11 + 5, 'No putative escape mutations', fontsize=TEXT_FONTSIZE, color='black', ha='center', va='center')
            continue
        for j in range(len(escape_loci)):
            x_pos = j * 11
            y_pos = index * 11
            rect = plt.Rectangle((x_pos, y_pos), 10, 10, linewidth=1, edgecolor='none', facecolor= row['color'], alpha=0.5)
            highlighter_ax.add_patch(rect)
    #Now label the sites on top of the rectangles
    for j in range(len(escape_loci)):
        x_pos = j * 11
        y_pos = len(par_hap_info_df) * 11 + 1
        if curr_par == '2E1' and sorted_escape_loci[j] in [462, 463]:
            highlighter_ax.text(x_pos+6, y_pos + 2, str(sorted_escape_loci[j]) + str('a'), fontsize=TEXT_FONTSIZE, color='black', ha='center', va='center')
        else:
            highlighter_ax.text(x_pos+6, y_pos + 2, sorted_escape_loci[j], fontsize=TEXT_FONTSIZE, color='black', ha='center', va='center')
    


    for index, row in par_hap_info_df.iterrows():
        if str(row['escape_alleles']) == 'Rare variants' or str(row['escape_alleles']) == 'nan':
            continue
        #Now I need to annotate the escape AAs
        curr_escape_muts= row['escape_alleles']
        curr_escape_muts = curr_escape_muts.split(',')
        for curr_mut in curr_escape_muts:
            only_nums = re.sub(r'[^\d]+', '', curr_mut) 
            only_letters = re.sub(r'[^\D]+', '', curr_mut)

            # Now I need to get the x position
            locus_pos = sorted_escape_loci.index(int(only_nums))
            x_pos = locus_pos * 11
            y_pos = index * 11
            curr_AA = only_letters

            if len(curr_AA) > 1:
                # If the AA is a string, I will take the first letter
                curr_AA = curr_AA[1]

            #Now I will plot the text
            highlighter_ax.text(x_pos+5, y_pos+5, curr_AA, fontsize=TEXT_FONTSIZE, color='black', ha='center', va='center')


    plt.ylim(-1, hap_count*11)
    plt.xlim(-1, len(escape_loci)*11)

    highlighter_ax.get_xaxis().set_visible(False)
    highlighter_ax.get_yaxis().set_visible(False)
    highlighter_ax.spines["top"].set_visible(False)
    highlighter_ax.spines["right"].set_visible(False)
    highlighter_ax.spines["bottom"].set_visible(False)
    highlighter_ax.spines["left"].set_visible(False)


    plt.savefig(currOutDir + curr_par +  '_3BNC117_hap_legend.png', dpi=300)
    plt.close()