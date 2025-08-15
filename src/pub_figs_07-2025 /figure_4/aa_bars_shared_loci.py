import os
import sys
sys.path.append('../../../bin/')
sys.path.append('../../../bin/wrappers/')
sys.path.append('../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams


import data_util as du
import dataset_metadata
from par_data_class import Pardata

#The code in this file makes figure 4 panel B

#Today I am going to investigate the examples in the 3BNC117 data set where
#participants have the same loci

#slide size parameters
params = {'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial', 'figure.titlesize': 6, 'axes.titlesize': 6}
rcParams.update(params)

inDir_3BNC = '../../../data/clyde_westfall_2024_final/3BNC117/'
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24']
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']

#Here is the output folder
outFolder = 'figure_4/'
outDir = '../../../results/pub_figs_07-2025/'
site_inpath = outDir + outFolder + '3BNC117_filtered_data.csv'

color_palette = sns.color_palette("tab20c")


# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0.01
# We will only use the two most frequent alleles at each site
MULTI_SEG = True

#I want to exclude alleles in the V5 region
V5_lower = 459
V5_upper = 470

RARE_FREQ_THRESH = 0.02
COLOR_PAL = 'Set2'
colors_seen = set()
################################################################################
#Load all of the sites shared across participants
site_par_df = pd.read_csv(site_inpath)
# site_par_df['Rounded_AA'] = site_par_df['hxb2_coord_AA'].round().astype(int)
# site_par_df = site_par_df.drop_duplicates(subset=['Rounded_AA', 'participant'])
site_par_df = site_par_df.drop_duplicates(subset=['hxb2_coord_AA', 'participant'])

#Get the sites that are shared across participants
shared_sites = site_par_df['hxb2_coord_AA'].value_counts()
shared_sites = shared_sites[shared_sites > 1]
shared_sites = shared_sites.index.tolist()
site_par_df = site_par_df[site_par_df['hxb2_coord_AA'].isin(shared_sites)]
site_par_df = site_par_df[~site_par_df['hxb2_coord_AA'].apply(lambda x: V5_lower <= x <= V5_upper)]


#I will get all sites with the AA rounded to the nearest integer
shared_sites = site_par_df['hxb2_coord_AA'].unique()

###############################################################################
#Get the list of participants

all_chosen_freqs = []

#I am going to go through each of the participants and get the allele
#frequencies at each of the sites of interest
for curr_par in par_list_3BNC:
    print(curr_par)

    #The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par  + \
                    '_NT.translated.fasta'
    
    #Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
    time_counts = seq_info_df['time_label'].value_counts()
    # time_counts = time_counts[time_counts > 20]
    seq_info_df = seq_info_df[seq_info_df['time_label'].isin(time_counts.index)]

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
    
    
    allele_freq_df['hxb2_AA_pos'] = allele_freq_df['position'].map(hxb2_nuc_coords_ath)

    print(allele_freq_df)

    #Now I will filter the allele frequencies to only include the sites of interest
    chosen_freqs = allele_freq_df[allele_freq_df['hxb2_AA_pos'].isin(shared_sites)].copy()
    #If there are no sites that varied, then skip this participant
    if len(chosen_freqs) == 0:
        print('No sites found for participant ' + curr_par)
        continue
    chosen_freqs['participant'] = curr_par
    all_chosen_freqs.append(chosen_freqs)

all_chosen_freqs = pd.concat(all_chosen_freqs, ignore_index=True)


num_colors = len(all_chosen_freqs['allele'].unique())
color_palette = sns.color_palette(COLOR_PAL, num_colors)
color_dict = {allele: color_palette[i] for i, allele in enumerate(all_chosen_freqs['allele'].unique())}

color_dict['-'] = 'black'
# color_dict['X'] = 'gray'

#Now I am going to make a plot for each site of interest
for curr_site in shared_sites:
    #First I will get the data from the site
    curr_site_df = all_chosen_freqs[all_chosen_freqs['hxb2_AA_pos'] == curr_site]
    site_participants = curr_site_df['participant'].unique()

    #Now make the plot
    print(curr_site)
    print(site_participants)
    print('------------')
    fig, ax = plt.subplots(1, len(site_participants), sharey=True,
                figsize = (0.75 * len(site_participants), 1.25))
    plot_freqs_df = []

    #Now go through each participant
    for index, curr_par in enumerate(site_participants):
        if len(site_participants) == 1:
            curr_ax = ax
        else:
            curr_ax = ax[index]
        curr_par_site_df = curr_site_df[curr_site_df['participant'] == curr_par]
        day_0_rare = curr_par_site_df[curr_par_site_df['time'] == 'D0']
        freqs_to_plot = curr_par_site_df.copy()
        


        #Check if there were gaps we filtered out
        gap_rows = []
        for name, group in freqs_to_plot.groupby('time'):
            if sum(group['freqs']) < 1:
                gap_freq = 1 - sum(group['freqs'])
                position = group['position'].values[0]
                time = name
                hxb2_AA_pos = group['hxb2_AA_pos'].values[0]
                participant = group['participant'].values[0]
                gap_rows.append([position, '-', time, gap_freq, hxb2_AA_pos, participant])
        if len(gap_rows) > 0:
            gap_df = pd.DataFrame(gap_rows, columns = ['position', 'allele', 'time', 'freqs', 'hxb2_AA_pos', 'participant'])
            for curr_time in freqs_to_plot['time'].unique():
                if curr_time not in gap_df['time'].values:
                    position = gap_df['position'].values[0]
                    gap_rows.append([position, '-', curr_time, 0, curr_site, curr_par])
            gap_df = pd.DataFrame(gap_rows, columns = ['position', 'allele', 'time', 'freqs', 'hxb2_AA_pos', 'participant'])
            freqs_to_plot = pd.concat([freqs_to_plot, gap_df], ignore_index=True)
        
        print(freqs_to_plot)
        freqs_to_plot['time_int'] = freqs_to_plot['time'].apply(lambda x: int(x[1:]) if x.startswith('W') else 0)
        freqs_to_plot.sort_values('time_int', inplace=True, ascending=True)

        #We will make a stacked bar plot of the frequencies
        bottom = np.zeros(len(freqs_to_plot['time_int'].unique()))

        for name, group in freqs_to_plot.groupby('allele'):
            if group['freqs'].sum() > 0.01:
                print(group)
                colors_seen.add(name)
            my_plot = curr_ax.bar(group['time'], group['freqs'], bottom = bottom, color = color_dict[name])
            
            bottom += group['freqs'].values
            print(bottom)
        
        curr_ax.set_title(curr_par)
        if index == 0:
            curr_ax.set_ylabel('Frequency')
        else:
            curr_ax.set_yticklabels([])
        curr_ax.set_xlabel('')
        curr_ax.set_ylim(0, 1)
        curr_ax.set_xticklabels(freqs_to_plot['time_int'].unique())


    plt.subplots_adjust(bottom=0.6, top=0.9, left=0.3, right=0.95) 
    dash_removed = False  
    if '-' in colors_seen:
        dash_removed = True
        colors_seen.remove('-')
    colors_seen_formatted = list(colors_seen)
    colors_seen_formatted = [name.replace(' ', '_') for name in colors_seen_formatted]
    if dash_removed:
        colors_seen_formatted.append('-')

    # colors_seen_formatted = set(colors_seen_formatted)
    # colors_seen_formatted = list(colors_seen_formatted)
    colors_seen_formatted.sort()
    plt.legend(handles = [mpl.patches.Patch(color = color_dict[name], label = name) for name in colors_seen_formatted],
                title = 'Allele', loc = 'lower center', bbox_to_anchor = (-1.25, -0.9), ncol = 7, columnspacing = 0.8)
    # plt.suptitle('Site ' + str(curr_site))
    # plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, left = 0.05, right = 0.95, top = 0.85, bottom = 0.4)
    print('Saving figure for site ' + str(curr_site))
    plt.savefig(outDir + outFolder + 'site_' + str(curr_site) + '.png', dpi = 300)
    plt.close()