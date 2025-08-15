import os
import sys
import numpy as np
import pandas as pd
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

import data_util as du
from par_data_class import Pardata

ALLELE_FREQ_THRESH = 0
MULTI_SEG = True
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2']
par_list = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']

TEXT_FONTSIZE = 6
params = {
          'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

#Today, I am going to be making allele frequency plots for all of the sites
#that we identified as selected.

focal_site_path =  '../../../../results/pub_figs_07-2025/supplemental_tables/3BNC117_all_data.csv'
identified_sites = pd.read_csv(focal_site_path)


alignment_file = '../../../../data/clyde_westfall_2024_final/3BNC117/%s/835_%s_NT.translated.fasta'
out_dir = '../../../../results/pub_figs_07-2025/supplemental_figs/allele_freq_traj/%s/'

#Get the focal sites we identified from our analysis


for curr_par in par_list:
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

    #Next, I will get only the focal sites that we are interested in
    focal_sites = identified_sites[identified_sites['participant'] == curr_par]

    #If MPL did not identify any sites for this participant, skip
    if focal_sites.empty:
        print(f"No identified sites for participant {curr_par}. Skipping...")
        continue



    #Get the allele frequencies for the focal sites
    focal_sites_arr = [hxb2_nuc_coords_hta[site] for site in focal_sites['hxb2_coord_AA'].values]

    allele_freq_df = allele_freq_df[allele_freq_df['position'].isin(focal_sites_arr)]
    unique_sites = allele_freq_df['position'].unique()
    unique_sites = [hxb2_nuc_coords_ath[x] for x in unique_sites]


    allele_freq_df['hxb2_coord_AA'] = allele_freq_df['position'].apply(
        lambda x: hxb2_nuc_coords_ath[x])
    allele_freq_df['time_label_int'] = allele_freq_df['time'].apply(
        lambda x: int(x[1:])*7)
    allele_freq_df.rename(columns={'allele': 'Allele'}, inplace=True)
    
    
    #Make a faceted plot for each participant of all of the allele frequencies
    num_plots = len(np.unique(allele_freq_df['hxb2_coord_AA']))

    fig, axs = plt.subplots(1, num_plots, figsize=(2 * num_plots, 2), sharex=True, sharey=True)

    for i, curr_site in enumerate(np.unique(allele_freq_df['hxb2_coord_AA'])):
        if num_plots == 1:
            axs = [axs]
        curr_out_dir = out_dir % curr_par
        if not os.path.exists(curr_out_dir):
            os.makedirs(curr_out_dir)
        
        curr_site_df = allele_freq_df[allele_freq_df['hxb2_coord_AA'] == curr_site]
        curr_site_df = curr_site_df.sort_values(by=['time_label_int', 'freqs'], ascending=False)
        sns.lineplot(data=curr_site_df,
                     x='time_label_int', y='freqs', ax=axs[i],
                     hue='Allele', palette='Set2', marker=None)
        axs[i].set_title(f"Site {curr_site}")
        axs[i].set_xlabel('Time (days since trial start)')
        axs[i].set_ylabel('Allele frequency')
        axs[i].set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(curr_out_dir + "allele_freq_traj_" + curr_par + ".png", dpi=300)
    plt.close()
