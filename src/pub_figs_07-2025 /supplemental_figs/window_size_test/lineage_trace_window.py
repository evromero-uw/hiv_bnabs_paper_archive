import sys
sys.path.append('../../../bin/')
sys.path.append('../../../bin/wrappers/')
sys.path.append('../../../data/clyde_westfall_2024_final/10-1074/')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

import matched_anc
import data_util as du
import dataset_metadata
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
from scipy import stats

#plot the estimates to show how accurate they are
MARKERSIZE = 15
LINEWIDTH = 0.5

params = {'figure.figsize': (2.5, 1.75), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

#I'll be sticking to the 60 nucleotide window, but in amino acids this is 20
WINDOW_SIZE = 300

# Here is the output folder
outDir = '../../../results/pub_figs_07-2025/figure_3/'

VARIABLE_LOOP_COORDS = [(384, 418), (294, 331), (155, 196), (129, 157)]
hxb2_var_set = [range(x[0], x[1] + 1) for x in VARIABLE_LOOP_COORDS]
hxb2_var_set = [item for sublist in hxb2_var_set for item in sublist]



# Here is the path to the data
inDir = '../../../data/clyde_westfall_2024_final/10-1074/'
par_list = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD4K', '1HD5K', '1HD6K', '1HD7K',
             '1HD9K', '1HD10K', '1HD11K']

time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2']

TIMEPOINT_DICT_OG = {'W1': ('D0'), 'W4' : ('D0', 'W1')}
REV_TIMEPOINT_DICT_OG = {('D0'): 'W1', ('D0', 'W1'): 'W4'}

TIMEPOINT_DICT_1HC2 = {'W4': ('D0')}
REV_TIMEPOINT_DICT_1HC2 = {('D0'): 'W4'}

TIMEPOINT_DICT_1HD6K = {'W8': ('D0')}
REV_TIMEPOINT_DICT_1HD6K = {('D0'): 'W8'}


###############################################################################
###############################################################################

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0
# We will only use the two most frequent alleles at each site
MULTI_SEG = True

# Initial frequency threshold
INITIAL_FREQ_THRESH = 0
use_from_file = False

if use_from_file:
    lineage_info_df = pd.read_csv(outDir + 'num_lineages_per_mutation_'+ \
                                  str(INITIAL_FREQ_THRESH) + '.csv')
else:
    #Here I am making a dataframe to store the number of lineages each mutation
    #is present on
    lineage_info_df = []

    for curr_par in par_list:
        print(curr_par)
        if curr_par == '1HC2':
            TIMEPOINT_DICT = TIMEPOINT_DICT_1HC2
            REV_TIMEPOINT_DICT = REV_TIMEPOINT_DICT_1HC2
        elif curr_par == '1HD6K':
            TIMEPOINT_DICT = TIMEPOINT_DICT_1HD6K
            REV_TIMEPOINT_DICT = REV_TIMEPOINT_DICT_1HD6K
        else:
            TIMEPOINT_DICT = TIMEPOINT_DICT_OG
            REV_TIMEPOINT_DICT = TIMEPOINT_DICT_OG

        #The path to the data and the resistance positions in hxb2 coordinates
        inFile = inDir + curr_par + '/885_' + curr_par  + \
                        '_NT_filtered.translated.fasta'
        hxb2_res_positions = dataset_metadata.RESISTANCE_POS_AA_HXB2


        #Now load the datastructures from the fasta
        seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
        hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, start_pos = 1)
        seq_info_df = seq_info_df[~seq_info_df['time_label'].isin(time_filter_out)]
        arr_res_set = set()
        for x in hxb2_res_positions:
            arr_res_set.add(x[0])
            arr_res_set.add(x[1])
        
        arr_res_set = [hxb2_nuc_coords_hta[x] for x in arr_res_set]
        arr_res_set = set(arr_res_set)
        arr_var_set = [hxb2_nuc_coords_hta[x] for x in hxb2_var_set]
        arr_var_set = set(arr_var_set)

        # I need to re calculate the segregating sites because I want to include
        # resistance sites
        all_seg_freq_dict = du.get_seg_sites(seqArr, arr_res_set,
                                allele_freq_thresh = ALLELE_FREQ_THRESH,
                                return_multi = MULTI_SEG,
                                include_drm_seg = True,
                                verbose = False)
        
        # Now I need to filter segregating sites to include only those with low
        # initial allele frequencies
        allele_freq_df, all_time_freqs, time_sample_sizes = du.seg_freq_by_timepoint(
                                                                all_seg_freq_dict,
                                                                seqArr,
                                                                seq_info_df)
        day_0_below_thresh = allele_freq_df[allele_freq_df['time'] == 'D0']
        day_0_below_thresh = day_0_below_thresh[day_0_below_thresh['freqs'] <= INITIAL_FREQ_THRESH]
        all_seg_freq_dict = {key: all_seg_freq_dict[key] for key in day_0_below_thresh['position']}


        # Now I will loop through every mutation and calculate the number of lineages
        # it stemmed from at day 0
        for curr_seg_site in all_seg_freq_dict.keys():

            # First, I will get the sequences in the window surrounding the site
            window_start = curr_seg_site - int(WINDOW_SIZE/2)
            if window_start < 0:
                arr_half = int(WINDOW_SIZE/2) + window_start
            else:
                arr_half = int(WINDOW_SIZE/2) 
            window_start = max(window_start, 0)
            window_end = curr_seg_site + int(WINDOW_SIZE/2)
            window_end = min(window_end, seqArr.shape[1])

            # I will only consider the sites that are within the window
            currArr_window_unmasked = seqArr[:, window_start:window_end]

            #I need to mask the focal loci for the ancestor matching
            currArr_window_first = seqArr[:, window_start:arr_half - 1]
            currArr_window_second = seqArr[:, arr_half + 2:window_end]
            currArr_window = np.concatenate((currArr_window_first, currArr_window_second), axis = 1)

            # Now, I need to collapse the identical sequences by timepoint
            seq_info_df['res_muts'] = 'None'
            collapsed_win_df = du.label_identical_seqs(currArr_window, seq_info_df)

            
            # Next, I need to separate out the day 0 sequences
            day_0_info = collapsed_win_df[collapsed_win_df['time_label'] == 'D0']

            # Loop through the timepoints
            for curr_time in TIMEPOINT_DICT.keys():
                #A dictionary where each key is an index at day 0 and each value is all the sequences
                #it is matched with as the closest match at the given timepoint
                closest_seqs = {}

                # Get the sequences for the current timepoint
                curr_time_info = collapsed_win_df[collapsed_win_df['time_label'] == curr_time]
                curr_time_inds = curr_time_info['seq_index'].to_numpy()
                curr_time_inds = curr_time_inds.astype(int)
                
                # Get the alleles at the site for the sequences at the current timepoint
                curr_time_alleles = currArr_window_unmasked[curr_time_inds, arr_half]
                curr_allele_set = set(curr_time_alleles) - set(['-'])

                # Loop through each allele at the site
                for curr_allele in curr_allele_set:
                    day_0_freq = allele_freq_df[(allele_freq_df['position'] == curr_seg_site) &\
                                (allele_freq_df['allele'] == curr_allele) &\
                                (allele_freq_df['time'] == 'D0')]['freqs'].values[0]
                    if day_0_freq > INITIAL_FREQ_THRESH:
                        continue
                    allele_inds = np.argwhere(curr_time_alleles == curr_allele)
                    allele_seq_inds = curr_time_inds[allele_inds].flatten()
                    
                    # A dictionary where each key is a sequence index at day 0 and
                    # each value is a number of the lineage that sequence is in
                    seq_lineage_dict = {}

                    # Find the ancestors for every sequence carrying that allele
                    for curr_allele_seq in allele_seq_inds:
                        curr_ancs, min_h_dist, percent_gap = matched_anc.find_closest_anc(
                                                                currArr_window,
                                                                curr_allele_seq,
                                                                day_0_info,
                                                                ignore_gap = True,
                                                                return_multi = True)
                        
                        # Now I need to count the number of unique day 0 lineages
                        # I did something similar in the 07-12-2024 script, and it would
                        # probably be best to move this functionality into the bin.
                        for curr_anc in curr_ancs:
                            if curr_anc in closest_seqs:
                                closest_seqs[curr_anc] += [curr_allele_seq]
                            else:
                                closest_seqs[curr_anc] = [curr_allele_seq]
    
                        seq_lineage_dict = matched_anc.count_lineages(closest_seqs)
                    
                    # Now I need to save the number of lineages for each allele
                    unique_lineages = [seq_lineage_dict.get(curr_ind) for curr_ind in allele_seq_inds]
                    unique_lineages = len(np.unique(np.array(unique_lineages)))



                    # I'll also get the frequency change for the allele
                    curr_freq = allele_freq_df[(allele_freq_df['position'] == curr_seg_site) &\
                                                (allele_freq_df['allele'] == curr_allele) &\
                                                (allele_freq_df['time'] == curr_time)]['freqs'].values[0]
                    freq_change = curr_freq - day_0_freq

                    if curr_seg_site in arr_res_set:
                        pos_label = 'Escape'
                    else:
                        pos_label = 'Other'
                    
                    
                    # if curr_seg_site in 
                    lineage_info_df.append([curr_par, curr_seg_site, curr_allele,
                                            curr_time, unique_lineages, freq_change,
                                            pos_label, hxb2_nuc_coords_ath[curr_seg_site]])




    # count the number of unique ones using the lineage counting code
    lineage_info_df = pd.DataFrame(lineage_info_df, columns = ['participant', 
                                                                'site', 'allele', 
                                                                'timepoint', 
                                                                'num_lineages',
                                                                'freq_change',
                                                                'pos_label', 'site_hxb2'])
    lineage_info_df.to_csv(outDir + 'num_lineages_per_mutation_'+ str(INITIAL_FREQ_THRESH) + '.csv')

###############################################################################
#################### Running the Regression Model #############################
###############################################################################
# Fit a linear regression model to the data
reg_results = smf.ols(formula = 'num_lineages ~ freq_change + freq_change:pos_label',
                      data = lineage_info_df).fit()
print(reg_results.summary())

# Conduct a likelihiood ratio test to compare the full model with a reduced model
reduced_model = smf.ols(formula = 'num_lineages ~ freq_change',
                        data = lineage_info_df).fit()
print(reduced_model.summary())
lr_test = reg_results.compare_lr_test(reduced_model)
print("Likelihood Ratio Test Results:")
print(lr_test)

print(f'Model p-values {reg_results.pvalues}')

#Get the t value for the interaction term
interaction_t_value = reg_results.tvalues['freq_change:pos_label[T.Other]']


# Get the p-value for the interaction term
interaction_p_value = stats.t.sf(np.abs(interaction_t_value),
                 df=reg_results.df_resid) * 2  # two-tailed test
print(f"Interaction p-value: {interaction_p_value}")

# Get the coefficients for the regression model
coefficients = reg_results.params
# Get the p-values for the regression model
p_values = reg_results.pvalues

# Save the model results to a text file
with open(outDir + 'regression_results.txt', 'w') as f:
    f.write('Regression Results:\n')
    f.write(reg_results.summary().as_text())


fig, ax = plt.subplots(1, 1)

sns.scatterplot(data = lineage_info_df, x = 'freq_change', y = 'num_lineages',
                hue = 'pos_label', palette=['black', 'red'],
                hue_order=['Other', 'Escape'], alpha = 0.5, s=MARKERSIZE,
                legend=False, marker = '.')

label_colors = {'Escape': 'red', 'Other': 'black'}
for curr_label in lineage_info_df['pos_label'].unique():
    x_pred = np.linspace(0, 1, 50)
    x_pred_df = pd.DataFrame({'freq_change': x_pred})
    x_pred_df['pos_label'] = curr_label

    # Get prediction results
    pred = reg_results.get_prediction(x_pred_df)
    # Get the confidence intervals for the predictions
    pred_summary = pred.summary_frame(alpha=0.05)
    print(pred_summary)
    ci_lower = pred_summary['mean_ci_lower']
    ci_upper = pred_summary['mean_ci_upper']
    pred_mean = pred_summary['mean']

    # Plot the predictions
    plt.plot(x_pred_df['freq_change'], pred_mean, color=label_colors[curr_label],
             label=curr_label, linewidth=LINEWIDTH)
    plt.fill_between(x_pred_df['freq_change'], ci_lower, ci_upper, color=label_colors[curr_label], alpha=0.3)

plt.xlim(0, max(lineage_info_df['freq_change']) + 0.05)
plt.xlabel('Frequency change')
plt.ylabel('Number of day 0 lineages')

Other_legend = plt.Line2D([0], [0], marker='.', color='w', markerfacecolor='black', markersize=MARKERSIZE/3, label='Other',
                          alpha = 0.5, markeredgewidth = 1, markeredgecolor = 'black')
Escape_legend = plt.Line2D([0], [0], marker='.', color='w', markerfacecolor='red', markersize=MARKERSIZE/3, label='Escape',
                           alpha = 0.5, markeredgewidth = 1, markeredgecolor = 'red')
handles = [Other_legend, Escape_legend]



plt.legend(handles= handles, title = 'Locus label', ncol = 1, bbox_to_anchor=(1, 0.7), loc = 'upper left',
           frameon=False, fontsize=6)
plt.subplots_adjust(bottom=0.2, left= 0.2, right = 0.7, top = 0.9)
plt.savefig(outDir + 'num_lineages_vs_freq_change_combined_' + str(INITIAL_FREQ_THRESH) + '.png', dpi = 300)

###############################################################################
