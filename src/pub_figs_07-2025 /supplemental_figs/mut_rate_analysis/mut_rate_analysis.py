import sys
import json
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/10-1074/')
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import linregress

import data_util as du
import dataset_metadata

#In this file I am making the supplemental plots of expected mutations (based on
# Zanini 2017 mutation rates) at HXB2 sites 332 and 334 in the 10-1074 cohort

TIMES_SEEN = 2

params = {'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial', 'figure.titlesize': 6,}
rcParams.update(params)

inDir = '../../../../results/pub_figs_03-2025/figure_2/'
inDir_DMS = '../../../../data/Radford2025/TRO11/'
color_dir = "../../../../data/no_green_scheme.json"
color_scheme2 = json.load(open(color_dir))
color_scheme2 = color_scheme2['colors']

# color_palette = sns.color_palette("tab20")
# color_scheme = {x: color_palette[i] for i, x in enumerate(color_scheme.keys())}
color_scheme = dataset_metadata.RESISTANCE_HUE_DICT
site_prefix_dict = {325: 'D/N325', 332: 'N332', 334:'S334'}

#Load the dataset
all_proportions = pd.read_csv(inDir + 'obs_vs_predicted_props.csv')
day0_resistance = pd.read_csv(inDir + 'day0_resistance.csv')
dms_entry_effects = pd.read_csv(inDir_DMS + 'TZM-bl_entry_func_effects.csv')


outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/mut_rate_analysis/'

######################## Plot the distributions ###############################

#Plot stacked bars of the distributions of the expected proportions

for curr_site in [332, 334]:
    # Filter to get the data for the current site
    curr_site_data = all_proportions[
                            all_proportions['site'] == curr_site]
    
    # Get the participants for the current site
    all_pars = curr_site_data['participant'].unique()
    fig, axs = plt.subplots(1, len(all_pars), 
                            figsize = (0.7 * len(all_pars), 2.5),
                            sharey = True)
    axs = axs.flatten()


    # Make a stacked bar of expected vs observed for each individual
    for i, curr_par in enumerate(all_pars):
        curr_ax = axs[i]
        curr_par_data = curr_site_data[curr_site_data['participant'] == curr_par]
        #It's very important to sort the data by AA so that the bars are in the same order
        curr_par_data = curr_par_data.sort_values(['AA'])

        curr_res_data = day0_resistance[day0_resistance['participant'] == curr_par]
        curr_res_data = curr_res_data[curr_res_data['site'] == curr_site]

        #stacked bar
        #This variable will keep track of the bottom of the bars, we need one entry for each bar
        bottom = np.array([0, 0])
        for curr_AA in sorted(curr_par_data['AA'].unique(), reverse= True):
            curr_AA_data = curr_par_data[curr_par_data['AA'] == curr_AA]
            exp_prop_weighted = curr_AA_data['exp_prop_weighted'].values[0]
            obs_prop = curr_AA_data['obs_prop'].values[0]

            #Add texture if resistant at day 0
            texture = ''
            if len(curr_res_data) > 0:
                if curr_AA in curr_res_data['AA'].values:
                    texture = '//'
            
            curr_AA_str = site_prefix_dict[curr_site] + curr_AA

            #Plot the bars
            curr_ax.bar(0, exp_prop_weighted, color = color_scheme2[curr_AA],
                    label = 'Expected', bottom = bottom[0], edgecolor = 'white',
                    hatch = texture)
            curr_ax.bar(1, obs_prop, color = color_scheme2[curr_AA],
                    label = 'Observed', bottom = bottom[1], edgecolor = 'white',
                    hatch = texture)
            bottom = bottom + np.array([exp_prop_weighted, obs_prop])
        
        curr_ax.set_title(curr_par)
        curr_ax.set_xticks([0, 1])
        curr_ax.set_xticklabels(['Expected', 'Observed'], rotation = 90)
    
    axs[0].set_ylabel('Proportion')

    #Make a legend
    legend_labels = curr_site_data['AA'].unique()
    legend_labels.sort()
    legend_colors = [color_scheme2[x] for x in legend_labels]
    legend_handlers = [plt.Rectangle((0,0),1,1, color = x) for x in legend_colors]
    
    plt.subplots_adjust(left = 0.05, wspace=0.2, top=0.85, bottom=0.2)
    plt.legend(legend_handlers, legend_labels, bbox_to_anchor=(1.05, 1), 
               loc='upper left', title = 'Amino acid', frameon = False)
    
    plt.suptitle('Site ' + str(curr_site))
    plt.savefig(outDir + 'mutation_rates_' + str(curr_site) + '.png', dpi = 300)
    plt.close()

######################## Plot proportion enriched vs dms ######################

#Plot the proportion enriched vs the DMS entry effects
dms_entry_effects = dms_entry_effects[dms_entry_effects['times_seen']>= TIMES_SEEN]
dms_entry_effects = dms_entry_effects.rename(columns = {'mutant': 'AA'})

all_proportions['obs_exp_ratio'] = all_proportions['obs_prop'] / all_proportions['exp_prop_weighted']
all_proportions['obs_exp_ratio'] = np.log(all_proportions['obs_exp_ratio'])
all_proportions['site'] = all_proportions['site'].astype(str)


for curr_site in all_proportions['site'].unique():
    # Filter to get the data for the current site
    curr_site_data = all_proportions[all_proportions['site'] == curr_site]
    curr_site_data = curr_site_data.sort_values(['AA'])
    zero_data = curr_site_data[curr_site_data['obs_exp_ratio'] == -np.inf]
    curr_site_data = curr_site_data[curr_site_data['obs_exp_ratio'] != -np.inf]

    #Merge the DMS entry effects with the proportions
    curr_site_data = curr_site_data.merge(
        dms_entry_effects[['site', 'AA', 'effect']], how = 'left',
                            on = ['site', 'AA'])
    
    corr_result = linregress(curr_site_data['obs_exp_ratio'],
                            curr_site_data['effect'], alternative = 'greater')
    regression_x = np.linspace(int(min(curr_site_data['obs_exp_ratio']) -1),
                            int(max(curr_site_data['obs_exp_ratio'])) +1, 100)
    regression_y = corr_result.intercept + corr_result.slope * regression_x



    plt.figure(figsize = (2.5, 2))
    for curr_AA in curr_site_data['AA'].unique():
        curr_AA_data = curr_site_data[curr_site_data['AA'] == curr_AA]
        curr_AA_str = site_prefix_dict[int(curr_site)] + curr_AA
        if len(curr_AA_data) > 0:
            plt.scatter(x = 'obs_exp_ratio', y = 'effect',
                    data = curr_AA_data, color = color_scheme[curr_AA_str],
                    marker = 'o', s = 2, label = curr_AA)

    plt.plot(regression_x, regression_y, color = 'black', linewidth = 0.5)
    plt.annotate(r'r = ' + str(round(corr_result.rvalue, 2)) +
                '\np = ' + '{:.2e}'.format(corr_result.pvalue),
                xy = (0.05, 0.95), xycoords='axes fraction',
                fontsize = 6, ha = 'left', va = 'top')
    plt.subplots_adjust(left = 0.21, right = 0.8, top = 0.90, bottom = 0.2)
    plt.legend(title = 'Amino acid', loc = 'upper left',
            bbox_to_anchor=(1, 1), fontsize = 6, frameon = False)
    plt.title('Site ' + str(curr_site))
    plt.ylabel('DMS entry effect')
    plt.xlabel('In vivo enrichment')
    plt.savefig(outDir + 'dms_entry_effects_' + str(curr_site) + '.png', dpi = 300)
    plt.close()