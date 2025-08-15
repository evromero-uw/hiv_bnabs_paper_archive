import os
import sys
sys.path.append('../../../bin/')
sys.path.append('../../../bin/wrappers/')
sys.path.append('../../../data/clyde_westfall_2024_final/3BNC117/')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

import data_util as du
import dataset_metadata
import data_processing_fig4 as fig4
from par_data_class import Pardata


from matplotlib import rcParams
from matplotlib.transforms import Affine2D
from matplotlib.markers import MarkerStyle

#The code in this file makes figure 4 panel a

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'
params = {'figure.figsize':(3.85, 1.86), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
MARKERSIZE = 4
MARKEREDGEWIDTH = 0.5

# For space purposes we will only plot between hxb2 250 and 500
PLOT_MIN_HXB2 = 250
PLOT_MAX_HXB2 = 500

rcParams.update(params)



inDir_3BNC = '../../../data/clyde_westfall_2024_final/3BNC117/'
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24']

#Here is the input folder
inFolder = 'figure_4/'
mpl_data_dir = '../../../results/pub_figs_07-2025/figure_4/MPL_analysis/'
enc_data_dir = '../../../results/pub_figs_07-2025/figure_4/Multi_encoded/'
ann_data_dir = '../../../results/pub_figs_07-2025/figure_4/data/'
supp_table_dir = '../../../results/pub_figs_07-2025/supplemental_tables/'

outFolder = '../../../results/pub_figs_07-2025/figure_4/'

dms_path = ann_data_dir + '032724_DMS_mean_sig_sites_cleaned.csv'
enc_path = enc_data_dir + 'multiply_encoded_3BNC117_filtered_input.csv'
add_path = ann_data_dir + 'additional_variation.csv'
var_path = mpl_data_dir + 'mpl_support_sites.csv'
adb_path = ann_data_dir + '3BNC117_antibody_database_rules.csv'

annotation_path = ann_data_dir + '3BNC117sigs.csv'

#combine the dms data, lanl sig sites, and antibody database sites which we
#will use for annotations
annotation_df = pd.read_csv(annotation_path)
adb_df = pd.read_csv(adb_path)
dms_df = fig4.format_dms_data_annotation(dms_path)
annotation_df = pd.concat([annotation_df, adb_df, dms_df], ignore_index=True)

#These are all of the regions of interest we will check
HXB2_RES_POS = dataset_metadata.RESISTANCE_POS_AA_HXB2 


region_coord_dict = {'Loop_D': (275, 283),
                    #  'V1': (131, 157),
                    #  'V2': (158, 196),
                    #  'V3': (296, 331),
                    #  'V4': (385, 418),
                     'V5': (460, 470),
                     'CD4_Binding_Loop': (364, 374),}

# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will only use all alleles at each site
MULTI_SEG = True
###############################################################################
#Plot parameters
marker_dict = {'shared_variability': 'o', 'frequency_increase': 'o',
               'multi_encodings': '|', 'additional_sites': 'o'}
marker_rotate_dict = {'multi_encodings': 45}
method_color_dict = {'shared_variability': 'red', 
                     'frequency_increase': 'gray',
                     'multi_encodings': 'black',
                     'additional_sites': 'mediumpurple'}
method_list = ['shared_variability', 'frequency_increase', 'multi_encodings',
               'additional_sites']

###############################################################################
# First I need to load the data and put it all in one big dataframe

collated_df = fig4.combine_all_data(enc_path = enc_path, 
                                    var_path = var_path,
                                    add_path = add_path)

collated_df = collated_df[~collated_df['method'].isin(['shared_variability'])]

# Next, I need to get all of the hxb2 positions across all individuals
hxb2_pos_unique = set()
for curr_par in par_list_3BNC:
    print(curr_par)

    # Here is the path to the data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + \
                '_NT.translated.fasta'
    
    # First, I need to load the data
    participant_dat = Pardata(inFile, 'clyde2024', curr_par)
    participant_dat.load_data_3BNC117(ALLELE_FREQ_THRESH, MULTI_SEG)

    # Next, I will get a couple of items out of the dataset
    hxb2_nuc_coords_hta = participant_dat.hxb2_nuc_coords_hta
    hxb2_nuc_coords_ath = participant_dat.hxb2_nuc_coords_ath

    hxb2_pos_unique.update(hxb2_nuc_coords_hta.keys())

    # Now I will get the sites which varied at timepoints with more than 
    # 10 sequences
    

    

#Last, I need to generate my dictionary which defines the mapping between
#participant coordinates and the shared heatmap coordinates
hxb2_pos_unique = list(hxb2_pos_unique)
hxb2_int_unique = range(0, len(hxb2_pos_unique))
hxb2_pos_unique.sort()
hxb2_pos_to_int_dict= dict(zip(hxb2_pos_unique, hxb2_int_unique))
hxb2_int_to_pos_dict = dict(zip(hxb2_int_unique, hxb2_pos_unique))

#Now I will map the hxb2 coordinates to integers
plot_df = []
for i, curr_par in enumerate(par_list_3BNC):
    #Get the data for the current participant
    curr_outs = collated_df[collated_df['participant'] == curr_par].copy()
    curr_out_sites = np.unique(curr_outs['hxb2_coord_AA'].values)
    curr_outs_int = []
    #Map the hxb2 coords to integer positions
    for curr_site in curr_out_sites:
        if curr_site in hxb2_int_to_pos_dict:
            curr_outs_int.append(hxb2_pos_to_int_dict[curr_site])
        else:
            matching_key = min(hxb2_int_to_pos_dict, key=lambda x: abs(x - curr_site))
            curr_outs_int.append(hxb2_pos_to_int_dict[matching_key])
    curr_outs_int = dict(zip(curr_out_sites, curr_outs_int))
    curr_outs['integer_site'] = curr_outs['hxb2_coord_AA'].map(curr_outs_int)
    curr_outs['par_y'] = np.full(len(curr_outs), i)
    plot_df.append(curr_outs)
plot_df = pd.concat(plot_df, ignore_index=True)


###############################################################################
#Trying to make a second plot where sites must be within 5AA of a binding site
#or previously annotated.
# Now I will make a line plot where I summarize all of the positions I'm 
# finding
fig, ax = plt.subplots(1, 1)


#First I need to filter the plot dataframe
#I'll start by getting a set of hxb2 resistance sites
int_res_set = set()
for curr_res in HXB2_RES_POS:
    start = hxb2_pos_to_int_dict[curr_res[0]]
    end = hxb2_pos_to_int_dict[curr_res[1]]
    int_res_set.update(range(start, end +1))

#I'll also get the previously annotated sites
annotation_df = annotation_df[annotation_df['Reference'] != 'Zhou2013a,Caskey2015,Scheid2016']
annotated_sites = annotation_df['HXB2pos'].values
annotated_sites_int = list(map(lambda x: hxb2_pos_to_int_dict[x], annotated_sites))
annotated_sites_int = set(annotated_sites_int)

all_sites_to_keep = int_res_set.union(annotated_sites_int)

plot_df.to_csv(supp_table_dir + '3BNC117_all_data.csv', index = False)

#Format the data from the supplemental table
###########################################################################
plot_df_save = plot_df.copy()
plot_df_save.sort_values(by=['hxb2_coord_AA', 'participant'], inplace=True)
plot_df_save['Putative Escape'] = plot_df_save['integer_site'].isin(all_sites_to_keep)
plot_df_save.drop(columns=['integer_site', 'par_y'], inplace=True)
plot_df_save.to_csv(supp_table_dir + '3BNC117_all_data.csv', index=False)
#Convert the hxb2 coordinates to strings for the supplemental table
new_AA_labels = []
for curr_label in plot_df_save['hxb2_coord_AA']:
    if curr_label % 1 == 0:
        new_AA_labels.append(str(int(curr_label)))
    else:
        remainder = curr_label % 1
        remainder = int(np.round(remainder * 1000))  # Convert to integer for display
        alphabet_remainder = ALPHABET[remainder]
        new_AA_labels.append(f"{int(curr_label)}{alphabet_remainder}")
plot_df_save['hxb2_coord_AA'] = new_AA_labels

plot_df_save['method'] = plot_df['method'].apply(lambda x: 'MPL' if x == 'frequency_increase' else x)
plot_df_save = plot_df_save[plot_df_save.method != 'additional_sites']
plot_df_save.rename(columns={'participant': 'Participant',
                             'method': 'Method'}, inplace=True)
plot_df_save.to_csv(supp_table_dir + '3BNC117_all_data_pretty.csv', index = False)
############################################################################


plot_df_filtered = plot_df[plot_df['integer_site'].isin(all_sites_to_keep)]

#Save the filtered dataframe
plot_df_filtered.to_csv(outFolder + '3BNC117_filtered_data.csv', index = False)
plot_df_filtered_save = plot_df_filtered.copy()
plot_df_filtered_save['method'] = plot_df_filtered['method'].apply(lambda x: 'individual_variability' if x == 'frequency_increase' else x)
plot_df_filtered_save.to_csv(supp_table_dir + '3BNC117_filtered_data.csv', index = False)

for i, curr_par in enumerate(par_list_3BNC):
    curr_outs = plot_df[plot_df['participant'] == curr_par].copy()

    #Represent each sequence in the window as a line
    plt.plot(hxb2_int_to_pos_dict.keys(), np.full(len(hxb2_int_to_pos_dict), i), color='black', linewidth=0.5, zorder = 2)


for i, curr_par in enumerate(par_list_3BNC):
    #Now, plot the variable loops and resistance positions
    for curr_res in HXB2_RES_POS:
        start = hxb2_pos_to_int_dict[curr_res[0]]
        end = hxb2_pos_to_int_dict[curr_res[1]]
        if i != 0:
            plt.fill_between(range(start, end), i-1, i+1, color='black', alpha=0.1, zorder=1)
        else:
            plt.fill_between(range(start, end), i-0.5, i+1, color='black', alpha=0.1, zorder=1)

#Reset the x-axis ticklabels to hxb2 coordinates
ax.set_ylim(-2, len(par_list_3BNC))
ax.set_yticks(range(len(par_list_3BNC)))
yticks = ax.get_yticks()
new_yticks = []
for curr_tick in yticks:
    if curr_tick in range(len(par_list_3BNC)):
        new_yticks.append(par_list_3BNC[int(curr_tick)])
    else:
        new_yticks.append('')
ax.set_yticklabels(new_yticks)

ax.set_xlim(hxb2_pos_to_int_dict[PLOT_MIN_HXB2], hxb2_pos_to_int_dict[PLOT_MAX_HXB2])
xticks = ax.get_xticks()
new_xticks = []
for curr_tick in xticks:
    if curr_tick in hxb2_int_to_pos_dict.keys():
        new_xticks.append(int(hxb2_int_to_pos_dict[curr_tick]))
    else:
        new_xticks.append('')
ax.set_xticklabels(new_xticks)


for step, curr_method in enumerate(method_list):
    print(curr_method)
    for i, curr_par in enumerate(par_list_3BNC):
        curr_outs = plot_df[plot_df['participant'] == curr_par].copy()

        if curr_method not in curr_outs['method'].values:
            continue
        curr_method_df = curr_outs[curr_outs['method'] == curr_method]
        curr_method = curr_method.strip()
        curr_marker = marker_dict[curr_method]
        curr_color = method_color_dict[curr_method]


        if curr_method in marker_rotate_dict:
            curr_rotate = marker_rotate_dict[curr_method]
            t = Affine2D().rotate_deg(curr_rotate)
            plt.plot(curr_method_df['integer_site'], curr_method_df['par_y'], zorder = 4, 
                        marker = MarkerStyle(curr_marker, 'left', t), linestyle = 'none', color = curr_color, markersize=MARKERSIZE,
                        markeredgecolor = 'black', label = curr_method, markeredgewidth = MARKEREDGEWIDTH)

        else:
            print(f"Plotting {curr_method} with marker {curr_marker} and color {curr_color}")
            plt.plot(curr_method_df['integer_site'], curr_method_df['par_y'], zorder = 3,
                        marker = curr_marker, linestyle='None', color = curr_color, markersize=MARKERSIZE,
                        markeredgecolor = 'none', markeredgewidth = 0, label = curr_method)



#Now get the sites that are annotated
for i, curr_par in enumerate(par_list_3BNC):
    curr_outs = plot_df[plot_df['participant'] == curr_par]
    all_sites = set(curr_outs['integer_site'].values)
    #Filter out caskey annotations so we are not doubledipping
    annotation_df = annotation_df[annotation_df['Reference'] != 'Zhou2013a,Caskey2015,Scheid2016']
    annotated_sites = []


    for curr_site in all_sites:
        if curr_site in all_sites_to_keep:
            annotated_sites.append(curr_site)
    annotated_sites_int = list(annotated_sites)
    plt.plot(annotated_sites_int, np.full(len(annotated_sites_int), i),  alpha = 1, markerfacecolor = 'none', markeredgecolor='black', markersize=MARKERSIZE, zorder = 3, marker = 'o',
                label = 'Resistance Sites', linestyle = 'none', markeredgewidth = MARKEREDGEWIDTH)


#Label the regions of interest
for region, (start, end) in region_coord_dict.items():
    if start < PLOT_MIN_HXB2 or end > PLOT_MAX_HXB2:
        print(f"Skipping region {region} as it is outside the inset range.")
        continue
    print(f"Plotting region: {region} from {start} to {end}")
    plot_start = hxb2_pos_to_int_dict[start]
    plot_end = hxb2_pos_to_int_dict[end]
    plt.fill_between(range(plot_start, plot_end), -1.5, -0.5, color='lightskyblue', zorder=4)


for curr_res in HXB2_RES_POS:
    start = hxb2_pos_to_int_dict[curr_res[0]]
    end = hxb2_pos_to_int_dict[curr_res[1]]
    plt.fill_between(range(start, end), -0.5, 0, color='black', alpha=0.1, zorder=3)



#make the legend
indiv_point = mpl.lines.Line2D([0], [0], marker= marker_dict['frequency_increase'],
                                color= 'none', label='Identified by MPL',
                              markerfacecolor= method_color_dict['frequency_increase'], 
                              markersize=MARKERSIZE, markeredgecolor='none')
enc_point = mpl.lines.Line2D([0], [0], marker=MarkerStyle(marker_dict['multi_encodings'], 'left', 
                                Affine2D().rotate_deg(marker_rotate_dict['multi_encodings'])), color='none', label='Multiply encoded variability',
                              markersize=MARKERSIZE, markeredgecolor='black', markeredgewidth=MARKEREDGEWIDTH)
annotation_point = mpl.lines.Line2D([0], [0], marker='o', color='none', markeredgecolor='black', label='Identified site is escape associated',
                              markerfacecolor='none', markersize=MARKERSIZE, markeredgewidth= MARKEREDGEWIDTH)
additional_point = mpl.lines.Line2D([0], [0], marker='o', color='none', markeredgecolor='none', label='Additional variation',
                              markerfacecolor=method_color_dict['additional_sites'], markersize=MARKERSIZE, markeredgewidth=MARKEREDGEWIDTH)

region_shade = mpl.patches.Patch(color='black', alpha=0.2, label='3BNC117 contact region')



legend = plt.legend(handles=[indiv_point, enc_point, additional_point, annotation_point, region_shade],
                    bbox_to_anchor=(0.46, 1), loc='lower center', ncol=3,
                    handletextpad = 0.3, frameon=False, columnspacing=0.5)
legend.get_frame().set_alpha(None)

plt.subplots_adjust(bottom=0.2, top=0.8, left=0.1, right=0.95)
plt.savefig(outFolder + 'summary_plot_final.png', dpi = 300)
plt.close()