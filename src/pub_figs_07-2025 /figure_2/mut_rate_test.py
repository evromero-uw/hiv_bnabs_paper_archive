import sys
import json
sys.path.append('../../../bin/')
sys.path.append('../../../bin/wrappers/')
sys.path.append('../../../data/clyde_westfall_2024_final/10-1074/')
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import binomtest

import data_util as du
import dataset_metadata

#The code in this file produces figure 2 panels e-g

params = {'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

#In this file I am going to try replicating megan's analysis
#Where she calculates the ancestral mutation rates for each amino acid
mut_rate_path = '../../../data/zanini2017/mutation_rates.csv'



outDir = '../../../results/pub_figs_07-2025/figure_2/'

#The path to the data for 10-1074
inDir_1074 = '../../../data/clyde_westfall_2024_final/10-1074/'
par_list_1074 = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD4K', '1HD5K', '1HD6K',
                 '1HD7K', '1HD10K', '1HD11K']
RES_1074_HXB2_NT = [973, 994, 1000]
SUSCEPTIBLE_AA = {325: ['D', 'N'], 332: ['N'], 334: ['S', 'T']}

WEEK_4_PARS = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD4K', '1HD5K', '1HD7K',
               '1HD10K', '1HD11K']
WEEK_8_PARS = ['1HD6K']


############################### Helper Functions ##############################

def enumerate_single_muts(anc_codon, mut_rates_df, site):
    """ Enumerate all single mutations for a given codon
    """
    nucleotides = ['A', 'C', 'G', 'T']
    #Store all single mutations and their rates
    all_single_muts = []

    #Mutate each position
    for pos in range(3):
        for nt in nucleotides:
            if nt != anc_codon[pos]:
                new_codon = anc_codon[:pos] + nt + anc_codon[pos+1:]

                #Only add it if it's a valid codon and it's not the ancestral AA
                check_AA = du.translate_seq(new_codon, error = False)
                anc_AA = du.translate_seq(anc_codon, error = False)
                if check_AA and check_AA[0] != anc_AA[0]:
                    if check_AA[0] in SUSCEPTIBLE_AA[site]:
                        continue
                    #If it is a valid codon, get the mutation rate
                    mut_string = f"{anc_codon[pos]}>{nt}".lower()
                    mut_rate = mut_rates_df.loc[
                        mut_rates_df['nucmut'] == mut_string,'rate'].values[0]
                    
                    #Add everything to the dataframe
                    all_single_muts.append({'codon': new_codon, 
                                            'rate': mut_rate, 
                                            'anc_codon': anc_codon,
                                            'AA' : check_AA[0]})

    
    #Return a dataframe with each possible codon mutation and its rate
    all_single_muts = pd.DataFrame(all_single_muts)                
    all_single_muts['proportion'] = all_single_muts['rate'] / all_single_muts['rate'].sum()

    return all_single_muts

###############################################################################
mut_rate_df = pd.read_csv(mut_rate_path)
time1_all_pars = []
time2_all_pars = []
expected_proportions = []
observed_proportions = []

day0_resistance = []

#We need to loop through every participant and get the ancestral codons and week
# 4/8 codon breakdowns

for curr_par in par_list_1074:
    # The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir_1074 + curr_par + '/885_' + curr_par  + \
                    '_NT_filtered.fasta'
    hxb2_res_positions = dataset_metadata.RESISTANCE_POS_NT_HXB2
    
    # Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    
    # Slice the array to get the codons
    res_arr_dict = {}
    for res_start in RES_1074_HXB2_NT:
        res_start_arr = hxb2_nuc_coords_hta[res_start]
        res_arr_dict[res_start_arr] = seqArr[:, res_start_arr:res_start_arr + 3]
    
    # Now I want to make a dataframe of every sequence and its codons
    seq_indices = np.array(seq_info_df['seq_index'].values).astype(int)

    for curr_site in res_arr_dict.keys():
        curr_site_arr = res_arr_dict[curr_site]
        curr_site_codons = [''.join(curr_site_arr[i]) for i in range(len(curr_site_arr))]
        curr_site_codons = np.array(curr_site_codons)
        curr_site_codons = curr_site_codons[seq_indices]

        #Put it in a dataframe column of seq info df
        curr_site_hxb2 = hxb2_nuc_coords_ath[curr_site]
        curr_site_hxb2 = int(np.ceil(curr_site_hxb2/3))
        seq_info_df['codon_' + str(curr_site_hxb2)] = curr_site_codons

        # Now I want to get the amino acids for each codon
        seq_info_df['AA_' + str(curr_site_hxb2)] = [du.translate_seq(seq, error = False)[0] for seq in\
                                  seq_info_df['codon_' + str(curr_site_hxb2)].values]
    
    seq_info_df['participant'] = curr_par
    
    # Next I need to filter to get only the relevant timepoints
    time1_data = seq_info_df[seq_info_df['time_label'] == 'D0'].copy()
    if curr_par in WEEK_4_PARS:
        time2_data = seq_info_df[seq_info_df['time_label'] == 'W4'].copy()
    elif curr_par in WEEK_8_PARS:
        time2_data = seq_info_df[seq_info_df['time_label'] == 'W8'].copy()

    # I want to unstack my dataframes so that now each site observation is a row
    time1_data = pd.wide_to_long(time1_data,
                                    stubnames= ['codon', 'AA'],
                                    i = ['orig_name', 'seq_index', 'time_label', 
                                        'participant'],
                                    j = 'site',
                                    sep = '_')
    time2_data = pd.wide_to_long(time2_data,
                                    stubnames= ['codon', 'AA'],
                                    i = ['orig_name', 'seq_index', 'time_label', 
                                        'participant'],
                                    j = 'site',
                                    sep = '_')

    #Now I want to get the ancestral codons for each site
    time1_data = time1_data.groupby(['participant', 'site', 'codon', 'AA']
                                            ).agg(func= 'size').reset_index()

    time1_data = time1_data.rename(columns = {0: 'count'})
    time1_data = time1_data[time1_data['codon'] != '---']

    # Next, for each site, I am going to get the expected proportion of each 
    # codon I just need to get the proportions expected from each ancestral
    # codon then average them across codons

    for curr_site in time1_data['site'].unique():
        all_codon_results = []
        curr_site_data = time1_data[time1_data['site'] == curr_site]

        #Mark the day 0 resistance
        curr_resistance = curr_site_data[~curr_site_data['AA'].isin(SUSCEPTIBLE_AA[curr_site])]
        if len(curr_resistance) > 0:
            day0_resistance.append(curr_resistance)

        
        #Filter out any codons that are resistant at day 0
        curr_site_data = curr_site_data[curr_site_data['AA'].isin(SUSCEPTIBLE_AA[curr_site])]

        

        #Get the frequencies of each codon at day 0
        codon_freqs = curr_site_data['count'] / curr_site_data['count'].sum()
        codon_freqs = dict(zip(curr_site_data['codon'], codon_freqs))

        #Get the mutation rates of each codon
        for curr_codon in curr_site_data['codon'].unique():
            all_single_muts = enumerate_single_muts(curr_codon, mut_rate_df, curr_site)
            
            #Multiply the proportions by the codon frequencies to weight them
            all_single_muts['exp_prop_weighted'] = all_single_muts['proportion'] *\
                                                     codon_freqs[curr_codon]
            all_codon_results.append(all_single_muts)

        all_codon_results = pd.concat(all_codon_results)
        all_codon_results.drop(columns = ['rate', 'proportion', 
                                          'anc_codon', 'codon'], inplace = True)

        all_codon_results = all_codon_results.groupby('AA').agg(func = 'sum').reset_index()
        all_codon_results['site'] = curr_site
        all_codon_results['participant'] = curr_par
        expected_proportions.append(all_codon_results)

    #Next, I want to get the actual observed proportions at the next timepoint
    time2_data = time2_data[time2_data['codon'] != '---']
    time2_data['count'] = 1
    time2_data = time2_data.groupby(['participant', 'site', 'AA']
                                            ).agg(func= 'sum').reset_index()
    time2_data = time2_data.rename(columns = {0: 'count'})

    

    #Get only the susceptible amino acids and calculate the proportions
    time2_data['susceptible'] = [SUSCEPTIBLE_AA[site] for site in time2_data['site']]
    for i, row in time2_data.iterrows():
        if row['AA'] in row['susceptible']:
            time2_data.drop(i, inplace = True)

    time2_data['obs_prop'] = time2_data.groupby(['participant', 'site'])['count'].transform(lambda x: x / x.sum())
    time2_data = time2_data.drop(columns = ['codon', 'susceptible'])

    
    observed_proportions.append(time2_data)

expected_proportions = pd.concat(expected_proportions)
observed_proportions = pd.concat(observed_proportions)

#Combine the expected and observed proportions
all_proportions = pd.merge(expected_proportions, observed_proportions,
                           how='outer', on = ['participant', 'site', 'AA'])
all_proportions = all_proportions.fillna(0)
for curr_par in all_proportions['participant'].unique():
    curr_par_data = all_proportions[all_proportions['participant'] == curr_par]

#Save the observed and expected proportions
all_proportions.to_csv(outDir + 'obs_vs_predicted_props.csv', index = False)

day0_resistance = pd.concat(day0_resistance)
day0_resistance.to_csv(outDir + 'day0_resistance.csv', index = False)


#################################### Test for enrichment ##################################

# I need to make a dataframe of the proportions vs the expected proportions for testing

testing_props_all = all_proportions.copy()
testing_props_all = testing_props_all.rename(columns = {'exp_prop_weighted': 'prop_expected',
                                                        'count': 'num_observed'})

#For each participant I want to calculate how many sequences were contributed
testing_props_all['participant_total'] = testing_props_all.groupby(['participant', 'site'])['num_observed'].transform(lambda x: x.sum())
testing_props_all['total_seqs'] = testing_props_all.groupby(['site'])['num_observed'].transform(lambda x: x.sum())
testing_props_all['participant_seq_prop'] = testing_props_all['participant_total'] / testing_props_all['total_seqs']

#Weight all of the expected proportions by the participant's contribution
testing_props_all['prop_expected'] = testing_props_all['prop_expected'] *\
                                      testing_props_all['participant_seq_prop']


testing_props_all = testing_props_all.drop(columns = ['obs_prop', 'participant',
                                                      'total_seqs'])
testing_props_all = testing_props_all[testing_props_all['site'] != 325]

#Aggregate over all participants
testing_props_all = testing_props_all.groupby(['site', 'AA']).agg(func = 'sum').reset_index()
testing_props_all['num_total'] = testing_props_all.groupby(['site'])['num_observed'].transform(lambda x: x.sum())


binom_test_results = []
for i, curr_row in testing_props_all.iterrows():
    test_result = binomtest(int(curr_row['num_observed']), int(curr_row['num_total']), curr_row['prop_expected'],
                       alternative = 'greater')
    binom_test_results.append(test_result.pvalue)


testing_props_all['p_val'] = binom_test_results
threshold = 0.05/ len(testing_props_all)
print('Threshold for significance:', threshold)
testing_props_all['significant'] = testing_props_all['p_val'] < threshold
testing_props_all['significant'] = testing_props_all['significant'].astype(bool)
testing_props_all['significant'] = testing_props_all['significant'].replace({True: 'Yes', False: 'No'})


testing_props_all.to_csv(outDir + 'mutation_rates_test_results.csv')
significance_df = testing_props_all[testing_props_all['significant'] == 'Yes']


##################### Test for enrichment excluding resistant amino acids #####
# Now I want to run the same test but exclude the resistant amino acids
testing_props_excl = all_proportions.copy()
testing_props_excl = testing_props_excl.rename(columns = {'exp_prop_weighted': 'prop_expected',
                                                        'count': 'num_observed'})

testing_props_excl = testing_props_excl[testing_props_excl['site'] != 325]

#Remove day 0 resistance make a column labeling it and then filter based
# on the column
day0_resistance['par_site_AA'] = [x[0] + '_' + str(x[1]) + '_' + \
                                  x[2] for x in zip(day0_resistance['participant'],
                                                    day0_resistance['site'], day0_resistance['AA'])]
testing_props_excl['par_site_AA'] = [x[0] + '_' + str(x[1]) + '_' +\
                                      x[2] for x in zip(testing_props_excl['participant'],
                                                        testing_props_excl['site'], testing_props_excl['AA'])]

testing_props_excl = testing_props_excl[~testing_props_excl['par_site_AA'].isin(day0_resistance['par_site_AA'])]

#For each participant I want to calculate how many sequences were contributed
testing_props_excl['participant_total'] = testing_props_excl.groupby(['participant', 'site'])['num_observed'].transform(lambda x: x.sum())
testing_props_excl['total_seqs'] = testing_props_excl.groupby(['site'])['num_observed'].transform(lambda x: x.sum())

testing_props_excl['participant_seq_prop'] = testing_props_excl['participant_total'] / testing_props_excl['total_seqs']
testing_props_excl['prop_expected'] = testing_props_excl['prop_expected'] *\
                                      testing_props_excl['participant_seq_prop']

testing_props_excl = testing_props_excl.drop(columns = ['obs_prop', 'participant',
                                                        'par_site_AA', 'total_seqs'])

#Now aggregate over all participants
testing_props_excl = testing_props_excl.groupby(['site', 'AA']).agg(func = 'sum').reset_index()
testing_props_excl['num_total'] = testing_props_excl.groupby(['site'])['num_observed'].transform(lambda x: x.sum())

binom_test_results = []
for i, curr_row in testing_props_excl.iterrows():
    test_result = binomtest(int(curr_row['num_observed']), int(curr_row['num_total']), curr_row['prop_expected'],
                       alternative = 'greater')
    binom_test_results.append(test_result.pvalue)

testing_props_excl['p_val'] = binom_test_results
threshold = 0.05/ len(testing_props_excl)
testing_props_excl['significant'] = testing_props_excl['p_val'] < threshold
testing_props_excl['significant'] = testing_props_excl['significant'].astype(bool)
testing_props_excl['significant'] = testing_props_excl['significant'].replace({True: 'Yes', False: 'No'})

testing_props_excl.to_csv(outDir + 'mutation_rates_test_results_no_preexist.csv')



############################# Plot the proportions vs expected ################

all_proportions['obs_exp_ratio'] = all_proportions['obs_prop'] / all_proportions['exp_prop_weighted']
all_proportions['obs_exp_ratio'] = np.log(all_proportions['obs_exp_ratio'])

fig, axs = plt.subplots(2, 1,
                        figsize = (2.5, 2.5))
axs = axs.flatten()
ax_334 = 1
ax_332 = 0

for curr_site in all_proportions['site'].unique():
    #Skip site 325 since we have less observations
    if curr_site == 325:
        continue
    elif curr_site == 334:
        curr_ax = axs[ax_334]
    else:
        curr_ax = axs[ax_332]

    curr_site_data = all_proportions[all_proportions['site'] == curr_site]
    curr_site_data = curr_site_data.sort_values(['AA'])
    zero_data = curr_site_data[curr_site_data['obs_exp_ratio'] == -np.inf]
    curr_site_data = curr_site_data[curr_site_data['obs_exp_ratio'] != -np.inf]

    # Filter out 1HD6K for site 334 since they had a different codon
    if curr_site == 334:
        curr_site_data = curr_site_data[curr_site_data['participant'] != '1HD6K']
        zero_data = zero_data[zero_data['participant'] != '1HD6K']
    
    curr_data_order = curr_site_data['AA'].unique()
    curr_data_order = sorted(curr_data_order)
    curr_site_data['AA'] = pd.Categorical(curr_site_data['AA'], categories = curr_data_order)
    

    sns.boxplot(data = curr_site_data, x = 'AA', y = 'obs_exp_ratio',
                fliersize = 0, fill= False, color = 'gray',
                width = 0.5, linewidth = 1, ax = curr_ax)
    sns.stripplot(data = curr_site_data, x = 'AA', y = 'obs_exp_ratio',
                color='black', ax = curr_ax, size = 2)

    sns.stripplot(x = zero_data['AA'], y = min(curr_site_data['obs_exp_ratio'])- 0.5,
                color= 'black', 
                jitter = 0.3, marker = 'x', linewidth= 0.5, size = 1.5, ax = curr_ax)

    curr_ax.annotate('Not observed', xy= (0.65, 0.07), xycoords = 'axes fraction',
            fontname = 'Arial', fontsize = 6, color = 'lightgrey', ha = 'left')
    
    y_max = max(curr_site_data['obs_exp_ratio']) + 0.6
    
    
    #Annotate the significant differences
    curr_sigs = significance_df[significance_df['site'] == curr_site]
    print(significance_df)
    print('ANNOTATING SITE', curr_site)
    if len(curr_sigs) > 0:
        print('curr_sigs', curr_sigs)
        for i, curr_row in curr_sigs.iterrows():
            #find the correct x position
            curr_AA = curr_row['AA']
            all_AAs = curr_site_data['AA'].unique()
            all_AAs = sorted(all_AAs)
            curr_AA_index = all_AAs.index(curr_AA)
            curr_x = curr_AA_index
            if curr_site == 334:
                curr_y = y_max - 0.6
            else:
                curr_y = y_max - 0.9
            print('curr_y', curr_y)
            curr_ax.annotate('*', xy = (curr_x, curr_y),
                        xycoords= 'data', fontsize = 10, color = 'red', ha = 'center')
            


    curr_ax.set_ylim(min(curr_site_data['obs_exp_ratio'])- 0.7, y_max)
    curr_ax.axhline(0, color = 'black', linestyle = '--', linewidth = 0.5)
    curr_ax.axhline(min(curr_site_data['obs_exp_ratio'])- 0.5, color = 'lightgrey',
                    linestyle = '--', linewidth = 0.5)
    #change the title padding
    curr_ax.set_title('Site ' + str(curr_site), y = 0.9)
    curr_ax.set_ylabel('Log(observed/expected)')
    curr_ax.set_xlabel('Amino acid')
    if curr_site == 332:
        curr_ax.set_xlabel('')


    [i.set_linewidth(0.5) for i in curr_ax.spines.values()]

plt.subplots_adjust(left = 0.2, right = 0.98, hspace = 0.5, top= 0.95, bottom = 0.15)
plt.savefig(outDir + 'mutation_rates_' + str(curr_site) + '_obs_exp_ratio.png', dpi = 300)
plt.close()
