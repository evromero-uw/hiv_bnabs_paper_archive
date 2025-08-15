import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import logomaker
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
from par_data_class import Pardata


params = {'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}

rcParams.update(params)

#The code in the file makes supplemental logo plots of the 3BNC117 contact
#regions which are comparable to the caskey et al 2015 logo plots.

CASKEY_HXB2_POS = [(270, 285), (360, 371), (455, 485)]
REGION_LABELS = ['Loop D', 'CD4 binding site', 'V5']

inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'
outDir_3BNC = '../../../../results/pub_figs_07-2025/supplemental_figs/logo_comparisons/'
par_list_3BNC = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']
time_filter_out = ['Rebound', 'screen', 'pre', 'Nadir', 'HXB2', 'W24']
time_list = ['D0', 'W4']
label_list = ['SGS', 'SMRT-UMI']
type_time = ['SGS D0', 'SMRT-UMI D0', 'SGS W4', 'SMRT-UMI W4']

#Get the putative escape loci
putative_escape_dir = '../../../../results/pub_figs_07-2025/figure_4/3BNC117_filtered_data.csv'
putative_escape_df = pd.read_csv(putative_escape_dir, index_col=None)


ALPHABET = '-abcdefghijklmnopqrstuvwxyz'


# We will only take the sites with a minor allele frequency above this 
# threshold
ALLELE_FREQ_THRESH = 0

# We will only use all alleles at each site
MULTI_SEG = True

###############################################################################

# Loop through each participant and create a logo plot
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

    # I need to get the sgs and porpid sequences separately
    seq_info_df['smrt_umi'] = ['KX' not in x for x in seq_info_df['orig_name']]
    seq_info_df['sample_type'] = ['SMRT-UMI' if x else 'SGS' for x in seq_info_df['smrt_umi']]

    # I will now format each of the datasets for making the logo plots
    # Create a dictionary to hold the dataframes for each time point, sample type, and region
    # The keys will be tuples of (time_label, sample_type, region_start, region_end)
    logo_df_dict = {}

    # Make a dictionary of positions to highlight, each key will be a tuple of (time_label, sample_type, start_pos, end_pos)
    highlight_positions = {}

    # Record the sample sizes for each time point and sample type
    sample_sizes = {}

    for time_label in time_list:
        for sample_type in label_list:
            # Filter the data for the current time point and sample type
            curr_df = seq_info_df[(seq_info_df['time_label'] == time_label) & 
                                  (seq_info_df['sample_type'] == sample_type)]
            if curr_df.empty:
                continue

            # Get the sequences for the current time point and sample type
            curr_seqArr = seqArr[curr_df['seq_index'], :]

            for curr_region in CASKEY_HXB2_POS:
                # Get the start and end positions for the current region
                start_pos, end_pos = curr_region
                start_pos_arr = hxb2_nuc_coords_hta[start_pos]
                end_pos_arr = hxb2_nuc_coords_hta[end_pos]
                
                # Get the sequences for the current region
                curr_region_seqArr = curr_seqArr[:, start_pos_arr:end_pos_arr + 1]

                # Get the unique characters in this region
                all_characters = curr_region_seqArr.flatten()
                all_characters = np.unique(all_characters)
                all_characters = [char if char != '-' else ' ' for char in all_characters]

                # Now I need to make a dataframe for each region where the rows are
                # the positions and the columns are the characters.
                # The entries will be the character heights.

                #iterate over the sequence array columns (positions) and add each position
                # as a row in the dataframe
                aa_freq_matrix = np.zeros((curr_region_seqArr.shape[1], len(all_characters)))
                for i in range(curr_region_seqArr.shape[1]):
                    curr_char_counts = pd.Series(curr_region_seqArr[:, i]).value_counts()
                    
                    #normalize the counts to get the character heights
                    curr_char_counts = curr_char_counts / curr_char_counts.sum()

                    for j, char in enumerate(all_characters):
                        aa_freq_matrix[i, j] = curr_char_counts.get(char, 0)
                
                # Create a DataFrame for the character frequencies
                aa_freq_df = pd.DataFrame(aa_freq_matrix,
                                          columns=all_characters,
                                          index=np.arange(start_pos_arr, end_pos_arr + 1))
                

                # Store the DataFrame in the dictionary
                logo_df_dict[(time_label, sample_type, start_pos, end_pos)] = aa_freq_df
            
            # Record the sample size for this time point and sample type
            sample_sizes[(time_label, sample_type)] = len(curr_df)
            
        #Now I need to find positions with new AAs after sequencing
        for curr_region in CASKEY_HXB2_POS:
            start_pos, end_pos = curr_region
            if (time_label, 'SGS', start_pos, end_pos) in logo_df_dict and \
               (time_label, 'SMRT-UMI', start_pos, end_pos) in logo_df_dict:
                sgs_df = logo_df_dict[(time_label, 'SGS', start_pos, end_pos)]
                smrt_umi_df = logo_df_dict[(time_label, 'SMRT-UMI', start_pos, end_pos)]

                # Find positions with new AAs in SMRT UMI compared to SGS
                # Loop through each position in the SMRT UMI DataFrame
                new_aa_positions_smrt = []
                new_aa_positions_sgs = []
                for pos in smrt_umi_df.index:
                    # Get the AAs at this position in both DataFrames
                    sgs_aa = sgs_df.loc[pos]
                    smrt_umi_aa = smrt_umi_df.loc[pos]

                    # Get all of the column names with non-zero values
                    sgs_aa_nonzero = sgs_aa[sgs_aa > 0].index.tolist()
                    smrt_umi_aa_nonzero = smrt_umi_aa[smrt_umi_aa > 0].index.tolist()

                    # Check if there are any AAs in SMRT UMI that are not in SGS
                    if any(aa not in sgs_aa_nonzero for aa in smrt_umi_aa_nonzero):
                        #Get the new AAs in SMRT UMI
                        new_aa = [aa for aa in smrt_umi_aa_nonzero if aa not in sgs_aa_nonzero]
                        pos_new_freq = smrt_umi_aa[new_aa].sum()
                        if pos_new_freq > 0.05:
                            new_aa_positions_smrt.append(pos)

                    if any(aa not in smrt_umi_aa_nonzero for aa in sgs_aa_nonzero):
                        #Get the new AAs in SGS
                        new_aa = [aa for aa in sgs_aa_nonzero if aa not in smrt_umi_aa_nonzero]
                        pos_new_freq = sgs_aa[new_aa].sum()
                        if pos_new_freq > 0.05:
                            new_aa_positions_sgs.append(pos)

                # Store the new AA positions in the highlight_positions dictionary
                highlight_positions[(time_label, 'SMRT-UMI', start_pos, end_pos)] = new_aa_positions_smrt
                highlight_positions[(time_label, 'SGS', start_pos, end_pos)] = new_aa_positions_sgs

    # Now, I need to figure out which regions have new variation after sequencing



    # Now I will plot the logos for each time point, sample type, and region
    # we have four sample type/ timepoint combos and three regions
    fig, axs = plt.subplots(4, 3, figsize=(7.5, 7.09), sharey=True, width_ratios=[0.75, 0.6, 2])
    untouched_rows = set(range(4))

    for curr_item in logo_df_dict.keys():
        time_label, sample_type, start_pos, end_pos = curr_item
        curr_df = logo_df_dict[curr_item]
        curr_samp_size = sample_sizes[(time_label, sample_type)]

        #get the correct axis for this item
        row_index = type_time.index(f"{sample_type} {time_label}")
        untouched_rows.discard(row_index)
        column_index = CASKEY_HXB2_POS.index((start_pos, end_pos))
        ax = axs[row_index, column_index]

        # Create the logo plot
        logo = logomaker.Logo(curr_df, ax=ax)
        curr_reg_label = REGION_LABELS[column_index]
        ax.set_title(f"{time_label} {sample_type}\n Region: {curr_reg_label}\n # of Sequences: {curr_samp_size}")


        # If there are putative escape mutations, highlight them
        curr_put_escape = putative_escape_df[putative_escape_df['participant'] == curr_par]
        if len(curr_put_escape) > 0:
            put_esc_pos = curr_put_escape['hxb2_coord_AA'].tolist()
            put_esc_pos = [hxb2_nuc_coords_hta[pos] for pos in put_esc_pos]
            for pos in put_esc_pos:
                if pos in curr_df.index:
                    logo.highlight_position(pos, color='lightskyblue', alpha=1)

        # Highlight positions with new AAs
        pos_to_highlight = highlight_positions.get(curr_item, [])
        for pos in pos_to_highlight:
            logo.highlight_position(pos, color='yellow', alpha=0.5)


        #make sure every x value has a tick label
        ax.set_xticks(curr_df.index)

        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels = [hxb2_nuc_coords_ath[int(x)] for x in labels]

        #label any nonint ticks with the 
        new_labels = []
        for curr_label in labels:
            if curr_label % 1 == 0:
                new_labels.append(str(curr_label))
            else:
                remainder = curr_label % 1
                remainder = int(np.round(remainder * 1000))  # Convert to integer for display
                alphabet_remainder = ALPHABET[remainder]
                new_labels.append(f"{int(curr_label)}{alphabet_remainder}")


        ax.set_xticklabels(new_labels)
        #rotate all of the tick labels
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
        if column_index == 0:
            ax.set_ylabel('Frequency')
        if row_index == 3:
            ax.set_xlabel('HXB2 Position')
    
    if untouched_rows:
        # If there are any untouched rows, remove them from the figure
        for row in sorted(untouched_rows, reverse=True):
            axs[row, 0].remove()
            axs[row, 1].remove()
            axs[row, 2].remove()

    # Now make a plot legend
    fig.legend(handles=[
        plt.Line2D([0], [0], color='lightskyblue', lw=4, label='Putative escape site', alpha = 1),
        plt.Line2D([0], [0], color='yellow', lw=4, label='AA specific to sequencing method at given time point', alpha = 0.5),
        plt.Line2D([0], [0], color='#C3E67C', lw=4, label='Both putative escape site and method specific AA', alpha = 1)
    ], loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=3, fontsize=6, frameon=False)

    # Adjust layout and show the plot
    plt.suptitle(f"3BNC117 Contact Regions - {curr_par}", fontsize=8)
    plt.subplots_adjust(hspace=0.7, wspace=0.1)
    plt.savefig(outDir_3BNC + f"{curr_par}_logo_comparisons.png", dpi=300, bbox_inches='tight')
    plt.close()





