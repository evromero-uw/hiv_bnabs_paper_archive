import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/3BNC117/')
import pandas as pd

import data_util as du

#In this file I am going to triage Alison's multiply encoded data,
#add the hxb2 positions, and then save the data to a file

poor_alignment_hand_checked = [139.006]

inDir = '../../../../results/pub_figs_07-2025/figure_4/Multi_encoded/'
inDir_3BNC = '../../../../data/clyde_westfall_2024_final/3BNC117/'

encoded_df = pd.read_csv(inDir + 'multiply_encoded_3BNC117_wide_updated_AFF.tsv', sep = ',')

#Drop rows where either encoding was only observed once
encoded_df = encoded_df[encoded_df['n_g1'] > 1]
encoded_df = encoded_df[encoded_df['n_g2'] > 1]


#Convert the alignment positions to hxb2 coordinates

#Make a dictionary to store all of the conversion dictionaries
all_par_hxb2_ath = {}
for curr_par in encoded_df['pt'].unique():
    #Get the path to the participant's data
    inFile = inDir_3BNC + curr_par + '/835_' + curr_par + '_NT.translated.fasta'

    #Now load the datastructures from the fasta
    seqArr, seq_info_df = du.fasta_to_dataStructs(inFile, clyde2024=True)
    hxb2_nuc_coords_hta, hxb2_nuc_coords_ath = du.hxb2_mapping_dict(seqArr, 
                                                                start_pos = 1)
    all_par_hxb2_ath[curr_par] = hxb2_nuc_coords_ath


encoded_df['hxb2_coord_AA'] = encoded_df.apply(
    lambda row: all_par_hxb2_ath[row['pt']][row['AA_pos']-1], axis=1)

#We handchecked alignment at all of the positions and now we will drop the
#positions that were poorly aligned
encoded_df['hxb2_coord_AA'] = encoded_df['hxb2_coord_AA'].astype(float)
#Round the hxb2 coordinates
encoded_df['hxb2_coord_AA'] = encoded_df['hxb2_coord_AA'].round(3)
encoded_df = encoded_df[~encoded_df['hxb2_coord_AA'].isin(poor_alignment_hand_checked)]

#Now we will save the data to a file
outFile = inDir + 'multiply_encoded_3BNC117_filtered.csv'
encoded_df.to_csv(outFile, index=False)

#Also save a version formatted for input into the summary plot
encoded_df = encoded_df[['hxb2_coord_AA', 'pt']]
encoded_df.rename(columns={'pt': 'participant'}, inplace=True)
encoded_df['method'] = 'multi_encodings'
encoded_df.to_csv(inDir + 'multiply_encoded_3BNC117_filtered_input.csv', index=False)
