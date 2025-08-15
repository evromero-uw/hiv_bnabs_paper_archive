import os
import sys
import ast
import numpy as np
import pandas as pd

#This file takes in the different data files and puts them into the same format for analysis.
#The data is then saved in a csv file for later use.

def format_dms_data(dms_path):
    """Format the dms data Elena G sent me in order to make the summary plot
    figure. Takes in the filepath and outputs a formatted dataframe which can
    be concatenated with the others produced by the other functions in this
    file.
    ---------------------------------------------------------------------------
    """

    #First read in the data
    dms_df = pd.read_csv(dms_path)

    #Get only the 3BNC117 data
    dms_df = dms_df[dms_df['Antibody'] == '3BNC117']
    dms_df = dms_df[dms_df['Confers'] == 'Resistance']

    dms_df = dms_df.rename(columns={'HXB2.pos': 'hxb2_coord_AA'})
    dms_df['method'] = 'DMS'
    dms_df = dms_df[['method', 'hxb2_coord_AA']]

    return dms_df


    # dms_df = dms_df.rename(columns={'DMS signatures found': 'sig_found',
    #                                 'Varying at DMS sites ': 'varying',
    #                                 'ID': 'participant'})
    
    # #
    
    # new_dms_df = []
    # for index, row in dms_df.iterrows():
    #     varying_sites = row['varying']
    #     varying_sites = varying_sites.split(',')
    #     varying_sites = [int(x.strip()) for x in varying_sites]
    #     for site in varying_sites:
    #         new_dms_df.append([row['participant'], site, row['sig_found']])
    # new_dms_df = pd.DataFrame(new_dms_df, columns=['participant', 'hxb2_coord_AA',
    #                                                'sig_found'])
    # new_dms_df['method'] = 'DMS'
    # return new_dms_df

def format_dms_data_annotation(dms_path):
    """Format the dms data Elena G sent me for incorporation into the
    annotation dataframe instead of the dataframe of varying sites
    """
    #First read in the data
    dms_df = pd.read_csv(dms_path)
    dms_df = format_dms_data(dms_path)
    dms_df = dms_df.rename(columns={'method': 'Reference',
                                    'hxb2_coord_AA': 'HXB2pos'})
    dms_df['Antibody'] = '3BNC117'
    dms_df['Confers'] = 'Resistance'
    dms_df['bNAbClass'] = 'CD4bs'
    return dms_df


def format_eeg_data(eeg_path):
    """Format the sites Elena G found in order to make the summary plot
    figure. Takes in the filepath and outputs a formatted dataframe which can
    be concatenated with the others produced by the other functions in this
    file.
    ---------------------------------------------------------------------------
    """
    #First read in the data
    eeg_df = pd.read_csv(eeg_path)
    eeg_df = eeg_df.rename(columns={'HXB2 pos': 'hxb2_coord',
                                    'Method' : 'method',
                                    'Region' : 'region',})
    eeg_df = eeg_df[eeg_df['Counts'] >= 3]
    eeg_df = eeg_df.drop(columns = ['Counts'])

    new_eeg_df = []
    for index, row in eeg_df.iterrows():
        par_ids = row['IDs']
        par_ids = par_ids.split(',')
        par_ids = [x.strip() for x in par_ids]
        curr_method = row['method']

        print(f"Processing row {index} with method: {curr_method}")

        if '(' in curr_method:
            curr_method = curr_method.split('(')
            curr_method = curr_method[0].strip()
        if ',' in curr_method:
            curr_method = curr_method.split(',')
        else: curr_method = [curr_method]

        for par_id in par_ids:
            for method in curr_method:
                if par_id == '':
                    continue
                if method.strip() == 'Variability':
                    method = 'shared_variability'
                #glycan variability sites are also shared variability sites
                if method.strip() == 'Glycan variability':
                    new_eeg_df.append([par_id,
                                    int(row['hxb2_coord']),
                                    'shared_variability',
                                    row['region']])
                if method.strip() == 'DMS':
                    new_eeg_df.append([par_id,
                                    int(row['hxb2_coord']),
                                    'shared_variability',
                                    row['region']])
                new_eeg_df.append([par_id, 
                                int(row['hxb2_coord']),
                                method.strip(),
                                row['region']])
    new_eeg_df = pd.DataFrame(new_eeg_df, columns=['participant', 'hxb2_coord_AA',
                                                    'method', 'region'])
    print(new_eeg_df)
    return new_eeg_df

def combine_all_data(enc_path, var_path, add_path):
    """Combine all the data into a single dataframe which can be used to make
    the summary plot figure. 
    ---------------------------------------------------------------------------
    """
    enc_df = pd.read_csv(enc_path)
    # enc_df['participant'] = [x[1:4] for x in enc_df['participant']]
    enc_df['method'] = 'multi_encodings'

    var_df = pd.read_csv(var_path)
    var_df['method'] = 'frequency_increase'

    add_df = pd.read_csv(add_path)
    add_df['method'] = 'additional_sites'

    all_data = pd.concat([enc_df, var_df, add_df],
                          ignore_index=True)
    cols_to_keep = ['participant', 'hxb2_coord_AA', 'method', 'region']
    all_data = all_data.drop(columns = [x for x in all_data.columns if x not in cols_to_keep])
    return all_data

