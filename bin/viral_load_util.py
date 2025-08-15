import math
import numpy as np


#This file will contain sequences for working with the viral load data

def get_rebound_closest(vl_df, sampling_times):
    """This function takes in the viral load filtered dataframe. It has columns
    for each timepoint and each row is the viral load for a patient at that 
    timepoint. The dataframe also has a column for rebound which indicates
    which sampling timepoint is closest to the rebound timepoint (since 
    we don't always have sequences from the rebound timepoint) 
    ---------------------------------------------------------------------------
    Params:
    -------
    vl_df: pd.DataFrame, a dataframe with viral load data for a single 
                    participant
    sampling_times: np.array, an array of the sampling timepoints, in string
                    format.

    Returns:
    --------
    rebound_tuple: a two-tuple where the zero element is the closest timepoint
                     to rebound and the one element is the actual rebound time
    """
    if len(np.unique(vl_df['Study_ID'])) > 1:
        raise ValueError('The dataframe should only contain data for one participant.')
    #Reformat our sampling times into integer days
    sampling_times = [int(x[1:]) if x[0] == 'D' else (7*int(x[1:]))\
                                      for x in sampling_times]
    sampling_times = np.unique(sampling_times)
    sampling_times = np.sort(sampling_times)
    sampling_times = np.asarray(sampling_times)

    #Get the rebound time
    rebound_time_str = vl_df['Nadir'].values[0]

    #If the rebound time is a string, we need to convert it to an integer
    if isinstance(rebound_time_str, str):
        rebound_time = np.nan
        if rebound_time_str[0] == 'D':
            rebound_time = int(rebound_time_str[1:])
        else:
            rebound_time = 7 * int(rebound_time_str[1:])

        #We need to find the closest sample to the rebound time
        closest = (np.abs(sampling_times - rebound_time)).argmin() + 1
        if closest == len(sampling_times):
            closest = np.nan
        elif closest == 0:
            closest = sampling_times[1]
        else:
            closest = sampling_times[closest]
        if not math.isnan(closest) and closest >= 7:
            closest = 'W' + str(closest // 7)
        else:
            closest = 'D' + str(closest)
    else:
        closest = np.nan

    return (closest, rebound_time_str)