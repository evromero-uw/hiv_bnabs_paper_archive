import numpy as np
import networkx as nx
from scipy.spatial.distance import hamming

#This file will contain utilities for performing the matching of sequences
#to ancestors at previous time points based on hamming distance.

def find_closest_anc(masked_seq_arr, query_ind, candidate_info, verbose=False,
                      return_multi=False, ignore_gap=False):
    """ This function takes in a masked sequence array where each row is
    a sequence and each column is a base at a specific position. Given a query
    index, the function will find the closest ancestor to the query sequence
    based on hamming distance. Only sequences with indices listed in 
    candidate_info will be considered as candidate ancestors. 
    ---------------------------------------------------------------------------
    Params:
    -------
    masked_seq_arr: np.array, shape (n, m), where n is the number of sequences
                    and m is the length of the sequences. Each row is a 
                    sequence and each column is a base at a specific position.
    query_ind: int, the row index of the query sequence in masked_seq_arr
    candidate_info: pd.DataFrame, contains the sequence info for each candidate
                    ancestor. Includes a column 'seq_index' which contains the
                    row index of the sequence in masked_seq_arr.
    verbose: bool, whether to print out additional information.
    ignore_gap: bool, whether to ignore positions with gaps in either sequence
                when calculating the hamming distance.
    return_multi: bool, whether to return multiple ancestors or just the closest
                ancestor.

    Returns:
    --------
    closest_ancs: np.array, an array of the row indices of the closest ancestors
                    to the query sequence.
    min_h_dist: float, the hamming distance between the query sequence and the
                    closest ancestor.
    percent_gap: float, the percentage of positions used for the match to the
                    closest ancestor.

    """

    #Get the query sequence
    query_seq = masked_seq_arr[query_ind, :]

    #Get the candidate sequences
    candidate_inds = candidate_info['seq_index'].values
    candidate_inds = candidate_inds.astype(int)
    candidate_seqs = masked_seq_arr[candidate_inds, :]
    if verbose:
        print(f"Query sequence: {query_seq}")
        print(f"Candidate sequences: {candidate_seqs}")

    #Calculate the hamming distance between the query sequence and each
    #candidate sequence
    hamming_dists = np.zeros(candidate_seqs.shape[0])
    percent_gap = np.zeros(candidate_seqs.shape[0])
    for i in range(candidate_seqs.shape[0]):
        possible_anc = candidate_seqs[i, :]

        if not ignore_gap:
            #Get indices where either sequence has a gap
            gap_indices = np.logical_or(query_seq == '-', possible_anc == '-')
            if verbose:
                print(f"Gap indices: {gap_indices}")

            #Don't consider all gap sequences as ancestors
            if np.sum(gap_indices) == len(possible_anc):
                if verbose:
                    print(f"Gap sequence: {possible_anc}")
                #Set the hamming distance to a large number
                hamming_dists[i] = 800
                percent_gap[i] = 1
                continue
            

            #Remove the gap indices from the sequences
            seq1 = query_seq[~gap_indices]
            seq2 = possible_anc[~gap_indices]
            if verbose:
                print(f"Seq1: {seq1}")
                print(f"Seq2: {seq2}")

            hamming_dists[i] = hamming(seq1, seq2)
            percent_gap[i] = np.sum(gap_indices) / len(gap_indices)

        else:
            seq1 = query_seq
            seq2 = possible_anc

            hamming_dists[i] = hamming(seq1, seq2)
            percent_gap[i] = None

    #Find the index of the sequence with the smallest hamming distance
    if verbose:
        print(f"Hamming distances: {hamming_dists}")
        print(f"Percent gap: {percent_gap}")
    min_h_dist = np.min(hamming_dists)
    if min_h_dist == 800:
        return None, None, None
    min_ind = np.argmin(hamming_dists)
    if return_multi:
        min_ind = list(np.where(hamming_dists == min_h_dist)[0])
    
    percent_gap_for_anc = percent_gap[min_ind]
    closest_anc_ind = candidate_inds[min_ind]

    return closest_anc_ind, min_h_dist, percent_gap_for_anc

def count_lineages(closest_seqs):
    """This function will take in a dictionary where each key is a sequence
    at day 0 and each value is the descendants of that sequence at the
    next timepoints. The function will then count the number of lineages
    present in the data (how many unique, unshared ancestors are present).
    ---------------------------------------------------------------------------
    Params:
    -------
    closest_seqs: dict, a dictionary where each key is a sequence at day 0 and
                    each value is the descendants of that sequence at the 
                    following timepoints.
    Returns:
    --------
    seq_lineage_dict: dict, a dictionary where each key is a sequence and each
                        value is the lineage it belongs to.
    """

    #Loop through day 0 ancestors to make sets of lineages
    anc_seqs = list(closest_seqs.keys())

    #First make sets of each of the keys ancestors have
    ancs_to_merge = []
    non_singleton_ancs = []

    #Compare all pairs of ancestors and see if they share any descendants
    for i in range(len(anc_seqs)):
        anc_i = anc_seqs[i]
        anc_i_descendants = set(closest_seqs[anc_i])

        for j in range(i+1, len(anc_seqs)):
            anc_j = anc_seqs[j]
            anc_j_descendants = set(closest_seqs[anc_j])

            if not anc_i_descendants.isdisjoint(anc_j_descendants):
                ancs_to_merge.append((anc_i, anc_j))
                non_singleton_ancs.append(anc_i)
                non_singleton_ancs.append(anc_j)

            elif anc_i in anc_j_descendants:
                ancs_to_merge.append((anc_i, anc_j))
                non_singleton_ancs.append(anc_i)
                non_singleton_ancs.append(anc_j)

            elif anc_j in anc_i_descendants:
                ancs_to_merge.append((anc_i, anc_j))
                non_singleton_ancs.append(anc_i)
                non_singleton_ancs.append(anc_j)
    

    #Now we'll merge the sets of lineages
    merge_graph = nx.Graph(ancs_to_merge)
    # merged_ancs = nx.connected_components(merge_graph)
    merged_lineages = list(map(tuple, nx.connected_components(merge_graph)))
    singleton_ancs = [tuple([anc]) for anc in anc_seqs if anc not in non_singleton_ancs]
    merged_lineages.extend(singleton_ancs)

    #A dictionary where each key is a sequence and each value is the lineage it
    #belongs to
    seq_lineage_dict = {}
    
    #Now, we need to add all of the descendants of the ancestors to the sets
    for i in range(len(merged_lineages)):
        curr_lineage = list(merged_lineages[i])
        descendants = []

        for anc in curr_lineage:
            descendants.extend(closest_seqs[anc])
        
        #Now update the lineage dictionary
        for seq in curr_lineage + descendants:
            seq_lineage_dict[seq] = i

    return seq_lineage_dict