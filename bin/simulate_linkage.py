import sys
sys.path.append('../bin/')
sys.path.append('../data/clyde2023/')
import data_util
import numpy as np
from functools import reduce


#This module contains the functions necessary to simulate scenarios of complete
#linkage disequilibrium and also equilibrium. These simulations can be used to
#form null expectations.

def initialize_simArr(seqArr, place_res_muts, hxb2_res_positions=None, 
                     arr_res_positions=None):
    """Takes in a sequence array, and the positions of the resistance sites in
    array and hxb2 coordinates. Then it generates an array that matches the 
    size of the sequence array but contains only consensus sequences. Lastly,
    if place_res_muts is true, it copies the original sequence array entries
    around the resistance sites into the new consensus simulation array. NOTE:
    THE FIRST SEQUENCE IN THE ARRAY MUST BE THE HXB2 REFERENCE SEQUENCE.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    seqArr: np.array, an array where each row is a sequence and each column is
                a position in that sequence. The sequences are aligned and the
                first row is the reference sequence (HXB2).
    place_res_muts: bool, True if the original sequences around the resistance
                sites should be copied into the simulation array. False if the
                simulation array should only contain consensus sequences.
    hxb2_res_positions: tuple, a tuple where each entry is itself a tuple 
                containing the start and end positions of the resistance 
                mutation in HXB2 coordinates ((start, end) ,(start2, end2))
    arr_res_positions: tuple, a tuple where each entry is itself a tuple
                containing the start and end positions of the resistance
                mutation in array coordinates.


    Returns:
    --------
    simArr: np.array, an array where each row is a sequence and each column is
                a position in that sequence. The array contains the same number
                of (nonreference) sequences as the input sequence array but the
                sequences are all consensus sequences except for at the
                specified resistance sites.
    """
    #Make a copy of the sequence array with the reference sequence and then
    #remove the reference sequence from the main array
    seqArr_with_ref = seqArr.copy()
    seqArr= seqArr[1:, :]

    #Get the consensus sequence for the array minus HXB2
    consensus_seq = data_util.generate_consensus_seq(seqArr)

    #Now make an array where every row is the consensus sequence which will be
    #and represents each initialized genome for the simulation
    simulationArr = np.repeat([consensus_seq], seqArr.shape[0], axis = 0)

    #Our work is done if we don't need to copy the resistance sequences
    if not place_res_muts:
        return simulationArr

    #Now copy the original sequences around the resistance sites
    masked_arr, res_pos_seqs = data_util.mask_hxb2_coords(seqArr_with_ref,
                        start_pos = hxb2_res_positions[0][0],
                        end_pos = hxb2_res_positions[0][1],
                        arr_start_pos = 1,
                        second_site = (hxb2_res_positions[1][0],
                                        hxb2_res_positions[1][1]))
    
    start_pos_1 = min(arr_res_positions[0][0], arr_res_positions[1][0])
    end_pos_1 = res_pos_seqs[0].shape[1] + start_pos_1

    start_pos_2 = max(arr_res_positions[0][0], arr_res_positions[1][0])
    end_pos_2 = res_pos_seqs[1].shape[1] + start_pos_2

    simulationArr[:, start_pos_1:end_pos_1] = res_pos_seqs[0][1:, :]
    simulationArr[:, start_pos_2:end_pos_2] = res_pos_seqs[1][1:, :]

    return simulationArr


def place_muts_LD(haplotype_list, seqArr, arr_res_pos, seg_freq_dict,
                 freq_tols = (0,0)):
    """Simulates a scenario of nearly complete linkage disequilibrium using
    the SNP frequency data from the in vivo sequence array.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    haplotype_list: list, a list of numpy arrays where each array consists of 
                a single haplotype where each identical row is a sequence and
                each column is a position in that sequence. The list is sorted
                from highest to lowest frequency.
    seqArr: np.array, an array where each row is a sequence and each column is
                a position in that sequence. The sequences are aligned and the
                first row is the reference sequence (HXB2).
    arr_res_pos: tuple, a tuple where each entry is itself a tuple containing
                the start and end positions of the resistance mutation in array
                coordinates ((start, end) ,(start2, end2)).
    seg_freq_dict: dict, a dictionary where each key is the array position of a
                segregating SNP and the value is a tuple containing the minor
                allele in position 0 and its frequency in position 1.
    freq_tols: the tolerance for how much lower the frequency of a haplotype
            can be than the frequency of the mutation that's being placed
            on it. The first entry is the tolerance for the frequency and 
            the second is the tolerance for the number of sequences mutated
            (which will be floored to the number of sequences in the
            haplotype if it is within the tolerance) otherwise an error 
            will be thrown.
    Returns:
    --------
    haplotype_list: list, the updated list of haplotypes with mutations placed
                on them.
    """
    #Sort the segregating sites by their frequency
    sorted_seg_freqs = data_util.get_closest_res(seg_freq_dict, arr_res_pos,
                                                        freq_sorted=True)

    freq_sorted_res_haps = get_res_mut_haps(seqArr[1:], arr_res_pos, 
                                                        hxb2_incl = False)

    #Get the highest frequency drug resistance mutation and add it to the list
    #that we will be looping through
    highest_freq_res_hap = freq_sorted_res_haps[0]
    highest_freq_res_hap = ('resistance', None, 0, highest_freq_res_hap[0],
                            highest_freq_res_hap[1])
    sorted_seg_freqs.append(highest_freq_res_hap)
    sorted_seg_freqs.sort(key = lambda x: x[4], reverse = True )

    #Loop through the segregating sites by their frequencies
    for curr_seg_site in sorted_seg_freqs:

        #If the segregating site is a resistance mutation, then we need to
        #place each of the individual haplotypes
        if curr_seg_site[0] == 'resistance':
            print('Placing resistance mutations')
            haplotype_list = sorted(haplotype_list, key = lambda x: x.shape[0],
                                    reverse = True)
            haplotype_list = place_res_muts_LD_recursive(
                                                freq_sorted_res_haps,
                                                arr_res_pos,
                                                haplotype_list, [])[0]

            #Now, merge any haplotypes that only differ at resistance sites
            haplotype_list = merge_res_haps(haplotype_list, arr_res_pos)
            continue

            

        #Get the minor allele's frequency
        minor_allele_freq = curr_seg_site[4]
        seg_site_arr_pos = curr_seg_site[0]
        seg_site_allele = curr_seg_site[3]

        #Calculate the haplotype frequencies
        total_haps = reduce(lambda x, y: x + y, 
                            [x.shape[0] for x in haplotype_list])
        #make a dictionary where the index is the key and the value is the freq
        hap_freqs = dict([(i, x.shape[0]/total_haps)\
                    for i,x in enumerate(haplotype_list)])

        #Get the haplotypes that are big enough to place the mutation on
        big_haps = [i for i in range(len(haplotype_list))\
                    if -freq_tols[0] < (hap_freqs[i] - minor_allele_freq)]

        
        
        #Otherwise place the mutation on a haplotype that can fit it
        #Next we need to choose a haplotype to mutate and place it on that
        #haplotype
        hap_to_mutate = np.random.choice(big_haps)
        hap_to_mutate = haplotype_list.pop(hap_to_mutate)
        
        num_mutated = int(np.round(minor_allele_freq * total_haps))

        if hap_to_mutate.shape[0] < num_mutated:
            print(num_mutated, hap_to_mutate.shape[0])
            if (num_mutated - hap_to_mutate.shape[0]) <= freq_tols[1]:
                num_mutated = hap_to_mutate.shape[0]
            else:
                print('Haplotype too small')
                quit()

        #Split the haplotype into two arrays
        mutated_hap = hap_to_mutate[:num_mutated, :]
        og_hap = hap_to_mutate[num_mutated:, :]


        #Change the sequence of the mutated haplotype
        mutated_hap[:, seg_site_arr_pos] = seg_site_allele

        #Add the new haplotypes to the list
        if mutated_hap.shape[0] > 0:
            haplotype_list.append(mutated_hap)
        if og_hap.shape[0] > 0:
            haplotype_list.append(og_hap)
    
    return haplotype_list

def place_muts_LE(simulationArr, arr_res_pos, seg_freq_dict):
    """Simulates a scenario of linkage equilibrium using the SNP frequency data
    from the in vivo sequence array.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    simulationArr: np.array, an array where each row is a sequence and each 
                column is a position in that sequence. All fo the positions 
                are consensus sequences except for the resistance sites.
    seg_freq_dict: dict, a dictionary where each key is the array position of a
                segregating SNP and the value is a tuple containing the minor
                allele in position 0 and its frequency in position 1.
    arr_res_pos: tuple, a tuple where each entry is itself a tuple containing
                the start and end positions of the resistance mutation in array
                coordinates ((start, end) ,(start2, end2)).
    
    Returns:
    --------
    haplotype_list: list, the updated list of haplotypes with mutations placed
                on them.
    """
    #Randomize the order of the sequences so that DRM's are not all grouped
    #together
    np.random.shuffle(simulationArr)

    #Get a list of segregating sites that aren't in the resistance region
    #Each element in the list is a three tuple (segregating site, closest 
    #resistance, distance to closest)
    seg_sites_list = data_util.get_closest_res(seg_freq_dict, arr_res_pos)

    #Sort by the distance to the closest resistance mutation
    seg_sites_list.sort(key = lambda x: x[2])

    #I need to save the simulation array as a list of haplotypes
    haplotype_list = [simulationArr]

    #Now loop through the segregating sites and place the mutations on the genomes
    for curr_site in seg_sites_list:
        site_pos = curr_site[0]
        site_allele = seg_freq_dict[site_pos][0]
        site_freq = seg_freq_dict[site_pos][1]

        #Make a list to store the updated haplotypes in
        updated_hap_list = []
        
        total_mutated = 0

        #Place the mutation on each haplotype in linkage equilibrium
        for curr_hap in haplotype_list:

            num_mutated = np.random.binomial(n = curr_hap.shape[0],
                                                    p = site_freq, size = 1)[0]
            total_mutated += num_mutated

            #If the haplotype is a single sequence
            if curr_hap.shape[0] == 1:
                if num_mutated == 1:
                    curr_hap[:, site_pos] = site_allele
                updated_hap_list.append(curr_hap)
                continue
            
            #Split the haplotype into two arrays
            mutated_hap = curr_hap[:num_mutated, :]
            og_hap = curr_hap[num_mutated:, :]

            #Change the sequence of the mutated haplotype
            mutated_hap[:, site_pos] = site_allele

            #Add the new haplotypes to the list
            if mutated_hap.shape[0] > 0:
                updated_hap_list.append(mutated_hap)
            if og_hap.shape[0] > 0:
                updated_hap_list.append(og_hap)

        #Update the haplotype list so the haplotypes are now split on this
        #segregating site
        haplotype_list = updated_hap_list

    return haplotype_list

###############################################################################
################################ Helper Functions #############################
###############################################################################
def get_res_mut_haps(seqArr, arr_res_positions, hxb2_incl):
    """Takes in a sequence array and a tuple containing the array coordinates
    of positions with resistance mutations. Returns a frequency sorted list of
    tuples where the zero element is the haplotype and the first element is the
    frequency of that haplotype. 
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    seqArr: np.array, an array where each row is a sequence and each column is
                a position in that sequence. 
    arr_res_pos: tuple, a tuple where each entry is itself a tuple containing
                the start and end positions of the resistance mutation in array
                coordinates Ex: ((start1, end1) ,(start2, end2)).
    hxb2_incl: bool, true if the hxb2 reference sequence is included as the
                first row of the sequence array. 

    Returns:
    --------
    hap_freq_list: list, a list of tuples where the zero element is the haplotype
                and the first element is the frequency of that haplotype. The 
                list is sorted from highest to lowest frequency.
    """
    #Remove the reference sequence if it's included
    if hxb2_incl:
        seqArr = seqArr[1:]

    #Make a list to store the resistance mutation slices in
    region_slices = []

    #Get the nucleotides in each resistance region
    for res_region in arr_res_positions:
        curr_slice = seqArr[:, res_region[0]:res_region[1] + 1]
        region_slices.append(curr_slice)

    #Next concatenate the regions together into one array
    res_arr = np.hstack(tuple(region_slices))

    #And make each element a string which we can recover the haplotypes from
    res_arr = np.apply_along_axis(lambda x: ''.join(x), 1, res_arr)

    #We want to get all of the haplotypes and group them by frequency
    res_haps, res_counts = np.unique(res_arr, return_counts = True)
    res_freqs = res_counts / sum(res_counts)

    #Now loop through the haplotypes and add them to the dictionary
    hap_freq_list = [(res_haps[i], res_freqs[i]) for i in range(len(res_haps))]
    hap_freq_list.sort(key = lambda x: x[1], reverse = True)

    return hap_freq_list

def place_res_muts_LD_recursive(freq_sorted_res_haps, arr_res_positions, 
                                unplaced_hap_list, placed_hap_list,
                                verbose = False):
    """Recursively places all of the resistance mutations on the haplotypes in 
    the placed haplotype list. This process is more difficult than placing other
    mutations because the additional minor alleles can necessitate that
    homoplasies are introduced into the array.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    freq_sorted_res_haps: list, a list of tuples produced by get_res_mut_haps
                where the zero element is the haplotype and the first element
                is the frequency of that haplotype. The list is sorted from
                highest to lowest frequency.
    arr_res_positions: tuple, a tuple where each entry is itself a tuple
                containing the start and end positions of the resistance
                mutation in array coordinates. NOTE: The array res positions
                must be in fed in the same order as was fed to get_res_mut_haps
                to produce the freq_sorted_res_haps list.
    placed_hap_list: list, a list of numpy arrays where each array consists of 
                a single haplotype where each row is a sequence and each column
                is a position in that sequence. All of these haplotypes have
                already had resistance mutations placed on them. This list is 
                sorted from highest to lowest frequency.
    unplaced_hap_list: list, the same format as the placed_hap_list except that
                these haplotypes have not had resistance mutations placed on
                them yet. This list is sorted from highest to lowest frequency.
    verbose: bool, True if the function should print out information for 
                debugging purposes. False if it should not.

    Returns:
    --------
    placed_hap_list: list, the updated list of haplotypes with resistance
                mutations placed on them.
    unplaced_hap_list: list, the updated list of haplotypes that still need
                resistance mutations placed on them.
    """
    if verbose:
        print('Placed haps: ' + str(len(placed_hap_list)))
        print('Unplaced haps: ' + str(len(unplaced_hap_list)))
        print('DRMs left: ' + str(len(freq_sorted_res_haps)))
        print('-------------------------------------------')

    #Check if there are any DRMs left to place
    if len(freq_sorted_res_haps) == 0:
        if verbose:
            print('Base case reached')
        return placed_hap_list, unplaced_hap_list
    
    #Get the next haplotype to place
    curr_drm_hap = freq_sorted_res_haps.pop(0)
    curr_drm_freq = curr_drm_hap[1]

    #Make a list of haplotypes to place the drm on
    haps_to_place = []
    haps_to_place_freq = 0

    #Get the frequencies of haplotypes without resistance mutations
    total_hap_count = reduce(lambda x, y: x + y, [x.shape[0] for x in \
                                unplaced_hap_list + placed_hap_list])
    hap_freqs_dict, sorted_hap_freqs, unplaced_hap_dict = \
        get_hap_freqs(unplaced_hap_list, total_hap_count)

    #Get rid of the unplaced hap list so we don't accidentally use it again
    del unplaced_hap_list

    #Loop through the unplaced haplotypes and merge them until we can place the
    #drm on them
    while haps_to_place_freq < curr_drm_freq:

        #Get the biggest haplotype
        big_hap_freq = sorted_hap_freqs.pop(0)
        dict_entry_haps = hap_freqs_dict[big_hap_freq]
        if len(dict_entry_haps) > 1:
            big_hap_ind = dict_entry_haps.pop(np.random.randint(0,
                                                    len(dict_entry_haps)-1))
        else:
            big_hap_ind = dict_entry_haps.pop(0)
        big_hap = unplaced_hap_dict[big_hap_ind]

        #Add it to the list of haplotypes to place the drm on
        haps_to_place.append(big_hap)
        haps_to_place_freq += big_hap_freq

        #Remove it from the dictionary of unplaced haplotypes
        if len(dict_entry_haps) == 0:
            del hap_freqs_dict[big_hap_freq]
        else:
            hap_freqs_dict[big_hap_freq] = dict_entry_haps
        
        del unplaced_hap_dict[big_hap_ind]

    #Count how many copies of the DRM to place
    left_to_place = int(curr_drm_freq * total_hap_count)

    #Shuffle the haplotypes so the haplotype with a DRM segregating is random
    np.random.shuffle(haps_to_place)

    #Make a place to store any haps that are left without a mutation
    new_unplaced = []

    #Now we need to place the drm on the haplotypes
    #I am using a greedy approach where I place as many copies of DRMs as
    #possible each time. So that only one haplotype will have the DRM 
    #segregating on it
    for curr_hap in haps_to_place:
        curr_hap_len = curr_hap.shape[0]

        #If we need to place fewer DRMs than the length of the haplotype
        if left_to_place < curr_hap_len:
            curr_hap_mut = curr_hap[:left_to_place, :]

            #get the unplaced copies and add them back into the unplaced dict
            curr_hap_no_mut = curr_hap[left_to_place:, :]
            new_unplaced.append(curr_hap_no_mut)
            left_to_place = 0

        #Otherwise place the DRMs on the full haplotype    
        else:
            curr_hap_mut = curr_hap
            left_to_place -= curr_hap_len

        #place the DRM
        curr_hap_mut = place_single_drm(arr_res_positions, curr_drm_hap[0],
                                curr_hap_mut)
        placed_hap_list.append(curr_hap_mut)

    #Now update and sort the haplotype lists
    placed_hap_list.sort(key = lambda x: x.shape[0], reverse = True)
    unplaced_hap_list = new_unplaced + list(unplaced_hap_dict.values())
    unplaced_hap_list.sort(key = lambda x: x.shape[0], reverse = True)

    #Recursively call the function to place the other DRMs
    return place_res_muts_LD_recursive(freq_sorted_res_haps, arr_res_positions,
                                        unplaced_hap_list, placed_hap_list)

def place_single_drm(arr_res_positions, drm_hap, hap_arr):
    """Given a array representing haplotype and a single drm, places the drm
    on the haplotype and returns the resulting array.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    arr_res_positions: tuple, a tuple where each entry is itself a tuple
                containing the start and end positions of the resistance
                mutation in array coordinates. NOTE: The array res positions
                must be in fed in the same order as was fed to get_res_mut_haps
                to produce the freq_sorted_res_haps list.
    drm_hap: str, a string representing the haplotype with the drm mutations
                placed on it.
    """
    #Check how long each of the resistance blocks is
    res_block_lengths = [x[1] - x[0] + 1 for x in arr_res_positions]

    #Separate out the different resistance blocks
    split_drms = []
    for i in range(len(res_block_lengths)):
        curr_res_block = drm_hap[:res_block_lengths[i]]
        split_drms.append(curr_res_block)
        drm_hap = drm_hap[res_block_lengths[i]:]
    
    #Now loop through the blocks and add them to the haplotype array
    for i in range(len(arr_res_positions)):
        curr_res_start = arr_res_positions[i][0]
        curr_res_end = arr_res_positions[i][1]

        #Get the nucleotides in the resistance block
        curr_res_block = split_drms[i]

        #Make an array holding the nucleotides in the resistance block
        curr_res_block_arr = np.array(list(curr_res_block))
        curr_res_block_arr = np.repeat([curr_res_block_arr], 
                                       hap_arr.shape[0], axis = 0)
        
        #Now add the resistance block to the haplotype array
        hap_arr[:, curr_res_start:curr_res_end + 1] = curr_res_block_arr
        
    return hap_arr

def get_hap_freqs(haplotype_list, num_sequences):
    """ Takes in a list of haplotypes where each element is a numpy array of 
    identical sequences and an integer indicating the total number of sequences
    in the analysis. Then, returns a set of data structures used for
    enumerating over the haplotypes by frequency.
    ---------------------------------------------------------------------------
    Returns:
    --------
    hap_freqs_dict: dict, a dictionary where the keys are the frequencies of
                the haplotypes and the values are the indices of the haplotypes
                in the haplotype list.
    sorted_hap_freqs: list, a list of the frequencies of the haplotypes sorted
                from highest to lowest.
    unplaced_hap_dict: dict, a dictionary where the keys are the indices of the
                haplotypes in the haplotype list and the values are the 
                haplotypes themselves.
    """
    #Make a dictionary where the key is the freq and the value is the index
    hap_freqs_dict = {}
    sorted_hap_freqs = []
    for i,x in enumerate(haplotype_list):
        curr_freq = x.shape[0]/num_sequences
        sorted_hap_freqs.append(curr_freq)
        if curr_freq in hap_freqs_dict.keys():
            hap_freqs_dict[curr_freq].append(i)
        else:
            hap_freqs_dict[curr_freq] = [i]
    
    #Make a dictionary where the index is the key and the value is the haplotype
    unplaced_hap_dict = dict([(i , x) for i,x in enumerate(haplotype_list)])

    sorted_hap_freqs = sorted(sorted_hap_freqs, reverse = True)

    return hap_freqs_dict, sorted_hap_freqs, unplaced_hap_dict

def merge_res_haps(haplotype_list, arr_res_pos):
    """Takes in a list of haplotypes where each element is a numpy array of 
    identical sequences and merges any haplotypes that only differ at the
    resistance sites.
    ---------------------------------------------------------------------------
    Params:
    -------
    haplotype_list: list, a list of haplotypes where each element is a numpy
                array of identical sequences.
    arr_res_pos: tuple, a tuple where each entry is itself a tuple containing
                the start and end positions of the resistance mutation in array
                coordinates ((start, end) ,(start2, end2)).

    Returns:
    --------
    merged_hap_list: list, the updated list of haplotypes with resistance
                mutations placed on them.
    """
    identical_dict = {}
    arr_res_pos = list(arr_res_pos)
    arr_res_pos.sort(key = lambda x: x[0])

    #Compare each set of haplotypes to each other
    for i in range(len(haplotype_list)):
        first_hap = haplotype_list[i]
        if first_hap.shape[0] == 0:
            continue
        first_before = first_hap[0, :arr_res_pos[0][0]]
        first_mid = first_hap[0, arr_res_pos[0][1]+1:arr_res_pos[1][0] + 1]
        first_after = first_hap[0, arr_res_pos[1][1]+1:]

        for j in range(i + 1, len(haplotype_list)):
            second_hap = haplotype_list[j]
            if second_hap.shape[0] == 0:
                continue
            second_before = second_hap[0, :arr_res_pos[0][0]]
            second_mid = second_hap[0, arr_res_pos[0][1]+1:arr_res_pos[1][0] + 1]
            second_after = second_hap[0, arr_res_pos[1][1]+1:]

            before_identical = np.array_equal(first_before, second_before)
            mid_identical = np.array_equal(first_mid, second_mid)
            after_identical = np.array_equal(first_after, second_after)

            if before_identical and mid_identical and after_identical:
                if i in identical_dict.keys():
                    identical_dict[i].append(j)
                else:
                    identical_dict[i] = [j]

    #Now merge the haplotypes that are identical
    merged_hap_list = []
    already_merged = set()
    for i in range(len(haplotype_list)):
        #Skip already merged haplotypes
        if i in already_merged:
            continue
        #If we have to merge
        elif i in identical_dict.keys():
            curr_hap = haplotype_list[i]
            identical_indices = identical_dict[i]
            for j in identical_indices:
                curr_hap = np.vstack((curr_hap, haplotype_list[j]))
                already_merged.add(j)
            merged_hap_list.append(curr_hap)
        #Otherwise just add the haplotype to the list
        else:
            merged_hap_list.append(haplotype_list[i])
    
    return merged_hap_list