import os
import sys
sys.path.append('../../bin/')
import numpy as np
import simulate_linkage


# This file runs the LD and LE null simulations for a directory of data
# from a single participant.

def run_one_LD_sim(participant_dat, freq_tols, hxb2_seq):
    """ This function runs one iteration of the LD null simulation and 
    generates a sequence array of simulated sequences in nearly complete
    linkage disequilibrium.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    participant_dat: Pardata object, the object containing the participant
                        data we will simulate from (matching allele frequencies
                        and segregating sites).
    freq_tols: the tolerance for how much lower the frequency of a haplotype
                        can be than the frequency of the mutation that's being 
                        placed on it. The first entry is the tolerance for the
                        frequency and the second is the tolerance for the
                        number of sequences mutated (which will be floored to
                        the number of sequences in the haplotype if it is 
                        within the tolerance) otherwise an error will be thrown
    hxb2_seq: np.array, the HXB2 sequence

    Returns:
    --------
    sim_seqArr: np.array, the simulated sequence array. The first row is the
                        hxb2 sequence and the rest are the simulated sequences
    """
    #Step 1: Initialize the linkage disequilibrium simulation array
    ###########################################################################
    simulationArr_LD = simulate_linkage.initialize_simArr(
                                            participant_dat.seq_arr,
                                            place_res_muts=False)
    
    #We will actually store the sequences as a list of haplotypes
    haplotype_list_LD = [simulationArr_LD]

    #Step 2: Simulate the evolution of the sequences under total linkage 
    # disequilibrium
    ############################################################################
    haplotype_list_LD = simulate_linkage.place_muts_LD(haplotype_list_LD, 
                                        participant_dat.seq_arr, 
                                        participant_dat.arr_res_positions,
                                        participant_dat.seg_freq_dict,
                                        freq_tols = freq_tols)
    simulationArr_LD = np.vstack(haplotype_list_LD)    
    simulationArr_LD = np.vstack([hxb2_seq, simulationArr_LD])
    return simulationArr_LD

def run_one_LE_sim(participant_dat, hxb2_res_positions, hxb2_seq):
    """ This function runs one iteration of the LE null simulation and 
    generates a sequence array of simulated sequences in linkage equilibrium.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    participant_dat: Pardata object, the object containing the participant
                        data we will simulate from (matching allele frequencies
                        and segregating sites).
    hxb2_res_positions: tuple, the positions of the resistance mutations in
                        HXB2 coordinates
    hxb2_seq: np.array, the HXB2 sequence

    Returns:
    --------
    sim_seqArr: np.array, the simulated sequence array. The first row is the
                        hxb2 sequence and the rest are the simulated sequences
    """
    #Step 1: Initialize the linkage equilibrium simulation array
    ###########################################################################
    simulationArr_LE = simulate_linkage.initialize_simArr(
                                    participant_dat.seq_arr,
                                    place_res_muts = True,
                                    hxb2_res_positions =
                                        hxb2_res_positions,
                                    arr_res_positions =
                                        participant_dat.arr_res_positions)

    #We will actually store the sequences as a list of haplotypes
    haplotype_list_LE = [simulationArr_LE]

    #Step 2: Simulate the evolution of the sequences under total linkage
    #equilibrium
    ###########################################################################
    haplotype_list_LE = simulate_linkage.place_muts_LE(simulationArr_LE,
                                            participant_dat.arr_res_positions,
                                            participant_dat.seg_freq_dict)
    simulationArr_LE = np.vstack(haplotype_list_LE)
    simulationArr_LE = np.vstack([hxb2_seq, simulationArr_LE])
    return simulationArr_LE
 


