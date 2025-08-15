import os
import sys
sys.path.append('../../../../bin/')
sys.path.append('../../../../bin/wrappers/')
sys.path.append('../../../../data/clyde_westfall_2024_final/10-1074/')
import math
import numpy as np
import pandas as pd
import seaborn as sns
import networkx
import matplotlib as mpl
import matplotlib.pyplot as plt

import matched_anc
import data_util as du
import dataset_metadata
from matplotlib import rcParams
from par_data_class import Pardata
from networkx.algorithms import bipartite

#The code in this file is for making the lineage tracing networks of escape 
#mutations in the 10-1074 cohort.

inDir = '../../../../data/clyde_westfall_2024_final/10-1074/'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/lineage_networks/'

#plot the estimates to show how accurate they are
params = {'figure.figsize':(3.55, 4.3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
FONTSIZE = 6

WT_GRAY = '#bebebe'
rcParams.update(params)


#I am not including 1HC2 and 1HD6K because they do not have both W1 and W4 timepoints
par_list = ['1HB3', '1HC2', '1HC3', '1HD1', '1HD4K', '1HD5K', '1HD6K', '1HD7K',
             '1HD9K', '1HD10K', '1HD11K']


#I will only use a 300bp window to avoid some of the effects of recombination
WINDOWSIZE = 300

TIMEPOINT_DICT_OG = {'W1': ('D0'), 'W4' : ('D0', 'W1')}
REV_TIMEPOINT_DICT_OG = {('D0'): 'W1', ('D0', 'W1'): 'W4'}

TIMEPOINT_DICT_1HC2 = {'W4': ('D0')}
REV_TIMEPOINT_DICT_1HC2 = {('D0'): 'W4'}

TIMEPOINT_DICT_1HD6K = {'W8': ('D0')}
REV_TIMEPOINT_DICT_1HD6K = {('D0'): 'W8'}


#Setting some aesthetic parameters
DESC_SPREAD = 50
ANC_SPREAD = 150

my_colors = mpl.colormaps['plasma'].resampled(4).colors
my_colors = {'A' : my_colors[0],
             'C' : my_colors[1],
             'G' : my_colors[2],
             'T' : my_colors[3]}

#Count all segregating sites
ALLELE_FREQ_THRESH = 0
MULTI_SEG = True

#Loop through each participant and make the network + highlighter plots
for SAMPLE_PAR in par_list:
    print(SAMPLE_PAR)
    if SAMPLE_PAR == '1HC2':
        TIMEPOINT_DICT = TIMEPOINT_DICT_1HC2
        REV_TIMEPOINT_DICT = REV_TIMEPOINT_DICT_1HC2
    elif SAMPLE_PAR == '1HD6K':
        TIMEPOINT_DICT = TIMEPOINT_DICT_1HD6K
        REV_TIMEPOINT_DICT = REV_TIMEPOINT_DICT_1HD6K
    # elif curr_par == '1HD5K':
    #     TIMEPOINT_DICT = TIMEPOINT_DICT_1HD5K
    #     REV_TIMEPOINT_DICT = REV_TIMEPOINT_DICT_1HD5K
    else:
        TIMEPOINT_DICT = TIMEPOINT_DICT_OG
        REV_TIMEPOINT_DICT = TIMEPOINT_DICT_OG

    #The path to the data and the resistance positions in hxb2 coordinates
    inFile = inDir + SAMPLE_PAR + '/885_' + SAMPLE_PAR  + '_NT_filtered.fasta'
    hxb2_res_positions = dataset_metadata.RESISTANCE_POS_NT_HXB2


    #Construct the data object and load the data
    participant_dat = Pardata(inFile, 'clyde2024', SAMPLE_PAR)
    participant_dat.load_data_10_1074(hxb2_res_positions, ALLELE_FREQ_THRESH, MULTI_SEG)
    seq_info_df = participant_dat.seq_info_df
    seq_arr = participant_dat.seq_arr
    seg_freq_dict = participant_dat.seg_freq_dict
    arr_res_set = participant_dat.arr_res_set


    #Slice Window/2 bp on either side of the resistance positions
    arr_res_positions = participant_dat.arr_res_positions
    arr_res_min = min(arr_res_positions[0][0], arr_res_positions[1][0])
    arr_res_max = max(arr_res_positions[0][1], arr_res_positions[1][1])
    arr_res_window_min = arr_res_min - int(WINDOWSIZE/2)
    arr_res_window_max = arr_res_max + int(WINDOWSIZE/2)
    seq_arr = seq_arr[:, arr_res_window_min:arr_res_window_max]


    #Now, for each sequence, we'll get a list of the segregating mutations 
    #Perhaps we should make a dictionary where each key is a sequence 
    #index and each value is a set of minor alleles that sequence contains
    seg_set_dict = du.make_seg_sets(seg_freq_dict, participant_dat.seq_arr, seq_info_df,
                                    arr_res_window_min, arr_res_window_max, arr_res_set)


    #Make a dictionary mapping sequence names to their array indices
    index_to_name = dict(zip(seq_info_df['seq_index'], seq_info_df['orig_name']))

    #Get the hxb2 position of the start of the array now that it's been sliced
    hxb2_res_window_min = participant_dat.hxb2_nuc_coords_ath[arr_res_window_min]
    hxb2_res_window_max = participant_dat.hxb2_nuc_coords_ath[arr_res_window_max]


    #A dictionary where each key is an index at day 0 and each value is all the sequences
    #it is matched with as the closest match at week 1
    closest_seqs = {}
    #The dictionary above with each key and value reversed
    closest_seqs_reversed = {}


    #A dictionary where each key is an index at day 0 and each value is the minimum
    #hamming distance between that sequence and its closest match at week 1
    min_dist_dict = {}

    #Subset the window and mask the resistance positions
    masked_seq_arr, removed_arrs = du.mask_hxb2_coords(seq_arr,
                                        start_pos = hxb2_res_positions[0][0],
                                        end_pos = hxb2_res_positions[0][1],
                                        arr_start_pos = hxb2_res_window_min,
                                        second_site = hxb2_res_positions[1])
    
    #Make a dataframe counting the number of identical sequences at each timepoint
    #This will be used for labeling identical sequences
    labeled_counts = du.label_identical_seqs(masked_seq_arr, seq_info_df)
    labeled_counts['num_identical'] = [len(x) + 1 for x in labeled_counts['identical_seqs']]
    labeled_counts = labeled_counts[['orig_name', 'seq_index', 'num_identical']]

    for curr_timepoint in TIMEPOINT_DICT.keys():

        #Get the sequences from the current timepoint
        curr_time_info = seq_info_df[seq_info_df['time_label'] == curr_timepoint].copy()


        #Get the sequences from the previous timepoint
        prev_timepoint = TIMEPOINT_DICT[curr_timepoint]
        if isinstance(prev_timepoint, str):
            prev_timepoint = [prev_timepoint]

        prev_time_info = seq_info_df[seq_info_df['time_label'].isin(prev_timepoint)].copy()
        prev_time_info = prev_time_info.reset_index(drop = True)
        prev_time_info = du.label_identical_seqs(masked_seq_arr, prev_time_info, ignore_time = True)

        
        curr_time_info = curr_time_info.reset_index(drop = True)
        curr_time_info = du.label_identical_seqs(masked_seq_arr, curr_time_info, ignore_time = True)
            
        #For each sequence, find it's closest match in the previous timepoint
        for index, row in curr_time_info.iterrows():
            curr_seq_ind = row['seq_index']

            #Find the closest ancestor
            closest_anc_ind, min_hdist, percent_gap = \
                matched_anc.find_closest_anc(masked_seq_arr, curr_seq_ind,
                                            prev_time_info, return_multi = True, ignore_gap = True)
            
            #Get the info for the closest ancestor and save it
            min_dist_dict[curr_seq_ind] = min_hdist

            if isinstance(closest_anc_ind, np.int64):
                closest_anc_ind = [closest_anc_ind]


            closest_seqs_reversed[curr_seq_ind] = closest_anc_ind


            for i in closest_anc_ind:
                if i not in closest_seqs.keys():
                    closest_seqs[i] = [curr_seq_ind]
                else:
                    closest_seqs[i].append(curr_seq_ind)
    

    ###############################################################################
    ##################### Setting the Y Positions of Sequences ####################
    anc_seqs = closest_seqs.keys()
    day0_seqs = set(seq_info_df[seq_info_df['time_label'] == 'D0']['seq_index'])
    week1_seqs = set(seq_info_df[seq_info_df['time_label'] == 'W1']['seq_index'])
    week4_seqs = set(seq_info_df[seq_info_df['time_label'] == 'W4']['seq_index'])

    seq_list = []
    for curr_seq in anc_seqs:
        if isinstance(curr_seq,tuple):
            for i in curr_seq:
                if i in day0_seqs:
                    seq_list.append(i)
        else:
            if curr_seq in day0_seqs:
                seq_list.append(curr_seq)
    anc_seqs = seq_list

    # #I am manually moving a couple ancestors here for the main text visualization
    # anc_to_move_4 = anc_seqs.pop(4)
    # anc_to_move_0 = anc_seqs.pop(0)
    # anc_seqs = [anc_to_move_0, anc_to_move_4] + anc_seqs

    #Get each of the ancestor positions
    num_anc_seqs = len(anc_seqs)
    anc_positions = range(0, num_anc_seqs* ANC_SPREAD, ANC_SPREAD)
    curr_anc_pos_dict = dict(zip(anc_seqs, anc_positions))


    #Now, we'll get the position of each of the descendants at week 1
    desc_positions = {}
    desc_pos_order = []
    ghost_dict = {}
    anc_index = 0

    for key in anc_seqs:
        desc_list = closest_seqs[key]
        #loop through each of the descendants and if it isn't in the dictionary
        #add it
        for i in range(len(desc_list)):
            if desc_list[i] in desc_positions.keys():
                continue
            # If the ancestor is week 1 and the descendant is week 4, we'll handle
            # it in the next loop
            elif desc_list[i] in week4_seqs and key in week1_seqs:
                continue
            # add a ghost ancestor node if the ancestor is day 0 and the 
            # descendant is week 4
            elif desc_list[i] in week4_seqs and key in day0_seqs:
                if str(key) + "_ghost" not in desc_positions.keys():
                    desc_positions[str(key) + "_ghost"] = anc_index
                    desc_pos_order.append(str(key) + "_ghost")
                    ghost_dict[str(key) + "_ghost"] = desc_list[i]
                else: 
                    continue
            else:
                desc_positions[desc_list[i]] = anc_index
                desc_pos_order.append(desc_list[i])
            anc_index += 1


    #Finally, we will get the position of each of the descendants at week 4
    #first we'll add all of the sequences with week 4 ancestors and then we'll
    #add sequences which track back to a week 1 ancestor
    desc_2_positions = {}
    anc_index = 0

    #loop through the intermediate timepoint sequences
    for key in desc_pos_order:
        
        if isinstance(key, str) and 'ghost' in key:
            curr_ind = int(key.split('_')[0])
            curr_desc = closest_seqs[curr_ind]
        
        else:
            curr_ind = key
            if curr_ind not in closest_seqs.keys():
                continue
            curr_desc = closest_seqs[curr_ind]
      
        for desc in curr_desc:
            if desc not in desc_2_positions.keys():
                desc_2_positions[desc] = anc_index
                anc_index += 1


    ###########################################################################
    ############################ Make the Networkx Graph ######################
    fig, ax = plt.subplots(1,2, sharey = True)
    ax = ax.flatten()
    curr_ax = ax[1]
    my_graph = networkx.Graph()

    #Finally we will plot the networkx graphs as tripartite graphs
    print(ghost_dict)
    
    #First we will plot the day 0 nodes
    for curr_day0 in curr_anc_pos_dict.keys():
        my_graph.add_node(curr_day0, bipartite = 0)
    
    for curr_week1 in desc_positions.keys():
        my_graph.add_node(curr_week1, bipartite = 1)
    
    for curr_week4 in desc_2_positions.keys():
        if curr_week4 not in desc_positions.keys():
            my_graph.add_node(curr_week4, bipartite = 2)
    
    #Now we will add the edges
    for curr_week1 in desc_positions.keys():
        print('HAAA')
        if isinstance(curr_week1, str) and 'ghost' in curr_week1:
            node_int = int(curr_week1.split('_')[0])
            desc_edges = [node_int]
            #get all of the descendants of the ghost ancestor
            min_dist_opts = []
            min_dists = []
            for i in closest_seqs[node_int]:
              min_dists.append(min_dist_dict[i])
            min_dist = min(min_dists) * masked_seq_arr.shape[1]
            
        else:
            desc_edges = closest_seqs_reversed[curr_week1]
            print(min_dist_dict[curr_week1])
            min_dist = min_dist_dict[curr_week1] * masked_seq_arr.shape[1]

        #Get the edge weight and color for the current sequence
        if curr_week1 in closest_seqs.keys():
            descendant_info = closest_seqs[curr_week1]
            descendant_info = seq_info_df[seq_info_df['seq_index'].isin(descendant_info)]
            all_descendant_res_info = descendant_info['res_muts'].values
            all_descendant_res_info = [x[0] for x in all_descendant_res_info]
            all_descendant_res_info = list(set(all_descendant_res_info))


        elif isinstance(curr_week1, str) and 'ghost' in curr_week1:
            descendants = closest_seqs[int(curr_week1.split('_')[0])]
            if curr_week1 in ghost_dict.keys():
                descendants.append(ghost_dict[curr_week1])
            descendant_info = seq_info_df[seq_info_df['seq_index'].isin(descendants)]
            all_descendant_res_info = descendant_info['res_muts'].values
            all_descendant_res_info = [x[0] for x in all_descendant_res_info]
            all_descendant_res_info = list(set(all_descendant_res_info))

        edge_weight = 0.5
        color = WT_GRAY

        #added aug 11 am
        print(min_dist)
        if min_dist == 0:
            color = 'black'
        

        #Add the edge to the graph
        for edge in desc_edges:
            my_graph.add_edge(edge, curr_week1, weight = edge_weight, color = color)
    
    identical_nodes = []

    for curr_week4 in desc_2_positions.keys():
        #Get whether the current sequence has a S334N mutation
        curr_node_info = seq_info_df[seq_info_df['seq_index'] == curr_week4]
        edge_weight = 0.5
        min_dist = min_dist_dict[curr_week4] * masked_seq_arr.shape[1]

        edge_color = WT_GRAY

        if min_dist == 0:
            edge_color = 'black'

        desc_edges = closest_seqs_reversed[curr_week4]
        for edge in desc_edges:
            #Check if it has a ghost ancestor
            if str(edge) + "_ghost" in desc_positions.keys():
                if curr_week4 in week4_seqs:
                    my_graph.add_edge(str(edge) + "_ghost", curr_week4, weight = edge_weight, color = edge_color)
                elif min_dist == 0:
                    identical_nodes.append((str(edge) + "_ghost", curr_week4))

            elif edge in desc_positions.keys():
                my_graph.add_edge(edge, curr_week4, weight = edge_weight, color = edge_color)
        
    #Remove ghost nodes that are identical to week 1 nodes this can happen
    #if a week 0 sequence is identical to a week 1 sequence except for
    #the resistance mutations
    for curr_node in identical_nodes:
        if curr_node[0] in my_graph.nodes():
            my_graph.remove_node(curr_node[0])
            desc_positions.pop(curr_node[0])
    
    #Now adjust the positions of week 4 nodes so they are not dodging the deleted
    #nodes
    new_desc_dict = {}
    desc_positions_keys = list(desc_positions.keys())
    desc_pos_values = list(desc_positions.values())
    desc_pos_values.sort()
    desc_reverse = {v: k for k, v in desc_positions.items()}
    desc_pos_order = range(len(desc_positions_keys))
    for i in desc_pos_order:
        new_desc_dict[desc_reverse[desc_pos_values[i]]] = i
    desc_positions = new_desc_dict

    #Now color and label the nodes
    node_sizes = []
    node_colors = []
    edge_colors = []
    label_dict = {}
    alpha_list = []
    
    for curr_node in my_graph.nodes():

        if isinstance(curr_node, str) and curr_node.split('_')[-1] == 'ghost':
            node_colors.append('white')
            #edge_colors.append('black')
            alpha_list.append(0)
            node_sizes.append(0)
            label_dict[curr_node] = ''
            continue
        node_info = seq_info_df[seq_info_df['seq_index'] == curr_node]
        node_size = labeled_counts[labeled_counts['seq_index'] == curr_node]['num_identical'].values[0] * 2 + 5

        muts = node_info['res_muts'].values[0]
        muts = muts[0]

        if not muts:
            my_color = WT_GRAY
        else:
            my_color = dataset_metadata.RESISTANCE_HUE_DICT[muts]
        node_colors.append(my_color)
        node_sizes.append(node_size)
        #edge_colors.append('white')
        alpha_list.append(0.75)

        if curr_node in anc_seqs:
            label_dict[curr_node] = ''
        else:
            label_dict[curr_node] = str(round(min_dist_dict[curr_node]* masked_seq_arr.shape[1]))


    #Now we will plot the graph
    nodes = my_graph.nodes()

    #Day 0
    nodes_0  = set([n for n in nodes if  my_graph.nodes()[n]['bipartite']==0])

    #Week 1
    nodes_1  = set([n for n in nodes if  my_graph.nodes()[n]['bipartite']==1])

    #Week 4
    nodes_2  = set([n for n in nodes if  my_graph.nodes()[n]['bipartite']==2])

    #We will plot the ancestors at x = 1 and the descendants at x = 2
    pos = dict()
    pos.update( (n, (1, curr_anc_pos_dict[n])) for i, n in enumerate(nodes_0) ) # put nodes from day 0 at x=1
    pos.update( (n, (2, desc_positions[n] * DESC_SPREAD * 2)) for i, n in enumerate(nodes_1) ) # put nodes from week 1 at x=2
    pos.update( (n, (3, desc_2_positions[n] * DESC_SPREAD)) for i, n in enumerate(nodes_2) ) # put nodes from week 4 at x=3

    #for some reason network x only allows transparent nodes if there are no node labels
    networkx.draw(my_graph, pos=pos, ax = curr_ax, node_color = node_colors, with_labels = False, font_size = FONTSIZE, 
                  font_color = 'black', alpha = alpha_list, node_size = node_sizes, width = 0)

    for edge in my_graph.edges():
        weight = my_graph.edges[edge[0], edge[1]]['weight']
        color = my_graph.edges[edge[0], edge[1]]['color']
        networkx.draw_networkx_edges(my_graph, pos, edgelist=[edge], width=weight, edge_color = color,
                                     ax = curr_ax)

    curr_ax.annotate(r'$\uparrow$' +  "\nDay 0 \nancestor", xy = (0.075, -0.01), xycoords = 'axes fraction', ha = 'center', fontsize = FONTSIZE)
    curr_ax.annotate(r'$\uparrow$' +  "\nWeek 1 \ndescendant", xy = (0.49, -0.01), xycoords = 'axes fraction', ha = 'center', fontsize = FONTSIZE)
    curr_ax.annotate(r'$\uparrow$' +  "\nWeek 4 \ndescendant", xy = (0.91, -0.01), xycoords = 'axes fraction', ha = 'center', fontsize = FONTSIZE)
    curr_ax.set_title('Sequence matching network', fontname = 'Arial', fontsize = FONTSIZE, x = 0.5, y = 0.9)
    y_min = curr_ax.get_ylim()[0]

    # curr_ax.set_ylim(y_min, max(anc_positions)+ 100)

    # Make legend
    # for n in range(150, np.max(node_sizes) + 100, 500):
    #     print(n)
    #     if n > 0 and n % 50 == 0:
    #         plt.plot([], [], 'o', markersize = int(np.sqrt(n)), label = f"{int((n-100)/50)}", color = 'black', alpha = 0.5)
    #         plt.legend(bbox_to_anchor=(0.5, 0.5), frameon = True,
    #                    title = 'Identical\nSequences')



    ###############################################################################
    ############################ Highlighter Plots ################################
    ###############################################################################



    #Try first just plotting all of the day 0 sequences on a line/ scatter plot
    #add all of the marker coordinates to lists
    Ax_coords = []
    Ay_coords = []
    Cx_coords = []
    Cy_coords = []
    Gx_coords = []
    Gy_coords = []
    Tx_coords = []
    Ty_coords = []

    for i in range(len(anc_seqs)):
        #Get the segregating mutations for that sequence which are already in a set
        curr_key = anc_seqs[i]
        curr_set = seg_set_dict[index_to_name[curr_key]]

        seq_y = anc_positions[i]

        for curr_mut in curr_set:
            if curr_mut[1] == 'A':
                Ax_coords.append(curr_mut[0])
                Ay_coords.append(seq_y)
            elif curr_mut[1] == 'C':
                Cx_coords.append(curr_mut[0])
                Cy_coords.append(seq_y)
            elif curr_mut[1] == 'G':
                Gx_coords.append(curr_mut[0])
                Gy_coords.append(seq_y)
            elif curr_mut[1] == 'T':
                Tx_coords.append(curr_mut[0])
                Ty_coords.append(seq_y)
            else:
                print('Invalid base found:' + curr_mut[1])
        
        ax[0].axhline(y = seq_y, color = WT_GRAY, zorder = 1)

    hxb2_res_min = min([x[0] for x in hxb2_res_positions])
    hxb2_res_max = max([x[1] for x in hxb2_res_positions])

    res_tup_hxb2 = (hxb2_res_window_min, hxb2_res_window_max)

    res_out_hxb2 = []
    for curr_locus in res_tup_hxb2:
        whole_part = math.floor(curr_locus)
        frac_part = curr_locus - whole_part
        whole_part = math.ceil(int(whole_part)/3)
        frac_part = frac_part/3
        aa_locus = whole_part + frac_part
        res_out_hxb2.append(int(aa_locus))

    hxb2_res_window_min_AA = res_out_hxb2[0]
    hxb2_res_window_max_AA = res_out_hxb2[1]


    ax[0].scatter(Ax_coords, Ay_coords, color = my_colors['A'], label = 'A', marker='v', s = 5, zorder = 2)
    ax[0].scatter(Cx_coords, Cy_coords, color = my_colors['C'], label = 'C', marker='v', s = 5, zorder = 2)
    ax[0].scatter(Gx_coords, Gy_coords, color = my_colors['G'], label = 'G', marker='v', s = 5, zorder = 2)
    ax[0].scatter(Tx_coords, Ty_coords, color = my_colors['T'], label = 'T', marker='v', s = 5, zorder = 2)
    # ax[0].legend(title = 'Minor Allele', title_fontsize = 12, fontsize = 12, ncol = 4, loc = 'lower left')
    ax[0].set_xlim(arr_res_window_min-3, arr_res_window_max)
    ax[0].set_axis_off()
    # ax[0].set_xlabel('HXB2 position')

    #Now annotate the resistance positions

        #now annotate hxb2 positions
    annotation_y = min(anc_positions) - 1
    ax[0].annotate(r'$\uparrow$'+ "\n"  + str(hxb2_res_window_min_AA), xy = (-0.05, -0.01), xycoords = 'axes fraction', fontsize = FONTSIZE, fontname = 'Arial')
    ax[0].annotate(r'$\uparrow$' + "\n" + str(hxb2_res_window_max_AA), xy = (0.95, -0.01), xycoords = 'axes fraction', fontsize = FONTSIZE, fontname = 'Arial')
    ax[0].annotate("HXB2 coordinate AA", xy = (0.5, -0.01), xycoords = 'axes fraction', fontsize = FONTSIZE, fontname = 'Arial', ha = 'center')
    #curr_ax.annotate('HXB2 Coordinate AA', xy = (arr_res_window_min + (arr_res_window_max - arr_res_window_min)/2 + 9, annotation_y - 0.2), xycoords = 'data', fontsize = FONTSIZE, ha = 'center')
    #curr_ax.annotate(r'$\uparrow$' + "\n" + 'Escape loci (325-334)', xy = (arr_res_window_min + (arr_res_window_max - arr_res_window_min)/2 + 9, annotation_y), xycoords = 'data', fontsize = FONTSIZE, ha = 'center', fontname = 'Arial')    


    # ax[0].annotate('HXB2 Position', xy = (0.5, 0.05), xycoords = 'axes fraction', ha = 'center', fontname = 'Arial', fontsize = FONTSIZE)
    # ax[0].annotate(r'$\uparrow$'+ "\n"  + str(hxb2_res_window_min), xy = (0.075, -0.01), xycoords = 'axes fraction', fontname = 'Arial', fontsize = FONTSIZE)
    # ax[0].annotate(r'$\uparrow$' + "\n" + str(hxb2_res_window_max), xy = (arr_res_window_max-2.5, -45), xycoords = 'data', fontname = 'Arial', fontsize = FONTSIZE)
    # break_string = str(hxb2_res_min) + r' // '+ str(hxb2_res_max)
    # ax[0].annotate(break_string, xy = (hxb2_res_min, -25), xycoords = 'data', fontname = 'Arial', fontsize = FONTSIZE)
    ax[0].set_title('Day 0 segregating sites', fontname = 'Arial', fontsize = FONTSIZE, x = 0.5, y = 0.9)

 
    fig.suptitle(SAMPLE_PAR, fontname = 'Arial', fontsize = FONTSIZE, x = 0.5, y = 0.95)
    plt.subplots_adjust(left=0.1, right=0.95, top=1, bottom=0.05, hspace=0.2, wspace=0.2)


    plt.savefig(outDir + SAMPLE_PAR + 'Full_Anc_Network_Test.png', dpi = 300)
    plt.close()

    #Now make a legend for the node sizes
    node_legend_list = [1, 5, 10, 25, 50]
    node_legend_list = [n * 2 + 5 for n in node_legend_list]
    for n in node_legend_list:
        print(n)
        plt.plot([], [], 'bo', markersize = np.sqrt(n), label = str(int((n-5)/2)), color = 'black', alpha = 0.5)
        plt.legend(bbox_to_anchor=(0, 1), frameon = False, loc = 'upper left',
                   title = 'Number of identical sequences', fontsize = FONTSIZE, ncol = len(node_legend_list))
    plt.savefig(outDir + SAMPLE_PAR + 'Dot_size_legend.png', dpi = 300)
    plt.close()