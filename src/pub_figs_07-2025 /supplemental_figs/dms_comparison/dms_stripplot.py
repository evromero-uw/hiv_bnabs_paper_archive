import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

#This code makes the strip plots of escape scores at sites identified 
#by our selection scans in 0,1 or 2 participants

params = {'figure.figsize':(2, 2), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial', 'figure.titlesize': 6, 'axes.titlesize': 6}
rcParams.update(params)

ALPHABET = '-abcdefghijklmnopqrstuvwxyz'

#Input folder with the selection scan data
shared_dir = "../../../../results/pub_figs_07-2025/supplemental_tables/"


#Input folders with the DMS data
dms_indir = "../../../../data/Radford2025/"
TIMES_SEEN_FILTER_BF520 = 3
TIMES_SEEN_FILTER_TRO11 = 2

#Here is the output folder
outFolder = 'dms_comparison/'
outDir = '../../../../results/pub_figs_07-2025/supplemental_figs/'

def get_hxb2_coord_string(hxb2_coord):
    """
    Convert a HXB2 coordinate to a string representation.
    """
    if hxb2_coord % 1 == 0:
        
        return str(int(hxb2_coord))
    else:
        remainder = hxb2_coord % 1
        remainder = int(np.round(remainder * 1000))  # Convert to integer for display
        alphabet_remainder = ALPHABET[remainder]

        return f"{int(hxb2_coord)}{alphabet_remainder}"

###############################################################################
####################### Load and filter the DMS data ##########################
bf520_3BNC = pd.read_csv(dms_indir + "BF520/3BNC117_mut_effect.csv")
bf520_3BNC['strain'] = 'BF520'
bf520_3BNC['antibody'] = '3BNC117'
bf520_3BNC = bf520_3BNC[bf520_3BNC['times_seen'] >= TIMES_SEEN_FILTER_BF520]


tro11_3BNC = pd.read_csv(dms_indir + "TRO11/3BNC117_mut_effect.csv")
tro11_3BNC['strain'] = 'TRO11'
tro11_3BNC['antibody'] = '3BNC117'
tro11_3BNC = tro11_3BNC[tro11_3BNC['times_seen'] >= TIMES_SEEN_FILTER_TRO11]

# Concatenate the dataframes to gather all the DMS data together
dms_info_3BNC = pd.concat([bf520_3BNC, tro11_3BNC])
dms_info_3BNC = dms_info_3BNC[['site', 'wildtype', 'mutant', 'mutation',
                               'escape_mean', 'escape_median', 'escape_std',
                               'strain', 'antibody', 'times_seen']]

# Group the DMS data by HXB2 site and get the max escape value for each site
dms_info_3BNC = dms_info_3BNC.groupby(['site', 'strain']).agg(
                                        {'escape_mean': 'max'}).reset_index()

# Take the max escape value for each site across strains
dms_info_3BNC = dms_info_3BNC.groupby('site').agg(
    {'escape_mean': 'max', 'strain': 'first'}).reset_index()

###############################################################################
############################ Load our shared sites ############################
shared_sites = pd.read_csv(shared_dir + '3BNC117_all_data.csv')

#label each site with the HXB2 coordinate string
shared_sites['hxb2_coord_string'] = shared_sites['hxb2_coord_AA'].apply(get_hxb2_coord_string)
shared_sites = shared_sites.drop_duplicates(subset=['hxb2_coord_AA', 'participant'])

#Count how many sites are shared by how many participants
shared_sites = shared_sites.groupby('hxb2_coord_string').agg(
    {'participant': 'count'}).reset_index()
# shared_sites = shared_sites[shared_sites['participant'] >= 2]

positive_sel = []
num_participants = []

for index, row in dms_info_3BNC.iterrows():
    # Get the escape score for the current site
    coord_str = row['site']
    curr_score = shared_sites[shared_sites['hxb2_coord_string'] == coord_str]
    
    if len(curr_score) > 0:
        par_count = curr_score['participant'].values[0]
        num_participants.append(par_count)
        positive_sel.append(True)
    else:
        num_participants.append(0)
        positive_sel.append(False)

dms_info_3BNC['num_participants'] = num_participants
dms_info_3BNC['positive_sel'] = positive_sel
print(dms_info_3BNC.head())
print(dms_info_3BNC[dms_info_3BNC['positive_sel'] == True].head())

dms_info_3BNC.rename(columns={'num_participants': '# of participants'},
                      inplace=True)

#Maybe I should just compare the shared sites to the non shared sites?
fig, ax = plt.subplots(1, 1)
sns.stripplot(x='# of participants', y='escape_mean', data=dms_info_3BNC, hue ='# of participants',
              order=[0, 1, 2], palette={0: 'lightgrey', 1: 'blue', 2: 'red'},
              alpha=0.5, jitter=True, edgecolor='black')

sns.boxplot(x='# of participants', y='escape_mean', data=dms_info_3BNC, hue='# of participants',
            order=[0, 1, 2], palette={0: 'black', 1: 'blue', 2: 'red'},
            fliersize=0, fill=False, linewidth=1, zorder = 4)
#Remove the legend
plt.legend([],[], frameon=False)
plt.subplots_adjust(left=0.2, right=0.95, bottom=0.2, top=0.9)
plt.xlabel('# of participants with given site\nidentified by selection scan')
plt.ylabel('Maximum escape score at given site')
plt.title('3BNC117')

plt.savefig(outDir + outFolder + 'dms_stripplot_3BNC117.png', dpi=300)


# # Now we can plot the escape scores for the shared sites
# fig, ax = plt.subplots(1, 2, sharey=True, sharex=True)

# # Plot for BF520
# bf520_df = dms_info_3BNC[dms_info_3BNC['strain'] == 'BF520'].copy()

# #For some reason, pandas kept thinking my data was wideform
# bf520_df = pd.DataFrame( {'escape_mean': bf520_df['escape_mean'],
#                           'num_participants': bf520_df['num_participants']})
# sns.histplot(x = 'escape_mean', hue = 'num_participants', data=bf520_df,
#              ax = ax[0], palette={0: 'lightgrey', 1: 'blue', 2: 'red'},
#              binwidth=0.1, multiple='stack')
# ax[0].set_title('BF520')

# # Plot for TRO11
# tro11_df = dms_info_3BNC[dms_info_3BNC['strain'] == 'TRO11']

# #For some reason, pandas kept thinking my data was wideform
# tro11_df = pd.DataFrame( {'escape_mean': tro11_df['escape_mean'],
#                           'num_participants': tro11_df['num_participants']})
# sns.histplot(x = 'escape_mean', hue = 'num_participants', data=tro11_df,
#              ax = ax[1], palette={0: 'lightgrey', 1: 'blue', 2: 'red'},
#              binwidth=0.1, multiple='stack')
# ax[1].set_title('TRO11')
# plt.yscale('log')

# plt.show()
