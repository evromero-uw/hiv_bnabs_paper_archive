import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from matplotlib import rcParams

#In this file I am making the supplemental plots of v5 length and glycosylation

inDir = '../../../../results/pub_figs_07-2025/supplemental_figs/var_char_summary/'
var_char_df = pd.read_csv(inDir + 'var_char_summary.csv')
var_char_df = var_char_df.drop(columns=['Unnamed: 13', 'Unnamed: 14'])

params = {'figure.figsize':(7.05, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)

par_colors = {
    '1HB3': '#A6CEE3', '1HC2': '#1F78B4',
    '1HC3': '#B2DF8A', '1HD1': '#33A02C',
    '1HD4K': '#FB9A99','1HD5K': '#E31A1C',
    '1HD6K': '#FDBF6F','1HD7K': '#FF7F00',
    '1HD9K': '#CAB2D6','1HD10K': '#6A3D9A',
    '1HD11K': '#B15928',
    '2C5': '#A6CEE3','2C1': '#1F78B4',
    '2D1': '#B2DF8A','2E1': '#33A02C',
    '2E2': '#FB9A99','2E3': '#E31A1C',
    '2E4': '#FDBF6F','2E5': '#FF7F00',
    '2E7': '#CAB2D6'}

plot_1_list = ['NGlycos.V1', 'NGlycos.V2', 'NGlycos.V3',
              'NGlycos.V4', 'NGlycos.V5', 'TotalPNGs']
plot_1_rename = {
    'NGlycos.V1': 'Average # of PNGs V1',
    'NGlycos.V2': 'Average # of PNGs V2',
    'NGlycos.V3': 'Average # of PNGs V3',
    'NGlycos.V4': 'Average # of PNGs V4',
    'NGlycos.V5': 'Average # of PNGs V5',
    'TotalPNGs': 'Average # of PNGs Total'}

plot_2_list = ['Length.V1', 'Length.V2', 
              'Length.V4', 'Length.V5']
plot_2_rename = {
    'Length.V1': 'Average length of V1',
    'Length.V2': 'Average length of V2',
    'Length.V4': 'Average length of V4',
    'Length.V5': 'Average length of V5'}

#I want to make plots of each characteristic
print(var_char_df['TimePt'].unique())
time_list = ['D0', 'W1', 'W4', 'W8', 'W12']
var_char_df = var_char_df[var_char_df['TimePt'].isin(time_list)]
var_char_df['time_label_int'] = [int(x[1:]) for x in var_char_df['TimePt']]
print(var_char_df.head())

for curr_bnab in var_char_df['bnAb'].unique():
    curr_df = var_char_df[var_char_df['bnAb'] == curr_bnab]
    
    # Make a set of plots
    fig, axs = plt.subplots(2, 3, sharex=True, figsize=(7, 4))
    axs = axs.flatten()
    fig.suptitle(curr_bnab, fontsize = 10)

    for i, curr_col in enumerate(plot_1_list):
        curr_ax = axs[i]
        sns.lineplot(curr_df, x='time_label_int', y=curr_col, hue='ID',
                     palette=par_colors, ax=curr_ax, legend=False)
        curr_ax.set_xlabel('Time (weeks)')
        curr_ax.set_ylabel(plot_1_rename[curr_col])
        #Make the x-axis ticks only integers
        curr_ax.xaxis.set_major_locator(mticker.MultipleLocator(1))

        if i == len(plot_1_list) - 1:
            legend_handles = []
            legend_labels = curr_df['ID'].unique().tolist()
            for participant in curr_df['ID'].unique():
                legend_handles.append(plt.Line2D([0], [0], marker=None, 
                                                 color= par_colors[participant],
                                                label=participant,
                                                markersize=5))
            curr_ax.legend(loc='upper left', bbox_to_anchor=(1, 1.5), title='Participant',
                                            handles=legend_handles, fontsize=6, title_fontsize=6)
            curr_ax.set_xticks(curr_df['time_label_int'].unique())

    plt.subplots_adjust(hspace=0.3, wspace=0.3, top=0.9)
    plt.savefig(inDir + f'{curr_bnab}_var_char_plot_glycos.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Make a set of plots
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(4.5, 3))
    axs = axs.flatten()
    fig.suptitle(curr_bnab, fontsize = 10)

    for i, curr_col in enumerate(plot_2_list):
        curr_ax = axs[i]
        sns.lineplot(curr_df, x='time_label_int', y=curr_col, hue='ID',
                     palette=par_colors, ax=curr_ax, legend=False)
        curr_ax.set_xlabel('Time (weeks)')
        curr_ax.set_ylabel(plot_2_rename[curr_col])

        if i == len(plot_2_list) - 1:
            legend_handles = []
            legend_labels = curr_df['ID'].unique().tolist()
            for participant in curr_df['ID'].unique():
                legend_handles.append(plt.Line2D([0], [0], marker=None, 
                                                 color= par_colors[participant],
                                                label=participant,
                                                markersize=5))
            curr_ax.legend(loc='upper left', bbox_to_anchor=(1, 1.5), title='Participant',
                                            handles=legend_handles, fontsize=6, title_fontsize=6)
            curr_ax.set_xticks(curr_df['time_label_int'].unique())
    
    plt.subplots_adjust(hspace=0.3, wspace=0.3, top=0.9)
    plt.savefig(inDir + f'{curr_bnab}_var_char_plot_length.png', dpi=300, bbox_inches='tight')
    plt.close()