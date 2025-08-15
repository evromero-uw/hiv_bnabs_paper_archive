import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# The code in this file is for making the supplemental histograms of MPL
# selection coefficients for all participants.

mpl_dir = '../../../../results/3BNC/MPL/%s/analysis/'
coeff_file = 'total-selection-%s.csv'
out_dir = '../../../../results/pub_figs_07-2025/supplemental_figs/all_coeff_hist/'

par_list = ['2C1', '2C5', '2D1', '2E1', '2E2', '2E3', '2E4', '2E5', '2E7']

SEL_THRESH = 0.044

TEXT_FONTSIZE = 6
params = {
          'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6, 'font.family': 'sans-serif',
          'font.sans-serif': 'Arial'}
rcParams.update(params)


fig, axs = plt.subplots(3, 3, figsize = (7, 7), sharex=True, sharey=True)
axs = axs.flatten()
# Now for each participant, I will get the selection coefficients at these
# sites and plot them.
for i, par in enumerate(par_list):
    curr_ax = axs[i]

    # Get the selection coefficients for this participant
    curr_mpl_dir = mpl_dir % par
    curr_coeff_file = coeff_file % par
    mpl_coeff = pd.read_csv(curr_mpl_dir + curr_coeff_file)

    # Label the selection coefficients higher than the threshold
    mpl_coeff['label'] = mpl_coeff['s_MPL'].apply(lambda x: 'Above threshold' if x > SEL_THRESH else 'Below threshold')

    #plot the in vs out of range values
    g = sns.histplot(x='s_MPL', data=mpl_coeff, stat='count', multiple='stack',
                     hue = 'label', palette=['lightgray', 'red'],
                binwidth = 0.005, ax=curr_ax)
    
    curr_ax.axvline(SEL_THRESH, color='black', linestyle='--', label='Threshold = %.3f' % SEL_THRESH,
                    linewidth=0.5)

    #make y axis log scale
    curr_ax.set_yscale('log')
    curr_ax.set_xlim(-0.04, 0.15)
    curr_ax.set_xlabel('MPL selection coefficient\n (For nucleotide mutation)')
    curr_ax.set_title(par)
    
    #turn off the axis legend
    if i == 0:
        curr_ax.legend(loc='upper right')
    else:
        curr_ax.get_legend().remove()   




plt.suptitle('MPL selection coefficients', fontsize=TEXT_FONTSIZE, y=0.97)
plt.subplots_adjust(top = 0.92)
plt.savefig(out_dir + 'all_mpl_coefficients_hist.png', dpi=300)


plt.close()
