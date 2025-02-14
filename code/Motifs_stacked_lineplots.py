#!/usr/bin/env python3

# Modifying previous code to stack all models lineplots together
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Func
# Modify from get_ave_dsvm_per_bp_around_motif
def stacked_lineplot(matrices, colors, file_out, xrange=[-30,30,5]):
    # set up plot
    plt.figure(figsize=(9,4))
    for file in matrices:
        #load file
        model = file.rsplit('_neg', 1)[0].rsplit('.',1)[1]
        if model == 'Macular_RPE':
            sub_model = model
        else:
            sub_model = model.rsplit('11mers_',1)[1].rsplit('_no',1)[0]
        motif_df = pd.read_csv(file, index_col=0)
        # process file
        cols = motif_df.columns
        motif_df[cols] = motif_df[cols].apply(pd.to_numeric, errors='coerce')
        motif_df_means = motif_df.mean(axis = 0, skipna = True)
        # add to plot, using color from colors dict
        color = colors.get(model, "gray")
        if color == "gray":
            print(f'No color found for {model}, using gray')
        sns.lineplot(data=motif_df_means, color=color, label=sub_model)
    # Fix x axis
    x_ticks = plt.gca().get_xticks()
    n = x_ticks[-1]
    custom_ticks_centered = list(range(xrange[0], (xrange[1] + 1), xrange[2]))
    custom_ticks_original = [tick + (n/2) for tick in custom_ticks_centered]

    #centered_labels = [round((x - (n/2)), 1) for x in x_ticks]
    #plt.gca().set_xticks(x_ticks)
    plt.gca().set_xticks(custom_ticks_original)
    #plt.gca().set_xticklabels(centered_labels)
    plt.gca().set_xticklabels(custom_ticks_centered)
    plt.tick_params(axis='x',labelsize=16)
    plt.tick_params(axis='y',labelsize=16)
    plt.xticks(rotation=45)
    # Add legend
    plt.legend(title="Models", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    # save plot
    plt.savefig(file_out)
    # close plot
    plt.close()

def main():
    working_directory = '/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/dsvm/motifs/'
    os.chdir(working_directory)
    motifs = ['SOX2','LHX2', 'ONECUT1', 'OTX2', 'ASCL1', 'CRX']
    models = ['all_11mers_Early_Progenitors_noImmIn','all_11mers_Late_Progenitors_noImmIn',
        'all_11mers_Mature_Mullers_noMatIn', 'all_11mers_Photoreceptor_Bipolar_Precursors_noImmIn',
        'all_11mers_Developing_Bipolars_noImmIn','all_11mers_Mature_Bipolars_noMatIn',
        'all_11mers_Developing_Rods_noImmIn','all_11mers_Mature_Rods_noMatIn',
        'all_11mers_Developing_Cones_noImmIn','all_11mers_Mature_Cones_noMatIn',
        'all_11mers_AC_HC_GC_precursors_noImmIn','all_11mers_Developing_Amacrines_noImmIn',
        'all_11mers_Mature_Amacrines_noMatIn','all_11mers_Developing_Horizontals_noImmIn',
        'all_11mers_Mature_Horizontals_noMatIn','all_11mers_Developing_Ganglions_noImmIn','Macular_RPE']
    colors = {'all_11mers_Early_Progenitors_noImmIn':'#5DADE2','all_11mers_Late_Progenitors_noImmIn':'#4CAF50',
        'all_11mers_Mature_Mullers_noMatIn':'#EC7063', 'all_11mers_Photoreceptor_Bipolar_Precursors_noImmIn':'#E67E22',
        'all_11mers_Developing_Bipolars_noImmIn':'#EC407A','all_11mers_Mature_Bipolars_noMatIn':'#AD1457',
        'all_11mers_Developing_Rods_noImmIn':'#7D3C98','all_11mers_Mature_Rods_noMatIn':'#7D3C98',
        'all_11mers_Developing_Cones_noImmIn':'#F44336','all_11mers_Mature_Cones_noMatIn':'#B71C1C',
        'all_11mers_AC_HC_GC_precursors_noImmIn':'#6E2C00','all_11mers_Developing_Amacrines_noImmIn':'#2ECC71',
        'all_11mers_Mature_Amacrines_noMatIn':'#145A32','all_11mers_Developing_Horizontals_noImmIn':'#21618C',
        'all_11mers_Mature_Horizontals_noMatIn':'#154360','all_11mers_Developing_Ganglions_noImmIn':'#F4D03F','Macular_RPE':'#A9A9A9'}
    for motif in motifs:
        files_in = []
        for model in models:
            new_in = f'dsvm.{model}_neg_scores.{motif}.ext25.csv'
            files_in.append(new_in)
        stacked_lineplot(files_in, colors, f'lineplots/{motif}_neg_scores_lineplot.pdf')

if __name__=='__main__':
    main()
