#!/usr/bin/python3

#Recreate corbo/cohen heatmaps from gnb3 data and dsvm

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



# Plotter function

def make_heatmap(df, hm_ID, values, levels, val_range, file_out, transform = False):
    # filter df by hm_ID
    plot_df = df[df['heatmap'] == hm_ID].copy()
    if transform:
        plot_df['log2_values'] = np.log2(plot_df[values])
        # reshape to matrix by cluster and label with the cells being average value
        plot_mx = plot_df.pivot(index='cluster', columns='label', values='log2_values')
    else:
        plot_mx = plot_df.pivot(index='cluster', columns='label', values=values)
    plot_mx = plot_mx.reindex(columns=levels)
    # plot rows = cell type, columns = variants/label
    palette = sns.color_palette('RdBu_r', 256)
    plt.figure(figsize=(25, 6))
    sns.heatmap(plot_mx, annot=False, cmap=palette, vmin = val_range[0], vmax = val_range[1])
    plt.xticks(rotation=45)
    plt.savefig(file_out, format='pdf', bbox_inches='tight')
    plt.close()


working_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/'

# read in files
    # dSVM
colnames = ['name','dSVM']
BPC_dsvm = pd.read_csv(f'{working_dir}deltaSVM_BPC_predict_GNB.txt', sep='\t', names=colnames)
BPC_dsvm['cluster'] = 'Bipolar'
Rod_dsvm = pd.read_csv(f'{working_dir}deltaSVM_Rods_predict_GNB.txt', sep='\t', names=colnames)
Rod_dsvm['cluster'] = 'Rod'
MG_dsvm = pd.read_csv(f'{working_dir}deltaSVM_MG_predict_GNB.txt', sep='\t', names=colnames)
MG_dsvm['cluster'] = 'Glia'
# Xn
Gnb3_xn = pd.read_csv(f'{working_dir}GNB_Variant_Expression.txt', sep='\t')
    # Table labels
Table_labels = pd.read_csv(f'{working_dir}GNB_table_labels.txt', sep='\t')
# merge all tables as appropriate
Gnb3_dsvm = pd.concat([BPC_dsvm, Rod_dsvm, MG_dsvm])
Gnb3_dsvm = Gnb3_dsvm.merge(Table_labels, on='name')


Gnb3_xn = Gnb3_xn.merge(Table_labels, on='name')
Gnb3_xn = Gnb3_xn[Gnb3_xn['cluster'] != 'Neuron']

# From Table labels get heatmap levels
Label_dict = {}
with open(f'{working_dir}GNB_table_labels.txt', 'r') as labels:
    for line in labels:
        if line.startswith('name'):
            continue
        line = line.strip()
        name, heatmap, label = line.split('\t')
        if heatmap in Label_dict:
            Label_dict[heatmap].append(label)
        else:
            Label_dict[heatmap] = [label]

# build plots
for key in Label_dict:
    keyname = key.replace(' ', '_')
    make_heatmap(Gnb3_xn, key, 'AVG', Label_dict[key], [-3,3],f'{working_dir}Gnb3_xn_{keyname}.pdf', transform=True)
    make_heatmap(Gnb3_dsvm, key, 'dSVM', Label_dict[key], [-25, 25], f'{working_dir}Gnb3_dsvm_{keyname}.pdf')
#make_heatmap(Gnb3_xn, 'E-Box mutants', f'{working_dir}Gnb3_xn_Ebox_mutants.pdf')
