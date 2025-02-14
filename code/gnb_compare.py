#!/usr/bin/python3

# Massage data and compare GNB3 predictions to sc expression
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Make plotter function
def plot_comparisons(df, file_out):
    corr_coef, p_value = pearsonr(df['dSVM'], df['AVG'])
    slope, intercept = np.polyfit(df['dSVM'], df['AVG'], 1)
    line_eq = f'y = {slope:.2f}x + {intercept:.2f}'
    plt.figure(figsize=(8,6))
    sns.regplot(x=df['dSVM'], y=df['AVG'], scatter_kws={'s':80})
    plt.text(0.05, 0.95, f'{line_eq}\nr={corr_coef:.3f}', transform=plt.gca().transAxes, fontsize=14, verticalalignment='top')
    plt.title('GNB3 dSVM vs Expression')
    plt.xlabel('dSVM score')
    plt.ylabel('Average Expression')
    plt.savefig(file_out, format='pdf')



# Read in deltasvm files
Rod_dSVM = pd.read_csv('/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/deltaSVM_Rods_predict_GNB.txt', sep='\t',names=['name','dSVM'])
BPC_dSVM = pd.read_csv('/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/deltaSVM_BPC_predict_GNB.txt', sep='\t', names=['name', 'dSVM'])
MG_dSVM = pd.read_csv('/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/deltaSVM_MG_predict_GNB.txt', sep='\t', names=['name', 'dSVM'])

# Read in expression data
all_expression = pd.read_csv('/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/GNB_Variant_Expression.txt', sep='\t', header=0, usecols=['name','cluster', 'AVG'])
# Split expression data by cell class
Rod_xn = all_expression.loc[all_expression['cluster'] == 'Rod']
BPC_xn = all_expression.loc[all_expression['cluster'] == 'Bipolar']
MG_xn = all_expression.loc[all_expression['cluster'] == 'Glia']
# Merge  dataframes by enhancer mod

Rod_data = pd.merge(Rod_xn, Rod_dSVM, how= 'left', on='name')
BPC_data = pd.merge(BPC_xn, BPC_dSVM, how='left', on='name')
MG_data = pd.merge(MG_xn, MG_dSVM, how='left', on='name')

# plot for each cell class
plot_comparisons(Rod_data, '/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/Rod_scatter.pdf')
plot_comparisons(BPC_data, '/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/BPC_scatter.pdf')
plot_comparisons(MG_data, '/data/hps/assoc/private/cherrylab/user/lvand1/corbo_gnb3/MG_scatter.pdf')
