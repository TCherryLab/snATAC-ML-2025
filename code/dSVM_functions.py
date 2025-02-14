#!/usr/bin/env python3


# A collection of useful functions intended to be genetally applicable

# dependencies
import os
import re
import pandas as pd
import dask as dd
import pybedtools
import subprocess
from Bio import motifs
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr
from functools import lru_cache

# Getting sequences (reference and variant)

def get_fasta(bed, ref_fa='/data/hps/assoc/public/bioinformatics/annotations/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'):
    position_seq = bed.sequence(fi=ref_fa)
    fasta_str = open(position_seq.seqfn).read()
    fasta_list = fasta_str.split('\n')[1:-1] # skip header and empty row
    return ''.join(fasta_list)

@lru_cache(maxsize=None)
def get_reference_sequence(chrom, position, ref, slop=11, ref_fa='/data/hps/assoc/public/bioinformatics/annotations/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'):
    start = int(position) - slop
    stop = int(position) + slop + len(ref) - 1
    bed = pybedtools.BedTool(f'{chrom}\t{start}\t{stop}', from_string=True)
    sequence_str = get_fasta(bed, ref_fa).upper()
    return sequence_str

def get_alt_sequence(sequence, ref, alt, slop=11):
    if len(ref) > 1:
        pos_end = slop + len(ref) - 1
        alt_seq = sequence[:slop-1] + alt.upper() + sequence[pos_end:]
    else:
        alt_seq = sequence[:slop-1] + alt.upper() + sequence[slop:]
    return alt_seq

## Optional wrapper for applying above functions to dataframe
## This may need to be built fresh depending on dataframe structure
def generate_sequences(row, slop=11, ref_fa='/data/hps/assoc/public/bioinformatics/annotations/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'):
    ref_seq = get_reference_sequence(row['Chr'], row['Pos'], row['Ref'], slop=slop, ref_fa=ref_fa)
    alt_seq = get_alt_sequence(ref_seq, row['Ref'], row['Alt'], slop=slop)
    return pd.Series([ref_seq, alt_seq])


# deltaSVM functions

def write_fa(df, sequence_col, file_out, header_col='header'):
    with open(file_out, 'w') as fa:
        for index, row in df.iterrows():
            fa.write(f'>{row[header_col]}\n')
            fa.write(f'{row[sequence_col]}\n')

def run_dsvm(vocab_file, ref_fa, alt_fa, out_file, deltasvm='/data/hps/assoc/private/cherrylab/user/lvand1/Resources/deltasvm_script/deltasvm.pl'):
    command = ['perl', deltasvm, ref_fa, alt_fa, vocab_file, out_file]
    subprocess.run(command, check=True)

def get_dsvm_dict(dsvm_table):
    dsvm_dict={}
    with open(dsvm_table, 'r') as dsvm:
        for line in dsvm:
            key, value = line.strip().split('\t')
            dsvm_dict[key] = value
    return dsvm_dict


# Motif characterization

def parse_JASPAR(jaspar_matrices_file):
        IUPAC_dict = {'R':'[AG]', 'Y':'[CT]', 'S':'[GC]',
        'W':'[AT]', 'K':'[GT]', 'M':'[AC]', 'B':'[CGT]',
        'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'} # For degenerate consensus regex
        JASPAR_dict = {}
        with open(jaspar_matrices_file, 'r') as matrices:
                for matrix in motifs.parse(matrices, 'jaspar'): # function from Bio.motifs
                        motif_counts = matrix.counts # Bio.motifs
                        consensus = str(motif_counts.degenerate_consensus)
                        consensus = consensus.replace('N','') # Ns in consensus are not useful downstream
                        motif_name = matrix.name
                        for code in IUPAC_dict:
                                consensus = consensus.replace(code, IUPAC_dict[code]) # make consensus regexable
                        JASPAR_dict[motif_name]=consensus
        return JASPAR_dict

def motif_search(sequence, JASPAR_dict):
        motifs_set = set()
        for motif in JASPAR_dict:
                if re.search(JASPAR_dict[motif], sequence):
                        motifs_set.add(motif)
        return motifs_set

## Another wrapper for applying to dataframe. If used in conjunction with the above functions for sequence grabbing, this should work fine on its own.
## Generates four new columns
def motif_wrapper(row, JASPAR_dict):
    ref_motifs = motif_search(row['ref_seq'], JASPAR_dict)
    alt_motifs = motif_search(row['alt_seq'], JASPAR_dict)
    motifs_down = ref_motifs - alt_motifs
    motifs_up = alt_motifs - ref_motifs
    return pd.Series([ref_motifs, alt_motifs, motifs_down, motifs_up])



# Plotters

def histo_plotter_dask(df, column_of_interest, nbins, axis_range, file_out):
    # Plot histo COI
    data = df[column_of_interest].repartition(npartitions=1).compute()
    # clear figure to fix axes
    plt.figure()
    # Make histo
    plt.hist(data, bins=nbins, range=axis_range, color='blue')
    # Plot labels
    plt.xlabel(f'{column_of_interest} distribution')
    plt.ylabel('Frequency')
    plt.savefig(file_out, format='pdf')
    # re-clear to prevent matplotlib opening too many figures
    plt.close()


def scatter_dSVM(x, y, file_out, transform_method='log10'):
    if transform_method == 'log10':
        y_old = y
        y = np.log10(y_old)
        ylab = 'Log10(Odds Ratio)'
    elif transform_method == 'log2':
        y_old = y
        y = np.log2(y_old)
        ylab = 'Log2(Odds Ratio)'
    elif transform_method == 'none':
        ylab = 'Odds Ratio'
    else:
        ylab = 'Odds Ratio'
        print('No accepted method of transformation selected.')
        print('For y axis transformation, set transform_method to log10 or log2.')
        print('Continuing with untransformed data.')
    corr_coef, pvalue = pearsonr(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    line_eq = f'y = {slope:.2f}x + {intercept:.2f}'
    plt.figure(figsize=(8,6))
    sns.regplot(x=x, y=y, scatter_kws={'s':80})
    plt.text(0.05 , 0.95, f'{line_eq}\nr={corr_coef:.3f}', transform=plt.gca().transAxes, fontsize=14, verticalalignment='top')
    plt.xlabel('deltaSVM')
    plt.ylabel(ylab)
    plt.savefig(file_out, format='pdf')
    # Close figure to clear output
    plt.close()





