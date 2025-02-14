#!/usr/bin/env python3


import os
import re
import pandas as pd
import pybedtools
import subprocess
from Bio import motifs


# To do:
    # write functions for dsvm, jaspar
    # test and apply functions
    # finish off adding dsvm scores and jaspar calls to table
    # run on real tables



def get_fasta(bed, ref_fa='/data/hps/assoc/public/bioinformatics/annotations/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'):
    position_seq = bed.sequence(fi=ref_fa)
    fasta_str = open(position_seq.seqfn).read()
    fasta_list = fasta_str.split('\n')[1:-1] # skip header and empty row
    return ''.join(fasta_list)


# build fresh variant_seq for individual sequences
    # account for indels
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
        #print(f'reference: {ref}\nsequence being replaced: {sequence[slop-1:pos_end]}')
    else:
        alt_seq = sequence[:slop-1] + alt.upper() + sequence[slop:]
    #print(f'ref_seq: {sequence}\nalt_seq: {alt_seq}')
    return alt_seq


# wrapper function to apply sequence getters to rows in a df
def generate_sequences(row, slop=11, ref_fa='/data/hps/assoc/public/bioinformatics/annotations/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'):
    ref_seq = get_reference_sequence(row['Chr'], row['Pos'], row['Ref'], slop=slop, ref_fa=ref_fa)
    alt_seq = get_alt_sequence(ref_seq, row['Ref'], row['Alt'], slop=slop)
    return pd.Series([ref_seq, alt_seq])

# build dsvm function(s)
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


# build jaspar function
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

def motif_wrapper(row, JASPAR_dict):
    ref_motifs = motif_search(row['ref_seq'], JASPAR_dict)
    alt_motifs = motif_search(row['alt_seq'], JASPAR_dict)
    motifs_down = ref_motifs - alt_motifs
    motifs_up = alt_motifs - ref_motifs
    return pd.Series([ref_motifs, alt_motifs, motifs_down, motifs_up])



# For simple code fixing
# test file has indel and snp examples
def tester():
    tester_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/GWAS_scoring/'
    tester_file = '/data/hps/assoc/private/cherrylab/user/lvand1/GWAS_scoring/test_file.txt'
    motif_dict = parse_JASPAR('/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/JASPAR_redundant_motifs.txt')
    my_df = pd.read_csv(tester_file, sep='\t')
    my_df[['ref_seq', 'alt_seq']] = my_df.apply(generate_sequences, axis=1)
    my_df['header'] = my_df.apply(lambda row: f"{row['Chr']}:{row['Pos']}_{row['Ref']}>{row['Alt']}", axis=1)
    my_df[['ref_motifs','alt_motifs','motifs_down','motifs_up']] = my_df.apply(lambda row: motif_wrapper(row, JASPAR_dict=motif_dict), axis=1)
    print(my_df)
    write_fa(my_df, 'ref_seq', f'{tester_dir}tester_ref.fa')
    write_fa(my_df, 'alt_seq', f'{tester_dir}tester_alt.fa')
    run_dsvm('/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/11mers/all_11mers_Macular_RPE_t25k.txt', f'{tester_dir}tester_ref.fa', f'{tester_dir}tester_alt.fa', f'{tester_dir}tester_dsvm.txt')
    dsvm_dict =  get_dsvm_dict(f'{tester_dir}tester_dsvm.txt')
    my_df['dsvm'] = my_df['header'].map(dsvm_dict)
    print(my_df)

# The real analysis, with loops
def main():
    # Base files
    working_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/GWAS_scoring/'
    GWAS_tables = ['AMD', 'Glaucoma', 'DR', 'MacTel']
    models_path = '/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/11mers/'
    models_list = ['Macular_RPE_t25k', 'Mature_Cones_noMatIn_t25k_', 'Mature_Mullers_noMatIn_t25k_', 'Mature_Rods_noMatIn_t25k_', 'Developing_Ganglions_noImmIn_t25k_']
    motif_dict = parse_JASPAR('/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/JASPAR_redundant_motifs.txt')
    # Loop all tables
    for table in GWAS_tables:
        full_path = f'{working_dir}{table}_GWAS.txt'
        my_df = pd.read_csv(full_path, sep='\t')
        # Get ref/alt sequences
        my_df[['ref_seq', 'alt_seq']] = my_df.apply(generate_sequences, axis=1)
        my_df['header'] = my_df.apply(lambda row: f"{row['Chr']}:{row['Pos']}_{row['Ref']}>{row['Alt']}", axis=1)
        # Get JASPAR motif calls for sequences
        my_df[['ref_motifs','alt_motifs','motifs_down','motifs_up']] = my_df.apply(lambda row: motif_wrapper(row, JASPAR_dict=motif_dict), axis=1)
        # write .fa files
        write_fa(my_df, 'ref_seq', f'{working_dir}{table}_ref.fa')
        write_fa(my_df, 'alt_seq', f'{working_dir}{table}_alt.fa')
        # get dSVM scores for all models desired (specified above)
        for model in models_list:
            vocab = f'{models_path}all_11mers_{model}.txt' 
            run_dsvm(vocab, f'{working_dir}{table}_ref.fa', f'{working_dir}{table}_alt.fa', f'{working_dir}{model}_dsvm.txt')
            dsvm_dict =  get_dsvm_dict(f'{working_dir}{model}_dsvm.txt') # parse scores back into python
            my_df[f'dsvm_{model}'] = my_df['header'].map(dsvm_dict) # map scores to table
        # write out table
        my_df.to_csv(f'{working_dir}{table}_GWAS_dSVM_scored.txt', sep='\t', index=False)


if __name__ == '__main__':
    #tester()
    main()
