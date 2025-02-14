#!/usr/bin/env python3

import csv
from dSVM_functions import *
from collections import defaultdict

# Build JASPAR function or use from dSVM_functions
# parse_JASPAR
motif_dict = parse_JASPAR('/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/JASPAR_redundant_motifs.txt')
# set working directory (11mers in scATAC)
working_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/11mers/'
# Set list of models
models = ['Early_Progenitors_noImmIn','Late_Progenitors_noImmIn',
		'Mature_Mullers_noMatIn', 'Photoreceptor_Bipolar_Precursors_noImmIn',
		'Developing_Bipolars_noImmIn','Mature_Bipolars_noMatIn',
		'Developing_Rods_noImmIn','Mature_Rods_noMatIn',
		'Developing_Cones_noImmIn','Mature_Cones_noMatIn',
		'AC_HC_GC_precursors_noImmIn','Developing_Amacrines_noImmIn',
		'Mature_Amacrines_noMatIn','Developing_Horizontals_noImmIn',
		'Mature_Horizontals_noMatIn','Ganglion_Precursors_noImmIn',
		'Developing_Ganglions_noImmIn','Macular_RPE']

# read in top 1% vocab to dict
# for each sequence, motif_search to a new level in dict
for model in models:
    file = f'{working_dir}top1_vocab_{model}_t25k_.txt'
    motif_counter = defaultdict(int)
    output = f'{working_dir}top1_vocab_{model}_motifcounts.txt'
    with open(file, 'r') as filein, open(output, 'w') as fileout:
        reader = csv.DictReader(filein, delimiter='\t', fieldnames=['sequence', 'score'])
        writer = csv.writer(fileout, delimiter='\t')
        for row in reader:
            for motif in motif_search(row['sequence'], motif_dict):
                motif_counter[motif] += 1
        sorted_motif_counts = sorted(motif_counter.items(), key=lambda x: x[1], reverse=True)
        writer.writerows(sorted_motif_counts)

