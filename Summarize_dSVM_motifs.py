#!/usr/bin/env python3

import sys
import subprocess
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Sasquatch recreation of ATimms motif summing script
# run in dSVM_comparisons mamba env
# find and add glob to env
# figure out bedtools in env

##parameters
delim = '\t'
##programs
#bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
##files etc
bt_genome = '/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/Hg38_chrsizes_chrM.txt'




######## Functions
# written by ATimms, edited as needed for updates to hpc


# Edited from original atimms function
# Takes raw bed with kluu satmut notation
def convert_dsvm_to_bed(infile, all_bed, neg_bed):
    score_dict = {}
    with open(infile, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip().split(delim)
            name = line[0]
            dsvm = float(line[1])
            extraction = r'(?P<chrom>chr[0-9XYM]+):(?P<start>\d+)-(?P<stop>\d+)\|Mut:(?P<ref>[atcgATCG])->(?P<mut>[atcgATCG])--(?P<pos_ref>\d+)/'
            match = re.match(extraction, name)
            if match:
                chrom = match.group('chrom')
                start = match.group('start')
                stop = match.group('stop')
                ref = match.group('ref')
                mut = match.group('mut')
                pos_ref = match.group('pos_ref')
            start = int(start) + int(pos_ref) - 1
            end = start + 1
            cse = '_'.join([chrom,str(start),str(end)])
            # print(cse, ref, mut, dvsm, line
            if ref != mut:
                if cse in score_dict:
                    score_dict[cse][0].append(dsvm)
                    if dsvm < 0:
                        score_dict[cse][1].append(dsvm)
                else:
                    if dsvm < 0:
                        score_dict[cse] = [[dsvm], [dsvm], ref]
                    else:
                        score_dict[cse] = [[dsvm], [], ref]
    # print(len(score_dict))
    all_bed_temp = all_bed.rsplit('.', 1)[0] + 'temp.bed'
    neg_bed_temp = neg_bed.rsplit('.', 1)[0] + 'temp.bed'
    with open(all_bed_temp, "w") as allt_fh, open(neg_bed_temp, "w") as negt_fh:
        for s in score_dict:
            # print(s, score_dict[s])
            all_ave = sum(score_dict[s][0]) / len(score_dict[s][0])
            if len(score_dict[s][1]) == 0:
                neg_ave = 'na'
            else:
                neg_ave = sum(score_dict[s][1]) / len(score_dict[s][1])
            allt_fh.write(delim.join(s.split('_') + [str(all_ave), score_dict[s][2]]) + '\n')
            negt_fh.write(delim.join(s.split('_') + [str(neg_ave), score_dict[s][2]]) + '\n')
    ##write final files
    with open(all_bed, "w") as all_fh:
        sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', all_bed_temp], stdout=all_fh)
        sort_file.wait()
    with open(neg_bed, "w") as neg_fh:
        sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', neg_bed_temp], stdout=neg_fh)
        sort_file.wait()


# motif_dict refers to simplification of motif names from homer ie 'CRX(Homeobox)': 'CRX'
# infile is a bed file from homer such as homer.KnownMotifs.hg38.191020.bed
# finds requested motif locations from homer database for hg38 (or other genome characterized by homer
def get_specific_motifs_from_homer_bed(infile, nomerge_out_prefix, merge_out_prefix, motif_dict):
    for motif in motif_dict:
        alt_motif = motif_dict[motif]
        no_merge_bed = nomerge_out_prefix + alt_motif + '.bed'
        final_bed = merge_out_prefix + alt_motif + '.bed'
        with open(infile, "r") as in_fh, open(no_merge_bed, "w") as tb_fh:
            for line in in_fh:
                line = line.rstrip().split(delim)
                chrom = line[0].replace('chr','')
                start = str(int(line[1]) - 1)
                end = line[2]
                homer_motif = line[3]
                if homer_motif == motif:
                    tb_fh.write(delim.join([chrom, start, end] + [motif]) + '\n')
        with open(final_bed, "w") as out_fh:
            sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', no_merge_bed], stdout=subprocess.PIPE)
            bt_merge = subprocess.Popen(['bedtools', 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
            bt_merge.wait()


# requires bedtools to be functional in env
# if bedtools works directly, change subprocess
# original version calls on local directory to personal bedtools download
# should be able to nab it from bioconda
def bedtools_intersect_dsvm_motif(dsvm_beds, motifs_beds, out_suffix):
    for dsvm_bed in dsvm_beds:
        for motifs_bed in motifs_beds:
            motif = motifs_bed.split('.')[3]
            out_bed = dsvm_bed.rsplit('.',1)[0] + '.' + motif + out_suffix
            with open(out_bed, "w") as naf_fh:
                hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', motifs_bed, '-b', dsvm_bed, '-wa', '-wb'], stdout=naf_fh)
                hom_bt_intersect.wait()


def get_counts_per_dsvm_and_motif(dsvm_beds, motifs_beds, out_file):
    dsvm_bed_count = 0
    header = ['dsvm_type']
    with open(out_file, "w") as out_fh:
        for dsvm_bed in dsvm_beds:
            dsvm_type = dsvm_bed.rsplit('.',1)[0]
            dsvm_bed_count += 1
            dsvm_scores = []
            for motifs_bed in motifs_beds:
                motif = motifs_bed.split('.')[3]
                header.append(motif)
                dvsm_motif_bed = dsvm_type + '.' + motif + '.bed'
                dvsm_list = []
                with open(dvsm_motif_bed, "r") as dmb_fh:
                    for line in dmb_fh:
                        line = line.rstrip().split(delim)
                        score = line[6]
                        if score != 'na':
                            dvsm_list.append(float(score))
                if len(dvsm_list) > 0:
                    average_dvsm = (sum(dvsm_list)/len(dvsm_list))
                else:
                    average_dsvm = np.nan
                    print(f'Divide by zero, no scores for {motif} in {dsvm_type}')
                dsvm_scores.append(str(round(average_dvsm,3)))
            all_bed = dsvm_type + '.bed'
            header.append('all')
            dvsm_list = []
            with open(all_bed, "r") as ab_fh:
                for line in ab_fh:
                    line = line.rstrip().split(delim)
                    score = line[3]
                    if score != 'na':
                        dvsm_list.append(float(score))
            average_dvsm = (sum(dvsm_list)/len(dvsm_list))
            dsvm_scores.append(str(round(average_dvsm,3)))
            if dsvm_bed_count == 1:
                out_fh.write(delim.join(header) + '\n')
            line_out = [dsvm_type] + dsvm_scores
            out_fh.write(delim.join(line_out) + '\n')


def bedtools_slop_motif_beds(in_beds):
    for in_bed in in_beds:
        ##bedtools slop
        out_bed = in_bed.rsplit('.', 1)[0] + '.ext25.bed'
        with open(out_bed, "w") as out_fh: 
            hom_bt_intersect = subprocess.Popen(['bedtools', 'slop', '-i', in_bed, '-g', bt_genome, '-b', '25'], stdout=out_fh)
            hom_bt_intersect.wait()


def get_ave_dsvm_per_bp_around_motif(dvsm_motif_beds):
    for dvsm_motif_bed in dvsm_motif_beds:
        ##make dict with all score per motif regions
        motif_region_svsm_dict = {}
        with open(dvsm_motif_bed, "r") as dmb_fh:
            for line in dmb_fh:
                line = line.rstrip().split(delim)
                motif_region = '_'.join(line[:3])
                motif_region_len = int(line[2]) - int(line[1])
                score = line[7]
                if motif_region in motif_region_svsm_dict:
                    motif_region_svsm_dict[motif_region][0].append(score)
                else:
                    motif_region_svsm_dict[motif_region] = [[score], motif_region_len]
        ##make new dict with just the regions that have dsvm score for all the reion
        new_motif_dvsm_dict = {}
        for mr in motif_region_svsm_dict:
            # print(dvsm_motif_bed, mr, motif_region_svsm_dict[mr][1], len(motif_region_svsm_dict[mr][0]))
            if motif_region_svsm_dict[mr][1] == len(motif_region_svsm_dict[mr][0]):
                new_motif_dvsm_dict[mr] = motif_region_svsm_dict[mr][0]
        ##make dict into a panda's dataframe and get
        motif_df = pd.DataFrame.from_dict(new_motif_dvsm_dict, orient='index')
        ##write file
        dvsm_motif_csv = dvsm_motif_bed.rsplit('.', 1)[0] + '.csv'
        motif_df.to_csv(dvsm_motif_csv)
        ##convert columns to floats
        cols = motif_df.columns
        motif_df[cols] = motif_df[cols].apply(pd.to_numeric, errors='coerce')
        ##get averages
        motif_df_means = motif_df.mean(axis = 0, skipna = True)
        # print(motif_df_means)
        # plt.figure(figsize=(6,4))
        #sns.lineplot(data=motif_df_means)
        #pdf_name = dvsm_motif_bed.rsplit('.', 1)[0] + '.pdf'
        #plt.savefig(pdf_name)
        #plt.close()

def get_ave_dsvm_per_bp_only_motif(dvsm_motif_beds):
    for dvsm_motif_bed in dvsm_motif_beds:
        ##make dict with all score per motif regions
        motif_region_svsm_dict = {}
        with open(dvsm_motif_bed, "r") as dmb_fh:
            for line in dmb_fh:
                line = line.rstrip().split(delim)
                motif_region = '_'.join(line[:3])
                motif_region_len = int(line[2]) - int(line[1])
                score = line[7]
                if motif_region in motif_region_svsm_dict:
                    motif_region_svsm_dict[motif_region][0].append(score)
                else:
                    motif_region_svsm_dict[motif_region] = [[score], motif_region_len]
        ##make new dict with just the regions that have dsvm score for all the reion
        new_motif_dvsm_dict = {}
        for mr in motif_region_svsm_dict:
            # print(dvsm_motif_bed, mr, motif_region_svsm_dict[mr][1], len(motif_region_svsm_dict[mr][0]))
            if motif_region_svsm_dict[mr][1] == len(motif_region_svsm_dict[mr][0]):
                new_motif_dvsm_dict[mr] = motif_region_svsm_dict[mr][0]
        ##make dict into a panda's dataframe and get
        motif_df = pd.DataFrame.from_dict(new_motif_dvsm_dict, orient='index')
        ##write file
        dvsm_motif_csv = dvsm_motif_bed.rsplit('.', 1)[0] + '.csv'
        motif_df.to_csv(dvsm_motif_csv)
        ##convert columns to floats
        cols = motif_df.columns
        motif_df[cols] = motif_df[cols].apply(pd.to_numeric, errors='coerce')
        ##get averages
        # motif_df_means = motif_df.mean(axis = 0, skipna = True)
        # print(motif_df_means)
        # plt.figure(figsize=(6,4))
        #sns.barplot(data=motif_df, color = 'blue', ci=None)
        #pdf_name = dvsm_motif_bed.rsplit('.', 1)[0] + '.pdf'
        #plt.savefig(pdf_name)
        #plt.close()


def main():
    ##run methods
    # edit all as required
    homer_bed = '/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/homer.KnownMotifs.hg38.191020.bed'
    working_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/dsvm'
    os.chdir(working_dir)
    ## motif bed from scanMotifGenomeWide.pl or from homer.KnownMotifs.hg38.191020.bed
    # motifs from Thomas et al heatmap fig 1
    # not all motifs are present in the homer database
    # it's possible to make novel motif bed files using homer but it's a pita
    homer_motif_dict = {'JunB(bZIP)': 'JUNB', 'Jun-AP1(bZIP)': 'JUN', 'Fra1(bZIP)':'FOS', 'Fosl2(bZIP)':'FOSL2',
            'Bach1(bZIP)':'BACH1', 'NF1-halfsite(CTF)':'NFIX', 'Maz(Zf)':'MAZ',
            'Lhx2(Homeobox)':'LHX2', 'RORa(NR)':'RORA', 'NFIL3(bZIP)':'NFIL3', 'HLF(bZIP)':'HLF', 'Mef2c(MADS)':'MEF2C',
            'Mef2d(MADS)':'MEF2D', 'Mef2a(MADS)':'MEF2A', 'Otx2(Homeobox)': 'OTX2', 'CRX(Homeobox)': 'CRX',
            'Meis1(Homeobox)':'MEIS1', 'HNF6(Homeobox)':'ONECUT1', 'Hnf6b(Homeobox)':'ONECUT2', 'Ascl2(bHLH)': 'ASCL2',
            'HEB(bHLH)':'TCF12', 'Ascl1(bHLH)':'ASCL1', 'CUX1(Homeobox)':'CUX1', 'NeuroD1(bHLH)':'NEUROD1', 
            'Ap4(bHLH)':'TFAP4', 'EBF1(EBF)': 'EBF1',  'Snail1(Zf)':'SNAI1', 'Erra(NR)':'ESRRA', 'ERRg(NR)':'ESRRG', 
            'Esrrb(NR)':'ESRRB', 'AP-2alpha(AP2)':'TFAP2A', 'HIF2a(bHLH)':'EPAS1'}
    homer_motif_prefix = 'homer.hg38.191020.'
    homer_motif_nomerge_prefix = 'homer.hg38.no_merge.'
    homer_motifs_beds = [homer_motif_prefix + i + '.bed' for i in homer_motif_dict.values()]
    homer_motifs_nomerge_beds = [homer_motif_nomerge_prefix + i + '.bed' for i in homer_motif_dict.values()]
    homer_motifs_ext25_beds = [homer_motif_nomerge_prefix + i + '.ext25.bed' for i in homer_motif_dict.values()]
    # Get motif specific bed files with and without merge
    #get_specific_motifs_from_homer_bed(homer_bed, homer_motif_nomerge_prefix, homer_motif_prefix, homer_motif_dict)
    predictions = ['all_11mers_AC_HC_GC_precursors_noImmIn', 'all_11mers_Late_Progenitors_noImmIn',
            'all_11mers_Developing_Amacrines_noImmIn', 'all_11mers_Mature_Amacrines_noMatIn',
            'all_11mers_Developing_Bipolars_noImmIn', 'all_11mers_Mature_Bipolars_noMatIn',
            'all_11mers_Developing_Cones_noImmIn', 'all_11mers_Mature_Cones_noMatIn',
            'all_11mers_Developing_Ganglions_noImmIn', 'all_11mers_Mature_Horizontals_noMatIn', 
            'all_11mers_Developing_Horizontals_noImmIn', 'all_11mers_Mature_Mullers_noMatIn', 
            'all_11mers_Developing_Rods_noImmIn', 'all_11mers_Mature_Rods_noMatIn',
            'all_11mers_Early_Progenitors_noImmIn', 'all_11mers_Photoreceptor_Bipolar_Precursors_noImmIn',
            'all_11mers_Ganglion_Precursors_noImmIn', 'Macular_RPE']
    for model in predictions:
        ## run methods by model
        dsvm_leah = 'dSVM_' + model + '_t25k_MegaOutgroup.txt'
        dsvm_all_bed = 'motifs/dsvm.' + model + '_all_scores.bed'
        dsvm_neg_bed = 'motifs/dsvm.' + model + '_neg_scores.bed'
        dsvm_motif_bed = 'motifs/dsvm.' + model + '_neg_scores.motifs.bed'
        dsvm_bed_files = [dsvm_all_bed, dsvm_neg_bed]
        motif_counts_file = 'motifs/' + model + '_motifs.counts.txt'
        dsvm_motifs_ext25_beds = glob.glob(f'motifs/*{model}*ext25.bed')
        dsvm_motif_only_beds = glob.glob(f'motifs/*{model}*.motif_only.bed')
        # Make dsvm bed files
        convert_dsvm_to_bed(dsvm_leah, dsvm_all_bed, dsvm_neg_bed)
        # bedtools intersect dsvm scores with homer beds
        print(homer_motifs_beds)
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_beds, '.bed')
        # get average dsvm per motif for the different score
        get_counts_per_dsvm_and_motif(dsvm_bed_files, homer_motifs_beds, motif_counts_file)
        ## Look at specific positions within motif
        # extend motif beds +-25
        bedtools_slop_motif_beds(homer_motifs_nomerge_beds)
        # intersect dsvm scores with homer beds extended
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_ext25_beds, '.ext25.bed')
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_nomerge_beds, '.motif_only.bed')
        ## Then count data and make graphs
        # 25 bp exts
        print(dsvm_motifs_ext25_beds)
        get_ave_dsvm_per_bp_around_motif(dsvm_motifs_ext25_beds)
        print(dsvm_motif_only_beds)
        get_ave_dsvm_per_bp_only_motif(dsvm_motif_only_beds)
        ## Plot?


def main2():
    ##run methods
    # edit all as required
    homer_bed = '/data/hps/assoc/private/cherrylab/user/lvand1/Resources/NGS_tools/SEQ_UTILS/homer.KnownMotifs.hg38.191020.bed'
    working_dir = '/data/hps/assoc/private/cherrylab/user/lvand1/scATAC_dsvm/dsvm'
    os.chdir(working_dir)
    ## motif bed from scanMotifGenomeWide.pl or from homer.KnownMotifs.hg38.191020.bed
    # motifs from Thomas et al heatmap fig 1
    # not all motifs are present in the homer database
    # it's possible to make novel motif bed files using homer but it's a pita
    homer_motif_dict_new = {'SOX2(HMG)':'SOX2'}
    homer_motif_prefix = 'homer.hg38.191020.'
    homer_motif_nomerge_prefix = 'homer.hg38.no_merge.'
    #get_specific_motifs_from_homer_bed(homer_bed, homer_motif_nomerge_prefix, homer_motif_prefix, homer_motif_dict_new)
    #homer_motif_dict = {'JunB(bZIP)': 'JUNB', 'Jun-AP1(bZIP)': 'JUN', 'Fra1(bZIP)':'FOS', 'Fosl2(bZIP)':'FOSL2',
        #'Bach1(bZIP)':'BACH1', 'NF1-halfsite(CTF)':'NFIX', 'Maz(Zf)':'MAZ', 'lhx4':'LHX4',
        #'Lhx2(Homeobox)':'LHX2', 'RORa(NR)':'RORA', 'NFIL3(bZIP)':'NFIL3', 'HLF(bZIP)':'HLF', 'Mef2c(MADS)':'MEF2C',
        #'Mef2d(MADS)':'MEF2D', 'Mef2a(MADS)':'MEF2A', 'Otx2(Homeobox)': 'OTX2', 'CRX(Homeobox)': 'CRX',
        #'Meis1(Homeobox)':'MEIS1', 'HNF6(Homeobox)':'ONECUT1', 'Hnf6b(Homeobox)':'ONECUT2', 'Ascl2(bHLH)': 'ASCL2',
        #'HEB(bHLH)':'TCF12', 'Ascl1(bHLH)':'ASCL1', 'CUX1(Homeobox)':'CUX1', 'NeuroD1(bHLH)':'NEUROD1',
        #'Ap4(bHLH)':'TFAP4', 'pou4f1':'POU4F1', 'pou4f2':'POU4F2', 'pou4f3':'POU4F3','EBF1(EBF)': 'EBF1',
        #'Snail1(Zf)':'SNAI1', 'rorb':'RORB', 'Erra(NR)':'ESRRA', 'ERRg(NR)':'ESRRG', 'Esrrb(NR)':'ESRRB',
        #'AP-2alpha(AP2)':'TFAP2A', 'HIF2a(bHLH)':'EPAS1', 'MITF(bHLH)':'MITF', 'PAX6(Paired,Homeobox)':'PAX6','SOX2(HMG)':'SOX2'} # avoiding making the big custom motifs table because I don't care for this one
    homer_motif_dict = {'sox2':'SOX2', 'lhx4':'LHX4'}
    homer_motifs_beds = [homer_motif_prefix + i + '.bed' for i in homer_motif_dict.values()]
    homer_motifs_nomerge_beds = [homer_motif_nomerge_prefix + i + '.bed' for i in homer_motif_dict.values()]
    homer_motifs_ext25_beds = [homer_motif_nomerge_prefix + i + '.ext25.bed' for i in homer_motif_dict.values()]
    # Get motif specific bed files with and without merge
    # In this version, this was done manually using custom motif files
    # Motif files modified from MEME motifs because Homer didn't have them
    # Get no_merge beds with scanMotifsGenomeWide.pl and then sort/merge
    predictions = ['all_11mers_Early_Progenitors_noImmIn','all_11mers_Late_Progenitors_noImmIn',
        'all_11mers_Mature_Mullers_noMatIn', 'all_11mers_Photoreceptor_Bipolar_Precursors_noImmIn',
        'all_11mers_Developing_Bipolars_noImmIn','all_11mers_Mature_Bipolars_noMatIn',
        'all_11mers_Developing_Rods_noImmIn','all_11mers_Mature_Rods_noMatIn',
        'all_11mers_Developing_Cones_noImmIn','all_11mers_Mature_Cones_noMatIn',
        'all_11mers_AC_HC_GC_precursors_noImmIn','all_11mers_Developing_Amacrines_noImmIn',
        'all_11mers_Mature_Amacrines_noMatIn','all_11mers_Developing_Horizontals_noImmIn',
        'all_11mers_Mature_Horizontals_noMatIn','all_11mers_Developing_Ganglions_noImmIn',
        'Macular_RPE','all_11mers_Ganglion_Precursors_noImmIn']
    for model in predictions:
        ## run methods by model
        dsvm_leah = 'dSVM_' + model + '_t25k_MegaOutgroup.txt'
        dsvm_all_bed = 'motifs/dsvm.' + model + '_all_scores.bed'
        dsvm_neg_bed = 'motifs/dsvm.' + model + '_neg_scores.bed'
        dsvm_motif_bed = 'motifs/dsvm.' + model + '_neg_scores.custom_motifs.bed'
        dsvm_bed_files = [dsvm_all_bed, dsvm_neg_bed]
        motif_counts_file = 'motifs/' + model + '_sox2_motifs.counts.txt'
        dsvm_motifs_ext25_beds = glob.glob(f'motifs/*{model}*ext25.bed')
        dsvm_motif_only_beds = glob.glob(f'motifs/*{model}*.motif_only.bed')
        # Make dsvm bed files
        convert_dsvm_to_bed(dsvm_leah, dsvm_all_bed, dsvm_neg_bed)
        # bedtools intersect dsvm scores with homer beds
        print(homer_motifs_beds)
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_beds, '.bed')
        # get average dsvm per motif for the different score
        get_counts_per_dsvm_and_motif(dsvm_bed_files, homer_motifs_beds, motif_counts_file)
        ## Look at specific positions within motif
        # extend motif beds +-25
        #bedtools_slop_motif_beds(homer_motifs_nomerge_beds)
        # intersect dsvm scores with homer beds extended
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_ext25_beds, '.ext25.bed')
        bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_nomerge_beds, '.motif_only.bed')
        ## Then count data and make graphs
        # 25 bp exts
        print(dsvm_motifs_ext25_beds)
        get_ave_dsvm_per_bp_around_motif(dsvm_motifs_ext25_beds)
        print(dsvm_motif_only_beds)
        get_ave_dsvm_per_bp_only_motif(dsvm_motif_only_beds)


if __name__=='__main__':
    main2()


