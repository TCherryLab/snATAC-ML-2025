#!/usr/bin/python
# coding: utf-8

from sys import argv
script, fasta, path, out = argv

mut_arr = ['A', 'T', 'G', 'C']

write_out_WT = f'{path}{out}_WT.fa'
write_out_Mut = f'{path}{out}_mut.fa'

with open(fasta, "r") as fi:
    fasta = fi.read().splitlines()


def mutateSeq(header, dna, wt_out, mut_out):
    for n in range(len(dna)):
        ref = dna[n]
        for alt in mut_arr:
            mut_header = header + '|Mut:' + ref + '->' + alt + '--' + str(n+1) + '/' + str(len(dna))     
            mut_seq = list(dna)
            mut_seq[n] = alt
            wt_out.extend([mut_header, dna])
            mut_out.extend([mut_header, ''.join(mut_seq)])


wt_out = []
mut_out = []
for read_ind in range(0, len(fasta), 2):
    header = fasta[read_ind]
    dna = fasta[read_ind+1]
    mutateSeq(header,dna,wt_out,mut_out)


fo = open(write_out_WT, "w")
fo.write("\n".join(wt_out))
fo.close()


fo = open(write_out_Mut, "w")
fo.write("\n".join(mut_out))
fo.close()
