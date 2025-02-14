#!/usr/bin/env python3
"""
dsvm2bedGraph.py
Written for python v3.?
Converts deltaSVM output score vector into 4 bedGraph outputs - one per alt allele. 
Attempt for user input of file names in command line
"""
from sys import argv
script, windows_in, dsvm_vec, out, dir = argv

A_bg_out = (f'{dir}A_{out}.bedGraph')
C_bg_out = (f'{dir}C_{out}.bedGraph')
G_bg_out = (f'{dir}G_{out}.bedGraph')
T_bg_out = (f'{dir}T_{out}.bedGraph')

#importing disease windows:
with open(windows_in, "r") as fi:
    fasta = fi.read().splitlines()

#breaking up disease windows into 1bp coords x 4muts in bed format:
nt_bed = []
for i in range(len(fasta)):
    bed = fasta[i].split("\t")
    for j in range(int(bed[1]),int(bed[2])):
        #repeat each bp coord x4 (once for each bp substitution):
        for k in range(4):
            nt_bed.append(bed[0] + '\t' + str(j) + '\t' + str(j+1))

#import delSVM score vector:
with open(dsvm_vec, "r") as fi:
    scores = fi.read().splitlines()

#append scores to tab delimited disease bases bedfile:
outA = []
outT = []
outG = []
outC = []
for i in range(len(nt_bed)):
    bg = nt_bed[i] + "\t" + scores[i]
    mutNum = i % 4
    if mutNum == 0:
        outA.append(bg)
    elif mutNum == 1:
        outT.append(bg)
    elif mutNum == 2:
        outG.append(bg)
    else:
        outC.append(bg)

#save to output files:
fo = open(A_bg_out, "w")
fo.write("\n".join(outA))
fo.write("\n")
fo.close()

fo = open(T_bg_out, "w")
fo.write("\n".join(outT))
fo.write("\n")
fo.close()

fo = open(G_bg_out, "w")
fo.write("\n".join(outG))
fo.write("\n")
fo.close()

fo = open(C_bg_out, "w")
fo.write("\n".join(outC))
fo.write("\n")
fo.close()
