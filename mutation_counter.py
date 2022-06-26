#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 05:17:55 2022

@author: hagar


***********not good yet - something about the alignment is wierd
This script was written to mesure the mutations in polio fasta
"""



#from sys import argv
from Bio.pairwise2 import format_alignment, align
from Bio import SeqIO
from os import listdir
from utils import get_sequences, hamming_distance
from math import floor
# seq_path = argv[1]
# ref_path = argv[2]
seq_path = "/mnt/project1/projects/POLIO/VDPV3_2022/tests_hagar/aligning_methods/de_novo/CNS/contig_based/"
ref_path = "/home/hagar/refs/polio/polio.fasta"

refs = get_sequences(ref_path)


mis_fruc={}

for seq_file in listdir(seq_path):
    seq = str(SeqIO.read(seq_path+seq_file, 'fasta').seq)
    seq_name = seq_file.split("-")[0]
    if "S1" in seq_file or "AY184219.1" in seq_file:
        mismatches = hamming_distance(list(seq), list(refs["AY184219.1"]))
        mis_fruc[seq_name+"_sabin1"] = floor(mismatches / (len(seq)) * 100)
    if "S2" in seq_file or "AY184220.1" in seq_file:
        mismatches = hamming_distance(list(seq), list(refs["AY184220.1"]))
        mis_fruc[seq_name+"_sabin2"] = floor(mismatches / (len(seq)) * 100)
    if "S3" in seq_file or "AY184221.1" in seq_file:
        mismatches = hamming_distance(list(seq), list(refs["AY184221.1"]))
        mis_fruc[seq_name+"_sabin3"] = floor(mismatches / (len(seq)) * 100)



alignments = align.globalxx("ACCGT", "ACG")
print(format_alignment(*alignments[0]))

