#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:18:33 2022

@author: hagar
"""

import pandas as pd
import utils
# from squars import dist_mat, squars
from threading import Lock

mutex = Lock()
utils.create_dirs(["alignments"])
regions_file = "/home/hagar/polio_regions/polio_regions.xlsx"
seq_file = "/home/hagar/refs/polio/polio.fasta"
regions = pd.ExcelFile(regions_file)
sabin1_regions = pd.read_excel(regions,"sabin1")
sabin2_regions = pd.read_excel(regions,"sabin2")
sabin3_regions = pd.read_excel(regions,"sabin3")
sabins={"AY184219.1": sabin1_regions, "AY184220.1": sabin2_regions, "AY184221.1": sabin3_regions}
seqs = utils.get_sequences(seq_file)
for index, row in sabin1_regions.iterrows():
    gene = sabin1_regions.at[index,"gene"]
    not_aligned = ""
    for sabin, df in sabins.items():
        start = df.at[index,"start"] - 1
        end = df.at[index,"end"] - 1 
        seq = seqs[sabin.split("_")[0]][start:end]
        not_aligned += ">"+sabin+'\n' + seq + '\n'
    f = open("alignments/" + gene + "_not_aligned.fasta", 'w')
    f.write(not_aligned)
    f.close()
    
    # mutex.acquire()
    # utils.mafft(0, "alignments/" + gene + "_not_aligned.fasta", "alignments/" + gene + "_aligned.fasta")
    # mutex.release()
    # aligned = utils.get_sequences("alignments/" + gene + "_aligned.fasta")
    # for sample, seq in aligned.items():
    #     aligned[sample] = list(seq)
    
    # dist_mut = dist_mat(aligned)
    # squars(dist_mut, gene + " Hamming Distance", aligned.keys())

    