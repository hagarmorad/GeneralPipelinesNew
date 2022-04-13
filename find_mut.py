#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:02:44 2022

@author: hagar
"""

import pandas as pd

f1 = open("/home/hagar/test/cmv1.fasta")
f2 = open("/home/hagar/test/cmv2.fasta")
fasta1 = ""
fasta2 = ""
for line in f1.readlines():
    if not line.startswith(">"):
        fasta1 += line.strip()
        
for line in f2.readlines():
    if not line.startswith(">"):
        fasta2 += line.strip()

fasta1 = list(fasta1)
fasta2 = list(fasta2)
noN = 0
N = 0
loc = [0]*len(fasta1)
for pos in range(len(fasta1)):
    if fasta1[pos] != fasta2[pos]:
        if fasta1[pos] == 'N' or fasta2[pos] == 'N' or fasta1[pos] == '-' or fasta2[pos] == '-':
            N += 1
            continue
        noN += 1
        loc[pos]=1


dict = {"cmv1": fasta1, "cmv2": fasta2, " ": loc}
df = pd.DataFrame(dict) 
 
df.to_csv('cmv.csv', index=False)