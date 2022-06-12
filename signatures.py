#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 14:06:34 2022

@author: hagar
"""
from sys import argv
from math import floor
import pandas as pd
from utils import translate_table, get_sequences, mutations_positions
from format_xl import save_format_xl
ambiguous_nucleotides = ["W", "Y", "R", "S", "D","K","M","V","H","B","X"]

"""
get a dictionary {sample : mutations list}
"""
def mutations_by_sample(mutations_position,sequences):
    mutations_by_sample = {}
    for sample, record in sequences.items():
        mutations = []
        for pos in mutations_position:
            mutations.append(record[pos])
        mutations_by_sample[sample] = mutations
    return mutations_by_sample


def get_regions(regions_csv):
    regions = {}
    f = open(regions_csv)
    for line in f.readlines():
        if line.upper().startswith("GENE"):
            continue
        line = line.split(",")
        regions[line[0]] = (int(line[1]),int(line[2]))
    f.close()
    return regions


def get_gene(mutations_positions_nt, regions):
    gene_names = []
    position_on_gene_nt = []
    position_on_gene_aa = []
    for mut in mutations_positions_nt:
        gene_name = "UTR"
        pos_gene_nt = ""
        pos_gene_aa = ""
        for gene, pos in regions.items():
            start = pos[0]
            end = pos[1]
            if mut in range(start,end):
                gene_name = gene
                pos_gene_nt = mut-start + 2
                pos_gene_aa = floor(pos_gene_nt/3) if pos_gene_nt%3 == 0 else  floor(pos_gene_nt/3) + 1
        gene_names.append(gene_name)
        position_on_gene_nt.append(pos_gene_nt)
        position_on_gene_aa.append(pos_gene_aa)

    return gene_names,position_on_gene_nt,position_on_gene_aa 

def aa_sum(df, sequences):
    aa_groups = pd.read_csv("AAproperties.txt", sep = '\t')
    for index, row in df.iloc[:,-(len(sequences)):].iterrows():
        row.reset_index(drop=True,inplace=True)
        no_x_aa = list(dict.fromkeys([x for x in row.to_list() if x != "X"]))
        if len(no_x_aa) > 1:
            df.at[index,"R/S"] =  "R"
            if len(no_x_aa) == 2:
                groups = aa_groups.loc[aa_groups['Abbv1'] == row.iloc[0], 'properties'].iloc[0] + "(" + row.iloc[0] +"), " 
                no_x_aa.remove(row.iloc[0])
                groups += aa_groups.loc[aa_groups['Abbv1'] == no_x_aa[0], 'properties'].iloc[0] + "(" + no_x_aa[0] +")"
                df.at[index, "aa_group"] = groups
        else:
            df.at[index,"R/S"] =  "S"


def get_single_aa(seq, position, start):
    pos_on_gene = position - start + 1

    mod = pos_on_gene % 3
    if mod == 0:
        codon_pos = (position - 2, position - 1, position)
    if mod == 1:
        codon_pos = (position, position + 1, position + 2)
    if mod == 2:
        codon_pos = (position - 1, position, position + 1)
        
    codon = seq[codon_pos[0]-1] + seq[codon_pos[1]-1] + seq[codon_pos[2]-1]
    aa = translate_table[codon] if not '-' in codon and not 'N' in codon else 'X'
    return aa


def get_all_aa(mutations_positions_nt, sequences, gene_names, regions):

    mutations_by_sample_aa = {}
    for sample, seq in sequences.items():
        sample_aa = []
        for i in range(len(mutations_positions_nt)):
            if gene_names[i] == 'UTR':
                sample_aa.append('X')
            else:
                pos = mutations_positions_nt[i]
                gene_start = regions[gene_names[i]][0] 
                aa = get_single_aa(seq, pos, gene_start)
                sample_aa.append(aa)
        mutations_by_sample_aa[sample] = sample_aa

    return mutations_by_sample_aa
    
def run(alignment_file,regions_csv,output):

    sequences = get_sequences(alignment_file)
    mutations_positions_nt = mutations_positions(sequences)
    mutations_by_sample_nt = mutations_by_sample(mutations_positions_nt,sequences)
    regions = get_regions(regions_csv)
    gene_names, position_on_gene_nt, position_on_gene_aa = get_gene(mutations_positions_nt, regions)
    get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    mutations_by_sample_aa = get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    
    df = pd.DataFrame()
    df["gene_name"] = gene_names
    df["nt_position_on_gene"] = position_on_gene_nt
    df["nt_position_on_genome"] = mutations_positions_nt
    df["nt_position_on_genome"] += 1
    for sample, mut in mutations_by_sample_nt.items():
        df[sample+"_NT"] = mut
    df["aa_position_on_gene"] = position_on_gene_aa
    for sample, mut in mutations_by_sample_aa.items():
        df[sample+"_AA"] = mut    
    
    aa_sum(df, sequences)
    save_format_xl(df, len(sequences)-1)

    
if __name__ == "__main__":
    run(argv[1], argv[2], argv[3])