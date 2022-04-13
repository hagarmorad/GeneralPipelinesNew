#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 14:06:34 2022

@author: hagar
"""

import sys
from math import floor
import pandas as pd
import utils
alignment_file = sys.argv[1]
regions_csv = sys.argv[2]
output = sys.argv[3]

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
        line = line.split(",")
        regions[line[0]] = (int(line[1]),int(line[2]))
    f.close()
    return regions


def translate(seq):
	table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
    
	seq = seq.replace("\n", "")
	seq = seq.replace("\r", "")
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			if "N" in codon or "-" in codon or any(x in codon for x in ambiguous_nucleotides):
				protein+= "X"
			else:
				protein+= table[codon]
	return protein

def translate_sequences(sequences,cds):
    aa_seqs = {}
    for sample, record in sequences.items():
        aa_seqs[sample] = translate(str(record[cds[0]-1:cds[1]]))
    return aa_seqs

def aa_mutations(sequences,cds):
    mutations_by_sample={}
    for sample, record in sequences.items():
        aa_muts = []     
        for mut_pos in mutations_positions_nt:
            if mut_pos < cds[0] or mut_pos > cds[1]:
                 aa_muts.append("X")
            else:
                aa_muts.append(record[int((mut_pos-cds[0]+1)/3)])
        mutations_by_sample[sample] = aa_muts
    return mutations_by_sample


def get_gene():
    gene_names = []
    position_on_gene_nt = []
    position_on_gene_aa = []
    for mut in mutations_positions_nt:
        for gene, pos in regions.items():
            start = pos[0]
            end = pos[1]
            if mut in range(start,end):
                gene_name = gene
                pos_gene_nt = mut-start + 2
                pos_gene_aa = floor((mut-start)/3)
        gene_names.append(gene_name)
        position_on_gene_nt.append(pos_gene_nt)
        position_on_gene_aa.append(pos_gene_aa)

    return gene_names,position_on_gene_nt,position_on_gene_aa 

def aa_sum(df):
    aa_groups = pd.read_csv("/home/hagar/useful/AAproperties.txt", sep = '\t')
    for index, row in df.iloc[:,-(len(sequences)):].iterrows():
        row.reset_index(drop=True,inplace=True)
        no_x_aa = list(dict.fromkeys([x for x in row.to_list() if x != "X"]))
        if len(no_x_aa) > 1:
            df.at[index,"R/S"] =  "R"
            if len(no_x_aa) == 2:
                groups = aa_groups.loc[aa_groups['Abbv1'] == row.iloc[0], 'properties'].iloc[0] + ", "
                no_x_aa.remove(row.iloc[0])
                groups += aa_groups.loc[aa_groups['Abbv1'] == no_x_aa[0], 'properties'].iloc[0]
                df.at[index, "aa_group"] = groups
        else:
            df.at[index,"R/S"] =  "S"
         
    
sequences = utils.get_sequences(alignment_file)
mutations_positions_nt = utils.mutations_positions(sequences)
mutations_by_sample_nt = mutations_by_sample(mutations_positions_nt,sequences)
regions = get_regions(regions_csv)
cds = regions.pop("CDS")
translated = translate_sequences(sequences,cds)
mutations_by_sample_aa = aa_mutations(translated,cds)
gene_names, position_on_gene_nt, position_on_gene_aa = get_gene()


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

aa_sum(df)
df.to_csv(output, index = False)
    