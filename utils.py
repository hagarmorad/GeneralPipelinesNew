#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:03:13 2022

@author: hagar
"""
import os
import sh
import subprocess
import shutil
import multiprocessing as mp
from Bio import SeqIO

#mafft commands
FIX = "eval \"$(conda shell.bash hook)\""
NEXT_ACT = "conda activate --stack nextstrain"
AUGUR = "  augur align --sequences alignment/all_not_aligned.fasta --reference-sequence %(reference)s --output alignment/all_aligned.fasta"
MAFFT = "mafft %(not_aligned)s > %(aligned)s"
DEACT ="conda deactivate"
#split bam file
SPLIT = "bamtools split -in %(bam)s -reference"

#return list of touple (R1,R2) file names
def get_r1r2_list(fastq_path):
    r1r2_list = []
    skip_files=["R2", "Undetermined", "unpaired", "singletons"]
    for r1 in os.listdir(fastq_path):
        if any(skip_file in r1 for skip_file in skip_files) or "fast" not in r1:
            continue
        r2 = r1.replace("R1","R2")
        sample = r1.split("_")[0].split(".fastq")[0] #sample short name
        r1r2_list.append([sample,r1,r2])
    return r1r2_list
 
#split all bam files by segments
def split_bam(dir):
    for bam_file in os.listdir(dir):        
        if "sorted" in bam_file and "bai" not in bam_file:
            subprocess.call(SPLIT % dict(bam=dir+bam_file), shell=True)
            os.remove(dir + bam_file)
  
#removes dir is exists and recreates it
def create_dirs(dirs):
        for dir in dirs:
            if os.path.exists(dir):
                shutil.rmtree(dir)
            os.makedirs(dir)
            
def remove_from_name(dir, to_remove):
    for file in os.listdir(dir):
        file = dir + file
        os.rename(file, file.replace(to_remove, ''))

#change first line in all files in a directory to ">sample_name"
def change_header(dir):
    for file in os.listdir(dir):
        new_header = ">" + file.split(".fa")[0]
        file = dir + file        
        sh.sed("-i", "1s/.*/" + new_header + "/", file)
        

def mafft(reference,not_aligned, aligned):
#    subprocess.call(FIX, shell=True)  
#    subprocess.call(NEXT_ACT, shell=True)   
#    subprocess.call(AUGUR % dict(reference=reference), shell=True)
#    subprocess.call(DEACT, shell=True)
    subprocess.call(MAFFT % dict(not_aligned=not_aligned, aligned=aligned), shell=True)

def get_sequences(alignment_file):
    sequences = {}
    alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
    for sample, record in alignment.items():
        sequences[sample] = str(record.seq).upper()

        sequences[sample.split("-")[0]] = sequences.pop(sample)

    return sequences




#get mutations positions list by comparing al sequences to each other.
#"-" and "N" are excluded.
def mutations_positions(sequences, no_n = 1):
    mutations_positons = []
    seq_length = len(next(iter(sequences.values())))
    for pos in range(seq_length-1):
        temp = ""
        for sample, record in sequences.items():
            if not temp:
                temp = record[pos]
            if no_n and record[pos] in ["-","N"]:
                break
            if not temp == record[pos] :
                mutations_positons.append(pos)
                break
        continue
    return mutations_positons

def run_mp(threads, func, arg):
    with mp.Pool(threads) as pool:
        pool.map(func,arg)
        pool.close()
        pool.join()