#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:26:50 2022

@author: hagar
"""
import os
import subprocess
import utils
fastq_path = "/mnt/data/projects/Merav/VDPV3_2022/fastq/raw/"
polio_refs = "/home/hagar/refs/polio/polio.fasta"

FILTER = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b -F 2 - | samtools fastq > %(output_path)s%(sample)s.fastq"
FILTER2 = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s| samtools view -@ 32 -b -f 4 -f 8 - | samtools fastq > %(output_path)s%(sample)s.fastq"
RUN_SPADES = "spades --12 %(fastq)s -o %(output_path)s --rna"
def filter():
    for fastq in os.listdir(fastq_path):
        if "R1" in fastq:
            r1 = fastq_path+fastq
            r2 = r1.replace("R1","R2")
            sample = fastq.split("_")[0]
            subprocess.call(FILTER % dict( reference=polio_refs, r1=r1, r2=r2, output_path=fastq_path+"test/", sample=sample), shell=True)
        
#run spades on paired end fastq files. the list must contain sample, r1 r2 file names 
def run_spades(fastq):
        sample = fastq.split(".fastq")[0]
        utils.create_dirs(["spades/spades_results/"+sample])
        subprocess.call(RUN_SPADES % dict( fastq=fastq_path+fastq, output_path="/mnt/data/projects/Merav/girl2/no_polio_de_novo_analysis/spades/spades_results/"+ sample + "/"), shell=True)
         
def filter_unmapped_with_no_mate(fastq):
    sample = fastq.split(".fastq")[0]
    r1 = fastq_path+fastq
    r2 = r1.replace("R1","R2")
    subprocess.call(FILTER2 % dict( reference=polio_refs, r1=r1, r2=r2, output_path=fastq_path+"no_polio_reads/unmapped_pairs/", sample=sample), shell=True)
        
        
#filter()
fastq_list=[]
for fastq in os.listdir("/mnt/data/projects/Merav/girl2/fastq/no_polio/"):
    fastq_list.append(fastq)
#utils.run_mp(5,run_spades,fastq_list)
#filter_unmapped_with_no_mate(fastq_list[0])
#utils.run_mp(5,filter_unmapped_with_no_mate,fastq_list)
filter()