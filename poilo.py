#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:16:08 2022

@author: hagar
"""


RUN_SPADES = "spades --12 %(sample)s -o %(output_path)s --rna"
FILTER = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b -f 2 - | samtools fastq > %(output_path)s%(sample)s.fastq"
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 -B 10000 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #strict mapping (high penalty)
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #filter out common reads
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #no filters
BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
import utils
from generalPipeline import  MAPPED_BAM, SORT, general_pipe
from deNovo import de_novo
import subprocess
import os
from utils import SPLIT

class polio(general_pipe):
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq) 
        utils.create_dirs([self.fastq+"polio_reads"]) #temp comment

    #filter the fastq files to contain only mapped reads 
    def filter_not_polio(self,sample_r1_r2):
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2]       
        subprocess.call(FILTER % dict( reference=self.reference, r1=self.fastq + r1, r2=self.fastq + r2, output_path=self.fastq+"polio_reads/", sample=sample), shell=True)
        print(1)
    
              
    #override
    #map each sample to its reference
    def bam(self,sample_r1_r2):
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2] 
        #subprocess.call(BWM_MEM_FASTQ % dict(reference=self.reference ,fastq=self.fastq+sample+".fastq", out_file=sample, output_path="BAM/"), shell=True) #map to reference
        subprocess.call(BWM_MEM_FASTQ % dict(reference=self.reference ,r1=self.fastq+r1, r2=self.fastq+r2, out_file=sample, output_path="BAM/"), shell=True) #map to reference
        subprocess.call(SPLIT % dict(bam="BAM/"+sample+".bam"), shell=True)
        os.remove("BAM/"+sample+".bam")

    def map_bam(self):

        for bam in os.listdir("BAM/"):
            sample_ref = bam.split(".bam")[0]
            subprocess.call(MAPPED_BAM % dict(sample=sample_ref, output_path="BAM/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample_ref, output_path="BAM/"), shell=True) 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    