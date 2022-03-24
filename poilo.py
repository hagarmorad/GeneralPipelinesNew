#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:16:08 2022

@author: hagar
"""


RUN_SPADES = "spades --12 %(sample)s -o %(output_path)s --rna"
FILTER = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b -f 2 - | samtools fastq > %(output_path)s%(sample)s.fastq"
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 -B 10000 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
#BWM_MEM_CONTIGS = "bwa mem -v1 -t 32  %(reference)s %(sample_fasta)s | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
BWM_MEM_CONTIGS = "bwa mem -v1 -t 32  %(reference)s %(sample_fasta)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
import utils
from generalPipeline import  MAPPED_BAM, SORT
from deNovo import de_novo
import subprocess
import os
from utils import SPLIT

class polio(de_novo):
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq) 
        #utils.create_dirs([self.fastq+"polio_reads","BAM/fastq_based", "BAM/contig_based"]) #temp comment
        utils.create_dirs(["BAM/fastq_based", "BAM/contig_based"])

    #filter the fastq files to contain only mapped reads 
    def filter_not_polio(self,sample_r1_r2):
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2]       
        subprocess.call(FILTER % dict( reference=self.reference, r1=self.fastq + r1, r2=self.fastq + r2, output_path=self.fastq+"polio_reads/", sample=sample), shell=True)
        print(1)
    
        
    def run_spades(self,sample_r1_r2):
        sample = sample_r1_r2[0]
        utils.create_dirs(["spades/spades_results/"+sample])
        subprocess.call(RUN_SPADES % dict( sample=self.fastq+sample+".fastq", output_path="spades/spades_results/"+ sample + "/"), shell=True)
              
    #override
    #map each sample to its reference
    def bam(self,sample_r1_r2):
        sample = sample_r1_r2[0]
        subprocess.call(BWM_MEM_CONTIGS % dict(reference=self.reference, sample_fasta="spades/spades_results/"+ sample + "/"+"transcripts.fasta", out_file=sample, output_path="BAM/contig_based/"), shell=True) #map to reference
        subprocess.call(BWM_MEM_FASTQ % dict(reference=self.reference ,fastq=self.fastq+sample+".fastq", out_file=sample, output_path="BAM/fastq_based/"), shell=True) #map to reference
        subprocess.call(SPLIT % dict(bam="BAM/fastq_based/"+sample+".bam"), shell=True)
        os.remove("BAM/fastq_based/"+sample+".bam")
        subprocess.call(SPLIT % dict(bam="BAM/contig_based/"+sample+".bam"), shell=True)
        os.remove("BAM/contig_based/"+sample+".bam")
    
    def map_bam(self):
        for bam in os.listdir("BAM/contig_based/"):
            sample_ref = bam.split(".bam")[0]
            subprocess.call(MAPPED_BAM % dict(sample=sample_ref, output_path="BAM/contig_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample_ref, output_path="BAM/contig_based/"), shell=True)         
        
        for bam in os.listdir("BAM/fastq_based/"):
            sample_ref = bam.split(".bam")[0]
            subprocess.call(MAPPED_BAM % dict(sample=sample_ref, output_path="BAM/fastq_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample_ref, output_path="BAM/fastq_based/"), shell=True) 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    