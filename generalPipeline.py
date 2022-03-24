#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:37:02 2022

@author: hagar
"""

import subprocess
import os
import pysam
from statistics import mean
import csv
import utils

#alignment conmmands
INDEX = "bwa index %(reference)s"
BWM_MEM = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b - > %(output_path)s%(sample)s.bam"
MAPPED_BAM = "samtools view -@ 32 -b -F 4 %(output_path)s%(sample)s.bam > %(output_path)s%(sample)s.mapped.bam"
SORT = "samtools sort -@ 32 %(output_path)s%(sample)s.mapped.bam -o %(output_path)s%(sample)s.mapped.sorted.bam"
DEPTH = "samtools depth -a %(bam_path)s%(bam_file)s > %(depth_path)s%(sample)s.txt"
CNS = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t 0.6 -m 1 -p %(cns_path)s%(sample)s.fa"
CNS5 = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t 0.6 -m 5 -p %(cns5_path)s%(sample)s.fa"
SAMTOOLS_INDEX = "samtools index -@ 32 %(bam_path)s%(bam_file)s"
#MAFFT commands
ALL_NOT_ALIGNED =   "cat %(dir)s > alignment/all_not_aligned.fasta"

#report commands
BREADTH_CNS5 = "$(cut -f3 QC/depth/%(sample)s.txt | awk '$1>5{c++} END{print c+0}')"
#directories for output

            
class general_pipe():
    
    def __init__(self, reference, fastq):
        self.reference = reference
        self.fastq = fastq
        self.r1r2_list = utils.get_r1r2_list(self.fastq)
    #TODO  
    def bam(self, sample_r1_r2):
        subprocess.call(INDEX % dict(reference=self.reference), shell=True)
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2]            
        subprocess.call(BWM_MEM % dict(reference=self.reference,r1=self.fastq + r1, r2=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
        subprocess.call(MAPPED_BAM % dict(sample=sample, output_path="BAM/"), shell=True) #keep mapped reads
        subprocess.call(SORT % dict(sample=sample, output_path="BAM/"), shell=True)         
    
    #find mapping depth and consensus sequence 
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path):
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file:
                print(bam_path + bam_file)
                sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                subprocess.call(DEPTH % dict(bam_path=bam_path, bam_file=bam_file, depth_path=depth_path, sample=sample), shell=True) 
                #consensus
                subprocess.call(CNS % dict(bam_path=bam_path, bam_file=bam_file, cns_path=cns_path, sample=sample), shell=True) 
                subprocess.call(CNS5 % dict(bam_path=bam_path, bam_file=bam_file, cns5_path=cns5_path, sample=sample), shell=True)
                #remove qual files
                os.remove(cns_path+sample+".qual.txt")
                os.remove(cns5_path+sample+".qual.txt")

    def mafft(self):
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS_5/*"), shell=True)
        utils.mafft(self.reference)
    
    #write report.csv - mapping analysis
    def results_report(self, bam_path, depth_path, output_report):
        
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
    
        writer.writerow(['sample', 'reference', 'mapped%','mapped_reads','total_reads','cov_bases','coverage%','coverage_CNS_5%','mean_depth','max_depth','min_depth'])
    
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    subprocess.call(SAMTOOLS_INDEX % dict(bam_path=bam_path, bam_file=bam_file), shell=True) 
                    total_reads = pysam.AlignmentFile(bam_path+bam_file).count(until_eof=True)
                    coverage_stats = pysam.coverage(bam_path+bam_file).split("\t")
                    mapped_reads = int(coverage_stats[11])
                    mapped_percentage = round(mapped_reads/total_reads*100,4) if total_reads else ''
                    cov_bases =  int(coverage_stats[12])
                    coverage = float(coverage_stats[13])
            
                    #depth 
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    depths = [int(x.split('\t')[2]) for x in open(depth_path+sample+".txt").readlines()]
                    depths = [i for i in depths if i != 0]
                    mean_depth = str(round(mean(depths),3)) if depths else ''
                    min_depth = min(depths) if depths else ''
                    max_depth = max(depths) if depths else ''
                    
                    breadth_cns5 = len([i for i in depths if i > 5])
                    genome_size = sum(1 for line in open(depth_path+sample+".txt"))
                    coverage_cns5 = round(breadth_cns5 / genome_size * 100,3)  if genome_size else ''
                    
                    #f=open("fasta/selected_references/"+sample+".fasta")
                    #reference = f.readline().split("man ")[1].split("strain")[0]
                    #f.close()
                    writer.writerow([sample, reference, mapped_percentage, mapped_reads, total_reads, cov_bases, coverage, coverage_cns5, mean_depth, max_depth, min_depth])
                
        f.close()

