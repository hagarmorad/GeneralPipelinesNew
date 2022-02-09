#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:52:43 2022

@author: hagar
"""

import argparse
from generalPipeline import general_pipe
from flu import flu
from deNovo import de_novo
import utils

dirs=['BAM','QC','CNS','CNS_5','alignment','QC/depth']

def parse_input():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input',required=True, help='fastq folder')
        parser.add_argument('-r','--reference' ,help='reference')
        parser.add_argument('-f', '--flu', action='store_true', help="influenza segements analysis") #store_true will store the argument as true
        parser.add_argument('-d', '--de_novo', action='store_true', help="de-novo analysis") #store_true will store the argument as true
        args = parser.parse_args()
        return args.reference, args.input, args.flu, args.de_novo        
        
if __name__ == "__main__":

    utils.create_dirs(dirs)
    reference, fastq, flu_flag, de_novo_flag= parse_input()
    if flu_flag:
        pipe = flu(reference,fastq)
    elif de_novo_flag:
        pipe = de_novo(reference,fastq)
        pipe.run_spades()     
        pipe.run_blast()
        sample_ref = pipe.choose_reference_filter_contigs()
        pipe.import_references(sample_ref)
    else:
        pipe = general_pipe(reference,fastq)
        
    pipe.bam("BAM/") 
    
    if flu_flag:    
        pipe.split_bam()
    
    pipe.cns_depth("BAM/","QC/depth/","CNS/","CNS_5/")
    
    if flu_flag:
        pipe.concat_samples()
        pipe.concat_segments()
    
    #pipe.mafft()
    pipe.results_report("BAM/", "QC/depth/", 'QC/report')
    