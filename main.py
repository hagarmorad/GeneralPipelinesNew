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
from poilo import polio
import utils
from threading import Lock

dirs=['BAM','QC','CNS','CNS_5','alignment','QC/depth']

def parse_input():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input',required=True, help='fastq folder')
        parser.add_argument('-r','--reference' ,help='reference')
        parser.add_argument('-p','--process' ,help='number of processes', default=1) #number of processes to run in perallel on supported tasks
        parser.add_argument('-f', '--flu', action='store_true', help="influenza segements analysis") #store_true will store the argument as true
        parser.add_argument('-d', '--de_novo', action='store_true', help="de-novo analysis") #store_true will store the argument as true
        parser.add_argument('--polio', action='store_true', help="PolioVirus analysis") #store_true will store the argument as true
        parser.add_argument('-c', '--cmv', action='store_true', help="cetomegalovirus (human herpesvirus 5) analysis") #store_true will store the argument as true
        args = parser.parse_args()
        return args.reference, args.input, args.flu, args.de_novo, args.polio, args.cmv, int(args.process)
        
if __name__ == "__main__":
    
        
    mutex = Lock()
 
    utils.create_dirs(dirs)
    reference, fastq, flu_flag, de_novo_flag, polio_flag, cmv_flag, process = parse_input()
    if not fastq.endswith("/"):
        fastq = fastq+"/"
    if flu_flag:
        pipe = flu(reference,fastq)
    elif de_novo_flag:
        pipe = de_novo(reference,fastq)
        #run spades multiprocessing
        utils.run_mp(process, pipe.run_spades, pipe.r1r2_list)
        pipe.run_blast()
        sample_ref = pipe.choose_reference_filter_contigs()
        pipe.import_references(sample_ref)
    elif polio_flag:
        pipe = polio(reference,fastq)
    else:
        pipe = general_pipe(reference,fastq)

    if polio_flag:
        #filterreads - keep only polio read 
        mutex.acquire()
        utils.run_mp(process, pipe.filter_not_polio, pipe.r1r2_list) #temp comment
        #pipe.filter_not_polio(pipe.r1r2_list[0]) #for debugging
        pipe.fastq = pipe.fastq + "polio_reads/"
        mutex.release()
        

    #mapping multiprocessing
    mutex.acquire()
    utils.run_mp(process, pipe.bam, pipe.r1r2_list)#temp comment
    mutex.release()
#   #pipe.bam(pipe.r1r2_list[0]) #for debuging
    
    if polio_flag:
        pipe.map_bam()
    
    if flu_flag:    
        utils.split_bam("BAM/")
    
    
    pipe.cns_depth("BAM/","QC/depth/","CNS/","CNS_5/") #temp comment
    
    if flu_flag:
        pipe.concat_samples()
        pipe.concat_segments()
    
    if flu_flag or cmv_flag: #need to fix fl
        utils.mafft(reference, "alignment/all_not_aligned.fasta", "alignment/all_aligned.fasta")
        
    pipe.results_report("BAM/", "QC/depth/", 'QC/report') #temp comment
    