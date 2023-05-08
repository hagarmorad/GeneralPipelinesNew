#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:52:43 2022

@author: hagar
"""

import argparse
from pipelines.generalPipeline import general_pipe
from viruses.flu import flu
from pipelines.deNovo import de_novo
from viruses.polio import polio
from utils import utils, parse_gb_file
from threading import Lock
from mutations import signatures

dirs=['BAM','QC','CNS','CNS_5','QC/depth']

def parse_input():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', help='fastq folder')
        parser.add_argument('-r','--reference' ,help='reference')
        parser.add_argument('-p','--process' ,help='number of processes', default=1) #number of processes to run in perallel on supported tasks
        parser.add_argument('-mr','--multi_ref', action='store_true', help='multi fasta reference') #
        parser.add_argument('-f', '--flu', action='store_true', help="influenza segements analysis") #store_true will store the argument as true
        parser.add_argument('-d', '--de_novo', action='store_true', help="de-novo analysis") #store_true will store the argument as true
        parser.add_argument('--polio', action='store_true', help="PolioVirus analysis") #store_true will store the argument as true
        parser.add_argument('-c', '--cmv', action='store_true', help="cetomegalovirus (human herpesvirus 5) analysis") #store_true will store the argument as true
        parser.add_argument('-gb','--gb_file' ,help='insert gb file and get reference regions report')
        parser.add_argument('-m','--mutations_table' , action='store_true',help='mutations table reprort. gb file flag is mandatory')
        parser.add_argument('-rg','--regions_file' ,help='insert gene regions file for mini')
        parser.add_argument('--mini' , action='store_true',help='run only mutation analysis. this flag requires --input flag as alignment file.')
        parser.add_argument('--skip_spades' , action='store_true',help='skip spades analysis used in polio and de novo classes. turn on this flag only if you already run spades once')
        parser.add_argument('-v','--vcf' , action='store_true',help='generates vcf files using gatk4. fill the identity column in the report file')
        parser.add_argument('--cnsThresh' ,help='Minimum frequency threshold(0 - 1) to call consensus. (Default: 0.6)', default=0.6)
        args = parser.parse_args()
        return args.reference, args.input, args.flu, args.de_novo, args.polio, args.cmv, \
                   int(args.process), args.gb_file, args.regions_file, args.mutations_table, \
                       args.mini, args.skip_spades, args.vcf, args.cnsThresh, args.multi_ref
        
if __name__ == "__main__":
    
    DEBUG = False
    mutex = Lock()
    
    reference, fastq, flu_flag, de_novo_flag, polio_flag, cmv_flag, \
        process, gb_file, regions_file, mutations_flag, mini, \
            skip_spades, vcf, cnsThresh, multi_ref_flag = parse_input()
    
    if not mini:      
        utils.create_dirs(dirs) 
        if not reference:
            raise ValueError("reference sequence is required.")
            
    if fastq and reference and not mini:
        if not fastq.endswith("/"):
            fastq = fastq+"/"
        
        if flu_flag:
            pipe = flu(reference,fastq)
        
        elif de_novo_flag:
            pipe = de_novo(reference,fastq) #temp comment
            if not skip_spades:
                #run spades multiprocessing
                mutex.acquire()
                utils.run_mp(process, pipe.run_spades, pipe.r1r2_list)
                mutex.release()
            pipe.run_blast()
            sample_ref = pipe.choose_reference_filter_contigs()
            pipe.import_references(sample_ref)
        
        elif polio_flag:
            pipe = polio(reference,fastq)
        
        else:
            pipe = general_pipe(reference,fastq)
        
        
       
        if polio_flag:
            #filter reads - keep only polio read 
            mutex.acquire()
            utils.run_mp(process, pipe.filter_not_polio, pipe.r1r2_list) if not DEBUG else pipe.filter_not_polio(pipe.r1r2_list[0])
            pipe.fastq = pipe.fastq + "polio_reads/"
            mutex.release()
            
            #run spades
            if not skip_spades:
                mutex.acquire()
                utils.run_mp(process, pipe.run_spades, pipe.r1r2_list)
                mutex.release()
            
    
        #mapping multiprocessing
        mutex.acquire()
        utils.run_mp(process, pipe.bam, pipe.r1r2_list) if not DEBUG else pipe.bam(pipe.r1r2_list[0])
        mutex.release()
        
        if polio_flag:
            pipe.map_bam()
        
        if flu_flag or multi_ref_flag:    
            utils.split_bam("BAM/")
        

        pipe.cns_depth("BAM/","QC/depth/","CNS/","CNS_5/", cnsThresh) #temp comment
        
        if flu_flag:
            pipe.concat_samples()
            pipe.concat_segments()
        
        if vcf:
            utils.create_dirs(["VCF"])
            mutex.acquire()
            if polio_flag or de_novo_flag:
                pipe.variant_calling("BAM/contig_based/", "VCF/contig_based/")
                pipe.variant_calling("BAM/fastq_based/", "VCF/fastq_based/")
            else:    
                pipe.variant_calling()
            mutex.release()
            pipe.results_report("BAM/", "QC/depth/", 'QC/report', vcf=1) #temp comment
            
        else:    
            pipe.results_report("BAM/", "QC/depth/", 'QC/report') #temp comment
        
    if cmv_flag or (mutations_flag and not mini):  # TODO - fix flu - im not sure this aligner fits
        utils.create_dirs(['alignment'])
        utils.mafft(reference, "alignment/all_not_aligned.fasta", "alignment/all_aligned.fasta")
    
    if flu_flag:
        utils.create_dirs(['alignment'])
        pipe.mafft()
                
        
    if gb_file:
        parse_gb_file.parse(gb_file)
        regions_file = gb_file.replace(".gb", "_regions.csv")

    if mutations_flag or mini:
        if not gb_file and not regions_file:
            raise ValueError("gene bank file or regions file is required.")
        utils.create_dirs(["reports"])
        aligned = fastq if mini else "alignment/all_aligned.fasta" # if mini flag is on, user must insert an alignment file instead of fastq.
        signatures.run(aligned, regions_file, "reports/mutations.xlsx")
        