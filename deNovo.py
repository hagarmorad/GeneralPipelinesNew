#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:29:54 2022

@author: hagar

this script gets fastq path of paired end read with no reference.
it runs rna SPAdes and saves the output for each sample in spades/spades_results.
in rnaSPAdes the contigs + scalffolds are in "transcript.fasta".

"""
import utils
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import shutil
from generalPipeline import general_pipe, INDEX, MAPPED_BAM, SORT, DEPTH, CNS, CNS5


#alignment conmmands
BWM_MEM_CONTIGS = "bwa mem -v1 -t 32 %(reference)s %(sample_fasta)s | samtools view -@ 32 -b - > %(output_path)s%(sample)s.bam"
BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b - > %(output_path)s%(sample)s.bam"


#database = "/mnt/data3/code/bard/references/virome/refseq_new_hagar/viral.1.1.genomic.fna"
RUN_SPADES = "spades -1 %(r1)s -2 %(r2)s -o %(output_path)s --rna"
RUN_BLAST = "blastn -db %(db)s -query %(query)s -out %(output_file)s -outfmt \"6 qseqid sseqid stitle qlen qcovs score \""


class de_novo(general_pipe):
      
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq) 
        
    #run spades on paired end fastq files. the path should contain R1 and R2 fastq files.
    def run_spades(self):
        utils.create_dirs(["spades/spades_results"])
        for sample, r1, r2 in utils.get_r1r2_list(self.fastq):
            utils.create_dirs(["spades/spades_results/"+sample])
            subprocess.call(RUN_SPADES % dict( r1=self.fastq + r1, r2=self.fastq + r2, output_path="spades/spades_results/"+ sample + "/"), shell=True)
    
    #run blastn on spades rna trascript.fasta  
    def run_blast(self):
        utils.create_dirs(["blast"])
        for spades_result in os.listdir("spades/spades_results/"):
            subprocess.call(RUN_BLAST % dict(db=self.reference, query="spades/spades_results/" + spades_result + "/transcripts.fasta", output_file="blast/" + spades_result + ".txt"), shell=True)
            
    #load blast output and choose the reference of the longest highest score contig. 
    def choose_reference_filter_contigs(self):
        utils.create_dirs(["fasta/","fasta/selected_contigs","fasta/all_contigs","fasta/selected_references"])
        sample_ref = {} #dict of the selected reference for each sample
        for sample, r1, r2 in utils.get_r1r2_list(self.fastq):   
           df = pd.read_csv('blast/'+sample+'.txt', sep="\t", header=None)
           df.columns = ["contig_seq-id", "reference_seq-id", "reference_title", "contig_length" ,"coverage(contig_to_ref)","score"] 
           #filter the highest score of each contig
           approved_fields = df.groupby("contig_seq-id")['score'].max()
           df = df[df['score'].isin(approved_fields)]
           df.to_csv('blast/filtered_'+sample+".csv", index=False)
           
           #create multi-fasta contains only the contigs blast found
           contigs_list = df['contig_seq-id'].tolist()
           #write fasta file with the selected contigs
           filtered_file = open("fasta/selected_contigs/"+sample+".fasta", 'a')
           #move spades contigs 
           shutil.copyfile("spades/spades_results/"+sample+"/transcripts.fasta","fasta/all_contigs/" +sample+".fasta")
           
           for record in SeqIO.parse("fasta/all_contigs/"+sample+".fasta", "fasta"):
               if record.description in contigs_list:
                   filtered_file.write(record.format("fasta"))
           filtered_file.close()
       #shutil.rmtree("spades/spades_results/")
            #select the reference of the longest contig
           reference = df.loc[df['contig_length'].idxmax()]['reference_seq-id'].split("|")[1].split("|")[0]
           sample_ref.update({sample:reference})
        return sample_ref
    
    #create reference fasta 
    def import_references(self, sample_ref):
        for sample, ref in sample_ref.items():
            for record in SeqIO.parse(self.reference, "fasta"):
                if ref in record.description:
                    ref_file = open("fasta/selected_references/"+sample+".fasta", 'a')
                    SeqIO.write(record,ref_file,"fasta")
                    ref_file.close()
                    continue
    #override
    #map each sample to its reference
    def bam(self,output_path):
        utils.create_dirs(["BAM/contig_based","BAM/fastq_based"])
        #align selected contig 
        for sample, r1, r2 in utils.get_r1r2_list(self.fastq):  
            fasta = sample + ".fasta"
            subprocess.call(INDEX % dict(reference="fasta/selected_references/"+fasta), shell=True)
            subprocess.call(BWM_MEM_CONTIGS % dict(reference="fasta/selected_references/"+fasta, sample_fasta="fasta/selected_contigs/"+fasta, sample=sample, output_path="BAM/contig_based/"), shell=True) #map to reference
            subprocess.call(BWM_MEM_FASTQ % dict(reference="fasta/selected_references/"+fasta ,r1=self.fastq+r1, r2=self.fastq+r2, sample=sample, output_path="BAM/fastq_based/"), shell=True) #map to reference
            subprocess.call(MAPPED_BAM % dict(sample=sample, output_path="BAM/fastq_based/"), shell=True) #keep mapped reads
            subprocess.call(MAPPED_BAM % dict(sample=sample, output_path="BAM/contig_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample, output_path="BAM/fastq_based/"), shell=True)         
            subprocess.call(SORT % dict(sample=sample, output_path="BAM/contig_based/"), shell=True)         
    
    #override
    #run general pipeline function twice (contigs based and fastq based)
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        utils.create_dirs([depth_path+contig_dir, depth_path+fastq_dir, cns_path+contig_dir, cns_path+fastq_dir, cns5_path+contig_dir,cns5_path+fastq_dir])
        super().cns_depth(bam_path+contig_dir, depth_path+contig_dir, cns_path+contig_dir, cns5_path+contig_dir)
        super().cns_depth(bam_path+fastq_dir, depth_path+fastq_dir, cns_path+fastq_dir, cns5_path+fastq_dir)
       
        
    #override
    #write report twice (contigs based and fastq based)
    def results_report(self, bam_path, depth_path, output_report):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"        
        super().results_report(bam_path+contig_dir, depth_path+contig_dir, output_report+"_contig_based")
        super().results_report(bam_path+fastq_dir, depth_path+fastq_dir, output_report+"_fastq_based")
        
        

        
        
        
        