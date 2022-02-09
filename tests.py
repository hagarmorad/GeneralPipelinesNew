#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:30:20 2022

@author: hagar
"""
import subprocess
INDEX = "bwa index %(reference)s"
subprocess.call(INDEX % dict(reference="/home/hagar/covid19/REF_NC_045512.2.fasta"), shell=True)

import pysam
samfile = pysam.AlignmentFile("BAM/14000.bam", "rb")
pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
for read in samfile.fetch():
    if read.is_paired:
        pairedreads.write(read)

pairedreads.close()
samfile.close()