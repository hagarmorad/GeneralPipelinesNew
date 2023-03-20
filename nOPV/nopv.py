#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 09:46:29 2023

@author: hagar
"""

import pandas as pd
from sys import argv

def run(muts_file):
    muts = pd.read_excel(muts_file)
    
    
    
if __name__ == "__main__":
    run(argv[1])