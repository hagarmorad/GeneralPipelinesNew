#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 09:49:42 2022

@author: hagar
"""

import os
import sys 

spades_path = sys.argv[1]

for file in os.listdir(spades_path):
    if "