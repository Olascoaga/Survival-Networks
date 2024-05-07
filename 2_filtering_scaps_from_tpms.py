# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 08:03:29 2024

@author: olask
"""

import pandas as pd

with open('SCAPs.txt', 'r') as file:
    targets = file.read().splitlines()
    
sample_attributes = pd.read_csv('SampleAttributes.csv')
df = pd.read_csv('TPMs.csv', index_col=1)
df = df[df.index.isin(targets)]
df.to_csv('SCAPs_TPMS.csv')
