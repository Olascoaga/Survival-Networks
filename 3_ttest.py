# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:42:48 2024

@author: olask
"""

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

sample_attributes = pd.read_csv('SampleAttributes.csv')
sample_attributes = sample_attributes[sample_attributes['GROUP'] != 'M']
data = pd.read_csv('SCAPS_TPM.csv', index_col=0)
data = data.sort_index()

with open("sample_exclusion.txt", "r") as exclusion_file:
    exclusion = exclusion_file.readlines()
    exclusion = list(map(lambda s: s.strip(), exclusion))

def filtering(data, age):
    group_samples = sample_attributes[sample_attributes['GROUP'] == age]['SAMPID']
    group_data = data[data.columns.intersection(group_samples)]
    tissues = sorted(sample_attributes['SMTSD'].unique().tolist())
    tissues = [x for x in tissues if x not in exclusion]
    expression = pd.DataFrame(columns=list(data.index), index=tissues)
    data = data.T
    expression = expression.T
    for tissue in tissues:
        samples = sample_attributes[sample_attributes['SMTSD'] == tissue]
        samples = samples[samples['GROUP'] == age]
        samples = samples['SAMPID'].tolist()
        df = data[data.index.isin(samples)]
        df = np.log10(df.mean(axis=0) + 1).tolist()
        expression[tissue] = df
    return expression

young = filtering(data, 'Y')
elderly = filtering(data, 'E')

t_test_results = {}
for tissue in young.columns:
    t_stat, p_value = ttest_ind(young[tissue], elderly[tissue])
    t_test_results[tissue] = {'t_statistic': t_stat, 'p_value': p_value}

report = pd.DataFrame.from_dict(t_test_results, orient='index')
significant = report[report['p_value'] < 0.05]

print("Resultados de la prueba t de Student entre Young y Elderly:\n")
print(report)
print("\nTejidos con diferencias significativas (p-value < 0.05):\n")
print(significant)
