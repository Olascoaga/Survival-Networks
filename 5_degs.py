# -*- coding: utf-8 -*-
"""
Created on Mon May  6 16:16:32 2024

@author: olask
"""

import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests


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
    return group_data, tissues

young, young_tissues = filtering(data, 'Y')
elderly, elderly_tissues = filtering(data, 'E')

pvals = pd.DataFrame(index=data.index, columns=young_tissues)

tissues = sample_attributes['SMTSD'].unique().tolist()
tissues = [x for x in tissues if x not in exclusion]

for gene in data.index:
    for tissue in tissues:
        samples = sample_attributes[sample_attributes['SMTSD'] == tissue]
        samples = samples['SAMPID'].tolist()

        young_samples = young.loc[[gene]]
        young_samples = young_samples.filter(items=samples)

        elderly_samples = elderly.loc[[gene]]
        elderly_samples = elderly_samples.filter(items=samples)

        elderly_array = elderly_samples.values.flatten()
        young_array = young_samples.values.flatten()

        statistic, p_value = mannwhitneyu(elderly_array, young_array, alternative='two-sided')

        pvals.loc[gene, tissue] = p_value


fc = pd.DataFrame(index=data.index, columns=young_tissues)
for gene in data.index:
    for tissue in tissues:
        samples = sample_attributes[sample_attributes['SMTSD'] == tissue]
        samples = samples['SAMPID'].tolist()

        young_samples = young.loc[[gene]]
        young_samples = young_samples.filter(items=samples)

        elderly_samples = elderly.loc[[gene]]
        elderly_samples = elderly_samples.filter(items=samples)

        young_mean = young_samples.loc[gene, :].mean()
        elderly_mean = elderly_samples.loc[gene, :].mean()

        fc.loc[gene, tissue] = (elderly_mean / young_mean)
        
fc = fc.astype(float)
fc = np.log2(fc + 0.01)
pvals = pvals.astype(float)

pvals_array = pvals.values.flatten()
reject_bh, pvals_bh_corrected, _, _ = multipletests(pvals_array, method='fdr_bh')
pvals_corrected = pd.DataFrame(pvals_bh_corrected.reshape(pvals.shape), index=pvals.index, columns=pvals.columns)
pvals = -np.log10(pvals_corrected)

x = np.arange(1, 40)
y = np.arange(1, 34)
X, Y = np.meshgrid(x, y)

tissues_labels = ["Adipose - Subcutaneous",
                  "Adipose - Visceral", 
                  "Adrenal Gland", 
                  "Artery - Aorta",
                  "Artery - Coronary",
                  "Artery - Tibial",
                  "Breast", 
                  "Cultured fibroblasts",
                  "EBV-transformed lymphocytes",
                  "Colon - Sigmoid",
                  "Colon - Transverse",
                  "Esophagus - Gastroesophageal",
                  "Esophagus - Mucosa",
                  "Esophagus - Muscularis",
                  "Heart - Atrial Appendage",
                  "Heart - Left Ventricle",
                  "Lung",
                  "Minor Salivary Gland",
                  "Muscle - Skeletal",
                  "Nerve - Tibial",
                  "Ovary",
                  "Pancreas",
                  "Prostate",
                  "Skin - Not Sun Exposed",
                  "Skin - Sun Exposed",
                  "Small Intestine",
                  "Spleen",
                  "Stomach",
                  "Testis",
                  "Thyroid",
                  "Uterus",
                  "Vagina",
                  "Whole Blood"]

fc.columns = tissues_labels
etiquetas_x = fc.index
etiquetas_y = fc.columns

fc_values_ordered = fc.values.ravel('F')
pvals_values_ordered = pvals.values.ravel('F')

# Graficar los datos
plt.figure(figsize=(12, 13), dpi=1200)  
scatter = sns.scatterplot(
    x=X.flatten(), 
    y=Y.flatten(), 
    hue=fc_values_ordered,  
    hue_norm=(-1, 1),
    palette='coolwarm',
    size=pvals_values_ordered,  
    size_norm=(0, 13),
    sizes=(20, 300), 
    legend=False, 
)  

plt.xticks(ticks=np.arange(1, 40), labels=etiquetas_x, fontsize=13, rotation='vertical')
plt.yticks(ticks=np.arange(1, 34), labels=etiquetas_y, fontsize=14)

plt.tight_layout()
plt.savefig('dges.png', dpi=1200) 
plt.show()