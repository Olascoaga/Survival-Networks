import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joypy


sns.set(font="Arial")
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=2)

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

tissues_labels.reverse()
young = filtering(data, 'Y',)
young.columns = tissues_labels
elderly = filtering(data, 'E')
elderly.columns = tissues_labels

young = young.reset_index()
young = young.rename(columns={'index': 'Tissue'})
elderly = elderly.reset_index()
elderly = elderly.rename(columns={'index': 'Tissue'})

young_melted = pd.melt(young, id_vars='Tissue', var_name='Gene', value_name='Young')
elderly_melted = pd.melt(elderly, id_vars='Tissue', var_name='Gene', value_name='Elderly')
result = pd.merge(young_melted, elderly_melted, on=['Gene', 'Tissue'], how='outer')
result = result[['Gene', 'Tissue', 'Young', 'Elderly']]

# Draw Plot
plt.figure(figsize=(20,10), dpi= 1200)
fig, axes = joypy.joyplot(result, column=['Young', 'Elderly'], by="Gene",
                          range_style='own',
                          linewidth=1,
                          fill=True,
                          overlap = 0.5,
                          legend=False,
                          grid='y',
                          fade=True,
                          figsize=(9,10))

plt.xlabel('log10(TPM + 1)')
plt.savefig('joyplot.png', dpi=1200) 
plt.show()

