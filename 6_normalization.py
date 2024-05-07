import os
import pandas as pd
import qnorm

def filter_genes(df):
    """
    Filtra genes con bajos conteos.
    """
    df_filtered = df.loc[df.iloc[:, 1:].mean(axis=1) >= 10]
    df_filtered = df_filtered.loc[df_filtered.iloc[:, 1:].apply(lambda x: (x == 0).sum() / len(x), axis=1) < 0.5]
    return df_filtered

samples = pd.read_csv('SampleAttributes.csv') 
def filter_samples(tissue, samples, group):
    filtered_samples = samples[(samples['SMTSD'] == tissue) & (samples['GROUP'] == group)]
    use_samples = filtered_samples['SAMPID'].tolist()
    use_samples.insert(0, 'Description')
    return use_samples

processed_tsv_files = [file for file in os.listdir() if file.endswith(".tsv")]

for file_name in processed_tsv_files:
    # Leer el archivo TSV
    df = pd.read_csv(file_name, sep='\t', index_col=0)
    df = filter_genes(df)
    df = qnorm.quantile_normalize(df, axis=1)
    df = df.round(3)
    # Partir en E y Y
    tissue = file_name.split(".")[0]
    young_samples = filter_samples(tissue, samples, 'Y')
    elderly_samples = filter_samples(tissue, samples, 'E')
    
    # Dejar solo las columnas correspondientes a las muestras jóvenes
    df_young = df.loc[:, young_samples[1:]]
    name = f"{tissue}_Y.tsv"
    df_young.to_csv(name, sep='\t', index=True)
    df_elderly = df.loc[:, elderly_samples[1:]]
    name = f"{tissue}_E.tsv"
    df_elderly.to_csv(name, sep='\t', index=True)
