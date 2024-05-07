import pandas as pd

# Leer SampleAttributes.csv y filtrar muestras una sola vez
samples = pd.read_csv('SampleAttributes.csv') 
with open("sample_inclusion.txt", "r") as inclusion_file:
    inclusion = [line.strip() for line in inclusion_file]

def filter_samples(tissue, group, samples):
    filtered_samples = samples[(samples['SMTSD'].isin(tissue)) & (samples['GROUP'] == group)]
    use_samples = filtered_samples['SAMPID'].tolist()
    use_samples.insert(0, 'Description')
    return use_samples
    
groups = ['Y', 'E']

# Leer GTEx.gct una vez fuera del bucle
gtex_data = pd.read_csv("GTEx.gct", sep='\t', skiprows=2)

for tissue in inclusion:
    for group in groups:
        # Filtrar los datos según sea necesario dentro del bucle
        df = gtex_data[filter_samples([tissue], group, samples)]
        name = f"{tissue}_{group}.tsv"
        df.to_csv(name, sep='\t', index=False)
        del df  # Liberar memoria del DataFrame utilizado en cada iteración
