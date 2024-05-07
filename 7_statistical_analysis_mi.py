import os
import pandas as pd
from scipy.stats import mannwhitneyu

"""
.sort files are in networks directory
"""

directorio = '.'
datos = []

for archivo in os.listdir(directorio):
    if archivo.endswith('.sort'):
        nombre_tejido, grupo = archivo.split('_')[0], archivo.split('_')[1]
        
        df = pd.read_csv(os.path.join(directorio, archivo), header=None, sep='\s+', names=['Gen1', 'Gen2', 'Valor'])
        df['Tejido'] = nombre_tejido
        df['Grupo'] = grupo
        datos.append(df)

data = pd.concat(datos, ignore_index=True)

def filtering(df, group):
    filtered_df = df[df['Grupo'] == group]
    filtered_df = filtered_df.drop(filtered_df.columns[:2], axis=1)
    return filtered_df

young = filtering(data, 'Y')
elderly = filtering(data, 'E')


mann_whitney_results = {}
for tejido in data['Tejido'].unique():
    tejido_elderly = elderly[elderly['Tejido'] == tejido]['Valor']
    tejido_young = young[young['Tejido'] == tejido]['Valor']
    _, p_value = mannwhitneyu(tejido_elderly, tejido_young)
    mann_whitney_results[tejido] = p_value

reporte_mann_whitney = pd.DataFrame.from_dict(mann_whitney_results, orient='index', columns=['Valor de probabilidad'])

reporte_mann_whitney.to_csv('reporte_mann_whitney.csv')

# Imprimir el reporte
print("Reporte del test de Mann-Whitney U para comparar los valores de 'elderly' y 'young' en cada tejido:")
print(reporte_mann_whitney)
