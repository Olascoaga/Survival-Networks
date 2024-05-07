import os
import pandas as pd

df_estadisticas = pd.read_csv("tissue_statistics.csv")
archivos_sort = [archivo for archivo in os.listdir() if archivo.endswith('.sort')]

if not os.path.exists('redes'):
    os.makedirs('redes')

def filtrar_archivo_sort(archivo_sort, columna_tejido, columna_umbral, tipo):
    nombre_tejido = archivo_sort.split('_')[0]
    umbral = df_estadisticas.loc[df_estadisticas['Tejido'] == nombre_tejido, columna_umbral].values[0]
    df_sort = pd.read_csv(archivo_sort, sep='\t', header=None)
    df_filtrado = df_sort[df_sort[2] > umbral]
    nombre_archivo_filtrado = f"redes/{nombre_tejido}_{tipo}_filtrado.sort"
    df_filtrado.to_csv(nombre_archivo_filtrado, sep='\t', index=False, header=False)
    print(f"Archivo {nombre_archivo_filtrado} guardado.")

for archivo_sort in archivos_sort:
    if "_Y_" in archivo_sort:
        filtrar_archivo_sort(archivo_sort, "Tejido", "Y", "Y")
    elif "_E_" in archivo_sort:
        filtrar_archivo_sort(archivo_sort, "Tejido", "E", "E")
