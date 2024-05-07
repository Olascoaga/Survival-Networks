from collections import defaultdict
import pandas as pd
from mpl_chord_diagram import chord_diagram
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'Arial'

def add_rank_column_to_variable(interactions):
    rank = len(interactions)  # Comenzar desde el rango más alto
    prev_interaction_strength = None
    interactions_with_rank = []
    for interaction in interactions:
        gene1, gene2, interaction_strength = interaction
        interaction_strength = float(interaction_strength)

        if interaction_strength != prev_interaction_strength:
            rank -= 1  # Disminuir el rango para la siguiente interacción más débil
            prev_interaction_strength = interaction_strength

        interactions_with_rank.append((gene1, gene2, interaction_strength, rank))
    
    return interactions_with_rank

def read_genes_from_file(genes_file):
    with open(genes_file, 'r') as f:
        genes = [line.strip() for line in f.readlines()]
    return genes

def extract_tissue_name(file_name):
    tissue_name = file_name.split('_')[0]
    return tissue_name

def generate_chord_diagram(sort_file_Y, sort_file_E, tissue_name):
    genes = read_genes_from_file("SCAPs.txt")

    with open(sort_file_Y, 'r') as f:
        interactions_Y = [line.strip().split('\t') for line in f.readlines()]

    interactions_with_rank_Y = add_rank_column_to_variable(interactions_Y)

    with open(sort_file_E, 'r') as f:
        interactions_E = [line.strip().split('\t') for line in f.readlines()]

    interactions_with_rank_E = add_rank_column_to_variable(interactions_E)
    interactions_matrix_Y = pd.DataFrame(index=genes, columns=genes)

    for interaction in interactions_with_rank_Y:
        gene1, gene2, interaction_strength, rank = interaction
        interactions_matrix_Y.at[gene1, gene2] = rank
        interactions_matrix_Y.at[gene2, gene1] = rank 

    interactions_matrix_Y = interactions_matrix_Y.apply(pd.to_numeric, errors='coerce')
    interactions_matrix_Y.fillna(0, inplace=True)
    interactions_matrix_Y = interactions_matrix_Y * 10
    interactions_matrix_E = pd.DataFrame(index=genes, columns=genes)

    for interaction in interactions_with_rank_E:
        gene1, gene2, interaction_strength, rank = interaction
        interactions_matrix_E.at[gene1, gene2] = rank
        interactions_matrix_E.at[gene2, gene1] = rank  # Debido a que es simétrica, se llena en ambas direcciones

    interactions_matrix_E = interactions_matrix_E.apply(pd.to_numeric, errors='coerce')
    interactions_matrix_E.fillna(0, inplace=True)

    interactions_matrix_E = interactions_matrix_E * 10

    fig, axs = plt.subplots(1, 2, figsize=(18, 9))
    chord_diagram(interactions_matrix_Y.values, 
                  names=interactions_matrix_Y.index, 
                  rotate_names=True, 
                  fontsize=17,
                  sort=None,
                  pad=5, 
                  alpha=0.6,
                  ax=axs[0])
    
    #axs[0].set_title(f"{tissue_name} - Y")
    chord_diagram(interactions_matrix_E.values, 
                  names=interactions_matrix_E.index, 
                  rotate_names=True, 
                  fontsize=17,
                  sort=None,
                  pad=5, 
                  alpha=0.6,
                  ax=axs[1])

    #axs[1].set_title(f"{tissue_name} - E")
    # Guardar la figura en formato PNG con 600 DPI
    plt.savefig(f'{tissue_name}_diagram.png', dpi=1200)
    # Mostrar la figura
    plt.show()

files_in_directory = os.listdir()

# Iterar sobre los archivos en el directorio
for file_name in files_in_directory:
    if file_name.endswith('_filtrado.sort'):
        tissue_name = file_name.split('_')[0]
        group_Y_file = f"{tissue_name}_Y_filtrado.sort"
        group_E_file = f"{tissue_name}_E_filtrado.sort"
        generate_chord_diagram(group_Y_file, group_E_file, tissue_name)
