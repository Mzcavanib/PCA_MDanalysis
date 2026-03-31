import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import networkx as nx
import pandas as pd

# --- Paso 1: Universo y selección ---
u = mda.Universe("final.gro", "final.xtc")
ca = u.select_atoms("protein and name CA")

G = nx.Graph()
cutoff = 4.5  # Å

# --- Paso 1: Construcción del grafo de contactos ---
for ts in u.trajectory:
    coords = ca.positions
    dist_matrix = distance_array(coords, coords)
    for i in range(len(ca)):
        for j in range(i+1, len(ca)):
            if dist_matrix[i, j] <= cutoff:
                res_i = ca[i].resid
                res_j = ca[j].resid
                if G.has_edge(res_i, res_j):
                    G[res_i][res_j]["weight"] += 1
                else:
                    G.add_edge(res_i, res_j, weight=1)

# --- Paso 2: Cargar contribuciones de PCA ---
df = pd.read_csv("residue_contributions.csv")

# --- Paso 3: Añadir atributos de PCA a los nodos ---
for _, row in df.iterrows():
    resid = row["resid"]
    if resid in G.nodes:
        G.nodes[resid]["resname"] = row["resname"]
        G.nodes[resid]["PC1_contrib"] = row["PC1_contrib"]
        G.nodes[resid]["PC2_contrib"] = row["PC2_contrib"]

# --- Exportar ---
nx.write_graphml(G, "residue_network_with_pca.graphml")
print("Archivo 'residue_network_with_pca.graphml' generado.")

