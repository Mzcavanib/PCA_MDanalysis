import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import pandas as pd

u = mda.Universe("final.gro", "final.xtc")
ca = u.select_atoms("protein and name CA")

pca_analysis = pca.PCA(u, select="protein and name CA").run()

try:
    p_components = pca_analysis.results.p_components
except AttributeError:
    p_components = pca_analysis.p_components

n_atoms = len(ca)
pc1 = p_components[:, 0].reshape(n_atoms, 3)
pc2 = p_components[:, 1].reshape(n_atoms, 3)

pc1_magnitudes = np.linalg.norm(pc1, axis=1)
pc2_magnitudes = np.linalg.norm(pc2, axis=1)

resids = ca.resids
resnames = ca.resnames
df = pd.DataFrame({
    "resid": resids,
    "resname": resnames,
    "PC1_contrib": pc1_magnitudes,
    "PC2_contrib": pc2_magnitudes
})

df.to_csv("residue_contributions.csv", index=False)
print("File 'residue_contributions.csv'.")

