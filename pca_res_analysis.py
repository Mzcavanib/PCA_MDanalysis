import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# --- Archives ---
file_a = "residue_network_with_pca_A.graphml"
file_b = "residue_network_with_pca_B.graphml"

# --- Load Networks ---
G_a = nx.read_graphml(file_a)
G_b = nx.read_graphml(file_b)

# --- Convert nodes to DataFrame ---
df_a = pd.DataFrame.from_dict(dict(G_cmp.nodes(data=True)), orient="index")
df_a["Condition"] = "a"

df_b = pd.DataFrame.from_dict(dict(G_nocmp.nodes(data=True)), orient="index")
df_b["Condition"] = "b"

# --- Link ---
df_all = pd.concat([df_a, df_b])

# --- Residue ID ---
if "resid" in df_all.columns:
    df_all["ResidueID"] = df_all["resid"].astype(int)
else:
    df_all["ResidueID"] = df_all.index.astype(int)

# --- Figure with two panels ---
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Panel A: PC1
for cond, df in df_all.groupby("Condition"):
    df_sorted = df.sort_values("ResidueID")
    axes[0].plot(df_sorted["ResidueID"], df_sorted["PC1_contrib"], label=cond)
axes[0].set_ylabel("PC1")
axes[0].text(0.01, 0.95, "A", transform=axes[0].transAxes,
             fontsize=14, fontweight="bold", va="top", ha="left")
axes[0].legend()

# Panel B: PC2
for cond, df in df_all.groupby("Condition"):
    df_sorted = df.sort_values("ResidueID")
    axes[1].plot(df_sorted["ResidueID"], df_sorted["PC2_contrib"], label=cond)
axes[1].set_xlabel("Residuo")
axes[1].set_ylabel("PC2")
axes[1].text(0.01, 0.95, "B", transform=axes[1].transAxes,
             fontsize=14, fontweight="bold", va="top", ha="left")
axes[1].legend()

plt.tight_layout()
plt.show()

