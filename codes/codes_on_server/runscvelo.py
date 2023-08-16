import scvelo as scv
import anndata as anndata
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys

scv.set_figure_params('scvelo')
scv.logging.print_version()
sc.logging.print_versions()

output_dir = "../data/output/scvelo/"

if not os.path.exists(output_dir):
  os.mkdir(output_dir)

# data from R
plot_dpi = 600
cells = pd.read_csv(output_dir + "cells.csv")
print(cells)
genes = pd.read_csv(output_dir + "genes.csv")
print(genes)
sisters = pd.read_csv(output_dir + "sisters.csv", index_col = "cell")
drugSens = pd.read_csv(output_dir + "DrugSens.csv", index_col = "cell")
geneMeta = pd.read_csv(output_dir + "gene_meta.csv")

# 1. Loading data
# The 'loom' dataset is the output of 'Velocyto' 
experiment = ["ctrl_1", "ctrl_2", "NK_1", "NK_2", "Olaparib_1", "Olaparib_2", "Carbo_1", "Carbo_2"]

adata = scv.read_loom("../../2021.11.26_KURAMOCHI_AML_lineage_tracing/data/raw/velocyto/scLT-2-pre-ctrl-1.loom", 
                        cleanup = True, sparse = True, var_names = "Accession", validate = False)
adata.obs_names = [x.replace("x", "") for x in adata.obs_names]
adata.var_names_make_unique()
for e in experiment[1:]:
  sample = "scLT-2-pre-" + e
  sample = sample.replace("_", "-")
  tmp = scv.read_loom(
    "../../2021.11.26_KURAMOCHI_AML_lineage_tracing/data/raw/velocyto/" +
    sample + ".loom", 
    cleanup = True, sparse = True, var_names = "Accession", validate = False)
  tmp.obs_names = [x.replace("x", "") for x in tmp.obs_names]
  tmp.var_names_make_unique()
  adata = adata.concatenate(tmp, join='outer') # scv.utils.merge(adata, tmp)

# 2. Filter cells and genes wich have been removed by QC in Seurate
adata.obs_names = ["".join(x.replace("scLT-2-pre-", "").replace("scLT-", "").replace("BT-", "").split("-")[:2]) for x in adata.obs_names]
print(adata.obs_names)
print(adata.var_names)
adata = adata[cells.x, genes.x]

# 3. Filter on gene level and calculate moments

# We filtered genes by setting the minimum number of counts (both unspliced and 
# spliced) as 20. There are a lot of other parameters for filtering data provided 
# by `filter_and_normalize` function. You can find more details 
# [here](https://scvelo.readthedocs.io/scvelo.pp.filter_and_normalize.html#scvelo.pp.filter_and_normalize).

scv.pp.filter_and_normalize(adata, min_shared_counts = 20) # change the parameters for filtering
scv.pp.moments(adata, n_neighbors = 30) # change the parameters for calculate moments

# 4. Calculating RNA velocity

## 1. Recover and calculate dynamics

# We use the "dynamical" model to calculate the RNA velocity. You can find the 
# discription of models used by scVelo in [this page](https://scvelo.readthedocs.io/about.html)

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')

## 2. RNA velocity confidence

scv.tl.velocity_graph(adata)
scv.tl.umap(adata)
scv.tl.velocity_confidence(adata)
scv.tl.tsne(adata)
scv.pl.scatter(adata, color='velocity_confidence', perc=[2,98], dpi = plot_dpi, 
  save = "all_BT_velocity_confidence.png")

# 3. Calculate latent time

scv.tl.recover_latent_time(adata)

# 5. Calculate cell-to-cell transition probability matrix

# Formula for cell-to-cell transition probabilities:
# $$ \widetilde{\pi}_{ij} = \frac{1}{z_i}exp(\pi_{ij}/\sigma) $$
# From the velocity graph $\pi_{ij}$, with row-normalization $z_i$ and kernel 
# width $\sigma$ (scale parameter $\lambda = \sigma^{-1}$).

cell_to_cell = scv.utils.get_transition_matrix(adata)

# 6. Label with sisters

# adata.obs["sister_in_BT"] = sisters['sisters']
# scv.pl.velocity_embedding_stream(adata, basis='umap', color = "sister_in_BT", 
#   dpi = plot_dpi, color_map = "gist_ncar_r", 
#   save = "all_BT_velocity_sister_umap.png")


# 7. lable with pre-resistant and pre-sensitive cells in BT

adata.obs["drug_sensitive"] = drugSens['drugSens']
adata.obs['drug_sensitive'] = adata.obs['drug_sensitive'].astype('category')
adata.obs["drug_sensitive"]

# scv.pl.velocity_embedding_stream(adata, basis='umap', 
#     groups = ["pre-resistant", "pre-sensitive"],
#     color = "drug_sensitive", 
#     dpi = plot_dpi, 
#     legend_loc = "right margin",
#     save = "all_BT_velocity_drug_sensitive_umap.png")

## 8. Save files

adata.write(output_dir + 'all_BT_velocity_object.h5ad')
set(adata.obs.drug_sensitive)

# adata.write_loom(output_dir + 'RNA_velocity_BT.loom') # ValueError: INF and NaN not allowed in loom matrix
obs = pd.DataFrame(data=adata.obs)
obs.to_csv(output_dir + 'RNA_velocity_obs_all_BT.csv')
velocity = pd.DataFrame(data=adata.layers['velocity'], index=adata.obs_names, 
  columns=adata.var_names)
velocity.dropna(axis=1, how='all').to_csv(output_dir + 'RNA_velocity_all_BT.csv')
unspliced = pd.DataFrame.sparse.from_spmatrix(adata.layers['unspliced'], 
  index=adata.obs_names, columns=adata.var_names)
unspliced.to_csv(output_dir + 'RNA_velocity_unspliced_all_BT.csv')
spliced = pd.DataFrame.sparse.from_spmatrix(adata.layers['spliced'], 
  index=adata.obs_names, columns=adata.var_names)
spliced.to_csv(output_dir + 'RNA_velocity_spliced_all_BT.csv')

transition = pd.DataFrame(adata.uns["velocity_graph"].todense(), 
  index=adata.obs_names, columns=adata.obs_names)
transition.to_csv(output_dir + 'RNA_velocity_transition_probability_all_BT.csv')

transition = pd.DataFrame(adata.var)
transition.to_csv(output_dir + 'RNA_velocity_var_all_BT.csv')

for e in experiment:
  e = e.replace("_", "")
  cells = [i for i in adata.obs_names if i.startswith(e)]
  pd.DataFrame.sparse.from_spmatrix(cell_to_cell[cells, cells], index=cells, 
    columns=cells).to_csv(output_dir + 'RNA_velocity_cell_to_cell_transition_probability_' + e + '.csv')

