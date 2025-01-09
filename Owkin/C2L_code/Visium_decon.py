import sys
import os 

# logging
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "w")

import scanpy as sc
import squidpy as sq
import cell2location as c2l
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import warnings
import torch

has_gpu = torch.cuda.is_available()

snakemake.params.ncell = 8

#### Read Visium sample:
print("First read in the data.")
_st_sample = sq.read.visium(snakemake.input.counts.split("/filtered_feature_bc_matrix")[0],
                            counts_file=snakemake.input.counts.split("/")[-1],
                            library_id=snakemake.wildcards.sample)

_st_sample.var_names_make_unique()

## Remove the genes that are not in scRNAseq data:
print("Now read in the signature matrix.")
inf_aver=pd.read_csv(snakemake.input.sig_mat, index_col=0, header=0)

gene_ids_sc = inf_aver.index.tolist()
gene_ids_visium = _st_sample.var.index.tolist()
intersect = list(set(gene_ids_sc).intersection(set(gene_ids_visium)))

print("Length of intersection of genes between Visium and scRNAseq")
print(len(intersect))
if len(intersect) == 0:
    raise ValueError("There is no intersection between Visium and scRNAseq gene lists.")
_st_sample.var

print("Filter AnnData object to get only intersection genes.")
st_sample = _st_sample[:, _st_sample.var.index.isin(intersect)]
print(st_sample)

print("Set up the cell state data frame.")
cell_state_df = inf_aver.loc[inf_aver.index.isin(intersect), :]

if not all(x.is_integer() for x in st_sample.X.toarray().flatten()):
    print("WARNING: Values in X will be round as SpotClean generated float values as raw counts, and is not currently handeled by cell2location.")
    st_sample.X = np.round(st_sample.X)
    
print("Run anndata setup function")
c2l.models.Cell2location.setup_anndata(st_sample)

print("Assemble and instantiate the cell2location model")
model = c2l.models.Cell2location(st_sample,
                                 cell_state_df=cell_state_df.loc[st_sample.var.index, :],
                                 N_cells_per_location=snakemake.params.ncell)
print("View the anndata setup")
model.view_anndata_setup()

print("Train the model.")
model.train(max_epochs=10_000, batch_size=None, train_size=1, use_gpu=has_gpu)

print("Save the loss plot")
fig, ax = plt.subplots()
model.plot_history(ax=ax)
fig.savefig(snakemake.output.loss_celltype_visium)

adata_st = model.export_posterior(st_sample,
                                  sample_kwargs={"num_samples": 2000,
                                                 "batch_size": model.adata.n_obs,
                                                 "use_gpu": has_gpu})
adata_st.obs[adata_st.uns["mod"]["factor_names"]] = adata_st.obsm["q05_cell_abundance_w_sf"]

adata_st.obs['patient'] = snakemake.wildcards.sample

# Visualize cell abundance mapped to the spatial coordinates:
slide = c2l.utils.select_slide(adata_st, batch_key='patient', s=snakemake.wildcards.sample)

with plt.rc_context({"figure.figsize": [4.5, 5]}):
    sc.pl.spatial(slide,
                  color=adata_st.uns["mod"]["factor_names"],
                  ncols=4,
                  size=1.3,
                  img_key="hires",
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0,
                  vmax="p99.2", show=False)
    plt.savefig(snakemake.output.cellTypeMap)


adata_st.write_h5ad(snakemake.output.post_c2l_anndata)