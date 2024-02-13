import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import os
from scvi.model.utils import mde
sc.set_figure_params(figsize=(12, 12))
import scrublet as scr
import matplotlib.colors as clr
cmap = clr.LinearSegmentedColormap.from_list('customcolors', ['#ffdee2','#5c3b82'], N=256)
from scipy.sparse import csr_matrix
from scipy.io import mmread

sc.settings.set_figure_params(
    dpi=100, facecolor="white", fontsize=13, figsize=(8, 8), frameon=False
)

# Read

atlas = {}

# Separated experiments were parsed to AnnData objects and stored within the atlas dict

# Concatenate

adata = ad.concat(atlas.values(), 
                  merge='same', 
                  join = 'outer', 
                  index_unique='_', 
                  keys = atlas.keys(), 
                  label = 'batch'
                 )
sc.pp.filter_cells(adata, min_genes = 15)
sc.pp.filter_genes(adata,min_cells = 15)

sc.pl.highest_expr_genes(adata)

# %% CALCULATE QC
adata.var["mit"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(
    adata,
    # layer='counts',
    qc_vars=["mit", "ribo",],
    percent_top=[20],
    log1p=True,
    inplace=True,
)

adata = adata[adata.obs['pct_counts_mit'] < 22]
adata = adata[adata.obs['total_counts'] < 55000]
adata = adata[adata.obs['total_counts'] > 500]

#%% DOUBLETS FOR MULTIPLE
exps = []
for sample in adata.obs["sample"].unique():
    tdata = adata[adata.obs["sample"] == sample, :].copy()
    print(sample)
    print(f"Total number of cells: {tdata.n_obs}")
    if (tdata.n_obs < 50): # don't look for doublets in case of small ds
        exps.append(tdata)
        print('PASS')
        print('_________')
        pass
    
    # find and remove doublets
    else: 
        scrub = scr.Scrublet(tdata.X, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=15,
        )
        tdata = tdata[~predicted_doublets, :]
        exps.append(tdata)


# concat
adata = ad.concat(exps, 
          merge='same', 
          join = 'outer', 
                 )
print(adata.shape)
del (exps, tdata)

# normalize
adata.raw = adata.copy()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers["norm"] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers["normlog"] = adata.X.copy()

# scoring
cell_cycle_genes = "MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, MCM6, CDCA7, DTL, PRIM1, UHRF1, MLF1IP, HELLS, RFC2, RPA2, NASP, RAD51AP1, GMNN, WDR76, SLBP, CCNE2, UBR7, POLD3, MSH2, ATAD2, RAD51, RRM2, CDC45, CDC6, EXO1, TIPIN, DSCC1, BLM, CASP8AP2, USP1, CLSPN, POLA1, CHAF1B, BRIP1, E2F8, HMGB2, CDK1, NUSAP1, UBE2C, BIRC5, TPX2, TOP2A, NDC80, CKS2, NUF2, CKS1B, MKI67, TMPO, CENPF, TACC3, FAM64A, SMC4, CCNB2, CKAP2L, CKAP2, AURKB, BUB1, KIF11, ANP32E, TUBB4B, GTSE1, KIF20B, HJURP, CDCA3, HN1, CDC20, TTK, CDC25C, KIF2C, RANGAP1, NCAPD2, DLGAP5, CDCA2, CDCA8, ECT2, KIF23, HMMR, AURKA, PSRC1, ANLN, LBR, CKAP5, CENPE, CTCF, NEK2, G2E3, GAS2L3, CBX5, CENPA".split(
    ", "
)
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

s_genes = list(set(s_genes) & set(adata[:, :].var.index))
g2m_genes = list(set(g2m_genes) & set(adata[:, :].var.index))

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# %% HVGENES
nhvg = 1330
adata.uns['log1p']['base'] = None

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=nhvg,
    batch_key="source",
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
)

# remove cell cycle from hivar
mito_genes = adata.var_names[adata.var_names.str.startswith("MT-")]
ribo_genes = adata.var_names.str.startswith(('RPL','RPS'))

for gl in (cell_cycle_genes, mito_genes, ribo_genes): 
    gl = list(set(gl)&set(adata.var_names))
    print(adata[:,gl].var['highly_variable'].sum())
    adata.var.loc[gl,'highly_variable'] = False

# %% SUBSET HVG
adata = adata[:, adata.var.highly_variable]
adata.shape

# %% SCALE AND PCA
sc.pp.scale(adata, zero_center=True, max_value=10)
sc.tl.pca(adata)
sc.pl.pca_overview(adata, color='source')

# HARMONY

import scanpy.external as sce

batch_var = "sample" 

sce.pp.harmony_integrate(
    adata,
    batch_var,
    plot_convergence=True,
    max_iter_harmony=35,
)

# NEIGHBORS

knn, pc = 30, 15
sc.pp.neighbors(adata, n_neighbors=knn, n_pcs=pc, use_rep="X_pca_harmony", method='umap')
sc.pl.scatter(adata, basis="pca_harmony", 
              color=batch_var, 
             )

# UMAP

sc.tl.umap(
    adata,
    min_dist=0.4,
)  

# CLUSTERS
res = .4
sc.tl.leiden(
    adata,
    resolution=res,
)
clusters = 'leiden'

sc.pl.umap(
    adata,
    color = ['leiden', 'CXCL12', 'LEPR',  'BGLAP', 'CDH5', 'MYH11', 'EBF3', 'source', ], 
    size=10,
    layer='norm',
    vmax='p99.8',
    cmap=cmap,
    ncols=3,
    legend_loc='on data',
)

# Heatmap
adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='mes', use_raw=False, 
                        layer='normlog', 
                        pts=True
                       )
sc.tl.dendrogram(adata, groupby='leiden')

sc.tl.filter_rank_genes_groups(adata, key='mes', 
                               min_in_group_fraction=0.5, max_out_group_fraction=0.45,
                               key_added='rg_filtered'
                              )

hcmap = clr.LinearSegmentedColormap.from_list('customcolors', ['lightgray','blue'], N=256)
sc.settings.set_figure_params(
    dpi=100, facecolor="white", fontsize=20, figsize=(4, 4), frameon=False
)

sc.pl.rank_genes_groups_matrixplot(adata, n_genes=3, use_raw=False, 
                                   cmap='bwr',
                                   key='mes',
                                   layer = 'normlog', 
                                   min_logfoldchange=4,
                                   dendrogram=True, 
                                   vmax=5,
                                  ) 


