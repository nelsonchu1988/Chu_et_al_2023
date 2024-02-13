## Start

import scvelo
import os
from scipy.io import mmread
from scipy.io import mmwrite
from scipy.stats import median_abs_deviation
import loompy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scrublet as scr
import random
random.seed(0)
sc.settings.verbosity = 3   
import matplotlib
import matplotlib.colors as clr
cmap = clr.LinearSegmentedColormap.from_list('customcolors', ['#ffdee2','#5c3b82'], N=256)
sc.settings.set_figure_params(
    dpi=100, facecolor="white", fontsize=13, figsize=(8, 8), frameon=False
)

## Read files

adata = sc.read('lund_bone_marrow.h5ad')
adata.var_names_make_unique()

adata.X = adata.layers['spliced'] + adata.layers['unspliced']
adata.var.rename(index={adata.var[adata.var['gene_ids']=='ENSG00000076706'].index[0]: 'MCAM'}, inplace=True)
sc.pp.filter_cells(adata, min_genes = 100)
sc.pp.filter_genes(adata, min_cells = 10)
adata

## Process

exps = []
for sample in adata.obs["sample"].unique():
    tdata = adata[adata.obs["sample"] == sample, :].copy()
    print(sample)
    print(f"Total number of cells: {tdata.n_obs}")
    if tdata.n_obs < 100:
        exps.append(tdata)
        pass
    # find and remove doublets
    scrub = scr.Scrublet(tdata.X, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=15,
    )
    tdata = tdata[~predicted_doublets, :]
    exps.append(tdata)

adata = exps[0].concatenate(exps[1:])
adata.obs_names = adata.obs_names.str[0:-2]

print(adata.shape)
del (exps, tdata)

### QC

adata.var["mit"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mit"],
    percent_top=[20],
    log1p=True,
    inplace=True,
)

adata = adata[adata.obs['n_genes_by_counts']>800]
adata = adata[adata.obs['total_counts']<50000]

### Step1

sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
adata.layers["norm"] = adata.X.copy()
adata.raw=adata.copy()
sc.pp.log1p(adata)
adata.layers["normlog"] = adata.X.copy()

# Scoring
cell_cycle_genes = 'MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, MCM6, CDCA7, DTL, PRIM1, UHRF1, MLF1IP, HELLS, RFC2, RPA2, NASP, RAD51AP1, GMNN, WDR76, SLBP, CCNE2, UBR7, POLD3, MSH2, ATAD2, RAD51, RRM2, CDC45, CDC6, EXO1, TIPIN, DSCC1, BLM, CASP8AP2, USP1, CLSPN, POLA1, CHAF1B, BRIP1, E2F8, HMGB2, CDK1, NUSAP1, UBE2C, BIRC5, TPX2, TOP2A, NDC80, CKS2, NUF2, CKS1B, MKI67, TMPO, CENPF, TACC3, FAM64A, SMC4, CCNB2, CKAP2L, CKAP2, AURKB, BUB1, KIF11, ANP32E, TUBB4B, GTSE1, KIF20B, HJURP, CDCA3, HN1, CDC20, TTK, CDC25C, KIF2C, RANGAP1, NCAPD2, DLGAP5, CDCA2, CDCA8, ECT2, KIF23, HMMR, AURKA, PSRC1, ANLN, LBR, CKAP5, CENPE, CTCF, NEK2, G2E3, GAS2L3, CBX5, CENPA'.split(', ')
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
s_genes = list(set(s_genes)&set(adata[:,:].var.index))
g2m_genes = list(set(g2m_genes)&set(adata[:,:].var.index))
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
len(cell_cycle_genes)

nhvg = 2000
sc.pp.highly_variable_genes(adata, n_top_genes=nhvg, batch_key='sample')
adata = adata[:, adata.var.highly_variable]

# Scaling and PCA

sc.pp.scale(
    adata, zero_center=True
)  
sc.tl.pca(adata)
sc.pl.pca_overview(adata)
adata.layers["scaled"] = adata.X.copy()

# Embedding
knn, pc = 40, 20
batch_var = 'sample'

sce.pp.harmony_integrate(
    adata,
    batch_var,
    plot_convergence=True,
    max_iter_harmony=35,
)

sc.pp.neighbors(adata, n_neighbors=knn, n_pcs=pc, use_rep="X_pca_harmony")
sc.pl.scatter(adata, basis="pca_harmony", color=batch_var)

sc.tl.umap(adata) 
res = .4
sc.tl.leiden(adata, resolution = res)

ds=50

for i in ['leiden', 'CXCL12', 'LEPR', 'NGFR']:
    sc.pl.umap(adata, color=i, size=30, 
               legend_loc='on data',
               palette=sc.pl.palettes.default_102, 
               ncols=1, 
               layer='norm',
               title='',
               cmap = cmap,
               vmax='p99.8',
          )

## Heatmaps

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='all', use_raw=False, 
                        layer='normlog', 
                        pts=True
                       )
sc.tl.dendrogram(adata, groupby='leiden')

sc.tl.filter_rank_genes_groups(adata, key='all', 
                               min_in_group_fraction=0.5, max_out_group_fraction=0.45,
                               key_added='rg_filtered'
                              )

sc.pl.rank_genes_groups_matrixplot(adata, n_genes=3, use_raw=False, 
                                   cmap='bwr',
                                   key='rg_filtered',
                                   layer = 'normlog', 
                                   min_logfoldchange=2,
                                   dendrogram=True, 
                                  ) 

## Subset

sel_clusters = ['0',]
adata.obs['sel'] = adata.obs['leiden'].isin(sel_clusters)
adata.obs['sel'] = adata.obs['sel'].astype('category')

sel_umis = adata[adata.obs['leiden'].isin(sel_clusters),].obs_names
sel_umis = pd.Series(sel_umis)
sel_umis.to_csv('sel_umis.csv', index=None)

sc.pl.umap(adata, color=['leiden', 'CXCL12', 'LEPR', 'sel'], size=10, legend_loc='on data',
           palette=sc.pl.palettes.default_102, ncols=4, use_raw=False,
          ) 

## Step 2

adata = S.copy()
adata = adata[sel_umis,].copy()

# Remove samples that hardly contribute
selected_samples = adata.obs.pivot_table(index='sample', values='n_genes', aggfunc='count')>20
selected_samples = list(selected_samples[selected_samples['n_genes']].index)
adata = adata[adata.obs['sample'].isin(selected_samples),].copy()

sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
adata.layers["norm"] = adata.X.copy()
adata.raw=adata.copy()
sc.pp.log1p(adata)
adata.layers["normlog"] = adata.X.copy()

# Scoring
cell_cycle_genes = 'MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, MCM6, CDCA7, DTL, PRIM1, UHRF1, MLF1IP, HELLS, RFC2, RPA2, NASP, RAD51AP1, GMNN, WDR76, SLBP, CCNE2, UBR7, POLD3, MSH2, ATAD2, RAD51, RRM2, CDC45, CDC6, EXO1, TIPIN, DSCC1, BLM, CASP8AP2, USP1, CLSPN, POLA1, CHAF1B, BRIP1, E2F8, HMGB2, CDK1, NUSAP1, UBE2C, BIRC5, TPX2, TOP2A, NDC80, CKS2, NUF2, CKS1B, MKI67, TMPO, CENPF, TACC3, FAM64A, SMC4, CCNB2, CKAP2L, CKAP2, AURKB, BUB1, KIF11, ANP32E, TUBB4B, GTSE1, KIF20B, HJURP, CDCA3, HN1, CDC20, TTK, CDC25C, KIF2C, RANGAP1, NCAPD2, DLGAP5, CDCA2, CDCA8, ECT2, KIF23, HMMR, AURKA, PSRC1, ANLN, LBR, CKAP5, CENPE, CTCF, NEK2, G2E3, GAS2L3, CBX5, CENPA'.split(', ')
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
s_genes = list(set(s_genes)&set(adata[:,:].var.index))
g2m_genes = list(set(g2m_genes)&set(adata[:,:].var.index))
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

nhvg = 2000
sc.pp.highly_variable_genes(adata, n_top_genes=nhvg,batch_key='sample') 
adata = adata[:, adata.var.highly_variable]

# REGRESSION

regr_parameters = ['n_genes_by_counts']
sc.pp.regress_out(adata, regr_parameters, n_jobs=24)
adata.layers["regr"] = adata.X.copy()

# Scaling and PCA

sc.pp.scale(
    adata, zero_center=True
)  
sc.tl.pca(adata)
adata.layers["scaled"] = adata.X.copy()

# Embedding
knn, pc = 35, 15
batch_var = 'sample'

sce.pp.harmony_integrate(
    adata,
    batch_var,
    plot_convergence=True,
    max_iter_harmony=35,
)

sc.pp.neighbors(adata, n_neighbors=knn, n_pcs=pc, use_rep="X_pca_harmony")
sc.pl.scatter(adata, basis="pca_harmony", color=batch_var)
sc.tl.umap(adata)
res = .2
sc.tl.leiden(adata, resolution = res,random_state=0, )
sc.pl.umap(adata, color='leiden')

# Remove the traces of other cell types
adata = adata[adata.obs['leiden'].isin(['0','1','2'])].copy()
sc.pl.umap(adata, color='leiden')

clusters = 'leiden'
sc.pl.umap(adata, color=[clusters], 
           size=30, 
           vmax='p99.5', 
           cmap=cmap, 
           legend_loc='on data',
           legend_fontsize=20,

### Heatmaps

adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='mes', use_raw=False, 
                        layer='normlog', 
                        pts=True #corr_method='bonferroni', rankby_abs=True,tie_correct=True, 
                       )
sc.tl.dendrogram(adata, groupby='leiden')


sc.tl.filter_rank_genes_groups(adata, key='mes', 
                               # min_fold_change=2, 
                               min_in_group_fraction=0.5, max_out_group_fraction=0.45,
                               key_added='rg_filtered'
                              )

hcmap = clr.LinearSegmentedColormap.from_list('customcolors', ['lightgray','blue'], N=256)
sc.settings.set_figure_params(
    dpi=100, facecolor="white", fontsize=20, figsize=(4, 4), frameon=False
)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, 
                                key='mes',
                                use_raw=False, 
                                swap_axes=True, 
                                vmin=0,
                                vmax=3, 
                                cmap='bwr', 
                                layer='scaled', 
                                figsize=(10,7), 
                               )

### Genes specific for different clusters

gl = ['SPARCL1', 'HFM1', 'FLRT2', 'FOSB', 'SOCS3', 'JUND']

for g in gl:
    sc.pl.umap(
        adata,
        color=g,
        size=60,
        use_raw=True,
        vmax = 'p99.8',
        cmap=cmap,
    )
           

# UNITvelo

import unitvelo as utv
import scvelo as scv

sc.settings.set_figure_params(dpi=100, facecolor="white", figsize=(5, 5), frameon=False)
velo = utv.config.Configuration()
velo.R2_ADJUST = True
velo.IROOT = None
velo.FIT_OPTION = "1"
velo.GPU = 0

from numba import cuda 
device = cuda.get_current_device()
device

label = clusters
exp_metrics = {}
adata = utv.run_model(adata, label, config_file=velo)
scv.pl.velocity_embedding_stream(
    adata, color='leiden', dpi=100, title="UNIT VELO", frameon=False,
    linewidth=1.2,
    alpha=0.4,
)

scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', title='')