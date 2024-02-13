## Import and initiate

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import os

from scvi.model.utils import mde
sc.settings.set_figure_params(
    dpi=100, facecolor="white", fontsize=13, figsize=(6, 6), frameon=False
)
import scrublet as scr
import matplotlib.colors as clr
cmap = clr.LinearSegmentedColormap.from_list('customcolors', ['#ffdee2','#5c3b82'], N=256)
cmap

import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread

## Read

adata = sc.read('nelson_main_umap.h5ad')
sc.pl.umap(adata, color='leiden')

gl = ['MCAM', 'NGFR']
for gene in gl:
    sc.pl.umap(adata, 
                color=gene, 
                frameon=False, 
                size=22,
                ncols=4,
                use_raw=True,
                # layer='norm',
                # legend_loc='on data',
                vmax='p99.8',
                cmap = cmap,
                # basis="X_mde",
               save='_full_%s.png'%(gene))

### Subset 

# keep the mes cluster
adata = adata[adata.obs['leiden']=='0']

# remove depleted experiments
samples_filtered = list((adata.obs['sample'].value_counts()[adata.obs['sample'].value_counts()>10]).index)
adata = adata[adata.obs['sample'].isin(samples_filtered)].copy()

### Stats

# count positive cells
cd_stats = pd.DataFrame(index=adata.obs_names,)
gl = ['NGFR', 'CXCL12', 'MCAM']

th = 1
for i in cd_stats.columns:
    print(i)
    print( str((cd_stats[i]>=th).sum()) + ' cells' )
    print(str(round(( ((cd_stats[i]>=th) & (cd_stats['CXCL12']>5)).sum() / cd_stats.shape[0]), 2)) + '%')
    print('-----')

for i in gl:
    sc.pl.umap(adata, 
               color=i,
               vmin=0, vmax=2, 
               size=10, 
               cmap=cmap)

### Table with counts

# count positive cells

cd_stats = pd.DataFrame(index=adata.obs_names,)
gl = ['NGFR', 'CXCL12', 'MCAM']

for gene in gl:
    cd_stats[gene] = adata.raw[:,gene].X.todense()
cd_stats = cd_stats.astype(int)

### NGFR

cxcl_threshold = 5
th = 1
gene = 'NGFR'

cd_stats['ngfr_dp'] = 'na'
cd_stats.loc[((cd_stats['CXCL12']>cxcl_threshold) & (cd_stats[gene]<th)), 'ngfr_dp'] = 'sp' # single positive
cd_stats.loc[((cd_stats['CXCL12']>cxcl_threshold) & (cd_stats[gene]>=th)), 'ngfr_dp'] = 'dp' # double positive
cd_stats.loc[((cd_stats['CXCL12']<=cxcl_threshold) & (cd_stats[gene]<th)), 'ngfr_dp'] = 'dn' # double negative

adata.obs['ngfr_dp'] = cd_stats['ngfr_dp'].copy() # add labels to adata

# rankgenes
genekey = 'ngfr'
adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby='ngfr_dp', groups=['dp'], reference='sp', 
                        method='wilcoxon', 
                        layer='normlog', 
                        use_raw=False, 
                        key_added=genekey
                       )

# Plot genes

sc.pl.rank_genes_groups(adata, 
                        key=genekey,
                        n_genes=20, 
                        fontsize =20, 
                       )

#### Volcano

# CREATE PVALS TABLE

keyz = list(adata.uns[genekey].keys())[1::]

df = pd.DataFrame()
for key in keyz:
    df[key] = [x[0] for x in adata.uns[genekey].get(key)]

df.index = df['names']
df = df[df['pvals_adj']<0.05]
df = df.sort_values('logfoldchanges', ascending=False)

#%% VOLCANO
# log2fold threshold
logft_ur = df['logfoldchanges'].sort_values(ascending=False)[20] # top 20 upreg to annotate
logft_dr = df['logfoldchanges'].sort_values(ascending=True)[20] # top 20 downreg to annotate
logft_dr = -0.5
logft = 0.5
y = 'pvals_adj'
x = 'logfoldchanges'

ax = df.plot(x=x, y=y, 
             kind='scatter',logx=False, logy=True,
             title = '%s'%(gene),
             )

# draw borders
plt.axhline(y = 0.05, color = 'r', linestyle = '-')
plt.axvline(x = logft, color = 'r', linestyle = '-')
plt.axvline(x = -logft, color = 'r', linestyle = '-')
plt.gca().invert_yaxis()
ax.grid(False)

# add annotations
for idx, row in df.iterrows():
    if (((row[x] < logft_dr) or (row[x] > logft_ur)) and (row[y] < 0.05)):
        ax.annotate(row['names'], 
                    (row[x], row[y]), 
                    xytext=((-15 + (np.sign(row[x])*20)),-5),
                    textcoords='offset points', 
                    family='sans-serif', 
                    fontsize=10
                   )

#### Pathways

import gseapy

res = df.copy()
res['Symbol'] = res.index.copy()

pval = 0.05
lfc = 0.5

# take up- and downregulated
upreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) > lfc)]['Symbol']
print(upreg.shape)
downreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) < -lfc)]['Symbol']
print(downreg.shape)

gene_set_names = pd.Series(gseapy.get_library_name(organism='Human'))
gene_set = ['KEGG_2021_Human']

gene_set_names

enr_res_up = gseapy.enrichr(gene_list=upreg,
                     organism='Human',
                     gene_sets=gene_set,
                        )
# store results
df_res_ur = pd.DataFrame(enr_res_up.results).set_index('Term')
# plot bars
gseapy.barplot(enr_res_up.res2d,title='Upreg %s'%(gene), top_term=20)

### MCAM

gene = 'MCAM'
gene_ph = 'mcam_dp'
genekey = gene


cd_stats[gene_ph] = 'na'
cd_stats.loc[((cd_stats['CXCL12']>5) & (cd_stats[gene]<th)), gene_ph] = 'sp'
cd_stats.loc[((cd_stats['CXCL12']>5) & (cd_stats[gene]>=th)), gene_ph] = 'dp'
cd_stats.loc[((cd_stats['CXCL12']<=cxcl_threshold) & (cd_stats[gene]<th)), gene_ph] = 'dn'
adata.obs[gene_ph] = cd_stats[gene_ph].copy()

adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby=gene_ph, groups=['dp'], reference='sp', method='wilcoxon', layer='normlog', use_raw=False, key_added=genekey)

# SQUARES

sc.pl.rank_genes_groups(adata, 
                        key=genekey,
                        n_genes=20, 
                        fontsize =20, 
                       )

### Volcano

#%% CREATE PVALS TABLE
keyz = list(adata.uns[genekey].keys())[1::]

df = pd.DataFrame()
for key in keyz:
    df[key] = [x[0] for x in adata.uns[genekey].get(key)]

df.index = df['names']
df = df[df['pvals_adj']<0.05]
df = df.sort_values('logfoldchanges', ascending=False)


#%% VOLCANO

# remove outliers
df = df[(df['logfoldchanges']<30) & (df['pvals_adj'] > 1e-100)]

# log2fold threshold
logft_ur = df['logfoldchanges'].sort_values(ascending=False)[20] # upreg
logft_dr = df['logfoldchanges'].sort_values(ascending=True)[20] #downreg
logft_dr = -0.5
logft = .5
y = 'pvals_adj'
x = 'logfoldchanges'


ax = df.plot(x=x, y=y, 
             kind='scatter',logx=False, logy=True,
             title = '%s'%(gene),
             )

# draw lines
plt.axhline(y = 0.05, color = 'r', linestyle = '-')
plt.axvline(x = logft, color = 'r', linestyle = '-')
plt.axvline(x = -logft, color = 'r', linestyle = '-')
plt.gca().invert_yaxis()
ax.grid(False)

# add annotations
for idx, row in df.iterrows():
    if (((row[x] < logft_dr) or (row[x] > logft_ur)) and (row[y] < 0.05)):
        ax.annotate(row['names'], (row[x], row[y]), 
                    xytext=((-15 + (np.sign(row[x])*20)),-5),
                    textcoords='offset points', 
                    family='sans-serif', fontsize=10)

### Pathways

res = df.copy()
res['Symbol'] = res.index.copy()

pval = 0.05
lfc = 0.5

upreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) > lfc)]['Symbol']
print(upreg.shape)
downreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) < -lfc)]['Symbol']
print(downreg.shape)

import gseapy
gene_set_names = pd.Series(gseapy.get_library_name(organism='Human'))
gene_set = ['KEGG_2021_Human']

enr_res_up = gseapy.enrichr(gene_list=upreg,
                     organism='Human',
                     gene_sets=gene_set,
                        )
# store results
df_res_ur = pd.DataFrame(enr_res_up.results).set_index('Term')
# plot bars
gseapy.barplot(enr_res_up.res2d,title='Upreg %s'%(gene), top_term=20)

## Three genes

cd_stats = pd.DataFrame(index=adata[adata.obs['leiden']=='0',].obs_names,)
for gene in ['NGFR', 'CXCL12', 'MCAM']:
    cd_stats[gene] = adata[adata.obs['leiden']=='0',].raw[:,gene].X.todense()
cd_stats = cd_stats.astype(int)
cd_stats

gene_ph = 'ppn'
genekey = 'ppn'
cxcl_threshold = 5
th = 1

# HIGHLIGHT DOUBLE POSITIVES
cd_stats[gene_ph] = 'na' 
cd_stats.loc[((cd_stats['CXCL12']>cxcl_threshold) & (cd_stats['MCAM']>=th) & (cd_stats['NGFR']<th)), gene_ph] = 'mcam_ppn' # MCAM+ NGFR- CXCL12+
cd_stats.loc[((cd_stats['CXCL12']>cxcl_threshold) & (cd_stats['NGFR']>=th) & (cd_stats['MCAM']<th)), gene_ph] = 'ngfr_ppn' # NGFR MCAM- CXCL12-
cd_stats.loc[((cd_stats['CXCL12']>cxcl_threshold) & (cd_stats['NGFR']>=th) & (cd_stats['MCAM']>=th)), gene_ph] = 'triple' # triple positive
adata.obs[gene_ph] = cd_stats[gene_ph].copy()


# RANK GENES
adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby=gene_ph, groups=['mcam_ppn'], reference='ngfr_ppn', method='wilcoxon', layer='normlog', use_raw=False, key_added=genekey)

cd_stats['ppn'].value_counts(normalize=False)

cd_stats['ppn'].value_counts(normalize=True)

# SQUARES

sc.pl.rank_genes_groups(adata, 
                        key=genekey,
                        n_genes=20, 
                        fontsize =20, 
                       )

### Volcano

#%% CREATE PVALS TABLE
# read the rest and join with expressed
keyz = list(adata.uns[genekey].keys())[1::]

df = pd.DataFrame()
for key in keyz:
    df[key] = [x[0] for x in adata.uns[genekey].get(key)]

df.index = df['names']
df = df[df['pvals_adj']<0.05]
df = df.sort_values('logfoldchanges', ascending=False)

#%% VOLCANO

# remove outliers
df = df[(df['logfoldchanges']<30) & (df['pvals_adj'] > 1e-100)]

# log2fold thr
logft_ur = df['logfoldchanges'].sort_values(ascending=False)[20] # upreg
print(logft_ur)
logft_dr = df['logfoldchanges'].sort_values(ascending=True)[20] #downreg
logft_dr = -0.5
logft = .5
y = 'pvals_adj'
x = 'logfoldchanges'


plt = df.plot(x=x, y=y, 
             kind='scatter',logx=False, logy=True,
             )
plt.axhline(y = 0.05, color = 'r', linestyle = '-')
plt.axvline(x = logft, color = 'r', linestyle = '-')
plt.axvline(x = -logft, color = 'r', linestyle = '-')
plt.invert_yaxis()
plt.grid(False)

for idx, row in df.iterrows():
    if (((row[x] < logft_dr) or (row[x] > logft_ur)) and (row[y] < 0.05)):
        plt.annotate(row['names'], (row[x], row[y]), 
                    xytext=((-15 + (np.sign(row[x])*20)),-5),
                    textcoords='offset points', 
                    family='sans-serif', fontsize=10)

### Pathways

import gseapy

res = df.copy()
res['Symbol'] = res.index.copy()
res.head()

pval = 0.05
lfc = 0.5

upreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) > lfc)]['Symbol']
downreg = res[(res.pvals_adj < pval) & ((res.logfoldchanges) < -lfc)]['Symbol']

all_sign = res[(res.pvals_adj < pval) & (abs(res.logfoldchanges) > lfc)]['Symbol']
gene_set_names = pd.Series(gseapy.get_library_name(organism='Human'))
gene_set = ['KEGG_2021_Human']

enr_res_up = gseapy.enrichr(gene_list=upreg,
                     organism='Human',
                     gene_sets=gene_set,
                        )
# store results
df_res_ur = pd.DataFrame(enr_res_up.results).set_index('Term')
# plot bars
gseapy.barplot(enr_res_up.res2d, top_term=20)

gene_set_names = pd.Series(gseapy.get_library_name(organism='Human'))
gene_set = ['KEGG_2021_Human']

enr_res_down = gseapy.enrichr(gene_list=downreg,
                     organism='Human',
                     gene_sets=gene_set,
                        )
# store results

# plot bars
gseapy.barplot(enr_res_down.res2d,title='downreg dp vs dp %s'%(gene), top_term=20, cutoff=0.1)

### Venn



# bool df with thresholds

df = cd_stats
df['NGFR'] = df['NGFR'].astype(bool)
df['MCAM'] = df['MCAM'].astype(bool)
df['CXCL12'] = df['CXCL12']>5

venn = []
for gene in df.columns:
    a = df.columns[df.columns!=gene][0]
    b = df.columns[df.columns!=gene][1]
    out = (df[gene] & ~df[a] & ~df[b]).sum()
    print(gene, out)
    venn.append(out)
    
    out = (df[gene] & df[a] & ~df[b]).sum()
    print(gene, a, out)
    venn.append(out)
    
    out = (df[gene] & ~df[a] & df[b]).sum()
    print(gene, b, out)
    venn.append(out)
    print()
    
out = (df[gene] & df[a] & df[b]).sum()
print(gene, b, a, out)
venn.append(out)

from matplotlib_venn import venn2  
from matplotlib import pyplot as plt 
  
# depict venn diagram 
venn2(subsets = (4956, 833, 456), set_labels = ('NGFR', 'MCAM')) 
plt.show()

# scVI embedding and clustering

samples_filtered = list((adata.obs['sample'].value_counts()[adata.obs['sample'].value_counts()>10]).index)
adata = adata[adata.obs['sample'].isin(samples_filtered)].copy()

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")

vae = scvi.model.SCVI(adata, n_layers=2, n_latent=15, dispersion='gene-batch',)
import torch
torch.set_float32_matmul_precision('medium')

vae.train(use_gpu=True)

adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"], device='cuda')

# remove outliers
adata = adata[adata.obsm['X_mde'][:,0] >-2]
adata = adata[adata.obsm['X_mde'][:,1] < 4]

# cluster by X_scVI
sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_scVI")
sc.tl.umap(adata,)
sc.tl.leiden(adata, key_added="leiden_scvi", resolution=.4)
clusters = 'leiden_scvi'

from matplotlib.colors import ListedColormap
col_dict = {0: "pink",
            1: "red",
            2: "orange",
            7: "green"}

# Custom color map
cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

clusters = 'leiden_scvi'
sc.pl.embedding(adata, 
                color=[clusters], 
                frameon=False, 
                size=50,
                ncols=4,
                vmax='p99.8', 
                color_map=cm, 
                basis="X_mde", 
               )