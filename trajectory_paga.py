import numpy as np
import scanpy as sc
import pandas as pd
import os

#Overall distribution of small molecule interference 

dge = pd.read_csv('../jiace/dge.csv',index_col=0)
cellinfo = pd.read_csv("./anno.csv",index_col=0)
cellinfo.index = cellinfo['cellname'] 
geneinfo = dge.index
gene = pd.DataFrame(geneinfo,columns='genename')
gene.index = geneinfo

adata = sc.AnnData(dge.values.T, obs=cellinfo, var = gene)
adata.obs['cluster']=cellinfo['cluster'].astype('str')

sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.draw_graph(adata)

sc.tl.paga(adata, groups='name')
sc.pl.paga(adata, color=['name'],edge_width_scale=0.5)
sc.settings.set_figure_params(dpi=400,dpi_save=200, facecolor='white')
sc.pl.paga(adata, color=['name'],edge_width_scale=0.1,frameon=False,fontsize=4)#,save='network.pdf')
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['cluster'], legend_loc='right margin',save='clusters1.pdf')
sc.pl.draw_graph(adata, color=['name'], legend_loc='right margin')#,save='paga.pdf')

#Trajectory analysis for several specific treatments
adata = sc.read('../jiace/pbmc_1.h5ad')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='name', legend_loc='right margin')
sc.tl.paga(adata, groups='name')
sc.pl.paga(adata, color=['name'],edge_width_scale=0.5)
sc.settings.set_figure_params(dpi=400,dpi_save=200, facecolor='white')
sc.pl.paga(adata, color=['name'],edge_width_scale=0.1,frameon=False,fontsize=4)#,save='network.pdf')
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['name'], legend_loc='right margin',save='pagaP0.pdf')
adata.uns['iroot'] = np.flatnonzero(adata.obs['name']  == 'P1')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['dpt_pseudotime'], legend_loc='right margin',save='pseudotime.pdf')
sc.pl.draw_graph(adata, color=['PRTG','PCAT14','MSX1','ZNF703','HOXA1','DSP','ARL4C','NR6A1','MT-CO3'], legend_loc='right margin',size=0.5,save='CH+RA+P1.pdf')