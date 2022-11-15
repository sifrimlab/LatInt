import pandas as pd
import scanpy as sc 
import anndata as ad
from typing import Tuple

def subsample(adata:ad.AnnData, cluster_key:str, target_cells:float=10000) -> ad.AnnData:
    """subsample.

    Parameters
    ----------
    adata : ad.AnnData
        adata
    cluster_key : str
        cluster_key
    target_cells : float
        target_cells

    Returns
    -------
    ad.AnnData

    """
    if not cluster_key:
        sc.pp.subsample(adata, fraction=0.2, n_obs=None, random_state=0, copy=False)
    
    adatas = [adata[adata.obs[cluster_key].isin([clust])] for clust in adata.obs[cluster_key].cat.categories]
    for dat in adatas:
        if dat.n_obs > target_cells:
            sc.pp.subsample(dat, n_obs=target_cells, random_state=0)
    adata_downsampled = adatas[0].concatenate(*adatas[1:])
    return adata_downsampled

def clusterLatentspace(adata:ad.AnnData,cluster_key=None,show_plot:bool=False):
    """clusterLatentspace.

    Parameters
    ----------
    adata : ad.AnnData
        adata
    cluster_key :
        cluster_key
        If no cluster key is given leiden clustering will be defaulted to.
    show_plot : bool
        show_plot
    """
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    if not cluster_key:
        sc.tl.leiden(adata,random_state=101,resolution=0.7)
        cluster_key = 'leiden'
    if show_plot:
        sc.pl.umap(adata,color=cluster_key,size=12)
        
def plotVariableLatent(adata:ad.AnnData,cluster_key='leiden',method='wilcoxon',num_features=5) -> Tuple:
    """plotVariableLatent.

    Parameters
    ----------
    adata : ad.AnnData
        adata
    cluster_key :
        cluster_key
    method :
        Different methods for ranking genes for characterizing groups. 't-test', 'wilcoxon' or 'logreg'.
    num_features :
        Number of highest variable features to return in dataframe

    Returns
    -------
    DotPlot of the different clusters
    DataFrame showing cluster_key and its most variable num_features 
    """
    P1 = sc.pl.dotplot(adata, var_names=adata.var_names, groupby=cluster_key)
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, key_added=method)
    var_features = pd.DataFrame(adata.uns[method]['names']).head(num_features)
    return P1,var_features 
