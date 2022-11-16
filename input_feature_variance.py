import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from typing import Optional, Union, Callable, Sequence
import anndata as ad

def feature_variance_per_latent_dimension(adata: ad.AnnData, 
                                          key: str, 
                                          n_features: int =  20, 
                                          plot_highly_variable_features: Optional[bool] = True, 
                                          plot_cov_matrices: Optional[bool] = False, 
                                          return_variable_features: Optional[bool] = True):
    
    """Calculate the n features that vary the most with each latent dimensions. Stores the results in adata.uns['variable_features_per_latent_dim']
    
    Parameters
    ----------
    adata : ad.AnnData
        annotated data matrix
    key : str
        key in adata.obsm[key] where the latent space matrix is stored (observations x features)
    n_features : int
        number of features to return. These are the features with the highest variances along the latent dimensions
    plot_highly_variable_features : bool
        whether to plot the count of highest scoring features per latent dimensions
    plot_cov_matrices : bool
        whether to plot the covariance matrices (features x latent_dims)
    return_variable_features : bool
        whether to return the (scaled) covariances and highest scoring genes or just store it in adata.uns['variable_features_per_latent_dim']
        
    """
    
    latent = adata.obsm["z_mvae"]
    original = adata.X # Get the normalized counts instead!

    cov = np.cov(original, latent, rowvar=False)
    sub_cov = cov[:original.shape[1], original.shape[1]:]
    scaled_covs = MinMaxScaler().fit_transform(sub_cov)
    
    if plot_cov_matrices:
        sns.heatmap(scaled_covs)
        sns.clustermap(scaled_covs)
        
    if plot_highly_variable_features:
        d = {}
        for dim in range(latent.shape[1]):
            feature_list = []
            for i, feature in enumerate(adata.var_names):
                if scaled_covs[i, dim] > 0.85:
                    feature_list.append(feature)

            d[f"LatentDim_{dim+1}"] = feature_list
            
        plt.figure(figsize = (20,8))
        ax = sns.barplot(x=list(d.keys()), y=[len(d[key]) for key in d.keys()])
        plt.ylabel("Count")
        for item in ax.get_xticklabels():
            item.set_rotation(90)
                
    variable_features_per_dimension = {}
    for dim in range(latent.shape[1]):
        recarray = np.recarray((n_features,), dtype=[("cov", float), ("scaled_cov", float), ("genes", object)])
        # Get 20 genes and their (scaled) variances
        highest_index = np.argsort(scaled_covs[:,dim])[::-1]
#         variable_features_per_dimension[f"LatentDim_{dim+1}"] = {
#             "cov": sub_cov[:,dim][highest_index][:n_features],
#             "scaled_cov": scaled_covs[:,dim][highest_index][:n_features],
#             "genes": np.array(adata.var_names[highest_index][:n_features], dtype=object)
#         }
        recarray["cov"] = sub_cov[:,dim][highest_index][:n_features]
        recarray["scaled_cov"] = scaled_covs[:,dim][highest_index][:n_features]
        recarray["genes"] = np.array(adata.var_names[highest_index][:n_features], dtype=object)
        variable_features_per_dimension[f"LatentDim_{dim+1}"] = recarray
        
    adata.uns["variable_features_per_latent_dim"] = variable_features_per_dimension
    
    if return_variable_features:
        return variable_features_per_dimension
    

def _savefig_or_show(
    writekey: str,
    show: Optional[bool] = None,
    dpi: Optional[int] = None,
    ext: str = None,
    save: Union[bool, str, None] = None,
):
    
    """
    Not actually used so ignore
    """
    if isinstance(save, str):
        # check whether `save` contains a figure extension
        if ext is None:
            for try_ext in ['.svg', '.pdf', '.png']:
                if save.endswith(try_ext):
                    ext = try_ext[1:]
                    save = save.replace(try_ext, '')
                    break
        # append it
        writekey += save
        save = True
    save = settings.autosave if save is None else save
    show = settings.autoshow if show is None else show
    if save:
        savefig(writekey, dpi=dpi, ext=ext)
    if show:
        pl.show()
    if save:
        pl.close()  # clear figure

def plot_feature_variance(
    adata: ad.AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: int = 20,
    gene_symbols: Optional[str] = None,
    key: Optional[str] = 'variable_features_per_latent_dim',
    fontsize: int = 8,
    ncols: int = 4,
    sharey: bool = True,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    ax: Optional[plt.Axes] = None,
    **kwds,
):
    
    """\
    Plot the n highest variable genes along the latent dim. Inspired by: https://github.com/scverse/scanpy/blob/d7e13025b931ad4afd03b4344ef5ff4a46f78b2b/scanpy/plotting/_tools/__init__.py
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix.
    groups : Union[str, Sequence[str]]
        The groups for which to show the gene ranking.
    n_genes : int
        Number of genes to show.
    fontsize : int
        Fontsize for gene names.
    ncols: int
        Number of panels shown per row.
    sharey : bool
        Controls if the y-axis of each panels should be shared. By passing
        `sharey=False`, each panel has its own y-axis range.
    """
    
    if 'n_panels_per_row' in kwds:
        n_panels_per_row = kwds['n_panels_per_row']
    else:
        n_panels_per_row = ncols
    if n_genes < 1:
        raise NotImplementedError(
            "Specifying a negative number for n_genes has not been implemented for "
            f"this plot. Received n_genes={n_genes}."
        )
        
    group_names = list(adata.uns[key].keys()) if groups is None else groups
    # one panel for each group
    # set up the figure
    n_panels_x = min(n_panels_per_row, len(group_names))
    n_panels_y = np.ceil(len(group_names) / n_panels_x).astype(int)

    from matplotlib import gridspec

    fig = plt.figure(
        figsize=(
            n_panels_x * plt.rcParams['figure.figsize'][0],
            n_panels_y * plt.rcParams['figure.figsize'][1],
        )
    )
    gs = gridspec.GridSpec(nrows=n_panels_y, ncols=n_panels_x, wspace=0.22, hspace=0.3)

    ax0 = None
    ymin = np.Inf
    ymax = -np.Inf
    
    for count, group_name in enumerate(group_names):
        gene_names = adata.uns[key][group_name]["genes"][:n_genes]
        scores = adata.uns[key][group_name]['cov'][:n_genes]

        # Setting up axis, calculating y bounds
        if sharey:
            ymin = min(ymin, np.min(scores))
            ymax = max(ymax, np.max(scores))

            if ax0 is None:
                ax = fig.add_subplot(gs[count])
                ax0 = ax
            else:
                ax = fig.add_subplot(gs[count], sharey=ax0)
        else:
            ymin = np.min(scores)
            ymax = np.max(scores)
            ymax += 0.3 * (ymax - ymin)

            ax = fig.add_subplot(gs[count])
            ax.set_ylim(ymin, ymax)

        ax.set_xlim(-0.9, n_genes - 0.1)

        # Mapping to gene_symbols
#         if gene_symbols is not None:
#             if adata.raw is not None and adata.uns[key]['params']['use_raw']:
#                 gene_names = adata.raw.var[gene_symbols][gene_names]
#             else:
#                 gene_names = adata.var[gene_symbols][gene_names]

        # Making labels
        for ig, gene_name in enumerate(gene_names):
            ax.text(
                ig,
                scores[ig],
                gene_name,
                rotation='vertical',
                verticalalignment='bottom',
                horizontalalignment='center',
                fontsize=fontsize,
            )

        ax.set_title('{}'.format(group_name))
        if count >= n_panels_x * (n_panels_y - 1):
            ax.set_xlabel('ranking')

        # print the 'score' label only on the first panel per row.
        if count % n_panels_x == 0:
            ax.set_ylabel('covariance')

    if sharey is True:
        ymax += 0.3 * (ymax - ymin)
        ax.set_ylim(ymin, ymax)

#     writekey = f"rank_genes_groups_{adata.uns[key]['params']['groupby']}"
#     savefig_or_show(writekey, show=show, save=save)
    plt.show()