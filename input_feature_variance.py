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

def feature_variance_per_latent_dimension(adata, key, n_features=20, plot_highly_variable_features=True, 
                                          plot_cov_matrices=False, return_variable_features=True):
    
    latent = adata.obsm["z_mvae"]
    original = adata.X

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
            ax.set_ylabel('score')

    if sharey is True:
        ymax += 0.3 * (ymax - ymin)
        ax.set_ylim(ymin, ymax)

#     writekey = f"rank_genes_groups_{adata.uns[key]['params']['groupby']}"
#     savefig_or_show(writekey, show=show, save=save)
    plt.show()