import scanpy as sc
import pandas as pd
import seaborn as sns

def plot_latent_distributions(adata:sc.AnnData, latent_key:str, subset:list=None) -> sc.AnnData:
    """ Plot latent distributions over the samples. Requires seaborn=0.12.0

    Parameters
    ----------
    adata : ad.AnnData
        Input scanpy object
    latent_key: str
        latent key stored in .obsm
    subset: list
        Subset of latent dimension to plot
    figsize: tuple
        Figure size

    Returns
    -------
    ad.Anndata

    """

    df = pd.DataFrame(ad.obsm[latent_key], index=ad.obs.index)
    if subset is not None:
        df = df.loc[:, subset]

    df = df.stack()
    df.name = 'val'
    df = df.reset_index()
   
    sns.displot(df,
        x='val',
        hue='level_1',
        kind = 'kde',
    )
