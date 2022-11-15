from typing import Union, Tuple
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import numpy as np



def plotLatentCorrelation(latent: Union[ad.AnnData, np.ndarray, str]) -> Tuple:
    if isinstance(latent, ad.AnnData):
        adata = latent
        latent = latent.X
    elif isinstance(latent, str):
        if latent.endswith("h5ad"):
            adata = ad.read_h5ad(latent)
            latent = adata.X
        elif latent.endswith("csv"):
            latent = np.genfromtxt(latent, delimiter=",")

    plt.matshow(np.corrcoef(latent, rowvar=False), cmap="coolwarm")
    plt.colorbar()
    plt.show()


        





































