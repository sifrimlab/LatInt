import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Union, Tuple
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype

def plotLatentCorrelation(latent: Union[ad.AnnData, np.ndarray, str], latent_key:str, meta_var = "") -> Tuple:
    if isinstance(latent, ad.AnnData):
        adata = latent
        latent = adata.obsm[latent_key]
        ticks = (range(len(adata.var_names)), adata.var_names)
        if not meta_var:
            corcoef = np.corrcoef(latent, rowvar=False)
        else:
            corcoef = np.corrcoef(latent, y=adata, rowvar=False)
            

    elif isinstance(latent, str):
        if latent.endswith("h5ad"):
            adata = ad.read_h5ad(latent)
            latent = adata.obsm[latent_key]
            corcoef = np.corrcoef(latent, rowvar=False)
        elif latent.endswith("csv"):
            latent = np.genfromtxt(latent, delimiter=",")
            corcoef = np.corrcoef(latent, rowvar=False)

    plt.matshow(corcoef, cmap="coolwarm")
    plt.colorbar()
    if "ticks" in locals():
        plt.xticks(*ticks, rotation=45)
        plt.yticks(*ticks)
    plt.show()

def plotLatentCorrelationWithMetaVariable(adata: ad.AnnData, latent_key:str, meta_var:str):

    # Sourced from: https://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays/30143754#30143754
    def corr2_coeff(A, B):
        # Rowwise mean of input arrays & subtract from input arrays themeselves
        A_mA = A - A.mean(1)[:, None]
        B_mB = B - B.mean(1)[:, None]

        # Sum of squares across rows
        ssA = (A_mA**2).sum(1)
        ssB = (B_mB**2).sum(1)

        # Finally get corr coeff
        return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

    latent = adata.obsm[latent_key]

    tmp_df = pd.DataFrame(adata.obs[meta_var])

    if not is_numeric_dtype(tmp_df[meta_var]):
        unique_clusters =sorted(tmp_df[meta_var].unique())

        unique_clusters_numbers = np.arange(1,len(unique_clusters)+1)
        cluster_mapping = {} # key = cluster string, value = cluster number
        for cluster, number in zip(unique_clusters, unique_clusters_numbers):
            cluster_mapping[cluster] = number
            tmp_df['meta_label'] = tmp_df[meta_var].map(cluster_mapping)
    else:
        # rename column to streamline code below
        tmp_df['meta_label'] = tmp_df[meta_var]


    tmp_array = np.array(tmp_df["meta_label"])
    tmp_array = np.expand_dims(tmp_array, 1)

    plt.bar(range(0,latent.shape[1]), corr2_coeff(latent.T,tmp_array.T).ravel(), color="#9467bd")
    plt.xticks(range(0,latent.shape[1]), range(0,latent.shape[1]))
    plt.xlabel("Latent dimensions")
    plt.ylabel("Correlation")
