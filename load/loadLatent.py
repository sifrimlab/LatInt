import os 
import numpy as np
import torch
from typing import Union
import scanpy as sc
import anndata as ad
from torch import nn
import torch.nn as nn
from utils import ImageDataset


def getLatentFromModel(model: Union[nn.Module, str], dataSet: torch.utils.data.Dataset, device = torch.device("cuda", 0), batch_size:int = 256, save:bool =False) -> np.ndarray:
    """Calculates latent space based on a trained pytorch model and its training dataset.

    Parameters
    ----------
    model : Union[nn.Module, str]
        Trained model, MUST implement an embed() function
    dataSet : torch.utils.data.Dataset
        Instantiation of the torch.Dataset class with the training data inside
    device :
        device to send model to for latent space calculation
    batch_size : int
        batch size, larger means faster calculation, but more memory usage
    save : bool
        Whether to save the latent space array to a file
        

    Returns
    -------
    np.ndarray where rows are observations and columns are their respective latent space values

    """

    if not isinstance(model, nn.Module):
        if isinstance(model, str):
            model = torch.load(model)
        else:
            raise TypeError("Input model is not of type str of nn.Module")

    if not isinstance(dataSet, torch.utils.data.Dataset):
        raise TypeError("Torch DataSet needed to load data and calculate latent space")

    if "embed" not in dir(model):
        raise KeyError("model does not implement embed function, latent space cannot be calculated")

    train_loader = dataLoader(dataSet, batch_size=batch_size, shuffle=False)

    with torch.no_grad():
        model.eval()
        model.to(device)
        
        latent_list = []
        for batch in tqdm(train_loader):
            latent, _, _ = model.embed(batch.to(device))
            latent_list.append(latent.to("cpu"))
        
        latent = torch.concat(latent_list, dim=0)
        latent = latent.numpy()
    if save:
        np.savetxt('./latent_space.csv', latent, sep=',')

    return latent

def getLatentFromAnnData(adata: Union[ad.AnnData, str]) -> np.ndarray:
    """ Gets latent space arrray from an anndata file

    Parameters
    ----------
    adata : Union[ad.AnnData, str]
        annData object or path to annData object

    Returns
    -------
    np.ndarray where rows are observations and columns are their respective latent space values

    """
    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)
    return adata.X

def getLatentFromCsv(path: str) -> np.ndarray: 
    """ Gets latent space from a csv file

    Parameters
    ----------
    path : str
        Path to csv file, where rows have to be observations and columns their respective latent space values

    Returns
    -------
    np.ndarray where rows are observations and columns are their respective latent space values

    """

    latent = np.genfromtxt(path, delimeter=',')

    return latent
