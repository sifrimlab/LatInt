import os 
import re
import torch
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from torch import nn
import torch.nn as nn
from typing import Union,List



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
        np.savetxt('./latent_space.csv', latent, delimiter=',')

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

    latent = np.genfromtxt(path, delimiter=',')

    return latent




def addMetadataFromPandas(adata: ad.AnnData, df_to_add: Union[pd.DataFrame, str], dim:str = "obs"):
    """ Add every column of a pandas as metadata of an adata object, either on the obs or var dimension

    Parameters
    ----------
    adata : ad.AnnData
        adata to add metadata to
    df_to_add : Union[pd.DataFrame, str]
        df of which the columns will be viewed as metadata
    dim : str
        flag that decides whether to add to the obs or var metadata (either 'obs' or 'var', default = 'obs')
    """
    if isinstance(df_to_add, str):
        df_to_add = pd.read_csv("df_to_add")

    if dim == "obs":
        for col in df_to_add:
            adata.obs[col] = list(df_to_add[col])
    elif dim == "var":
        for col in df_to_add:
            adata.var[col] = list(df_to_add[col])
    else:
        raise ValueError("Adding metadata to given dimension is not supported")



def exportLatent(adata:ad.AnnData,latent_key:str) -> ad.AnnData:
    """Takes anndata and latent key as input and returns a new anndata object of the specified latent space with the metadata of the original object

    Parameters
    ----------
    adata : ad.AnnData
        adata
    latent_key : str
        latent_key

    Returns
    -------
    ad.AnnData

    """
    bdata = ad.AnnData(adata.obsm[latent_key])
    addMetadataFromPandas(bdata,adata.obs)
    return bdata

def addMetadataFromFileList(adata:ad.AnnData, file_list: List, delim:str = "_"):
    """ Add metadata extracted for the basenames of filepaths to AnnData object

    the metadata is expected to be delimeted by a specific character, and the metata itself to be consisting of a pure string, and its value a pure int

    Parameters
    ----------
    adata : ad.AnnData
        adata to add metadata to
    file_list : List
        list of files where the order and length is equal to the observations of the adata object.
    delim :  str
        delimiter of the metadata elements in the filepaths
    """

    if adata.shape[0] != len(file_list):
        raise ValueError("File list is not equally as long as the length of the observations")

    variable_dict = {}

    # first we initiate all keys with the first one
    first_file =  os.path.splitext(os.path.basename(file_list[0]))[0].split(delim)

    for el in first_file:
        try:
            string = re.findall('[a-zA-Z]+', el)[0]
            number =  re.findall('\d+',el)[0]
            variable_dict[string] = []
        except IndexError:
            pass

    for file in file_list:
        file_split = os.path.splitext(os.path.basename(file))[0].split(delim)
        for el in file_split:
            try:
                string = re.findall('[a-zA-Z]+', el)[0]
                number =  re.findall('\d+', el)[0]
                variable_dict[string].append(number)
            except IndexError:
                pass
    for key, value in variable_dict.items():
        adata.obs[key] = value
