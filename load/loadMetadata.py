import re
import os
import numpy as np
import pandas as pd
import anndata as ad
from typing import Union,List



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
            adata.obs[col] = df_to_add[col]
    elif dim == "var":
        for col in df_to_add:
            adata.var[col] = df_to_add[col]
    else:
        raise ValueError("Adding metadata to given dimension is not supported")

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
