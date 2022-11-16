import torch 
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable

def plotTensor(tensor):
    plt.imshow(torch.squeeze(tensor.cpu().detach(), dim=0).permute(1,2,0))

def interactLatentSpaceImages(model:nn.Module, n_latent_dims:int ,dim:int ,value:float , device:str, plotting_function:Callable = plotTensor):
    """interactLatentSpaceImages.

    Parameters
    ----------
    model : nn.Module
        Model to interpret, must implement a decode() function
    n_latent_dims : int
        Number of latent dimensions in the model
    dim : int
        dimension to change
    value : float
        Value to change the dimension to
    device : str
        Device to compute on
    plotting_function : Callable
        Function that plot the decoded tensor. This will depend on the output shape of your decoder. Default is for our own usecase, where the output shape = [1, 1, 100, 100].
        In short, it should: remove empty dimension to make the tensor 2D and move the tensor to cpu



    Usage:

    ```
    from ipywidgets import interact, fixed 

    interact(
        interactLatentSpaceImages,
        model=fixed(model),
        device=fixed(device),
        n_latent_dims=fixed(n_latent_dims),
        dim = (0,n_latent_dims - 1,1),
        value=(-2,2,0.1)
    );
    ```

    """
    empty_nda = np.zeros((1,n_latent_dims))
    empty_nda[0,dim] = value
    changed_tensor = torch.from_numpy(empty_nda).float()

    decoded_tensor = model.decode(changed_tensor.to(device))

    plotTensor(decoded_tensor)

