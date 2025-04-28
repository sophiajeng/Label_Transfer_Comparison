import os
os.chdir('../')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import scipy
from scipy import sparse, io
import pandas as pd


sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
trvae_epochs = 500
surgery_epochs = 500

early_stopping_kwargs = {
    "early_stopping_metric": "val_unweighted_loss",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
adata_ref=sc.read("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/refmatrix.txt")
adata_query=sc.read("/home/exacloud/gscratch/mcweeney_lab/jengs/scarches/refmatrix.txt")
#sparsematrix = io.mmread('sparsematrix.txt')
#m_dense = sparsematrix.toarray()
ref_dense=adata_ref.toarray()
var_names = np.genfromtxt('rownames.txt', dtype=str)
col_names = np.genfromtxt('colnames.txt', dtype=str)
