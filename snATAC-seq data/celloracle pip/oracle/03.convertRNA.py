import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
from anndata import AnnData
from scipy.io import mmread
from scipy.sparse import csr_matrix
import pyarrow.parquet as pq
import pyprojroot
import matplotlib
import matplotlib.pyplot as plt
import importlib
import mycelloracle

mtx_dir='/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/RNA'
X = mmread(os.path.join(mtx_dir,"/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/RNA/SST.logCPM.gbyc.mtx")).astype("float32")
X = csr_matrix(X)
X_cbyg = csr_matrix.transpose(X)
mm = AnnData(X_cbyg)
cell_ids = pd.read_csv(os.path.join(mtx_dir, "barcode.txt"), header = None).values[:,0]
gene_ids = pd.read_csv(os.path.join(mtx_dir, "gene.txt"), header = None).values[:, 0]
meta_info = pd.read_csv(os.path.join(mtx_dir, "meta.data.csv"),index_col=0)

mm = AnnData(X = X_cbyg,
             obs = meta_info.reindex(cell_ids),
             var = pd.DataFrame(index = gene_ids))

mm.write_h5ad(filename = "/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/RNA/SST-layer.logCPM.h5ad")

out_dir = '/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/RNA/layer'
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok = True)
subclasses = mm.obs['layer'].unique()
for sc in subclasses:
    print(f"Get subclass {sc}")
    adata_local = mm[mm.obs['layer'].isin([sc])]
    adata_local.write(os.path.join(out_dir, f"SST.layer.{sc}.ann.hdf5"))
