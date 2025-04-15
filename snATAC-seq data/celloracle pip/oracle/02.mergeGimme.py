import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
import pyarrow.parquet as pq
import mycelloracle

def load_tfidf(f) -> pd.DataFrame:
    import pyarrow.parquet as pq
    r = pq.read_table(f)
    return r.to_pandas()

tfscan_dir = '/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle'

group_file = '/cluster/share/atac_group/mafas5/chen_ws/SST/rerun.txt'
with open(group_file) as f:
    lines = f.readlines()
    groups = [l.strip() for l in lines if len(l.strip()) > 1]

tfidfs = [load_tfidf(f"{tfscan_dir}/{i}.pdc.baseGRN.df.parquet")
          for i in groups]
r = pd.concat(tfidfs, axis = 0, ignore_index= True)
r = r.fillna(int(0))
del tfidfs

prefix='SST-layer'
r.to_parquet(f"{tfscan_dir}/{prefix}.all.baseGRN.df.parquet")

tfCount = r.sum(axis = 0, numeric_only = True)
peakCount = r.sum(axis = 1, numeric_only = True)
