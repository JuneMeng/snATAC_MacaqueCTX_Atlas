import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from gimmemotifs.motif import MotifConfig, read_motifs
from genomepy import Genome
import scanpy as sc
import mycelloracle


parser = argparse.ArgumentParser(description = "run gimmemotifs")
parser.add_argument('--pdcbedpe', type = str)
parser.add_argument('--outdir', type = str, default = "tfscan")
parser.add_argument('--threshold', type = int, default = 10)
parser.add_argument('--refgenome', type = str, default = "macFas5")

args = parser.parse_args()

# * main
# to celloracle supported file
print(f"TFscan for {args.pdcbedpe}")
pdc_bedpe = pd.read_table(args.pdcbedpe, sep = "\t", header = None)
pdc_bedpe.columns = ["pchr", "pstart", "pend", "dchr",
                     "dstart", "dend", "pair", "pearson", "n1", "n2"]

peak_id = pdc_bedpe['dchr']  + '_' + pdc_bedpe['dstart'].map(str) + '_'  + pdc_bedpe['dend'].map(str)
gene_short_name = pdc_bedpe['pair'].apply(lambda x: x.split('|')[0])
pdc = pd.concat([peak_id, gene_short_name], axis = 1)
pdc.columns = ["peak_id", 'gene_short_name']
mycelloracle.check_peak_format(peaks_df = pdc, ref_genome = args.refgenome)

# * tfcan
tfi = ma.TFinfo(peak_data_frame = pdc, ref_genome = args.refgenome)

# load motif
config = MotifConfig()
motif_dir = config.get_motif_dir()
path = os.path.join(motif_dir, "gimme.vertebrate.v5.0.pfm")
motifs = read_motifs(path)

tfi.scan(fpr = 0.02, motifs = motifs, verbose = True)
# save result
prefix = os.path.basename(args.pdcbedpe).replace(".bedpe", "")
tfi.to_hdf5(file_path=f"{args.outdir}/{prefix}.celloracle.tfinfo")

# * filter
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=args.threshold)
tfi.make_TFinfo_dataframe_and_dictionary(verbose = True)
df_tfi = tfi.to_dataframe()
# save result
df_tfi.to_parquet(f"{args.outdir}/{prefix}.baseGRN.df.parquet")
print("Done")
