import argparse
import sys
from collections import defaultdict
import numpy as np
import pandas as pd

def pandas_merge(tsv_files, outname):
	df_final = None
	for tsv in tsv_files:
		df = pd.read_csv(tsv, sep='\t')
		if df_final is None:
			df_final = df
		else:
			df_final = pd.merge(df_final, df, on='variant_id', how='outer')
	df_final.to_csv(outname, sep='\t', index=False, na_rep='nan')

def pandas_merge_fast(tsv_files, outname):
	dfs = (pd.read_csv(tsv, sep='\t').set_index("variant_id", drop=True) for tsv in tsv_files)
	df_final = pd.concat(dfs, axis=1, join='outer', copy=False)
	df_final.reset_index(drop=False, inplace=True)
	df_final.to_csv(outname, sep='\t', index=False, na_rep='nan')

parser = argparse.ArgumentParser(prog='merge_vcfs.py', description=__doc__)
parser.add_argument('tsvfiles', metavar='TSVFILES', nargs='+', help='List of TSV files to be combined into one table.')
parser.add_argument('outfile', metavar='OUTFILE', help='Name of the output TSV file.')
args = parser.parse_args()

tsv_files = args.tsvfiles
pandas_merge(tsv_files, args.outfile)
