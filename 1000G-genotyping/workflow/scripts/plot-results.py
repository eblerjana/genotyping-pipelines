import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse


def plot_concordances_all(files, outname, sources, variants):
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080']
	type_to_file = {}
	n_rows = 4
	n_cols = 5
	plt.figure(figsize=(19,20))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	plot_index = 1
	for var in variants:
		plt.subplot(n_rows, n_cols, plot_index)
		plot_index += 1
		all_samples = []
		x_values = []
		is_first = True
		for i,source in enumerate(sources):
			print('source', source)
			samples = []
			concordances = []
			if not (source, var) in type_to_file:
				continue
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				samples.append(fields[0])
				concordances.append(float(fields[1]))
			if is_first:
				all_samples = samples
			else:
				assert all_samples == samples
			is_first = False
			x_values = [i*5 for i in range(len(samples))]
			plt.plot(x_values, concordances, label=source, color=colors[i], marker='o')
		plt.title(var_to_name[var])
		plt.xticks(x_values, all_samples, rotation='vertical')
		plt.ylabel('weighted genotype concordance [%]')
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.tight_layout()
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)


def plot_fscores_all(files, outname, sources, variants):	
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080']
	type_to_file = {}
	n_rows = 4
	n_cols = 5
	plt.figure(figsize=(19,20))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	plot_index = 1
	for var in variants:
		plt.subplot(n_rows, n_cols, plot_index)
		plot_index += 1
		all_samples = []
		is_first = True
		for i,source in enumerate(sources):
			print('source', source)
			samples = []
			fscores = []
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				samples.append(fields[0])
				fscores.append(float(fields[3]))
			if is_first:
				all_samples = samples
			else:
				assert all_samples == samples
			is_first = False
			x_values = [i*5 for i in range(len(samples))]
			plt.plot(x_values, fscores, label=source, color=colors[i], marker='o')
		plt.title(var_to_name[var])
		plt.xticks(x_values, all_samples, rotation='vertical')
		plt.ylabel('adjusted F-score [%]')
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.tight_layout()
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)

parser = argparse.ArgumentParser(prog='plot-results.py', description="Plot concordances or precision/recall statistics.")
parser.add_argument('-files', metavar='FILES', nargs='+', help='files with results per sample.')
parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output file.')
parser.add_argument('-sources',  metavar='SOURCES', nargs='+', help='regions concordance_all')
parser.add_argument('-metric', metavar='METRIC', required=True, choices=['concordance', 'precision-recall-typable'], help='Which metric to use for plotting.')
parser.add_argument('-variants', metavar='VARIANTS', nargs='+', default=[], choices=['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex'])
args = parser.parse_args()

variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
if args.variants != []:
	variants = [v for v in args.variants]

if args.metric == 'concordance':
	plot_concordances_all(args.files, args.outname, args.sources, variants)
else:
	plot_fscores_all(args.files, args.outname, args.sources, variants)
