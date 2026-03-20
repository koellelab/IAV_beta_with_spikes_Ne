import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import glob
import ast
import argparse
try:
	from utils import map_arr, read_fasta, import_seqs
except:
	from scripts.utils import map_arr, read_fasta, import_seqs


def get_minor_allele(k):
	# minor allele from first file
	minor = np.argmax(np.where(
		(k[['A_x', 'C_x', 'G_x', 'T_x']] == k[['A_x', 'C_x', 'G_x', 'T_x']].max(axis=1).values[:,None]),
		-np.inf,
		k[['A_x', 'C_x', 'G_x', 'T_x']]), axis=1)
	# minor allele from second file
	minor_y = np.argmax(
		np.ma.masked_array(k[['A_y', 'C_y', 'G_y', 'T_y']].values, 
				mask=(np.arange(4) == np.argmax(k[['A_x', 'C_x', 'G_x', 'T_x']], axis=1)[:, None]) | 
					(k[['A_y', 'C_y', 'G_y', 'T_y']].values == 0)).\
			filled(-np.inf),
		axis=1)
	# consolidate
	# replace any -inf or 0 values with nan
	minor = np.where(
		k[['A_x', 'C_x', 'G_x', 'T_x']].values[np.arange(k.shape[0]),minor] == 0.0,
		np.where(minor_y <= 0, np.nan, minor_y),
		minor)
	return(minor)


def expand_seq_arr(paths):
	# ambiguity codes in McCrone references throws things off
	# R = A or G, so calling A, to revisit
	return(pd.concat([pd.DataFrame(list(import_seqs(open(i, 'r'))[1][0].upper().replace('R', 'A')), columns=['nuc']).\
			reset_index(names=['pos']).\
			assign(
				pos = lambda k: k.pos + 1,
				segment = i.split('/')[-1].split('.')[0].split("_")[1]) for
		i in paths]).\
		assign(freq = 1).\
		pivot(index=['segment', 'pos'],
			columns='nuc',
			values='freq').\
		reset_index())


def merge_var_ref_arr(var, ref):
	return(pd.concat([
			var,
			ref[['segment', 'pos']].\
				merge(var[['segment', 'pos']],
					on=['segment', 'pos'],
					how='left', indicator=True).\
				query('_merge == "left_only"').\
				drop('_merge', axis=1).\
				merge(ref,
					on=['segment', 'pos'],
					how='left')]).\
		sort_values(by=['segment', 'pos']))


def run():
	parser = argparse.ArgumentParser()
	# input files
	# processed vcf arrs
	parser.add_argument('--varFreq1', default=None)
	parser.add_argument('--varFreq2', default=None)
	parser.add_argument('--useIDs')
	parser.add_argument('--labels', nargs=3)
	parser.add_argument('--posMap')
	parser.add_argument('--minPos', type=int, default=200)
	args = parser.parse_args()
	#args.varFreq1 = 'data/processed_vcf_arr.csv'
	#args.varFreq2 = 'data/mccroneIAV_analysis_files/no_freq_cut.qual.snv_arr.csv'
	#args.posMap = 'data/genomes/all_aln_map.tsv'
	#args.useIDs = 'data/all_dates_mccrone.tsv'
	'''
	args.labels = [
		'Human IAV iSNV frequencies called by McCrone et al. (2018)',
		'Human IAV iSNV frequencies called by our pipeline',
		'jt_mm_compare']
	'''
	
	use_ids = pd.read_csv(args.useIDs, sep='\t',  keep_default_na=False).\
		assign(ID = lambda k: k.ID.astype(str),
			clinic = lambda k: k.clinic.astype(float))

	pos_map = pd.read_csv('data/genomes/all_aln_map.tsv', sep='\t').\
		assign(segment = lambda k: np.where(k.seg.isnull(), 'NA', k.seg),
			ref = lambda k: 
				k['data'].astype(str) + 
					'_' + k['ID'].astype(str) + 
					'_' + k['date'].astype(str) + 
					'_' + k['acc'].astype(str))\
		[['ref', 'segment', 'pos', 'upos']].\
		assign(segment_len = lambda k: k.groupby(['ref', 'segment']).pos.transform('max'))

	
	var_freq_x = pd.read_csv(args.varFreq1,  keep_default_na=False).\
		assign(
			freq1 = True,
			ID = lambda k: k.ID.astype(str),
			clinic = lambda k: k.clinic.astype(float))

	var_freq_y = pd.read_csv(args.varFreq2,  keep_default_na=False).\
		assign(
			freq2 = True,
			ID = lambda k: k.ID.astype(str),
			clinic = lambda k: k.clinic.astype(float))

	# filter based on use_ids 
	var_freq_x = var_freq_x.merge(use_ids[['data', 'ID', 'date', 'clinic']],
		on=['data', 'ID', 'date', 'clinic'],
		how='inner')
	var_freq_y = var_freq_y.merge(use_ids[['data', 'ID', 'date', 'clinic']],
		on=['data', 'ID', 'date', 'clinic'],
		how='inner')

	# add reference columns to use_ids
	use_ids_full = use_ids.\
		merge(
			pd.read_csv(args.varFreq1.replace('.csv', '_refs.csv')).\
				assign(ID = lambda k: k.ID.astype(str),
					segment = lambda k: k.segment.fillna('NA'))\
				[['data', 'ID', 'date', 'clinic', 'segment', 'ref_path', 'ref']].\
				drop_duplicates().\
				rename(columns={'ref_path': 'ref_path_x', 'ref': 'ref_x'}),
			how='left').\
		merge(
			pd.read_csv(args.varFreq2.replace('.csv', '_refs.csv')).\
				assign(ID = lambda k: k.ID.astype(str),
					segment = lambda k: k.segment.fillna('NA'))\
				[['data', 'ID', 'date', 'clinic', 'segment', 'ref_path', 'ref']].\
				drop_duplicates().\
				rename(columns={'ref_path': 'ref_path_y', 'ref': 'ref_y'}),
			how='left')

	if use_ids_full.shape[0]/8 != use_ids.shape[0]:
		raise Exception('full use IDs wrong shape')


	# merge
	# remerge in reference information
	# re-add coordinates from pos map
	merged_var_freq = \
		var_freq_x[['data', 'ID', 'date', 'clinic', 'duplicate', 'segment', 'upos', 'A', 'C', 'G', 'T']].\
			merge(
				var_freq_y[['data', 'ID', 'date', 'clinic', 'duplicate', 'segment', 'upos', 'A', 'C', 'G', 'T']],
				on=['data', 'ID', 'date', 'clinic', 'duplicate', 'segment', 'upos'],
				how='outer').\
			merge(use_ids_full[['data', 'ID', 'date', 'clinic', 'duplicate', 'segment', 'ref_path_x', 'ref_x', 'ref_path_y', 'ref_y']],
				on=['data', 'ID', 'date', 'clinic', 'duplicate', 'segment'],
				how='left').\
			merge(pos_map.drop_duplicates().rename(columns={'pos': 'pos_x', 'ref': 'ref_x', 'segment_len': 'segment_len_x'}),
				how='left',
				on=['ref_x', 'segment', 'upos']).\
			merge(pos_map.drop_duplicates().rename(columns={'pos': 'pos_y', 'ref': 'ref_y', 'segment_len': 'segment_len_y'}),
				how='left',
				on=['ref_y', 'segment', 'upos'])


	if merged_var_freq[['ref_path_x', 'ref_path_y']].isnull().any().any() | \
			merged_var_freq[['pos_x', 'pos_y']].isnull().any().any() :
		raise Exception('missing reference paths')


	filled_var_freq = []
	for idx in ['x', 'y']:
		for g, g_dat in merged_var_freq[merged_var_freq['A_' + idx].isnull()].groupby('ref_path_' + idx):
			#if (g_dat.ID.iloc[0] == "50993") & (g_dat.segment.iloc[0] == "NA") & (g_dat.upos.iloc[0] == 934):
			#	print('here')
			#	break
			ref_arr = expand_seq_arr([g_dat['ref_path_'+idx].iloc[0]]).\
							rename(columns={
								'pos': 'pos_'+idx, 
								'A': 'A_'+idx, 
								'C': 'C_'+idx, 
								'G': 'G_'+idx, 
								'T': 'T_'+idx}).\
							fillna(0)
			filled_var_freq.append(
				g_dat.drop(['A_'+idx, 'C_'+idx, 'G_'+idx, 'T_'+idx], axis=1).\
					merge(
						ref_arr,
						how='left',
						on=['segment', 'pos_'+idx]))
			if filled_var_freq[-1].shape[0] != g_dat.shape[0]:
				raise Exception('merging changed shape')
			

	# concatenate
	filled_var_freq = pd.concat(filled_var_freq + \
		[merged_var_freq[~merged_var_freq[['A_x', 'A_y']].isnull().any(axis=1)]])

	if filled_var_freq.shape[0] != merged_var_freq.shape[0]:
		raise Exception('filled var freq wrong shape')

	# add summary columns
	filled_var_freq = filled_var_freq.assign(
		major_allele_x = lambda k:
			np.array(['A', 'C', 'G', 'T'])\
				[k[['A_x', 'C_x', 'G_x', 'T_x']].values.argmax(axis=1)],
		minor_allele_x = lambda k: np.where(
			(k[['A_x', 'C_x', 'G_x', 'T_x']].values ==1).any(axis=1),
			'X',
			np.array(['A', 'C', 'G', 'T'])[np.argsort(k[['A_x', 'C_x', 'G_x', 'T_x']].values, axis=1)[:,-2]]),
		major_allele_y = lambda k:
			np.array(['A', 'C', 'G', 'T'])\
				[k[['A_y', 'C_y', 'G_y', 'T_y']].values.argmax(axis=1)],
		minor_allele_y = lambda k: np.where(
			(k[['A_y', 'C_y', 'G_y', 'T_y']].values ==1).any(axis=1),
			'X',
			np.array(['A', 'C', 'G', 'T'])[np.argsort(k[['A_y', 'C_y', 'G_y', 'T_y']].values, axis=1)[:,-2]]))

	
	# add depth
	# read in all depth files from my analysis
	# DP always merges in on x
	dp_path = "data/depth/*/*.tsv"
	dp_files = glob.glob(dp_path)
	# flag only do split once
	dp = pd.concat([pd.read_csv(i, sep='\t', header=None).\
				assign(i = i)\
			for i in dp_files]).\
		rename(columns={0: 'segment', 1: 'pos', 2:'depth'}).\
		assign(
			type = lambda k: [i.split('_')[0] for i in k.segment],
			segment = lambda k: [i.split('_')[1] for i in k.segment],
			data = lambda k: [i.split('/')[-2].split('_')[0] for i in k.i],
			ID = lambda k: [i.split('/')[-2].split('_')[1] for i in k.i],
			date = lambda k: [i.split('/')[-2].split('_')[2] for i in k.i],
			acc = lambda k: [i.split('/')[-2].split('_')[3] for i in k.i]).drop(['i'], axis=1)

	filled_var_freq = filled_var_freq.\
		merge(
			dp.merge(use_ids[['data', 'ID', 'acc', 'clinic']].drop_duplicates(),
				how='inner',
				on=['data','ID','acc']).\
				rename(columns={'pos': 'pos_x'}).\
				drop(['acc', 'type', 'data'], axis=1),
			how='left',
			on=['ID', 'date', 'clinic', 'segment', 'pos_x'])

	# filter on variant position
	filtered_var_freq = filled_var_freq.\
		query('(pos_x > @args.minPos) & ((segment_len_x - pos_x) > @args.minPos)').\
		query('(pos_y > @args.minPos) & ((segment_len_y - pos_y) > @args.minPos)')

	if filtered_var_freq[['A_x', 'A_y']].isnull().any().any():
		raise Exeption('null freqs')

	
	import matplotlib.pyplot as plt
	compare_snps = filtered_var_freq[
		(filtered_var_freq[['A_x', 'C_x', 'G_x', 'T_x']].values != \
			filtered_var_freq[['A_y', 'C_y', 'G_y', 'T_y']].values).any(axis=1)]

	compare_snps = compare_snps.\
		assign(minor_allele_x_freq_x = lambda k: 
				np.column_stack(
					[k[['A_x', 'C_x', 'G_x', 'T_x']].values,
						np.repeat(0, k.shape[0])])\
				[np.arange(k.shape[0]),
					map_arr(
						np.where(
							k.minor_allele_x.values == 'X',
							k.minor_allele_y.values,
							k.minor_allele_x.values), 
						{'A': 0, 'C': 1, 'G': 2, 'T': 3, 'X': 4})]).\
		assign(minor_allele_x_freq_y = lambda k: 
				np.column_stack(
					[k[['A_y', 'C_y', 'G_y', 'T_y']].values,
						np.repeat(0, k.shape[0])])\
				[np.arange(k.shape[0]),
					map_arr(
						np.where(
							k.minor_allele_x.values == 'X',
							k.minor_allele_y.values,
							k.minor_allele_x.values), 
						{'A': 0, 'C': 1, 'G': 2, 'T': 3, 'X': 4})])

	compare_snps.to_csv(f'data/{args.labels[2]}_snps.csv', index=None)
	compare_snps.query('depth >= 500').to_csv(f'data/{args.labels[2]}_snps.csv', index=None)
	
	os.makedirs('figures', exist_ok=True)
	fig, axs = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
	tmp = compare_snps
	axs.scatter(tmp['minor_allele_x_freq_y'], tmp['minor_allele_x_freq_x'],
		facecolor='none',
		edgecolor='#333333',
		alpha=0.5,
		zorder=2)

	axs.set_xlabel(args.labels[0])
	axs.set_ylabel(args.labels[1])
	axs.set_xlim(-0.025, 1.025)
	axs.set_ylim(-0.025, 1.025)
	axs.grid(axis='both', color='#eaeaea', zorder=0)
	axs.set_aspect('equal')
	fig.savefig(f'figures/{args.labels[2]}.pdf')
	plt.close()

	fig, axs = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
	tmp = compare_snps.query('(depth >= 500)')
	axs.scatter(tmp['minor_allele_x_freq_y'], tmp['minor_allele_x_freq_x'],
		facecolor='none',
		edgecolor='#333333',
		alpha=0.5,
		zorder=2)

	axs.set_xlabel(args.labels[0])
	axs.set_ylabel(args.labels[1])
	axs.set_xlim(-0.025, 1.025)
	axs.set_ylim(-0.025, 1.025)
	axs.grid(axis='both', color='#eaeaea', zorder=0)
	axs.set_aspect('equal')
	axs.set_title('>= 500 reads')
	fig.savefig(f'figures/{args.labels[2]}_500.pdf')
	plt.close()


	fig, axs = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
	tmp = compare_snps.query('(depth >= 500)')
	axs.scatter(tmp['minor_allele_x_freq_y'], tmp['minor_allele_x_freq_x'],
		facecolor='none',
		edgecolor=np.array(['#333333', 'indianred'])[1*tmp.duplicate.values],
		alpha=0.5,
		zorder=2)

	axs.set_xlabel(args.labels[0])
	axs.set_ylabel(args.labels[1])
	axs.set_xlim(-0.025, 1.025)
	axs.set_ylim(-0.025, 1.025)
	axs.grid(axis='both', color='#eaeaea', zorder=0)
	axs.set_aspect('equal')
	axs.set_title('>= 500 reads')
	fig.savefig(f'figures/{args.labels[2]}_500_labelled.pdf')
	plt.close()



if __name__ == "__main__":
    run()
