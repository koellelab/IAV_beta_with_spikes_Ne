import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob

# read in list of all data that were used in my analysis
all_dates = pd.read_csv('data/all_dates.tsv', sep='\t', keep_default_na=False).\
	assign(ID = lambda k: k.ID.astype(str)).\
	query('data == "mccrone"')

# read in sra run table and subset to what we want
sra = pd.read_csv('data/mccroneIAV_analysis_files/mccrone_SraRunTable.csv')\
		[['Run', 'Isolate', 'SPECID', 'Collection_date', 'Library Name']].\
	assign(Isolate = lambda k: k.Isolate.astype(str)).\
	merge(all_dates[['ID']].drop_duplicates().rename(columns={'ID': 'Isolate'}),
		how='inner',
		on='Isolate')

# get which are technical replicates from run table
dupes = sra[['Isolate', 'SPECID', 'Run']].drop_duplicates().\
	groupby('SPECID').size().reset_index(name='size').\
	query('size > 1').SPECID.drop_duplicates().values


if sra.Isolate.drop_duplicates().shape[0] != all_dates.ID.drop_duplicates().shape[0]:
	raise Exception('sra dimensions do not match all dates dimensions')

# read in raw mccrone data,
# remove variants where ref == alt
# subset to just the ones we want
jm = pd.read_csv('data/mccroneIAV_analysis_files/no_freq_cut.qual.snv.csv').\
	query('ref != var')\
	[['ENROLLID', 'SPECID', 'collect', 'home_collected', 'run', 'mutation',  'pos', 'ref', 'var', 'freq.var']].\
	assign(
		ENROLLID = lambda k: k['ENROLLID'].astype(str),
		segment = lambda k: [i.split('_')[0] for i in k.mutation],
		library = lambda k: k.run.str.replace('_', '').str.replace('victoria', 'vic'),
		date = lambda k: k['collect'],
		clinic = lambda k: k.home_collected == 0,
		duplicate = lambda k: np.isin(k.SPECID, dupes)).\
	assign(
		segment = lambda k: np.where(k.segment == 'M', 'MP', np.where(k.segment == "NR", "NA", k.segment))).\
	drop(['run'], axis=1).\
	merge(all_dates[['ID', 'date', 'clinic']].drop_duplicates().rename(columns={'ID': 'ENROLLID'}),
		how='inner',
		on=['ENROLLID', 'date', 'clinic'])


if jm.ENROLLID.drop_duplicates().shape[0] != all_dates.ID.drop_duplicates().shape[0]:
	raise Exception('sra dimensions do not match all dates dimensions')


if jm[['ENROLLID', 'date', 'clinic', 'duplicate']].drop_duplicates().duplicate.sum() != dupes.shape[0]:
	raise Exception('number of duplicates do not match')


# read in mccrone duplicate data and subset to just ones we want
# remove variants where ref == alt
jm_dupes = pd.read_csv('data/mccroneIAV_analysis_files/duplicate_sequences.csv').\
	query('ref != var').\
	merge(jm[['SPECID']].drop_duplicates(),
		left_on='SPECID_original',
		right_on='SPECID',
		how='inner')

# demonstrate that for those sequenced in duplicate
# not all variants in no_freq_cut.qual.snv.csv are in duplicate_sequences file
# thus, this file is useless to us
print(jm[np.isin(jm.SPECID, dupes)][['SPECID', 'mutation', 'freq.var']].\
	assign(jm = True).\
	merge(jm_dupes[['SPECID', 'mutation', 'freq1', 'freq2', 'SPECID1', 'SPECID2']].drop_duplicates().\
		assign(dupes = True),
		on=['SPECID', 'mutation'],
		how='outer').\
	query('dupes.isnull()'))

# consequently, we get the "called" variant frequency for each site
# ignoring which technical replicate was used
jm_arr = jm[['ENROLLID', 'date', 'clinic', 'duplicate', 'segment', 'pos', 'ref', 'var', 'freq.var']].\
	rename(columns={'ref': 'ref_nt', 'ENROLLID': 'ID'}).\
	pivot(index=['ID', 'date', 'clinic', 'duplicate', 'segment', 'pos', 'ref_nt'],
		columns='var', values='freq.var').\
	reset_index().fillna(0).\
	assign(
		A = lambda k: np.where(k.ref_nt == 'A', 1 - k[['C', 'G', 'T']].sum(axis=1), k.A),
		C = lambda k: np.where(k.ref_nt == 'C', 1 - k[['A', 'G', 'T']].sum(axis=1), k.C),
		G = lambda k: np.where(k.ref_nt == 'G', 1 - k[['A', 'C', 'T']].sum(axis=1), k.G),
		T = lambda k: np.where(k.ref_nt == 'T', 1 - k[['A', 'C', 'G']].sum(axis=1), k['T']),
		data = 'mccrone',
		ref = lambda k: 'mccrone_HK_' + k.segment + '_nan_nan')

# merge in universal positions file
pos_map = pd.read_csv('data/genomes/all_aln_map.tsv', sep='\t').\
			query('(data == "mccrone") & date.isnull()').\
			assign(ref = lambda k: 'mccrone_' + k.ID + '_nan_nan')

jm_arr = jm_arr.merge(pos_map[['ref', 'pos', 'upos']],
	how='left',
	on=['ref', 'pos']).\
	assign(upos = lambda k: k.upos.astype(int))

if jm_arr.upos.isnull().sum() > 0:
	raise Exception('null upos')

jm_arr.drop(['ref'], axis=1).to_csv('data/mccroneIAV_analysis_files/no_freq_cut.qual.snv_arr.csv', index=None)

# save list of included samples, ref, and reference path
all_dates[jm.segment.drop_duplicates().values]  = True
all_dates.drop(['acc'], axis=1).melt(id_vars = ['data', 'ID', 'date', 'duplicate', 'd', 'clinic'],
	var_name = 'segment').\
	drop(['value'], axis=1).\
	assign(
		ref = lambda k: 'mccrone_HK_' + k.segment + '_nan_nan',
		ref_path = lambda k: 'data/mccroneIAV_analysis_files/Hong_Kong_PC1A/HongKongPC1A_' + k.segment + '.fasta').\
	to_csv('data/mccroneIAV_analysis_files/no_freq_cut.qual.snv_arr_refs.csv', index=None)
