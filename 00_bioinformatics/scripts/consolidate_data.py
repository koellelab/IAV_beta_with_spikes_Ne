import glob
import pandas as pd
import numpy as np


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def import_seqs(fh):
    s_names = []
    s_list = []
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            s_list.append(s)
    fh.close()    
    return(np.array(s_names), s_list)


def process_vcf(vcf_path):
	with open(vcf_path, 'r') as fp:
		for line in fp:
			if line[0:2] != '##':
				header=line.split('\n')[0].split('\t')
				break
	try:
		vcf = \
			pd.read_csv(vcf_path, sep='\t', header=None, comment='#')
		vcf.columns = header
		# todo make this better
		info_split = np.array([i.split(';') for i in vcf.INFO])
		info_cols = {i.split("=")[0]: idx for idx,i in enumerate(info_split[0])}
		info_split = np.array([i.split('=')[1] for i in info_split.\
				flatten()]).\
			reshape((info_split.shape[0], -1))
		dp4_split = np.array([i.split(',') for i in info_split[:,info_cols['DP4']]]).astype(int).cumsum(axis=1)[:,[1,-1]]
		dp4_split[:,1] -= dp4_split[:,0]
		info_cols[info_split.shape[1]] = 'DP4_ref'
		info_cols[info_split.shape[1]+1] = 'DP4_alt'
		vcf[['AF', 'SB']] = \
			info_split[:,[info_cols['AF'], info_cols['SB']]].astype(float)
		vcf[['DP4_ref', 'DP4_alt']] = \
			dp4_split
		vcf['DP4'] = dp4_split.sum(axis=1)
		return(vcf.rename(
			columns={
				'#CHROM': 'segment',
				'REF': 'ref_nt',
				'ALT': 'iSNV_nt',
				'POS': 'iSNV_loc'}))
	except:
		return(pd.DataFrame(columns=
			['segment', 'iSNV_loc', 'ID', 'ref_nt', 'iSNV_nt', 
				'QUAL', 'FILTER', 'INFO', 'AF', 'DP4_ref', 
				'DP4_alt', 'DP4', 'SB']))


def process_vcf_af(vcf_path):
	with open(vcf_path, 'r') as fp:
		for line in fp:
			if line[0:2] != '##':
				header=line.split('\n')[0].split('\t')
				break
	try:
		vcf = \
			pd.read_csv(vcf_path, sep='\t', header=None, comment='#')
		vcf.columns = header
		# todo make this better
		info_split = np.array([i.split(';') for i in vcf.INFO])
		info_cols = {i.split("=")[0]: idx for idx,i in enumerate(info_split[0])}
		vcf['AF'] = np.array([i[info_cols['AF']].replace('AF=', '') for i in info_split]).astype(float)
		return(vcf.rename(
			columns={
				'#CHROM': 'segment',
				'REF': 'ref_nt',
				'ALT': 'iSNV_nt',
				'POS': 'iSNV_loc'})\
			[['segment', 'iSNV_loc', 'ID', 'ref_nt', 'iSNV_nt',
					'QUAL', 'AF']])
	except:
		return(pd.DataFrame(columns=
			['segment', 'iSNV_loc', 'ID', 'ref_nt', 'iSNV_nt', 
				'QUAL', 'FILTER', 'AF']))


# read in vcf files:
def get_subtype(x):
	x_split = x.split('_')
	if len(x_split) > 2:
		return(x_split[2])
	else:
		return('')


# need to read in McCrone metadata and get list of IDs
# for samples taken at clinic on the same day as household sample
# will add 0.5 days to the date for this sample
mccrone_runtable = pd.read_csv('data/mccroneIAV_analysis_files/mccrone_SraRunTable.csv')
homochronous_clinic_accs = mccrone_runtable[['Run', 'Collection_date', 'Isolate', 'SPECID']].\
	drop_duplicates().\
	assign(
		clinic = lambda k: 1 * (k.SPECID.str[0:2] == "MH"))\
	[['Run', 'clinic']].\
	rename(columns={'Run': 'acc'})


dp_path = "data/depth/*/*.tsv"
vcf_path = "data/variants_noBAQ/*/*.vcf"
vcf_raw_path = "data/variants_raw/*/*.vcf"

 
dp = pd.concat([pd.read_csv(i, sep='\t', header=None).\
			assign(i = i)\
		for i in glob.glob(dp_path)]).\
	rename(columns={0: 'segment', 1: 'pos', 2:'depth'}).\
	assign(
		type = lambda k: [i.split('_')[0] for i in k.segment],
		segment = lambda k: [i.split('_')[1] for i in k.segment],
		data = lambda k: [i.split('/')[-2].split('_')[0] for i in k.i],
		ID = lambda k: [i.split('/')[-2].split('_')[1] for i in k.i],
		date = lambda k: [i.split('/')[-2].split('_')[2] for i in k.i],
		acc = lambda k: [i.split('/')[-2].split('_')[3] for i in k.i]).drop(['i'], axis=1).\
	assign(
		segment_len = lambda k: k.groupby(['data', 'ID', 'date', 'acc', 'segment'])['pos'].transform('max'))


vcf = pd.concat([process_vcf(i).\
			assign(i = i)\
		for i in glob.glob(vcf_path)]).\
	assign(
		type = lambda k: [i.split('_')[0] for i in k.segment],
		subtype = lambda k: [get_subtype(i) for i in k.segment],
		segment = lambda k: [i.split('_')[1] for i in k.segment],
		data = lambda k: [i.split('/')[-2].split('_')[0] for i in k.i],
		ID = lambda k: [i.split('/')[-2].split('_')[1] for i in k.i],
		date = lambda k: [i.split('/')[-2].split('_')[2] for i in k.i],
		acc = lambda k: [i.split('/')[-2].split('_')[3] for i in k.i])\
	[['data', 'ID', 'acc', 'date', 'segment', 'type', 'subtype', 'iSNV_loc', 
		'ref_nt', 'iSNV_nt', 'AF', 'DP4_ref', 'DP4_alt', 'DP4', 'SB']]


vcf_raw = pd.concat([process_vcf_af(i).\
		assign(i = i)\
	for i in glob.glob(vcf_raw_path)]).\
	assign(
		type = lambda k: [i.split('_')[0] for i in k.segment],
		subtype = lambda k: [get_subtype(i) for i in k.segment],
		segment = lambda k: [i.split('_')[1] for i in k.segment],
		data = lambda k: [i.split('/')[-2].split('_')[0] for i in k.i],
		ID = lambda k: [i.split('/')[-2].split('_')[1] for i in k.i],
		date = lambda k: [i.split('/')[-2].split('_')[2] for i in k.i],
		acc = lambda k: [i.split('/')[-2].split('_')[3] for i in k.i])\
	[['data', 'ID', 'acc', 'date', 'segment', 'type', 'subtype', 'iSNV_loc', 
		'ref_nt', 'iSNV_nt', 'AF']]


# get only IDs with > 1 sampling timepoint or with home and clinic sample on same date
# add a boolean indicating which replicate we use in the case of replicates
# based on which had maximum depth
# merge in mccrone indicator for clinic sample
use = dp.\
		merge(homochronous_clinic_accs[['acc', 'clinic']],
			how='left',
			on='acc').\
		assign(clinic = lambda k: k.clinic.fillna(0)).\
		groupby(['type', 'data', 'ID']).\
		filter(lambda x: (x['date'].drop_duplicates().shape[0] > 1) | (x['clinic'].drop_duplicates().shape[0] > 1)).\
		groupby(['type', 'data', 'ID', 'clinic', 'date', 'acc']).\
		agg({'depth': 'mean'}).reset_index().\
		rename(columns={'depth': 'genome_depth'}).\
		sort_values(by='genome_depth', ascending=False).\
		assign(use = lambda k: k.groupby(['ID', 'clinic', 'date']).cumcount() == 0,
			duplicate = lambda k: k.groupby(['ID', 'clinic', 'date']).transform('size') > 1)

use_vcf = vcf.merge(use,
	on = ['data', 'ID', 'date', 'acc', 'type']).\
	query('use == True')
use_dp = dp.merge(use,
	on = ['data', 'ID', 'date', 'acc', 'type']).\
	query('use == True')
use_vcf_raw = vcf_raw.merge(use,
	on = ['data', 'ID', 'date', 'acc', 'type']).\
	query('use == True')

# get a list of sample IDs and add a day indicator that matches original study
# with exception of 0.5 adjustment for same-day samples (bio-replicates)
all_dates = pd.concat([
		use_dp[['data', 'ID', 'date', 'acc', 'duplicate']].drop_duplicates().query('data == "swine"').\
			merge(pd.read_csv('data/vaninsberghe_analysis_files/swine_SraRunTable.csv')[['Run', 'Sample Name']].\
					assign(d = lambda k: [int(i.split('_')[1].replace('day', '')) for i in k['Sample Name']])\
						[['Run', 'd']],
				left_on='acc',
				right_on='Run').\
			drop(['Run'], axis=1),
		use_dp[['data', 'ID', 'date', 'acc', 'duplicate']].drop_duplicates().query('data == "mccrone"').\
			merge(pd.read_csv('data/mccroneIAV_analysis_files/mccrone_all_meta.csv')[['ENROLLID', 'onset', 'collect']].\
					drop_duplicates().\
					assign(d = lambda k: (pd.to_datetime(k.collect) - pd.to_datetime(k.onset)).dt.days),
				left_on=['ID', 'date'],
				right_on=['ENROLLID', 'collect'],
				how='left').\
			drop(['onset', 'collect', 'ENROLLID'], axis=1).\
			assign(d = lambda k: k.d.astype(int))]).drop_duplicates().\
		merge(homochronous_clinic_accs,
			on='acc',
			how='left').\
		assign(clinic = lambda k: k.clinic.fillna(0)).\
		assign(d = lambda k: k.d + 0.5*k.clinic*k.groupby(['data', 'ID']).d.transform(lambda g: g.drop_duplicates().shape[0] == 1)).\
		groupby(['data', 'ID']).filter(lambda x: x.d.drop_duplicates().size > 1)

# remove excluded putatively coinfected samples
all_dates = all_dates[~((all_dates['data'].values == 'swine') & \
	(np.isin(all_dates['ID'].values, np.array(['20645', '37345', '20392', '37187', '37349']))))]

all_dates.to_csv('data/all_dates.tsv', sep='\t', index=None)

# make wide vcf, columns are days
vcf_wide_backup = use_vcf[['data', 'ID', 'date', 'clinic', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt', 'AF']].\
	merge(
		all_dates[['data', 'ID', 'date', 'clinic', 'd']].\
			drop_duplicates(),
		on=['data', 'ID', 'date', 'clinic']).\
	drop(['date'], axis=1).\
	pivot(index=['data', 'ID', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt'],
		columns='d',
		values='AF').\
	reset_index().\
	melt(id_vars = ['data', 'ID', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt'],
		value_name = 'AF').\
	merge(
		use_vcf_raw[['data', 'ID', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt', 'date', 'clinic', 'AF']].\
			rename(columns={'AF': 'AF_raw'}).\
			merge(all_dates[['data', 'ID', 'date', 'd', 'clinic']],
				on=['data', 'ID', 'date', 'clinic']).\
			drop(['date', 'clinic'], axis=1),
		on=['data', 'ID', 'd', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt'],
		how='outer').\
	merge(
		all_dates[['data', 'ID', 'd', 'date']],
		on=['data', 'ID', 'd'],
		how='left').\
	query('~date.isnull()').\
	assign(
		AF = lambda k: np.where(k.AF.isnull(), 0, k.AF),
		AF_raw = lambda k: np.where(k.AF_raw.isnull(), 0, k.AF_raw),
		d = lambda k: 'd' + k.d.astype(str))

vcf_wide = vcf_wide_backup.copy()

if vcf_wide.query('(AF > 0) & (AF_raw == 0)').shape[0] > 0:
	raise Exception('all variants in filtered AF files should be in raw AF files')

vcf_wide = vcf_wide.\
	melt(id_vars = ['data', 'ID', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt', 'd', 'date']).\
	pivot(index=['data', 'ID', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt'], columns=['d', 'variable'], values='value').\
	reset_index()

vcf_wide.columns = ['_'.join(i) if i[-1] != '' else i[0] for i in vcf_wide.columns]


# subset just to variants observed at least once in filtered file 
vcf_wide = vcf_wide[
	(vcf_wide[[i for i in vcf_wide.columns if 
		(('d' in i) and (i[-3:] == '_AF'))]] > 0).sum(axis=1) > 0]


if (~vcf_wide[[i for i in vcf_wide.columns if (('d' in i) and (i[-3:] == '_AF'))]].isnull().values).\
		sum(axis=1).min() < 2:
	raise Exception('all variants should have at least two time points of observation')

vcf_wide = vcf_wide.\
	merge(
		use_dp[['data', 'ID', 'date', 'clinic', 'segment', 'pos', 'depth']].\
			rename(columns={'pos': 'iSNV_loc'}).\
			merge(
				all_dates[['data', 'ID', 'date', 'd', 'clinic']].\
					drop_duplicates(),
				on=['data', 'ID', 'date', 'clinic']).\
			assign(d = lambda k: 'depth_d' + k.d.astype(str)).\
			pivot(
				index=['data', 'ID', 'segment', 'iSNV_loc'],
				columns='d', values='depth').\
			reset_index(),
		on=['data', 'ID', 'segment', 'iSNV_loc'],
		how='left').\
	merge(
		use_dp[['data', 'ID', 'segment', 'segment_len']].\
			groupby(['data', 'ID', 'segment'])['segment_len'].max().\
			reset_index(),
		how='left',
		on=['data', 'ID', 'segment'])

# write wide data
vcf_wide.to_csv('data/processed_vcf_wide.csv', index=False)


# format and write long data, drop non-observed iSNVs that are retained for pivoting to wide 
use_vcf = use_vcf.\
	merge(use_vcf_raw[['acc', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt', 'AF']].\
		rename(columns={'AF': 'AF_raw'}),
		how='left',
		on=['acc', 'segment', 'iSNV_loc', 'ref_nt', 'iSNV_nt']).\
	assign(AF = lambda k: np.where(k.AF.isnull(), 0, k.AF)).\
	merge(all_dates,
		on=['data', 'ID', 'date', 'clinic', 'duplicate', 'acc'],
		how='left').\
	merge(use_dp[['data', 'ID', 'date', 'acc', 'segment', 'pos', 'depth']].\
			rename(columns={'pos': 'iSNV_loc'}).drop_duplicates(),
		on=['data', 'ID', 'date', 'acc', 'segment', 'iSNV_loc'],
		how='left').\
	merge(
		use_dp[['data', 'ID', 'segment', 'segment_len']].\
			groupby(['data', 'ID', 'segment'])['segment_len'].max().\
			reset_index(),
		how='left',
		on=['data', 'ID', 'segment'])

use_vcf.to_csv('data/processed_vcf.csv', index=None)

# add universal coordinates to long-form variant data and convert to
# nuc freq matrix
pos_map = pd.read_csv('data/genomes/all_aln_map.tsv', sep='\t').\
	assign(
		seg = lambda k: np.where(k.seg.isna(), "NA", k.seg),
		ID = lambda k: k.ID.astype(str)).\
	assign(seg = lambda k: np.where(k.seg == 'M', 'MP', k.seg)).\
	rename(columns={'seg': 'segment'}).\
	assign(ref = lambda k: k['data'] + '_' + k['ID'] + '_' + k['date'] + '_' + k['acc'])

use_vcf_arr = use_vcf.\
	assign(AF = lambda k: k.AF.astype(float)).\
	pivot(index=
		['data', 'ID', 'acc', 'date', 'd', 'clinic', 'duplicate', 'type', 'subtype', 'segment', 'iSNV_loc',
			'ref_nt', 'depth'],
		columns='iSNV_nt', values='AF').reset_index().\
	assign(
		A = lambda k: np.where(k.ref_nt == 'A', 1 - k[['C', 'G', 'T']].sum(axis=1), k.A),
		C = lambda k: np.where(k.ref_nt == 'C', 1 - k[['A', 'G', 'T']].sum(axis=1), k.C),
		G = lambda k: np.where(k.ref_nt == 'G', 1 - k[['A', 'C', 'T']].sum(axis=1), k.G),
		T = lambda k: np.where(k.ref_nt == 'T', 1 - k[['A', 'C', 'G']].sum(axis=1), k['T'])).\
	fillna(0).\
	rename(columns={'iSNV_loc': 'pos'})
	
# need to merge in which reference was used
# merge in universal position
use_vcf_arr  = use_vcf_arr.\
	assign(ID = lambda k: k.ID.astype(str)).\
	merge(
		pd.read_csv('data/references.tsv', sep='\t', header=None,
				names=['data_ID', 'date', 'ref', 'reads']).\
			assign(data = lambda k: [i.split('_')[0] for i in k.data_ID],
				ID = lambda k: [i.split('_')[1] for i in k.data_ID]).\
			drop(['data_ID', 'reads', 'date'], axis=1),
		how='left',
		on=['data', 'ID']).\
	merge(pos_map[['ref', 'segment', 'upos', 'pos']],
		how='left',
		on=['ref', 'segment', 'pos'])


if use_vcf_arr.query('upos.isnull()').shape[0] > 0:
	raise Exception('some samples missing upos')

# save
use_vcf_arr.drop(['ref'], axis=1).to_csv('data/processed_vcf_arr.csv', index=None)
# save reference file
use_vcf_arr[['data', 'ID', 'date', 'duplicate', 'd', 'clinic', 'type', 'subtype', 'segment', 'ref']].drop_duplicates().\
	assign(
		ref_path = lambda k: 'data/IRMA/' + k.ref + '/'+k.type+'_'+k.segment+'_' + k.subtype + '.fasta').\
	assign(ref_path = lambda k:k.ref_path.str.replace('_.', '.')).\
	to_csv('data/processed_vcf_arr_refs.csv', index=None)




