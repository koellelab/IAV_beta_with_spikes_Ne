import pandas as pd
import numpy as np
import argparse
try:
	from utils import import_aln
except:
	from scripts.utils import import_aln


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--seqNameSep', default='|')
	args = parser.parse_args()
	#args.seqs = 'data/genomes/PA_aln copy.fasta
	#args.seqs = 'data/genomes/NS_aln.fasta'
	s_names, s_arr = import_aln(open(args.seqs, 'r'))
	# format names
	#s_names = np.array([i.split(args.seqNameSep)[0] for i in s_names])
	out = np.full(s_arr.shape, np.nan)
	out[np.where(s_arr != "-")] = \
		np.hstack([np.arange(i) for i in (s_arr != "-").sum(axis=1)])
	# adding 1 because variants coordinates are 1-indexed 
	df = pd.DataFrame(
		out+1, 
		index=s_names, 
		columns=np.arange(out.shape[1])+1).T.\
	reset_index().melt(id_vars='index')
	# everything below here specific code for this project
	df[['data', 'ID', 'date', 'acc']] = df['variable'].str.split('_', expand=True)
	df = df.drop(['variable'], axis=1).\
		rename(columns={'index': 'upos', 'value': 'pos'}).\
		query('~pos.isnull()').\
		assign(
			ID = lambda k: [i.replace('|', '_') for i in k['ID']],
			pos = lambda k: k.pos.astype(int))
		
	df.to_csv(args.seqs.replace('.fasta', '_map.tsv'), sep='\t', index=None)


if __name__ == "__main__":
    run()



