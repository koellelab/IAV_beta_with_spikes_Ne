library(tidyverse)

d = read_csv('data/processed_vcf_wide.csv',
	show_col_types=FALSE)

# set filtering params
min_af = 0.02
min_dp = 500
min_pos = 200

af_cols = colnames(d)[sapply(colnames(d), function(x){substr(x, nchar(x)-2, nchar(x))}) == '_AF']
raw_af_cols = colnames(d)[sapply(colnames(d), function(x){substr(x, nchar(x)-5, nchar(x))}) == 'AF_raw']
depth_cols = colnames(d)[substr(colnames(d), 1, 6) == 'depth_']

# identify which variants to keep
# which sites have a filtered AF >= min_af at at least one time point
keep = rowSums(d[af_cols] >= min_af, na.rm=TRUE) > 0


# which sites do not have an unfiltered AF >= min_af but are removed by filtering
keep = keep & (rowSums(
	(d[af_cols] == 0) & 
		(d[raw_af_cols] >= min_af),
	na.rm=TRUE) == 0)

# which sites have depth >= min_dp at all sampled time points
keep = keep & (rowSums(
	d[depth_cols] < min_dp,
	na.rm=TRUE) == 0)

# which sites are not located in the first or last min_pos of asegment
keep = keep & (d$iSNV_loc > min_pos) & ((d$segment_len - d$iSNV_loc) > min_pos)

# finally, filter wide variant file
keep_d = d[keep,1:ncol(d)]

write_csv(keep_d, 
	paste(c('data/processed_vcf_wide_filtered_minAF', min_af, '_minDP', min_dp, '_minPOS', min_pos, '.csv'), collapse=''))
