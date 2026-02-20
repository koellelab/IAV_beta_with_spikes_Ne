library(tidyverse)

d = read_csv("03_data/processed_vcf_wide_rep.csv",
	show_col_types=FALSE)

# set filtering params
min_af = 0.005
min_dp = 500
min_pos = 200

af_cols = colnames(d)[sapply(colnames(d), function(x){substr(x, nchar(x)-2, nchar(x))}) == "_AF"]
raw_af_cols = colnames(d)[sapply(colnames(d), function(x){substr(x, nchar(x)-5, nchar(x))}) == "AF_raw"]
depth_cols = colnames(d)[substr(colnames(d), 1, 6) == "depth_"]

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

# rearrange columns
# identify AF columns only
af_cols <- names(keep_d)[str_detect(names(keep_d), "^d-?\\d+\\.?\\d*_AF$")]

# sort AF columns by numeric day
af_cols_sorted <- af_cols[
  order(
    as.numeric(str_extract(af_cols, "(?<=d)-?\\d+\\.?\\d*"))
  )
]

# rebuild dataframe: non-AF cols stay put, AF cols reordered
keep_d <- keep_d %>%
  select(
    everything(),
    -all_of(af_cols)
  ) %>%
  bind_cols(keep_d %>% select(all_of(af_cols_sorted)))

# rename columns
keep_d <- keep_d %>%
  dplyr::rename(
    `d-2` = `d-2.0_AF`,
    d0    = d0.0_AF,
    d0.25 = d0.5_AF,
    d1    = d1.0_AF,
    d1.25    = d1.5_AF,
    d2    = d2.0_AF,
    d3    = d3.0_AF,
    d4    = d4.0_AF,
    d5    = d5.0_AF,
    d6    = d6.0_AF
  )

write.csv(keep_d, 
	paste0("03_data/processed_vcf_wide_minAF", min_af, ".csv"), row.names = FALSE, quote = FALSE, na = "")
