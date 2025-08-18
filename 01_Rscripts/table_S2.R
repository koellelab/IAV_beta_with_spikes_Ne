# Table S2
set.seed(27)
flumina_0.02_origin <- read_csv("03_data/processed_vcf_wide_0.02.csv", na = c(""))
flumina_0.02_origin$index <- 1:nrow(flumina_0.02_origin)

swine_flumina <- flumina_0.02_origin %>%
  filter(data == "swine") %>%
  mutate(iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt)) %>%
  rowwise()

# exclude rows with only one observation
swine_flumina <- swine_flumina %>%
  data.frame() %>%
  filter(rowSums(!is.na(.[, c(8:15)])) > 1)

# adjust dataset to be the same format as the McCrone data
results_list_flumina <- lapply(1:nrow(swine_flumina), function(i) {
  res <- find_overlapping_pairs(swine_flumina[i, (9:15)])
  if (is.null(res)) return(NULL)
  cbind("Individual ID" = swine_flumina$ID[i], 
        iSNV = swine_flumina$iSNV[i],
        res)
})

swine_flumina_alt <- results_list_flumina %>%
  do.call(rbind, .) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(across(everything(), unlist)) %>%
  setNames(c("ID", "iSNV", "DPS1", "DPS2", "freq1", "freq2")) %>%
  mutate(gen = (DPS2 - DPS1) * (24 / hourPerGen))

swine_flumina_alt$index <- c(1:nrow(swine_flumina_alt))

# exclude rows with both observation equal to 0
swine_flumina_alt <- filter(swine_flumina_alt, freq1 != 0 | freq2 != 0)

# swine IAV Data Subset 2
swine_flumina_alt$sample_time <- paste0(swine_flumina_alt$DPS1, "-", swine_flumina_alt$DPS2)

swine_flumina_sub2 <- swine_flumina_alt %>%
  filter(gen > 0, freq1 == 0) %>%
  group_by(ID, sample_time) %>% # select one sample per individual pig per time interval
  slice_sample(n = 1) %>%
  ungroup()

# Update the Subset column for matching rows
swine_flumina_alt <- swine_flumina_alt %>%
  mutate(Subset = if_else(index %in% swine_flumina_sub2$index, 
                          2, NA)) %>% # assign subsets 2
  group_by(ID, sample_time) %>%
  mutate(Subset = if_else(freq1 != 0 & gen > 0 & row_number() == which.min(abs(freq1 - 0.5)), 
                          1, Subset)) %>% # assign subsets 1
  ungroup()

# export table S2
table_S2 <- swine_flumina_alt[, c("ID", "iSNV", "DPS1", "DPS2", "freq1", "freq2", "Subset")]
colnames(table_S2) <- c("Individual ID", "iSNV", "DPS1", "DPS2", "DPS1 Frequency", "DPS2 Frequency", "Subset")

write.csv(table_S2, "03_data/table_S2.csv", row.names = FALSE, quote = FALSE)
