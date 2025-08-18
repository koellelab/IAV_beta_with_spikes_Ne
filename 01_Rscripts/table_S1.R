# Table S1 
set.seed(27)
flumina_0.02_origin <- read_csv("03_data/processed_vcf_wide_0.02.csv", na = c(""))
flumina_0.02_origin$index <- 1:nrow(flumina_0.02_origin)

dps_values <- c(-2, 0:6)

McCrone_IAV_flumina_data <- flumina_0.02_origin %>%
  filter(data == "mccrone") %>%
  rowwise() %>%
  mutate(DPS1 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][1],
         DPS2 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][2], 
         gen = (DPS2 - DPS1)*(24/hourPerGen),
         iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt),
         non_na_values = list(na.omit(c_across(c(`d-2`, d0:d6)))),
         freq1 = non_na_values[1], 
         freq2 = non_na_values[2]) %>% 
  select(-non_na_values) %>%
  ungroup()

# Human IAV Data Subset 2
McCrone_IAV_flumina_data_sub2 <- McCrone_IAV_flumina_data %>%
  filter(gen > 0, freq1 == 0) %>%
  group_by(ID) %>% 
  slice_sample(n = 1) %>%
  ungroup()

# Update the Subset column for matching rows
McCrone_IAV_flumina_data <- McCrone_IAV_flumina_data %>%
  mutate(Subset = if_else(index %in% McCrone_IAV_flumina_data_sub2$index, 
                          2, NA)) %>% # assign subsets 2
  group_by(ID) %>%
  mutate(Subset = if_else(freq1 != 0 & gen > 0 & row_number() == which.min(abs(freq1 - 0.5)), 
                          1, Subset)) %>% # assign subsets 1
  ungroup()

# export table S1
table_S1 <- McCrone_IAV_flumina_data[, c("ID", "iSNV", "DPS1", "DPS2", "freq1", "freq2", "Subset")]
colnames(table_S1) <- c("Individual ID", "iSNV", "DPS1", "DPS2", "DPS1 Frequency", "DPS2 Frequency", "Subset")

write.csv(table_S1, "03_data/table_S1.csv", row.names = FALSE, quote = FALSE)
