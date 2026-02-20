# Table S2
set.seed(27)
threshold <- 0.02
flumina_0.02_origin <- read_csv("03_data/processed_vcf_wide_minAF0.02.csv")[, c(1:6, 28:37)]
flumina_0.02_origin$index <- 1:nrow(flumina_0.02_origin)

dps_values <- c(-2, 0, 0.25, 1, 1.25, 2:6)

McCrone_IAV_flumina_rep <- flumina_0.02_origin %>%
  filter(data == "mccrone") %>%
  rowwise() %>%
  mutate(DPS1 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][1],
         DPS2 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][2], 
         gen = (DPS2 - DPS1)*(24/hourPerGen),
         iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt),
         non_na_values = list(na.omit(c_across(c(`d-2`, d0:d6)))),
         freq1 = non_na_values[1], 
         freq2 = non_na_values[2]) %>% 
  select(-c(non_na_values, DPS2)) %>%
  filter(gen == 1)

# export table S1
table_S2 <- McCrone_IAV_flumina_rep[, c("ID", "iSNV", "DPS1", "freq1", "freq2")]
colnames(table_S2) <- c("Individual ID", "iSNV", "DPS", "First Frequency", "Second Frequency", "Subset")

write.csv(table_S2, "03_data/table_S2.csv", row.names = FALSE, quote = FALSE)
