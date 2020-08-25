setwd("/Users/Desktop/melonoma")
kirc_02=readRDS("analysis_size.rds")



kirc_02@selection_results$gene=apply(data.frame(kirc_02@selection_results[,1]),1,FUN = function(x){strsplit(x,"_")[[1]][1]})
summary_multivar <- function(data, group_name) {
  data_clean <- data %>% 
    filter(group == group_name) %>%
    arrange(desc(selection_intensity)) %>%
    filter(selection_intensity > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(selection_intensity), # change sum to mean and sd
              mean_si = mean(selection_intensity),
              sd = sd(selection_intensity),
              max_si = max(selection_intensity),
              n_variant = n_distinct(variant)) %>%
    filter(n_variant > 1) 
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  return(merge_info)
}

stage_data <- 
  data.frame(variant = kirc_02@selection_results$variant,
             gene = kirc_02@selection_results$gene,
             selection_intensity = kirc_02@selection_results$selection_intensity,
             group = kirc_02@selection_results$progression)

stage1_info <- summary_multivar(stage_data, group_name = "Primary")
stage234_info <- summary_multivar(stage_data, group_name = "Metastasis")

fill_na <- function(x, fill = 0) {
  x = ifelse(is.na(x), fill, x)
  return(x)
}


sd_info <- stage_data %>% 
  filter(selection_intensity > 1) %>%
  group_by(gene) %>% summarise(sd = sd(selection_intensity))

stage_merge <- merge(stage1_info, stage234_info, by = "gene", all = T, 
                     suffixes = c(".1", ".234")) %>%
  mutate_at(c("cum_si.1", "cum_si.234", 
              "mean_si.1", "mean_si.234", 
              "sd.1", "sd.234",
              "n_variant.1", "n_variant.234"), fill_na) %>%
  mutate(n_variant_total = n_variant.1 + n_variant.234, 
         mean_si_total = (cum_si.1 + cum_si.234) / n_variant_total) %>%
  arrange(desc(n_variant_total))

stage_merge <- stage_merge %>%
  mutate(n_variant = paste(as.character(n_variant.1), 
                           as.character(n_variant.234), sep = "|"),
         topvar.1 = ifelse(is.na(variant.1), "NA", 
                           paste0(variant.1, "(", 
                                  round(max_si.1, 1),"x10^0)")),
         topvar.234 = ifelse(is.na(variant.234), "NA", 
                             paste0(variant.234, "(", 
                                    round(max_si.234, 1),"x10^0)")),
         topvar = paste(topvar.1, topvar.234, sep = "|"),
         mean_si = paste(round(mean_si.1, 2), round(mean_si.234, 2), sep = "|"),
         mean_si_total = round(mean_si_total, 2),
         sd = paste(round(sd.1, 2), round(sd.234, 2), sep = "|")
  )
stage_merge_sub <- stage_merge %>%
  select("gene", "n_variant", "mean_si", "sd", "topvar") %>%
  # arrange(desc(cum_si_total)) %>%
  mutate(format = "Primary | Metastasis")

names(stage_merge_sub) <- 
  c("gene", "#variant", "mean SI", "sd", "TopVariant", "format")

write.csv(stage_merge_sub, file = "tcga_melonoma_gene_summary.csv", quote = F, row.names = F)






