library(cancereffectsizeR)
library(cowplot)
library(dplyr)
library(ggplot2)
#setwd("/gpfs/ysm/project/townsend/zd72/melanoma")
setwd("/Users/dongzihan/Desktop/melanoma")

tcga_data=read.csv("tcga_melanoma_new_var.csv",header = T,stringsAsFactors = F)
ces_anal = CESAnalysis(genome = "hg19",progression_order = c("WM X Primary","WM X Metastasis","non-WM X Primary","non-WM X Metastasis"))
#ces_anal = CESAnalysis(genome = "hg19",progression_order = c("Primary","Metastasis"))
analysis_size <- load_maf(ces_anal, maf = tcga_data, progression_col = "condition",chain_file = "hg38ToHg19.over.chain")

analysis_size = trinuc_mutation_rates(analysis_size)
analysis_size = gene_mutation_rates(analysis_size,covariate_file ="skcm_pca")
analysis_size = annotate_variants(analysis_size)
analysis_size = ces_snv(analysis_size)
saveRDS(analysis_size,"analysis_size_TCGA_new_condition.rds")

kirc_02=analysis_size=readRDS("analysis_size_TCGA_new_condition.rds")
a=analysis_size@samples
#stage1=a[a$progression_name=="WM X Primary",]$Unique_Patient_Identifier
stage1=a[a$progression_name=="non-WM X Primary",]$Unique_Patient_Identifier

#stage234=a[a$progression_name=="non-WM X Primary",]$Unique_Patient_Identifier
stage234=a[a$progression_name=="non-WM X Metastasis",]$Unique_Patient_Identifier


KIRC_trinuc_stage <- 
  kirc_02@trinucleotide_mutation_weights$trinuc_proportion_matrix
tumor_stage1 <- KIRC_trinuc_stage[stage1, ]
tumor_stage234 <- KIRC_trinuc_stage[stage234, ]

tumor_stage1_avg <- 
  apply(tumor_stage1, 2, mean)
tumor_stage234_avg <- 
  apply(tumor_stage234, 2, mean)

tumor_stage1_avg_ordered <- 
  data.frame(average = tumor_stage1_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

tumor_stage234_avg_ordered <- 
  data.frame(average = tumor_stage234_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)
trinuc.mutation_data <- tumor_stage1_avg_ordered[, -1]

KIRC_trinuc_stage_heatmap_data <- 
  data.frame(deconstrucSig=trinuc.mutation_data$mutation,
             Upstream=trinuc.mutation_data$upstream,
             Downstream=trinuc.mutation_data$downstream,
             mutated_from=trinuc.mutation_data$mut_from,
             mutated_to=trinuc.mutation_data$mut_to,
             trinuc_Size_3_and_less=tumor_stage1_avg_ordered$average,
             trinuc_Size_3.1_and_up=tumor_stage234_avg_ordered$average) %>%
  mutate(mutation = paste0(mutated_from, "\u2192", mutated_to))

KIRC_trinuc_stage1_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3_and_less*100), color="white") + 
  scale_fill_gradient(low="white", high="orange", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3_and_less, 4)*100), size=2) +
  labs(title="non-WM X Primary", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))


KIRC_trinuc_stage234_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3.1_and_up*100), color="white") + 
  scale_fill_gradient(low="white", high="orange", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3.1_and_up, 4)*100), size=2) +
  labs(title="non-WM X Metastasis", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))

KIRC_trinuc_stage1_heatmap
ggsave("tcga_melanoma_trinuc_heatmap_WM_primary.png", width=8, height=2.25)

KIRC_trinuc_stage234_heatmap
ggsave("tcga_melanoma_trinuc_heatmap_non_WM_primary.png", width=8, height=2.25)

combined_trinuc_heatmap_stage <- 
  plot_grid(KIRC_trinuc_stage1_heatmap, KIRC_trinuc_stage234_heatmap, 
            labels = c("A", "B"), align="h", ncol=1)
combined_trinuc_heatmap_stage
ggsave("tcga_melanoma_trinuc_heatmap_non-WM_total.png", width=8, height=4.5)






