library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/dongzihan/Desktop/melanoma")

#kirc_02=analysis_size=readRDS("analysis_size_TCGA.rds")
kirc_02=analysis_size=readRDS("analysis_size_TCGA_new_condition.rds")
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
             group = kirc_02@selection_results$progression) %>%
  filter(selection_intensity > 1)

stage_data=stage_data[stage_data$group %in% c("non-WM X Primary","non-WM X Metastasis"),]
stage_data$group=droplevels(stage_data$group)

SR=data.frame(kirc_02@selection_results)
stage1_info <- summary_multivar(stage_data, group_name = "non-WM X Primary")
stage234_info <- summary_multivar(stage_data, group_name = "non-WM X Metastasis")

genes_Pri=as.character(stage1_info[order(stage1_info$mean_si,decreasing = T),]$gene[1:10])
genes_Meta=as.character(stage234_info[order(stage234_info$mean_si,decreasing = T),]$gene[1:10])

top_gene_stage=genes_Pri


top_gene_stage=genes_Meta


stage_plot_data <- stage_data %>% 
  filter(gene %in% top_gene_stage) %>%
  mutate(gene = factor(gene, levels = top_gene_stage))



si_boxplot <- function(data, group_names, genes, colormap1,colormap2, color_num1,color_num2,yticks,main) {
  palette <- c(brewer.pal(6, colormap1)[color_num1],brewer.pal(6, colormap2)[color_num2])
  myplt <-
    boxplot(selection_intensity ~ group*gene, data = data, boxwex=0.4,
            col = palette, xlab = "", ylab = "", 
            xaxt="n", yaxt="n")
  title(ylab = expression(paste("Selection intensity /", "10"^"4")), 
        mgp = c(2, 0, 0))
  title(xlab = "Gene", mgp = c(2, 0, 0))
  title(main=main)
  
  axis(1, mgp = c(0, 0.2, 0), 
       at = seq(1.5 , 20 , 2), 
       labels = genes, 
       tick=FALSE , cex.axis=0.5)
  
  axis(2, at = yticks * 1e4, las=2,
       labels = yticks, cex.axis=0.6)
  
  # Add the grey vertical lines
  for(i in seq(0.5 , 21 , 2)){ 
    abline(v=i, lty=1, col="grey")
  }
  
  # Add a legend
  legend("topright", legend = group_names, 
         col=palette, 
         pch = 15, bty = "n", pt.cex = 2, cex = 1,  horiz = F)
}
stage_plot_data$group=relevel(stage_plot_data$group,"non-WM X Primary")


#2.078740e+06  is the max of gene BRAF
#stage_plot_data=stage_plot_data[stage_plot_data$selection_intensity<max(stage_plot_data$selection_intensity),]
pdf("non-WM_non-WM_Primary_top10_genes.pdf", width = 8, height = 6)
si_boxplot(stage_plot_data, c( "non-WM X Primary","non-WM X Metastasis"), 
           top_gene_stage, "Oranges","Oranges",2,6, yticks = seq(0, 40, 2.5),main="Top10_non-WM_Primary_genes")
dev.off()
  


pdf("non-WM_non-WM_Met_top10_genes.pdf", width = 8, height = 6)
si_boxplot(stage_plot_data, c( "non-WM X Primary","non-WM X Metastasis"), 
           top_gene_stage, "Oranges","Oranges",2,6,yticks = seq(0, 40, 2.5),main="Top10_non-WM_Metastasis_genes")
dev.off()




