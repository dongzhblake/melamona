
setwd("/Users/Desktop/melanoma")
library(ggplot2)
library(ggrepel)
library(crayon)
mut_rate=list()
MAF_files=c("TCGA-SKCM-maf.txt","TCGA-SKCM-maf.txt","Halaban-SKCM-maf.txt","Halaban-SKCM-maf.txt","Yale-Gilead-maf.txt","Yale-Gilead-maf.txt")
covariate_file = c("skcm_pca","default","skcm_pca","default","skcm_pca","default")
for(i in 1:6){
  ces_anal = CESAnalysis(genome = "hg19", 
                         progression_order = c("Primary","Metastasis"))
  MAF=read.delim(MAF_files[i],stringsAsFactors = F)
if(i %in% c(1,2)){
analysis_size <- load_maf(ces_anal, maf = MAF, progression_col = "stage",chain_file = "hg38ToHg19.over.chain")
}
  else{
    analysis_size <- load_maf(ces_anal, maf = MAF, progression_col = "stage")
    
  }
analysis_size = trinuc_mutation_rates(analysis_size)
analysis_size = gene_mutation_rates(analysis_size,covariate_file =covariate_file[i])
mut_rate[[i]] <- data.frame(gene = analysis_size@mutrates[[1]],
                       Primary = analysis_size@mutrates[[2]],
                       Metastasis = analysis_size@mutrates[[3]])
}
i=
scatterplot_gene_mutrate <- ggplot(data = mut_rate[[i]], aes(x=Primary, y=Metastasis))+
  geom_point(size=1.5, color="steelblue", alpha=0.4) +
  ggtitle(paste0(MAF_files[i],"_Gene-level mutation rates")) +
  labs(subtitle =paste0("covariate_file:",covariate_file[i])) +
  xlab('Primary') +
  ylab('Metastasis') +
  geom_smooth(method="lm", color="navyblue") +
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  theme_bw() 
scatterplot_gene_mutrate








