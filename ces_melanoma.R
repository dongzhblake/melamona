library(cancereffectsizeR)
setwd("/Users/Desktop/melanoma")
tcga_data=read.delim("TCGA-SKCM-maf.txt",stringsAsFactors = F)

ces_anal = CESAnalysis(genome = "hg19",progression_order = c("Metastasis","Primary"))
analysis_size <- load_maf(ces_anal, maf = tcga_data,chain_file = "hg38ToHg19.over.chain",progression_col = "stage")
MAF=read.delim("Yale-Gilead-maf.txt",stringsAsFactors = F)
analysis_size <- load_maf(analysis_size, maf = MAF,progression_col = "stage")
MAF=read.delim("Halaban-SKCM-maf.txt",stringsAsFactors = F)
analysis_size <- load_maf(analysis_size, maf = MAF,progression_col = "stage")

analysis_size = trinuc_mutation_rates(analysis_size)
analysis_size = gene_mutation_rates(analysis_size,covariate_file ="skcm_pca")
analysis_size = annotate_variants(analysis_size)
analysis_size = ces_snv(analysis_size)

saveRDS(analysis_size,"analysis_size_all_old_condition.rds")




