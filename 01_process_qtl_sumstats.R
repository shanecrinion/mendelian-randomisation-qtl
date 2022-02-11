# #!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@nuigalway.ie
# 10-2-2022
# This script processes ugly non-standard QTL summary stats and formats for TwoSampleMR

# PROCESSES & OUTPUTS
# - clean and tidy QTL sumstats
# ├── import summary stats, lookup and header files
# ├── extract significant SNPs
# ├── extract annotations from lookup file
# ├── format data
# ├── clump data
# output: processed sumstats for QTL data
# - create metadata file
# ├── collect data name, description, source & SNP and gene counts
# output: metadata file for datasets used

# DATASETS
# -BRAIN
# - PSYCHencode ✓
# ├── Prefrontal cortex ✓
# - GTEx 
# ├── Cortex ✓
# ├── Hypothalamus
# ├── Full Brain 


# import libraries
suppressPackageStartupMessages({
  library(data.table) # import data
  library(TwoSampleMR)})

setwd('~/Desktop/files/data/eqtl/')

# set up metadata file
header <- paste("data", "region", "individuals","individuals.descr","count.variants","count.5e-8","count.duplicated","snp.gene.pairs","snps.unique", "genes.unique", "source", "data.description", sep = ',')
write(header,file="metadata.csv", append=TRUE)

write_metadata <- function(){
  line <- paste(name, region, individuals,  individuals_descr, count.variants, count.sig, count.duplicated, count.snps, count.snps.unique, count.genes.unique, source, data_descr, sep=",")
  write(line,file="metadata.csv", append=TRUE)}

## 1. PSYCHencode 
# info for metadata file
name <- "PSYCHencode"
region <- "prefrontal cortex"
individuals <- 1387
individuals_descr <- "Filtered adult samples with matching gene expression and genotypes (679 healthy controls + 497 schizophrenia + 172 bipolar disorder + 31 autism spectrum disorder and 8 affective disorder patients)" 
data_descr <-  "Full set of cis-eQTLs with no p-value or FDR filtering"
source <- "http://resource.psychencode.org/Datasets/Derived/QTLs/Full_hg19_cis-eQTL.txt.gz"


# import data
eQTL <- fread('psychencode/psychencode_Full_hg19_cis-eQTL.txt.gz',sep = ' ')
colnames(eQTL) <- names(fread("psychencode/DER-08b_hg19_eQTL.bonferroni.txt", sep="\t"))[1:14]
eQTL.lookup <- fread('~/Desktop/files/data/eqtl/psychencode/SNP_Information_Table_with_Alleles.txt', sep="\t")
eQTL.sig <- subset(eQTL, nominal_pval < 5e-8)

# get RSID and alleles for remaining SNPs
eQTL.sig$rsid <- eQTL.lookup$Rsid[match(eQTL.sig$SNP_id, eQTL.lookup$PEC_id)] 
eQTL.sig$ref <- eQTL.lookup$REF[match(eQTL.sig$SNP_id, eQTL.lookup$PEC_id)]
eQTL.sig$alt <- eQTL.lookup$ALT[match(eQTL.sig$SNP_id, eQTL.lookup$PEC_id)]

# calculate SE
tail <- 2
se <- abs(eQTL.sig$regression_slope/ qnorm(eQTL.sig$nominal_pval/tail))
eQTL.sig$se <- se

# save counts for metadata file
count.variants <- dim(eQTL)[1] 
count.sig <- dim(eQTL.sig)[1]
count.duplicated <- sum(duplicated(eQTL.sig$rsid))

# create a unique locator for each gene-SNP pair (many SNPs associated with > 1 gene)
eQTL.sig$locator <- paste(eQTL.sig$gene_id, eQTL.sig$rsid, eQTL.sig$SNP_id, sep="-")


# format for TwoSampleMR (locator col used for SNPs to preserve all unique gene-SNP pairs)
exposure_data <- format_data(dat = eQTL.sig, type = "exposure",snp_col ="locator",
                             beta_col = "regression_slope",se_col =  "se",
                             effect_allele_col = "alt",other_allele_col = "ref",pval_col = "nominal_pval")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_PSYCHencode_eQTL.csv')
rm(eQTL, eQTL.lookup, eQTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_PSYCHencode_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

# write metadata to file
write_metadata()

