#install.packages("googlesheets4")
library(googlesheets4)
library(data.table)
library(dplyr)
library(rtracklayer)
library(coloc)

## import significant results from MR analysis 
mr_sig <- 
  read_sheet("https://docs.google.com/spreadsheets/d/1bHrFBF3M1tLwSBX6R3h_iLKeu6NIvbzsFY1QddC2REg/edit?usp=sharing",
             sheet = "PsychENCODE (prefrontal cortex)")
mr_sig$entrez <- toupper(unlist(lapply(strsplit(mr_sig$`SNP (entrez_rsid_pos)`,'-'), '[[', 1)))
mr_sig$rsid <- unlist(lapply(strsplit(mr_sig$`SNP (entrez_rsid_pos)`,'-'), '[[', 2))
mr_sig$pos <- toupper(unlist(lapply(strsplit(mr_sig$`SNP (entrez_rsid_pos)`,'-'), '[[', 3)))
mr_sig$bp <- unlist(lapply(strsplit(mr_sig$pos,':'), '[[', 2))
mr_sig$chr <- paste0("chr",unlist(lapply(strsplit(mr_sig$pos,':'), '[[', 1)))  
mr_sig


# get headers 
expression_cols <- names(fread("~/Desktop/files/data/eqtl/psychencode/DER-08b_hg19_eQTL.bonferroni.txt", sep="\t"))
coloc_regions <- fread('~/Desktop/git/mendelian-randomisation-qtl/region_rs7460093.txt')
colnames(coloc_regions) <- expression_cols[1:14]

# calculate SE
tail <- 2
se <- abs(coloc_regions$regression_slope/ qnorm(coloc_regions$nominal_pval/tail))
coloc_regions$se <- se

# get lookup file for rsid of surrounding SNPs
#lookup_file <- fread('~/Desktop/files/data/eqtl/psychencode/SNP_Information_Table_with_Alleles.txt', sep="\t")

coloc_regions$rsid <- lookup_file$Rsid[match(coloc_regions$SNP_id, lookup_file$PEC_id)] 
coloc_regions$ref <- lookup_file$REF[match(coloc_regions$SNP_id, lookup_file$PEC_id)]
coloc_regions$alt <- lookup_file$ALT[match(coloc_regions$SNP_id, lookup_file$PEC_id)]

## get gwas data for these SNPs
#library(TwoSampleMR)
outcome_data <- fread('~/Desktop/files/data/PGC3-supp/EU/daner_PGC_SCZ_w3_90_0518d_eur.edited.gz')
outcome_data.subset <- subset(outcome_data, SNP %in% coloc_regions$rsid)
outcome_data.subset
outcome_data <- outcome_data.subset[,c(1:6,9:11)]
outcome_data$beta <- log(outcome_data$beta)
rm(outcome_data.subset)
dim(outcome_data) ; dim(coloc_regions)

# format for coloc
# need snp, beta, varbeta, snp, position, type
coloc_regions <- coloc_regions[,c(1,2,8,9,10,12,13,15,16,17,18)]
coloc_regions$SNP_id <- unlist(lapply(strsplit(coloc_regions$SNP_id, ":"), '[[', 2))
colnames(coloc_regions) <- c('gene_id', 'gene_chr', 'pos', "chr", "snp_pos", "pvalues", "beta", "se", "rsid", "ref", 'alt') 
coloc_regions$varbeta <- coloc_regions$se^2
coloc_regions$type <- 'quant'
colnames(coloc_regions)[6] <- "pvalues"
coloc_regions$sdY <- 1
names(outcome_data) <- c('chr', "rsid", "pos", "alt", "ref", "MAF", "beta", "se", "pvalues") 

outcome_data$N <- 130644
outcome_data$type <- "cc"
outcome_data$varbeta <- outcome_data$se^2

library(dplyr)

coloc_regions.arranged<-as.data.frame(coloc_regions) %>% group_by(rsid)  %>% filter(pvalues==min(pvalues)) %>% slice_min(order_by = pvalues)
coloc_regions.subset <- subset(coloc_regions.arranged, rsid %in%outcome_data$rsid)
coloc_regions.subset$pos <- as.numeric(coloc_regions.subset$pos)
coloc_regions.subset$chr <- as.numeric(sub("chr", "", coloc_regions.subset$chr))

coloc_list <- list(
  beta=coloc_regions.subset$beta,
  varbeta=coloc_regions.subset$varbeta,
  snp=coloc_regions.subset$rsid,
  pos=coloc_regions.subset$pos,
  type="quant", sdY=1)

check_dataset(coloc_list)

outcome_list <- list(
  beta=outcome_data$beta,
  varbeta=outcome_data$varbeta,
  snp=outcome_data$rsid,
  pos=outcome_data$pos,
  type="cc")

check_dataset(outcome_list)

results <- coloc.abf(coloc_list, outcome_list)
sink('coloc_rs7460093.txt')
results
sink()


