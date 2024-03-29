library(data.table)
library(googlesheets4)
library(dplyr)
library(rtracklayer)
library(coloc)
# Shane Crinion / shanecrinion@gmail.com / s.crinion1@nuigalway.ie
# 11-2-2022
# This script processes ugly non-standard GWAS summary stats and formats for TwoSampleMR

# PROCESSES & OUTPUTS
# - 1. Process outcomes
# ├── Import and format GWAS data
# ├── Convert OR to logOR where necessary
# ├── Create var variable
# output: processed GWAS data   
# - 2. Process exposures
# ├── Import and format eQTL data
# ├── Get RSID, ref and alt allele and gene names
# ├── Filter to one SNP per gene
# output: processed eQTL data
# - 3. format datasets
# ├── Format for coloc analysis
# ├── Create list of all required data
# output: data formatted for coloc.abf 
# - 4. Run coloc
# ├── Import required files for formatting
# ├── Perform coloc SNP by SNP
# output: Results files 


# output: metadata file for datasets used


# list location for each GWAS file
outcome_files <- 
  list(`Chronotype`=
         'chronotype/morning_person_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3_logORs.mrbaseformat.txt.gz',
       ASD='ASD/iPSYCH-PGC_ASD_Nov2017.mrbaseformat.tsv.gz', 
       ADHD='ADHD/adhd_eur_jun2017.mrbaseformat.tsv', 
       `BD` = "BIP/Mullins/bip_pgc3_nm_FRQA41917_FRQU371549.mrbaseformat.tsv.gz", 
       Insomnia='insomnia/Insomnia_sumstats_Jansenetal.mrformat.txt.gz', 
       `MDD`='MDD/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.mrbaseformat.tsv', 
       SZ = 'PGC3-supp/EU/daner_PGC_SCZ_w3_90_0518d_eur.edited.gz')
outcome_files <- 
  lapply(outcome_files, function(x) paste0('~/Desktop/files/data/', x))

psychencode_files <- list(headers=c('gene_id', 'gene_chr', 'gene_start', 'gene_end', 'strand', 'number_of_SNPs_tested', 'SNP_distance_to_TSS',
                                    "SNP_id", "SNP_chr", "SNP_start", "SNP_end", "pval_nominal", "beta", "top_SNP"),
                          gene_mapping='~/Desktop/files/current-projects/eQTL-MR/gencode.v19.genes.v7.patched_contigs.gtf',
                          snp_info = '~/Desktop/files/data/eqtl/psychencode/SNP_Information_Table_with_Alleles.txt',
                          coloc_regions='~/Desktop/files/data/eqtl/coloc/top_snps_500kb/psychencode/region_')

gtex_files <- list(headers=c('gene_id', 'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'pval_nominal', 'beta', 'se'),
                   snp_info="~/Desktop/files/data/gtexv8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                   gene_mapping="~/Desktop/files/data/gtexv8/gencode.v26.GRCh38.genes.gtf",
                   cortex="~/Desktop/files/data/eqtl/coloc/top_snps_500kb/gtex_cortex/region_gtexv8cortex_",
                   blood="~/Desktop/files/data/eqtl/coloc/top_snps_500kb/gtex_wholeblood/region_gtexv8wholeblood_",
                   hippocampus='~/Desktop/files/data/eqtl/coloc/top_snps_500kb/gtex_hippocampus/region_gtexv8hippocampus_',
                   hypothalamus='~/Desktop/files/data/eqtl/coloc/top_snps_500kb/gtex_hypothalamus/region_gtexv8hypothalamus_')


# function to format GWAS data, depending on each file
process_outcome <- function(results_file){
  phenotype <- results_file[1,]$outcome
  if(length(unique(results_file$outcome)) > 1){
    stop("Error! More than one type of outcome in results file")
    message(unique(results_file$outcome))}
  else if(phenotype %in% names(outcome_files[-1])){ # -1 refers to chronotype as logOR is in orig data
    print(as.character(outcome_files[phenotype]))
    outcome_data <- fread(as.character(outcome_files[phenotype]))
    print(str(outcome_data))
    outcome_data$beta <- log(outcome_data$beta)
    outcome_data$var <- outcome_data$se^2}
  else if(results_file[1,]$outcome=="Chronotype"){
    outcome_data <- fread(as.character(outcome_files[phenotype]))
    outcome_data$var <- outcome_data$se^2}
  else {outcome_data <- "Detected no matching outcome"}
  return(list(data=outcome_data, phenotype=phenotype))}


process_eqtl <- function(snp, eqtl, tissue, eqtl_source, filtering=2, lookup_file, gtf){
  if(eqtl_source=="psychencode" & tissue=="prefrontal cortex"){
    coloc_regions <- fread(paste0(unlist(psychencode_files['coloc_regions']), snp, '.txt'))
    colnames(coloc_regions) <- psychencode_files$headers
    coloc_regions$se <- abs(coloc_regions$beta/qnorm(coloc_regions$pval_nominal/2))
    # extract desired data from files
    coloc_regions$rsid <- lookup_file$Rsid[match(coloc_regions$SNP_id, lookup_file$PEC_id)] 
    coloc_regions$ref <- lookup_file$REF[match(coloc_regions$SNP_id, lookup_file$PEC_id)]
    coloc_regions$alt <- lookup_file$ALT[match(coloc_regions$SNP_id, lookup_file$PEC_id)]
  } else if(eqtl_source=="gtex"){
    coloc_regions<-fread(paste0(unlist(gtex_files[tissue]) ,snp, '.txt'))
    colnames(coloc_regions)<- gtex_files$headers
    coloc_regions$rsid <- lookup_file$rs_id_dbSNP151_GRCh38p7[match(coloc_regions$variant_id, lookup_file$variant_id)] 
    coloc_regions$ref <- lookup_file$ref[match(coloc_regions$variant_id, lookup_file$variant_id)]
    coloc_regions$alt <- lookup_file$alt[match(coloc_regions$variant_id, lookup_file$variant_id)]
    coloc_regions$gene_name <- gtf$gene_name[match(coloc_regions$gene_id, gtf$gene_id)] }
  if(filtering==1){
    coloc_regions<- as.data.frame(na.omit(coloc_regions)) %>% group_by(rsid)  %>% filter(pval_nominal==min(pval_nominal))
  } else if (filtering==2){
    desired_pair <- coloc_regions[coloc_regions$gene_id==eqtl & coloc_regions$rsid==snp,]
    remaining_snps <- as.data.frame(na.omit(coloc_regions[coloc_regions$rsid!=snp,])) %>% group_by(rsid)  %>% filter(pval_nominal==min(pval_nominal))
    coloc_regions <- rbind(as.data.frame(desired_pair), as.data.frame(remaining_snps))}
  return(coloc_regions)}


format_data <-function(coloc_regions, outcome_data, snp, eqtl, tissue, eqtl_source){
  # get subset of SNPs in both datasets
  outcome_data <- subset(outcome_data, SNP %in% coloc_regions$rsid)
  coloc_regions <- subset(coloc_regions, rsid %in% outcome_data$SNP)
  
  # get a joined dataset (in case other methods require harmonisation) - harmonisation doesn't change these results
  joined_data <- full_join(coloc_regions, outcome_data[!duplicated(outcome_data$SNP),], by=c("rsid"="SNP"))
  write.csv(joined_data, paste0('coloc/', outcome_data$outcome,'/', snp, '_', eqtl, '_ ', tissue, '_', eqtl_source, '.csv'))
  
  # make the colocalisation format lists
  coloc_list <- list(
    beta=coloc_regions$beta, # all the reference are different in the expression 
    varbeta=coloc_regions$se^2,
    snp=coloc_regions$rsid,
    #pos=coloc_regions$pos,
    type="quant", sdY=1)
  str(coloc_list)
  check_dataset(coloc_list)
  outcome_list <- list(
    beta=-1*outcome_data[!duplicated(outcome_data$SNP),]$beta,
    varbeta=outcome_data[!duplicated(outcome_data$SNP),]$var,
    snp=outcome_data[!duplicated(outcome_data$SNP),]$SNP,
    pos=outcome_data[!duplicated(outcome_data$SNP),]$pos,
   # N=outcome_data[!duplicated(outcome_data$SNP),]$N,
    #pvalues=outcome_data[!duplicated(outcome_data$SNP),]$pval,
    type="cc")
  check_dataset(outcome_list)
  str(outcome_list)
  formatted_data = list(eqtl=coloc_list, outcome=outcome_list)
  return(formatted_data)
}
#lookup_file <- fread(gtex_files$snp_info, sep="\t")
#gtf <- import(gtex_files$gene_mapping)

## create function to perform colocalisation
run_coloc <- function(sig_data, tissue, eqtl_source, filtering){
  for(trait in unique(sig_data$outcome)){
  outcome_data <-  process_outcome(subset(sig_data, outcome==trait))
  sprintf('outcome data is processed')
  if(eqtl_source=="gtex"){
    print('importing eqtl info files')
    expression_cols <- gtex_files$headers
    lookup_file <- fread(gtex_files$snp_info, sep="\t")
    gtf <- import(gtex_files$gene_mapping)
  } else if(eqtl_source=="psychencode"){
    #expression_cols <- psychencode_files$headers
    #gtf <- import(unlist(psychencode_files['gene_mapping'])) 
    lookup_file <- fread(unlist(psychencode_files['snp_info']), sep="\t")}
    for(i in 1:nrow(sig_data)){
      print(strsplit(as.character(sig_data[i,"SNP (entrez-rsid-pos)"]), '-') [[1]][2])
      top_snp <- strsplit(as.character(sig_data[i,"SNP (entrez-rsid-pos)"]), '-') [[1]][2]
      top_entrez <- toupper(strsplit(as.character(sig_data[i,"SNP (entrez-rsid-pos)"]), '-') [[1]][1])
      coloc_regions <- process_eqtl(snp= top_snp, eqtl=top_entrez, tissue=tissue, eqtl_source=eqtl_source, 
                                    filtering=filtering, lookup_file=lookup_file, gtf=gtf)
      data <- format_data(coloc_regions, outcome_data$data, snp=top_snp, eqtl=top_entrez, tissue=tissue, eqtl_source=eqtl_source)
      results <- coloc.abf(data$eqtl, data$outcome)
      phenotype <- stringr::str_replace_all(outcome_data$phenotype, "[[:punct:]]", "")
      phenotype <- stringr::str_replace_all(phenotype, " ", "")
      prefix <- paste0(top_snp, '_', top_entrez, '_', tissue, '_', eqtl_source, '_', phenotype)
      write.csv(x=results$summary, file=paste0('coloc/', phenotype, '/', prefix, '_coloc.beta_change.summary.csv'))
      write.csv(x=results$results, file=paste0('coloc/', phenotype, '/', prefix, '_coloc.beta_change.results.csv'))
      write.csv(x=results$priors, file=paste0('coloc/', phenotype, '/',prefix, '_coloc.beta_change.priors.csv'))}}}


# hypothalamus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hypothalamus)")
unique(data$outcome)
run_coloc(subset(data, outcome=="SZ"), tissue="hypothalamus", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="hypothalamus", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="hypothalamus", eqtl_source = 'gtex', filtering=2)

# hippocampus
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hippocampus)")
unique(data$outcome)
run_coloc(subset(data, outcome=="SZ"), tissue="hippocampus", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="hippocampus", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="hippocampus", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="ASD"), tissue="hippocampus", eqtl_source = 'gtex', filtering=2)

data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (whole blood)")
run_coloc(subset(data, outcome=="SZ"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="ADHD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="MDD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Insomnia"), tissue="blood", eqtl_source = 'gtex', filtering=2)


data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (cortex)")
run_coloc(subset(data, outcome=="SZ"), tissue="cortex", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="cortex", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="cortex", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="MDD"), tissue="cortex", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Insomnia"), tissue="cortex", eqtl_source = 'gtex', filtering=2)



data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (whole blood)")
unique(data$outcome)
run_coloc(subset(data, outcome=="SZ"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="MDD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="ADHD"), tissue="blood", eqtl_source = 'gtex', filtering=2)
run_coloc(subset(data, outcome=="Insomnia"), tissue="blood", eqtl_source = 'gtex', filtering=2)



data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "PsychENCODE (prefrontal cortex)")
unique(data$outcome)
run_coloc(subset(data, outcome=="SZ"), tissue="prefrontal cortex", eqtl_source = 'psychencode', filtering=2)
run_coloc(subset(data, outcome=="BD"), tissue="prefrontal cortex", eqtl_source = 'psychencode', filtering=2)
run_coloc(subset(data, outcome=="Chronotype"), tissue="prefrontal cortex", eqtl_source = 'psychencode', filtering=2)
run_coloc(subset(data, outcome=="MDD"), tissue="prefrontal cortex", eqtl_source = 'psychencode', filtering=2)
run_coloc(subset(data, outcome=="Insomnia"), tissue="prefrontal cortex", eqtl_source = 'psychencode', filtering=2)

data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1bHrFBF3M1tLwSBX6R3h_iLKeu6NIvbzsFY1QddC2REg/edit?usp=sharing",sheet = "GTEx v8 (whole blood)")


run_coloc(subset(data, outcome %in% c("Major Depressive Disorder")), 
          tissue="blood", eqtl_source = 'gtex', filtering=2)
       