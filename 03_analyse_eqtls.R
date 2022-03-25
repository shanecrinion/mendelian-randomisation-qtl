#cortex
library(TwoSampleMR)
setwd('~/Desktop/files/data/eqtl/')

View(available_outcomes())
exposure_data <- read.csv('processed_GTExCortex_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956", "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)
harm_data$SNP <- harm_data$pos.y

library(dplyr)
mr_wr <- mr_singlesnp(dat = harm_data)
singlesnp<- harm_data[4,]
example <- mr_singlesnp(singlesnp)

mr_wald_ratio(singlesnp$beta.exposure, singlesnp$beta.outcome, singlesnp$se.exposure, singlesnp$se.outcome)

dim(mr_wr)[1]
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))
dim(mr_wr)
library(rtracklayer)
my_gtf <- import('~/Desktop/files/current-projects/eQTL-MR/gencode.v19.genes.v7.patched_contigs.gtf')
my_gtf.df<- my_gtf 
mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]
dim(mr_wr)
names(mr_wr)
dim(subset(mr_wr, p<(0.05/dim(mr_wr)[1])))

write.csv(subset(mr_wr.sig, p<(0.05/4778)) %>% arrange(p),
          file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExCortex.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExCortex.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExCortex.byoutcome.alloutcomes.csv')

rm(exposure_data, harm_data, mr_wr, mr_wr.sig, outcome_data)

#blood
exposure_data <- read.csv('processed_GTExWholeBlood_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956", "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)



harm_data$SNP <- harm_data$pos.y
library(dplyr)
mr_wr <- mr_singlesnp(dat = harm_data)
dim(mr_wr)
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))


#library(rtracklayer)
#my_gtf <- import('~/Desktop/files/current-projects/eQTL-MR/gencode.v19.genes.v7.patched_contigs.gtf')
#my_gtf.df<- my_gtf 

mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]
dim(subset(mr_wr, p<0.05/dim(mr_wr)[1]))
write.csv(subset(mr_wr.sig, p<0.05/dim(mr_wr)[1]) %>% arrange(p),
          file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExBlood.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExWholeBlood.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExWholeBlood.byoutcome.alloutcomes.csv')
dim(mr_wr)

hypothalamus 
exposure_data <- read.csv('processed_GTExHypothalamus_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956", "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)
harm_data$SNP <- harm_data$pos.y
library(dplyr)
mr_wr <- mr_singlesnp(dat = harm_data)
dim(mr_wr)
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))

#library(rtracklayer)
#my_gtf <- import('~/Desktop/files/current-projects/eQTL-MR/gencode.v19.genes.v7.patched_contigs.gtf')
#my_gtf.df<- my_gtf 

mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]
dim(subset(mr_wr, p<0.05/dim(mr_wr)[1]))
dim(mr_wr)
write.csv(subset(mr_wr.sig, p<0.05/dim(mr_wr)[1]) %>% arrange(p),
          file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExHypothalamus.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExHypothalamus.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = '~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExHypothalamus.byoutcome.alloutcomes.csv')
dim(subset(mr_wr.sig, p<0.05/dim(mr_wr)[1]))
dim(mr_wr)[1]


# PSYCHencode
pnorm(abs(singlesnp$beta.outcome/singlesnp$beta.exposure)/
        (singlesnp$se.outcome/singlesnp$beta.exposure), 
      lower.tail = FALSE) * 2

((singlesnp$beta.outcome/singlesnp$beta.exposure)/ (singlesnp$se.outcome/singlesnp$beta.exposure) )^2

mr_wald_ratio(singlesnp$beta.exposure, singlesnp$beta.outcome, singlesnp$se.exposure, singlesnp$se.outcome)

chisq.test(((singlesnp$beta.outcome/singlesnp$beta.exposure)/ (singlesnp$se.outcome/singlesnp$beta.exposure) )^2, y=0)


# check intersection between
psychencode <- read.csv('~/Desktop/files/data/eqtl/processed_PSYCHencode_eQTL.clumped.csv')
cortex <- read.csv('processed_GTExWholeBlood_eQTL.clumped.csv')
length(intersect(cortex$SNP, psychencode$SNP))

psychencode.res <- read.csv('~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_PSYCHENCODE.alloutcomes.p005.csv')
cortex.res <- read.csv('~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExCortex.alloutcomes.p005.csv')
length(intersect(cortex.res$gene_name, psychencode.res$gene))

blood <- read.csv('~/Desktop/files/data/eqtl/processed_GTExWholeBlood_eQTL.clumped.csv')
blood.res <- read.csv('~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExWholeBlood.alloutcomes.p005.csv')
length(intersect(blood$SNP, psychencode$SNP))
length(intersect(blood.res$SNP, psychencode.res$SNP))
dim(blood)

hypot <- read.csv('~/Desktop/files/data/eqtl/processed_GTExHypothalamus_eQTL.clumped.csv')
hypot.res <- read.csv('~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExHypothalamus.alloutcomes.p005.csv')
hypot.bonf <- read.csv('~/Desktop/files/current-projects/eQTL-MR/results/mr_wr_GTExHypothalamus.alloutcomes.Bonf.csv')

pnorm(abs(-0.1176287)/0.05006685, lower.tail = FALSE) * 2

mr_wald_ratio
