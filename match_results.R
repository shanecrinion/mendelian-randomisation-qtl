## function
match_coloc_results <- function(results_files, data){
  for(i in results_files){
    info <- unlist(strsplit(i, '_'))
    #print(unlist(strsplit(info[1], '/')))
    snp_i <- unlist(strsplit(info[1], '/'))[2]; print(snp_i)
    entrez_i <- info[2] ; print(entrez_i)
    outcome_i <- info[5]; print(outcome_i)
    #print(i)
    snp_result <- subset(read.csv(paste0('~/Desktop/files/coloc/', i)), snp==snp_i) [['SNP.PP.H4']]
    print(snp_result)
    sum_result <- read.csv(stringr::str_replace(paste0('~/Desktop/files/coloc/', i), "results", "summary"), row.names = 1)[6,]
    print(str(data))
    data[data$rsid==snp_i & data$entrez==entrez_i  & data$outcome==outcome_i,][['snp_pp4']] <- snp_result
    data[data$rsid==snp_i & data$entrez==entrez_i & data$outcome==outcome_i,][['sum_pp4']] <- sum_result
  }
  message("need to manually edit the following rows:")
  results=list(coloc=data, editing_needed=distinct(data[match(duplicated(data$rsid), duplicated(data$entrez)),]))
  return(results)}

hypothalamus_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(hypothalamus_results$coloc, '~/Desktop/files/coloc/results_GTExHypothalamusv8_eQTL_GWAS.Bonf.withcoloc.csv')



# do this by dataset
# gtex cortex
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (cortex)")

data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv') # make recursive
results_files <- results_files[grepl("^rs.*_cortex_", results_files)]
cortex_results <- match_coloc_results(results_files = results_files, data = data)

write.csv(cortex_results$coloc, '~/Desktop/files/coloc/results_GTExCortexv8_eQTL_GWAS.Bonf.withcoloc.csv')

# do this by dataset
# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hippocampus)")
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv')
results_files <- results_files[grepl("hippocampus", results_files)]
hippocampus_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(hippocampus_results$coloc, '~/Desktop/files/coloc/results_GTExHippocampusv8_eQTL_GWAS.Bonf.withcoloc.csv')


# do this by dataset
# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (whole blood)")
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv', recursive = T)
results_files <- results_files[grepl("blood", results_files)]
results_files
blood_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(blood_results$coloc, '~/Desktop/files/coloc/results_GTExWholeBloodv8_eQTL_GWAS.Bonf.withcoloc.csv')


# do this by dataset
# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "PsychENCODE (prefrontal cortex)")

data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv', recursive = T)
results_files <- results_files[grepl("prefrontal", results_files)]
results_files
prefrontalcortex_results <- match_coloc_results(results_files = results_files, data = data)
View(prefrontalcortex_results$coloc)
write.csv(prefrontalcortex_results$coloc, '~/Desktop/files/coloc/results_PSYCHENCODE_eQTL_GWAS.Bonf.withcoloc.csv')


# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (cortex)")
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv', recursive = T)
results_files <- results_files[grepl("_cortex_", results_files)]
results_files
cortex_results <- match_coloc_results(results_files = results_files, data = data)
View(cortex_results$coloc)
write.csv(cortex_results$coloc, '~/Desktop/files/coloc/results_GTExCortexv8_eQTL_GWAS.Bonf.withcoloc.csv')


# do this by dataset
# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hypothalamus)")

data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'beta_change.results.csv',recursive = T)
results_files <- results_files[grepl("_hypothalamus_", results_files)]
hypothalamus_results <- match_coloc_results(results_files = results_files, data = data)

write.csv(hypothalamus_results$coloc, '~/Desktop/files/coloc/results_GTExHypothalamusv8_eQTL_GWAS.Bonf.flip.withcoloc.csv')

View(read.csv('~/Desktop/files/coloc/results_GTExHypothalamusv8_eQTL_GWAS.Bonf.flip.withcoloc.csv'))
View(read.csv('~/Desktop/files/coloc/results_GTExHypothalamusv8_eQTL_GWAS.Bonf.withcoloc.csv'))


# do this by dataset
# gtex hippocampus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hippocampus)")

data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv',recursive = T)
results_files <- results_files[grepl("_hippocampus_", results_files)]



data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",
  sheet = "GTEx (whole blood)")
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))
data$snp_pp4 <- 0
data$sum_pp4 <- 0
results_files <- list.files('~/Desktop/files', pattern = 'results.csv')
results_files <- results_files[grepl("^rs.*_blood_", results_files)]
blood_results <-  match_coloc_results(results_files = results_files, data = data)
hippocampus_results <- match_coloc_results(results_files = results_files, data = data)
View(hippocampus_results$coloc)
write.csv(hippocampus_results$coloc, '~/Desktop/files/coloc/results_GTExHippocampusv8_eQTL_GWAS.Bonf.withcoloc.csv')
write.csv(blood_results$coloc, '~/Desktop/files/coloc/results_GTEx_wholeblood_GWAS.Bonf.withcoloc.csv')

