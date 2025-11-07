setwd("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/HOS6362_Molc_Mark/Linkage Lounge Top Secret Work/Project_2/filtered_datasets_by_genotype_and_chromosome")

install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

library(GAPIT3)
library(data.table)
#install.packages("irlba")
library(irlba)
library(ggplot2)
#install.packages("ggrepel")  # only once
library(ggrepel)

pheno <- read.table("phenotype_filtered.tsv", head = TRUE)
geno <-  fread("geno_filtered.tsv", sep = "\t", header = TRUE, data.table = FALSE)
snp_map <-  read.table("snp_map_filtered.tsv", head = TRUE)
covariates <- read.table("covariates_filtered.tsv", head = TRUE)

#farm CPU and MLM (Q+k)
#GWAS
GAPIT_myFarmCPU = GAPIT(
  Y=pheno[,c(1,2,3,4,5)], #fist column is ID, 4 phenotypes
  GD=geno,
  GM=snp_map,
  PCA.total=3,
  model=c("FarmCPU", "MLM"),
  Multiple_analysis=TRUE)
