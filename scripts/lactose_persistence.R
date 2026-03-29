#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries
library(vcfR)
library(dplyr)
library(ggplot2)

#Need to merge VCF (Variant call format) data with metadata, want to align sampleIDs

#Read in VCF file
#Should output processed variant: 6776
vcf <- read.vcfR("../data/lactose_persistence.vcf", verbose = TRUE)

#Read metadata
metadata <- read.table("../data/20130606_g1k_3202_samples_ped_population.txt", header = TRUE, sep = " ") #Columns separated by spaces
names(metadata)

#Extract genotypes and sample IDs
#Extract genotype matrix (Converts genotypes into number of alternate alleles)
genotypes <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
nrow(genotypes) #6776, each row is a variant

#Transpose so samples are rows and SNPs are columns
genotypes_t <- t(genotypes)
ncol(genotypes_t) #6776

#Convert data to a dataframe to make sampleIDs a column
geno_df <- as.data.frame(genotypes_t) #nrow  = 
geno_df$SampleID <- rownames(geno_df)
geno_df <- geno_df %>%
  filter(SampleID %in% metadata$SampleID) %>% 

#Merge metadata with genotypes
geno_data <- geno_df %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))
nrow(geno_data) #3202 (3202 samples in 30X 1000 Genomes project)
ncol(geno_data) #6783 (6776 SampleIDs + SampleID + Sex + FamilyID + Population + Superpopulation + MotherID + FatherID )

#Check populations
print(table(geno_data$Population)) #Three letter codes for locations
print(table(geno_data$Superpopulation)) #Africa, America, East Asia, Europe, South Asia


#Write out CSV with variants and data
#write.csv(geno_data, "../data/genotype_data.csv")

# EDA ====
#Exploratory data analysis to identify which subpopulations to extract before doing clustering by genotypes

#Separate genotype columns from the metadata columns
genotype_columns <- geno_data %>% #ncol = 6776 variants
  select(-SampleID, -FamilyID, -FatherID, -MotherID, -Sex, -Population, -Superpopulation, )

#Calculate the sum of alleles for each SNP
# A sum of 0 means the variant (minor allele) is not present in the 3,202 samples
snp_sums <- colSums(genotype_columns, na.rm = TRUE)

#Identify SNPs with at least one variant (sum > 0)
#Remove SNPs where everyone has the variant (sum == 2 * number of samples (3202))
#Remove SNPs that are rare in the population (<0.001)
keep_snps <- snp_sums > 0 & snp_sums < 2*nrow(geno_data) & snp_sums > 0.001*nrow(geno_data)

#Subset your data to keep only informative SNPs
geno_data_filt <- geno_data[, c(names(which(keep_snps)), "SampleID", "Population", "Superpopulation", "Sex")]

# Check how many SNPs you have left
print(paste("Original SNPs:", ncol(genotype_columns)))
print(paste("Informative SNPs remaining:", sum(keep_snps)))

