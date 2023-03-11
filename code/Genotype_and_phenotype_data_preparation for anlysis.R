#########################################
# preparing the genotypic and phenotypic data for the spatial analysis
####################################
# gen = which(rownames(snps_dosage) %in%pyt_functional_1$Genotype ) #
# Subset the pheno data that has snp data


#################################
# prepare the genotypic data - read from the vcf
################################
library(vcfR)
snps = read.vcfR(file = "Imputed_data/Imputed/ABBmid_2014_Blacksburg_Warsaw_immputed.vcf")
head(snps[,1:20])

snps_num <- vcfR::extract.gt(snps,
                             element = "GT",
                             IDtoRowNames  = T, # to list the row name or not
                             as.numeric = T,
                             convertNA = T
)


snps_num_t <- t(snps_num)
snps_num_df <- data.frame(snps_num_t)
head(snps_num_df[,1:20])

###################################
# Read the phenotype data
##################################
pheno = read.csv("Phenotype_ABBmid_2014_Blacksburg_Warsaw.csv")
head(pheno)

#####################################################
# select the genotypes with SNPs and phenotype data
#######################################################

# identify the genotypes in the phenotype data found in the SNP
gt = which(pheno$germplasmName %in% rownames(snps_num_t) ) # identify the row number/genotype name those gentypes from pheno to snps

# Subset the phenotype data based on the genotypes with the SNPs data
phen1 = droplevels(pheno[gt,]) # subset the pheno data that found in the snp data

head(phen1)
dim(phen1)
length(rownames(snps_num_t)) # length row name of the snps
length(levels(as.factor(phen1$germplasmName))) # the length of the selected pheno data genotype levels

# if the legnths are equal no need subseting the snp data
length(rownames(snps_num_t)) == length(levels(as.factor(phen1$germplasmName)))
# Since all the genotype in the SNPs data found in the phenotype data
# but if the length is different we should subset the snp data in relation to the selected pheno genotypes

gsnp = which(rownames(snps_dosage) %in% levels(phen1$Genotype))
head(snps_dosage[,1:20])
snp2 = snps_dosage[rownames(snps_dosage)[gsnp],] # subset the snps based on teh genotypes with snp and also phenotype
dim(snp2)
length(levels(phen1$Genotype))
rownames(snp2) %in% levels(phen1$Genotype)
A = A.mat(snp2)
dim(A)

mix1 <- mmer(WAC~ year+location,
             random=~vsr(Genotype, Gu=A) ,
             rcov=~units,
             data=phen1)
summary(mix1)$varcomp
#### calculate heritability
vpredict(mix1, h1 ~ V1/(V1+V3) )

