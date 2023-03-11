#########################################################
# spatial analysis with sommer including Genomic data
#########################################################
#Import the phenotypic and genotypic data

phen = read.csv(file = "phenotype_YLDqtVal2014_2015.csv") # phenotypic data
head(phen[1:10,1:10])
dim(phen)
snp = read.csv("Genotype_YLdqt_2014_15Edited.csv") # read snp data
head(snp[1:10, 1:10])
rownames(snp) <- snp$X
snp = snp[,-1]
# factor variables
phen$location = as.factor(phen$locationName)
phen$year = as.factor(phen$studyYear)
phen$geno = as.factor(phen$germplasmName)
phen$rep = as.factor(phen$replicate)
phen$rowf = as.factor(phen$rowNumber)
phen$colf = as.factor(phen$colNumber)
phen$col = phen$colNumber
phen$row = phen$rowNumber
phen$env = phen$location:phen$year
phen$env = gsub(pattern = ":", replacement = "_", phen$env)
phen$env = as.factor(phen$env)
table(phen$env, phen$rep)
# Re-naming the traits

colnames(phen)[34] = "test_weight"
colnames(phen)[35] = "Yield"
colnames(phen)[33] = "Grain_protein"
str(phen)
phen$Yield = as.numeric(phen$Yield)
phen$test_weight = as.numeric(phen$test_weight)
phen$Grain_protein = as.numeric(phen$Grain_protein)
str(phen) # to check if the change is made
#################################################
#### estimating the addative relationship matrix
##################################################
library(rrBLUP)
library(sommer)
dim(snp)
head(snp[1:10,1:20])
A = rrBLUP::A.mat(snp)
dim(A)

####################################################
# select the genotypes with snp and phenotype data
####################################################
##Subset the phenotyp data set with genotype (snp)

library(dplyr)
gt = which(rownames(A) %in% levels(phen$geno)) # the list of the rows of pheno data that has the snp
phen1 = droplevels(phen[gn,]) # the phenotype data substed
gn = rownames(A) %in% levels(phen1$geno) # Checking it the geno name at Amatrix and levels of genotype in the phenotypic data


#################################
## Heritability and variance component estimation using sommer package
####################################################################
## 1. Signle location with marker
####################################
#Heritability
Result_sommer_Marker.h2 <- matrix(nrow = ne, ncol = nt) # the matrix to put the h2
colnames(Result_sommer_Marker.h2) <- Traits # naming the columns
rownames(Result_sommer_Marker.h2) <- Env # naming the row of the matrix

# # Define the list to put the BLUP mean value for each environment
#
# Result.Mean_sommer_sp_BLUP <- vector(mode = "list", length = length(Env)) #length(levels(met.fbn$env)))
# names(Result.Mean_sommer_sp_BLUP) <- Env

# Package used for analysis
library(sommer)

system.time(for (i in 1:ne){
  SL <- droplevels(subset(x = phen1, subset = env == Env[i]) )# Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait without NA
  ntt = length((TraitN))
  gen = length(levels(SL$geno))
  MAT <- matrix(nrow = gen, ncol = ntt) # create matrix based on the number of traits
  colnames(MAT)<- TraitN
  rownames(MAT) <- sort(levels(SL$geno))

  for (Trait in TraitN){
    eval(parse(text = paste("mixspM = mmer(",Trait," ~ 1,
                            random = ~ vsr(geno,Gu = A),
                            rcov = ~units,
                            data = SL,
                            tolParConv = 1e-6, verbose = FALSE)")))

    svsp1 = summary(mixspM)$varcomp
    rownames(svsp1) <- c("geno", "Residual")
    v.e1 <-svsp1["Residual", "VarComp"]
    v.g1 <-svsp1["geno", "VarComp"]
    Result_sommer_Marker.h2[i,Trait]<- round(v.g1/(v.g1 + v.e1/r),3)
    #########################################
    # # Predict mean using fitted value
    # #########################################
    #
    # fit = fitted(mixsp)
    # head(fitval)
    # fitval = fit$dataWithFitted # the data with fitted value
    #
    # fitval$AdjMean = fitval[,paste0(Trait,".fitted")]- fitval[,"C.fitted"] -
    #   fitval[,"R.fitted"] - fitval[,"rep.fitted"]- fitval[,"A:all.fitted"]
    #
    # adjMean = fitval %>% group_by(geno) %>% summarise_at(.vars = "AdjMean",
    #                                                      .funs = mean)
    # adjMean = as.data.frame(adjMean)
    # rownames(adjMean) = adjMean$geno

    ####################################################
    # To be used in case- it required more processing time
    #######################################################
    # prsp = predict(mixsp,classify = "geno")
    # pmsp = prsp$pvals[,1:3]
    # rownames(pmsp) = pmsp$geno

    ##Putting the mean data of the genotyeps in the matrix
    # for(name in seq(levels(SL$geno))){
    #   MAT[name,Trait] = adjMean[name,2]
    # }
    #
  }

  # Result.Mean_sommer_sp_BLUP[[i]] = MAT
}
)




##Note to estimate heritability with standard error
#vpredict(mix1, h1 ~ V1/(V1+V3/2) ) # heritability with stanard error

#####################################################
# Single location with marker and spatial analysis\
######################################################
#Heritability
Traits = c("test_weight","Grain_protein", "Yield")

phen1$env = phen1$location:phen1$year
phen1$env = gsub(pattern = ":", replacement = "_",x = phen1$env)
phen1$env = as.factor(phen1$env)
Env = levels(phen1$env)
ne = length(levels(phen1$env))
nt = length(Traits)
r = length(levels(phen1$rep))

Result_sommer_marker_sp.h2 <- matrix(nrow = ne, ncol = nt) # the matrix to put the h2
colnames(Result_sommer_marker_sp.h2) <- Traits # naming the columns
rownames(Result_sommer_marker_sp.h2) <- Env # naming the row of the matrix

# Define the list to put the BLUP mean value for each environment

# Result.Mean_sommer_sp_BLUP <- vector(mode = "list", length = length(Env)) #length(levels(met.fbn$env)))
# names(Result.Mean_sommer_sp_BLUP) <- Env

# Package used for analysis
library(sommer)

system.time(for (i in 1:ne){
  SL <- droplevels(subset(x = phen1, subset = env == Env[i]) )# Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait without NA
  ntt = length((TraitN))
  gen = length(levels(SL$geno))
  MAT <- matrix(nrow = gen, ncol = ntt) # create matrix based on the number of traits
  colnames(MAT)<- TraitN
  rownames(MAT) <- sort(levels(SL$geno))
  table(phen1$geno,phen1$rep)
  for (Trait in TraitN){
    eval(parse(text = paste("mixsp = mmer(",Trait," ~ 1,
                            random = ~ vsr(geno,Gu = A) + rowf+ colf + spl2Da(col,row),
                            rcov = ~units,
                            data = SL,
                            tolParConv = 1e-6, verbose = FALSE)")))

    svsp = summary(mixsp)$varcomp
    rownames(svsp) <- c("geno", "row", "col", "A:all", "Residual")
    v.e1 <-svsp["Residual", "VarComp"]
    v.g1 <-svsp["geno", "VarComp"]
    Result_sommer_marker_sp.h2[i,Trait]<- round(v.g1/(v.g1 + v.e1/r),3)

    #########################################
    # Predict mean using fitted value
    #########################################
    #
    # fit = fitted(mixsp)
    # head(fitval)
    # fitval = fit$dataWithFitted # the data with fitted value
    #
    # fitval$AdjMean = fitval[,paste0(Trait,".fitted")]- fitval[,"C.fitted"] -
    #   fitval[,"R.fitted"] - fitval[,"rep.fitted"]- fitval[,"A:all.fitted"]
    #
    # adjMean = fitval %>% group_by(geno) %>% summarise_at(.vars = "AdjMean",
    #                                                      .funs = mean)
    # adjMean = as.data.frame(adjMean)
    # rownames(adjMean) = adjMean$geno
    #
    # ####################################################
    # # To be used in case- it required more processing time
    # #######################################################
    # # prsp = predict(mixsp,classify = "geno")
    # # pmsp = prsp$pvals[,1:3]
    # # rownames(pmsp) = pmsp$geno
    #
    # ##Putting the mean data of the genotyeps in the matrix
    # for(name in seq(levels(SL$geno))){
    #   MAT[name,Trait] = adjMean[name,2]
    # }

  }

  # Result.Mean_sommer_sp_BLUP[[i]] = MAT
}
)





