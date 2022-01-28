###########################################################################
## 
## File Name: 'GS_1st_stage'
##
## Authors:    SB Fernandes < samuelfernandes@agronomo.eng.br >
##             KOG Dias < kaio.o.dias@ufv.br >  
##
## Date:      Jan 28th, 2022
##
## Contents:  Asreml models for obtaining adjusted means by location
##            
## input:  "pheno.csv"  
##
## output: "means.csv"
##
############################################################################
#setwd("...")
library(asreml)
library(reshape2)

#### loading phenotypes and subseting it by location #####
pheno <- read.csv("pheno.csv", stringsAsFactors = F)
pheno <- transform(pheno, LOC = factor(LOC),  GENO = factor(GENO),	range = factor(range),	row = factor(row),	BLOCK = factor(BLOCK))
pheno <- pheno[order(pheno$LOC, pheno$range, pheno$row, pheno$BLOCK, pheno$GENO),]
pheno12 <- droplevels(subset(pheno, LOC == "12EF"))
pheno13 <- droplevels(subset(pheno, LOC == "13EF"))
pheno14 <- droplevels(subset(pheno, LOC == "14EF"))
rm(pheno)

###### Phenotypic analysis by trait and by location #####
###### Biomass yield ######
Y_12 <- asreml(fixed = Y ~ GENO,
             random = ~ BLOCK,
             residual = ~ ar1(range):ar1(row),
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100, data = pheno12)

Y_13 <- asreml(fixed = Y ~ GENO,
             random = ~ BLOCK,
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100, data = pheno13)

Y_14 <- asreml(fixed = Y ~ GENO,
             random = ~ BLOCK,
             residual = ~ ar1(range):ar1(row),
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100, 
             data= pheno14)

###### Moisture  ######
M_12 <- asreml(fixed = M ~ GENO,
             random = ~ BLOCK,
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100,  data = pheno12)
M_12 <- update.asreml(M_12)

M_13 <- asreml(fixed = M ~ GENO,
             random = ~ BLOCK,
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100, data = pheno13)
# Update
new_model <- M_13
for(k in 1:2){
  temp_model <- update(new_model)
  new_model <- temp_model
}
M_13 <- new_model
rm(temp_model, new_model, k)

M_14 <- asreml(fixed = M ~ GENO,
             random = ~ BLOCK,
             na.action = na.method("include"),
             workspace = 64e6, maxit = 100, data = pheno14)
M_14 <- update(M_14)


###### Height 1 (30 DAP) ######
h1_12 <- asreml(fixed = h1 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno12)
h1_12 <- update(h1_12)

h1_13 <- asreml(fixed = h1 ~ GENO,
              random = ~ BLOCK,
              residual = ~ ar1(range):ar1(row),
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno13)

h1_14 <- asreml(fixed = h1 ~ GENO,
              random = ~ BLOCK,
              residual = ~ ar1(range):ar1(row),
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno14)
h1_14 <- update.asreml(h1_14)

###### Height 2 (60 DAP) ######
h2_12 <- asreml(fixed = h2 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno12)
h2_12 <- update.asreml(h2_12)

h2_13 <- asreml(fixed = h2 ~ GENO,
              random = ~ BLOCK,              
              residual = ~ ar1(range):ar1(row),
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno13)
# Update
new_model <- h2_13
for(k in 1:3){
  temp_model <- update(new_model)
  new_model <- temp_model
}
h2_13 <- new_model
rm(temp_model, new_model, k)

h2_14 <- asreml(fixed = h2 ~ GENO,
              random = ~ BLOCK,
              residual = ~ ar1(range):ar1(row),
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno14)
# Update
new_model <- h2_14
for(k in 1:2){
  temp_model <- update(new_model)
  new_model <- temp_model
}
h2_14 <- new_model
rm(temp_model, new_model, k)

###### Height 3 (90 DAP) ######
h3_12 <- asreml(fixed = h3 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno12)

h3_13 <- asreml(fixed = h3 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno13)

h3_14 <- asreml(fixed = h3 ~ GENO,
              random = ~ BLOCK,
              residual = ~ ar1(range):ar1(row),
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno14)
# Update
new_model <- h3_14
for(k in 1:6){
  temp_model <- update(new_model)
  new_model <- temp_model
}
h3_14 <- new_model
rm(temp_model, new_model, k)

###### Height 4 (120 DAP) ######
h4_12 <- asreml(fixed = h4 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno12)

h4_13 <- asreml(fixed = h4 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno13)

h4_14 <- asreml(fixed = h4 ~ GENO,
              random = ~ BLOCK,
              na.action = na.method("include"),
              workspace = 64e6, maxit = 100, data = pheno14)


####### AIC comparison #######
# KM1 = length(model1$gammas)
# lM1 = model1$logl
# KM2 = length(model2$gammas)
# lM2 = model2$logl
# (AIC1 = -2*lM1 + 2*KM1 ) #aic for model 1
# (AIC2 = -2*lM2 + 2*KM2 ) #aic for model 2
############################

####### Using asremlPlus #######
library(asremlPlus)
infoCriteria.asreml(Y_12)

###### predicted values ######
vars <- ls()
vars <- vars[grepl("._", vars)]
results <- data.frame(trait = character(), LOC = character(), 
                      GENO = character(), predicted = numeric())
#Adjusting means for all models at once
for(i in 1:length(vars)){
  model <- eval(parse(text = vars[i]))
  x <- predict(model, classify = "GENO",
               sed = TRUE, pworkspace = 128e6, maxit = 100)$pvals[,c(1,2)]
  results <-
    rbind(results, data.frame(
      trait = gsub("_..", "", vars[i]),
      LOC = gsub(paste(gsub("_..", "", vars[i]), "_", sep = ""), "", vars[i]),
      x
    ))
  rm(model)
}
colnames(results) <- c("trait", "LOC",  "GENO",   "pred")
#data in a wide format table
means <- dcast(results, LOC + GENO ~ trait, value.var = "pred")
#saving adjusted means for 2nd stage analysis
#write.csv(means, "means.csv")
