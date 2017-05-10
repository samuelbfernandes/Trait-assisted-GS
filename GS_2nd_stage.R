###########################################################################
## 
## File Name: 'GS_2nd_stage'
##
## Authors:    SB Fernandes < samuelfernandes@agronomo.eng.br >
##             KOG Dias < kaioolimpio@hotmail.com >  
## Date:      May 10th, 2017
##
## Contents:  Asreml-r codes for obtaining cross-validation results 
##            for single trait and multi-trait GBLUP
##
## input:  "snps.csv", "means.csv"  
##
## output: "ginv.csv", "means2.csv", "accuracies.csv", "gebv.csv"(one for each model)
##
## source: "CI.R"
############################################################################
#setwd("")
library(asreml)
require(cvTools)
require(plyr)
library(rrBLUP)
library(MASS)
library(reshape)
snps<-read.csv("snps.csv",row.names = 1)
means<-read.csv("means.csv",row.names = 1)
#using only phenotypes with snp information
means<-droplevels(means[means$GENO%in%rownames(snps),])
#changing names for numbers
names<-as.numeric(as.factor(rownames(snps)))
means$geno<-NA
for(i in 1:453){ 
  means$geno[which(rownames(snps)[i]==means$GENO)]<-names[i]
}
means<-means[,c("geno","LOC","h1", "h2", "h3", "h4", "M", "Y")]
means<-transform(means, LOC=factor(LOC),geno=factor(geno))


#### Creating a Relationship Matrix from SNPs information #####

kmatrix=A.mat(snps,return.imputed=F)
rm(snps,i)
# kmatrix<-as.matrix(read.table("kmatrix.txt", h=T)) ##### if already created  
#inverting relationship matrix
A<-ginv(kmatrix)
colnames(A)<-names
rownames(A)<-names
rm(names,kmatrix)

#changing inverted A matrix format to use in asreml
A[lower.tri(A)]<-NA
A<-na.omit(melt(A))
rownames(A)<-NULL
#ginv <- read.csv("ginv.csv")##### if saved 
ginv<-data.frame(A[,2],A[,1],A[,3])
colnames(ginv)<- c("Row", "Column", "GINV")
attr(ginv,"rowNames")<-1:453
rm(A)
#write.csv(ginv, "ginv.csv", row.names = FALSE)

###### Means over location ##########
means2<-ddply(means,.(geno), summarize, Y=mean(Y,na.rm=T), M=mean(M,na.rm=T),  h1=mean(h1,na.rm=T), h2=mean(h2,na.rm=T), h3=mean(h3,na.rm=T), h4=mean(h4,na.rm=T)  )
##### Area Under Growth Progress Curve #####
means2$A<-(((means2$h1+means2$h2)/2)*(30)+((means2$h2+means2$h3)/2)*(30)+((means2$h3+means2$h4)/2)*(30))
#write.csv(means2, "means2.csv")

#### sorts for cross-validation ####
sort<-list()
for(a in 1:30){
  for(j in 1:5){
    folds <- cvFolds(nlevels(means$geno),type ="random", K=5, R = 1)
    Sample<-cbind(folds$which,folds$subsets)
    cv<-split(Sample[,2], f=Sample[,1])
  }
  sort[[a]]<-cv  
}
rm(a, folds, j, cv, Sample)

################################
##### GS Cross-validation ######
################################

##### Standard  GS Y ####
acy<-c() ## store accuracies
gebvy<-1:453  ## store GEBV's
for(j in 1:length(sort)){
  ry<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    model <- asreml(Y ~1, random = ~giv(geno), 
                    ginverse=list(geno=ginv), 
                    maxiter=100,data = test) 
    x<-as.matrix(coef(model, pattern = 'geno')[1:453,])
    ry[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ray<-rbind(ry[[1]],ry[[2]], ry[[3]], ry[[4]], ry[[5]])
  ray<-ray[order(ray[,1]),]
  acy[j]<-cor(ray[,2], means2$Y)
  gebvy<-cbind(gebvy,ray[,2])
}
rm(test, x, i, j, model, ry, ray)

##### Standard  GS M and Indirect GS M ####
acm<-c()
acm_y<-c()
gebvm<-1:453
for(j in 1:30){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"M"]<-NA
    
    
    modelm <- asreml(M ~1, random = ~giv(geno), 
                     ginverse=list(geno=ginv), 
                     maxiter=100,data = test) 
    x<-as.matrix(coef(modelm, pattern = 'geno'))
    rm[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acm[j]<-cor(ram[,2], means2$M)
  acm_y[j]<-cor(-ram[,2], means2$Y)
  gebvm<-cbind(gebvm,ram[,2])
}
rm(test, x, i, j, modelm, ram, rm)

##### Standard  GS H1  and Indirect GS H1 ####
ach1<-c()
ach1_y<-c()
gebvh1<-1:453
for(j in 1:30){
  r1<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h1"]<-NA
    
    
    modelh1 <- asreml(h1 ~1, random = ~giv(geno), 
                      ginverse=list(geno=ginv), 
                      maxiter=100,data = test) 
    x<-as.matrix(coef(modelh1, pattern = 'geno'))
    r1[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra1<-rbind(r1[[1]],r1[[2]], r1[[3]], r1[[4]], r1[[5]])
  ra1<-ra1[order(ra1[,1]),]
  ach1[j]<-cor(ra1[,2], means2$h1)
  ach1_y[j]<-cor(ra1[,2], means2$Y)
  gebvh1<-cbind(gebvh1,ra1[,2])
}
rm(test, x, i, j, modelh1, ra1, r1)

##### Standard  GS H2  and Indirect GS H2 ####
ach2<-c()
ach2_y<-c()
gebvh2<-1:453
for(j in 1:30){
  r2<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h2"]<-NA
    
    
    modelh2 <- asreml(h2 ~1, random = ~giv(geno), 
                      ginverse=list(geno=ginv), 
                      maxiter=100,data = test) 
    x<-as.matrix(coef(modelh2, pattern = 'geno'))
    r2[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra2<-rbind(r2[[1]],r2[[2]], r2[[3]], r2[[4]], r2[[5]])
  ra2<-ra2[order(ra2[,1]),]
  ach2[j]<-cor(ra2[,2], means2$h2)
  ach2_y[j]<-cor(ra2[,2], means2$Y)
  gebvh2<-cbind(gebvh2,ra2[,2])
}
rm(test, x, i, j, modelh2, ra2, r2)

##### Standard  GS H3  and Indirect GS H3 #####
ach3<-c()
ach3_y<-c()
gebvh3<-1:453
for(j in 1:30){
  r3<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h3"]<-NA
    
    modelh3 <- asreml(h3 ~1, random = ~giv(geno), 
                      ginverse=list(geno=ginv), 
                      maxiter=100,data = test) 
    x<-as.matrix(coef(modelh3, pattern = 'geno'))
    r3[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra3<-rbind(r3[[1]],r3[[2]], r3[[3]], r3[[4]], r3[[5]])
  ra3<-ra3[order(ra3[,1]),]
  ach3[j]<-cor(ra3[,2], means2$h3)
  ach3_y[j]<-cor(ra3[,2], means2$Y)
  gebvh3<-cbind(gebvh3,ra3[,2])  
}
rm(test, x, i, j, modelh3, ra3, r3)

##### Standard  GS H4  and Indirect GS H4#####
ach4<-c()
ach4_y<-c()
gebvh4<-1:453
for(j in 1:30){
  r4<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h4"]<-NA
    
    modelh4 <- asreml(h4 ~1, random = ~giv(geno), 
                      ginverse=list(geno=ginv), 
                      maxiter=100,data = test) 
    x<-as.matrix(coef(modelh4, pattern = 'geno'))
    r4[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra4<-rbind(r4[[1]],r4[[2]], r4[[3]], r4[[4]], r4[[5]])
  ra4<-ra4[order(ra4[,1]),]
  ach4[j]<-cor(ra4[,2], means2$h4)
  ach4_y[j]<-cor(ra4[,2], means2$Y)
  gebvh4<-cbind(gebvh4,ra4[,2])   
}
rm(test, x, i, j, modelh4, ra4, r4)

##### Standard  GS A  and Indirect GS A ####
aca<-c()
aca_y<-c()
gebva<-1:453
for(j in 1:30){
  ra<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"A"]<-NA
    
    modela <- asreml(A ~1, random = ~giv(geno), 
                     ginverse=list(geno=ginv), 
                     maxiter=100,data = test) 
    x<-as.matrix(coef(modela, pattern = 'geno')[1:453,])
    ra[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  raa<-rbind(ra[[1]],ra[[2]], ra[[3]], ra[[4]], ra[[5]])
  raa<-raa[order(raa[,1]),]
  aca[j]<-cor(raa[,2], means2$A)
  aca_y[j]<-cor(raa[,2], means2$Y)
  gebva<-cbind(gebva,raa[,2])     
}
rm(test, x, i, j, modela, raa, ra)

################################
#####    Multitrait GS    ######
################################

### Multitrait GS YM ####
acymm<-c()
gebvymm<-1:453
for(j in 1:30){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"M"]<-NA
    
    multitraitm<-asreml(fixed=cbind(Y,M)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitm, pattern = 'geno')[1:453,])
    rm[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acymm[j]<-cor(ram[,2], means2$Y)
  gebvymm<-cbind(gebvymm,ram[,2])
}
rm(test, x, i, j, multitraitm, rm, ram)


### Multitrait GS YH1 ####
acyh1m<-c()
gebvyh1m<-1:453
for(j in 1:30){
  r1<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h1"]<-NA
    
    multitrait1<-asreml(fixed=cbind(Y,h1)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait1, pattern = 'geno')[1:453,])
    r1[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra1<-rbind(r1[[1]],r1[[2]], r1[[3]], r1[[4]], r1[[5]])
  ra1<-ra1[order(ra1[,1]),]
  acyh1m[j]<-cor(ra1[,2], means2$Y)
  gebvyh1m<-cbind(gebvyh1m,ra1[,2])
}
rm(test, x, i, j, multitrait1, r1, ra1)


### Multitrait GS YH2 ####
acyh2m<-c()
gebvyh2m<-1:453
for(j in 1:30){
  r2<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h2"]<-NA
    
    multitrait2<-asreml(fixed=cbind(Y,h2)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait2, pattern = 'geno')[1:453,])
    r2[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra2<-rbind(r2[[1]],r2[[2]], r2[[3]], r2[[4]], r2[[5]])
  ra2<-ra2[order(ra2[,1]),]
  acyh2m[j]<-cor(ra2[,2], means2$Y)
  gebvyh2m<-cbind(gebvyh2m,ra2[,2])
}
rm(test, x, i, j, multitrait2, r2, ra2)


### Multitrait GS YH3 ####
acyh3m<-c()
gebvyh3m<-1:453
for(j in 1:30){
  r3<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h3"]<-NA
    
    multitrait3<-asreml(fixed=cbind(Y,h3)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait3, pattern = 'geno')[1:453,])
    r3[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra3<-rbind(r3[[1]],r3[[2]], r3[[3]], r3[[4]], r3[[5]])
  ra3<-ra3[order(ra3[,1]),]
  acyh3m[j]<-cor(ra3[,2], means2$Y)
  gebvyh3m<-cbind(gebvyh3m,ra3[,2])
}
rm(test, x, i, j, multitrait3, r3, ra3)


### Multitrait GS YH4 ####
acyh4m<-c()
gebvyh4m<-1:453
for(j in 1:30){
  r4<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h4"]<-NA
    
    multitrait4<-asreml(fixed=cbind(Y,h4)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait4, pattern = 'geno')[1:453,])
    r4[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra4<-rbind(r4[[1]],r4[[2]], r4[[3]], r4[[4]], r4[[5]])
  ra4<-ra4[order(ra4[,1]),]
  acyh4m[j]<-cor(ra4[,2], means2$Y)
  gebvyh4m<-cbind(gebvyh4m,ra4[,2])
}
rm(test, x, i, j, multitrait4, r4, ra4)


### Multitrait GS YA ####
acyam<-c()
gebvyam<-1:453
for(j in 1:30){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"A"]<-NA
    
    multitraitm<-asreml(fixed=cbind(Y, A)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitm, pattern = 'geno')[1:453,])
    rm[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acyam[j]<-cor(ram[,2], means2$Y)
  gebvyam<-cbind(gebvyam,ram[,2])
}
rm(test, x, i, j, multitraita, ra, ram)

################################
#####   Trait-assisted GS ######
################################

#### Trait-assisted GS YM #####
acym<-c()
gebvym<-1:453
for(j in 1:length(sort)){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    multitraitm<-asreml(fixed=cbind(Y,M)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitm, pattern = 'geno')[1:453,])
    rm[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acym[j]<-cor(ram[,2], means2$Y)
  gebvym<-cbind(gebvym,ram[,2])
}
rm(test, x, i, j, multitraitm, rm, ram)

#### Trait-assisted GS YH1 #####
acyh1<-c()
gebvyh1<-1:453
for(j in 1:length(sort)){
  r1<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    
    multitrait1<-asreml(fixed=cbind(Y,h1)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait1, pattern = 'geno')[1:453,])
    r1[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra1<-rbind(r1[[1]],r1[[2]], r1[[3]], r1[[4]], r1[[5]])
  ra1<-ra1[order(ra1[,1]),]
  acyh1[j]<-cor(ra1[,2], means2$Y)
  gebvyh1<-cbind(gebvyh1,ra1[,2])
}
rm(test, x, i, j, multitrait1, r1, ra1)

#### Trait-assisted GS YH2 #####
acyh2<-c()
gebvyh2<-1:453
for(j in 1:length(sort)){
  r2<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    
    multitrait2<-asreml(fixed=cbind(Y,h2)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait2, pattern = 'geno')[1:453,])
    r2[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra2<-rbind(r2[[1]],r2[[2]], r2[[3]], r2[[4]], r2[[5]])
  ra2<-ra2[order(ra2[,1]),]
  acyh2[j]<-cor(ra2[,2], means2$Y)
  gebvyh2<-cbind(gebvyh2,ra2[,2])
}
rm(test, x, i, j, multitrait2, r2, ra2)

#### Trait-assisted GS YH3 #####
acyh3<-c()
gebvyh3<-1:453
for(j in 1:length(sort)){
  r3<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    
    multitrait3<-asreml(fixed=cbind(Y,h3)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait3, pattern = 'geno')[1:453,])
    r3[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra3<-rbind(r3[[1]],r3[[2]], r3[[3]], r3[[4]], r3[[5]])
  ra3<-ra3[order(ra3[,1]),]
  acyh3[j]<-cor(ra3[,2], means2$Y)
  gebvyh3<-cbind(gebvyh3,ra3[,2])
}
rm(test, x, i, j, multitrait3, r3, ra3)

#### Trait-assisted GS YH4 #####
acyh4<-c()
gebvyh4<-1:453
for(j in 1:length(sort)){
  r4<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    
    multitrait4<-asreml(fixed=cbind(Y,h4)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitrait4, pattern = 'geno')[1:453,])
    r4[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  ra4<-rbind(r4[[1]],r4[[2]], r4[[3]], r4[[4]], r4[[5]])
  ra4<-ra4[order(ra4[,1]),]
  acyh4[j]<-cor(ra4[,2], means2$Y)
  gebvyh4<-cbind(gebvyh4,ra4[,2])
}
rm(test, x, i, j, multitrait4, r4, ra4)

#### Trait-assisted GS YA #####
acya<-c()
gebvya<-1:453
for(j in 1:length(sort)){
  ra<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"Y"]<-NA
    
    multitraita<-asreml(fixed=cbind(Y, A)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraita, pattern = 'geno')[1:453,])
    ra[[i]]<-cbind(sort[[j]][[i]],x[1:453,1][sort[[j]][[i]]]) 
  }
  
  raa<-rbind(ra[[1]],ra[[2]], ra[[3]], ra[[4]], ra[[5]])
  raa<-raa[order(raa[,1]),]
  acya[j]<-cor(raa[,2], means2$Y)
  gebvya<-cbind(gebvya,raa[,2])
}
rm(test, x, i, j, multitraita, ra, raa)


################################
###  Indirect Multi-trait GS ###
################################

#### Indirect Multi-trait GS MH3 #####
means2$h3<-scale(means2$h3)
means2$M<-scale(means2$M)

gebvmh3<-1:453
acmh3<-c()
for(j in 1:length(sort)){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"M"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h3"]<-NA
    
    multitraitmh3<-asreml(fixed=cbind(M,h3)~trait,                    
                          random=~us(trait):giv(geno),
                          rcov=~units:us(trait),
                          ginverse=list(geno=ginv), 
                          workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitmh3, pattern = 'geno')[1:453,])
    x1<-as.matrix(coef(multitraitmh3, pattern = 'geno')[454:906,])
    x2<-cbind(x, x1)
    
    x2[,1]<-x2[,1]*-0.3925930726
    x2[,2]<-x2[,2]*0.8088956
    ind<-as.matrix(apply(x2, 1, sum))
    rm[[i]]<-cbind(sort[[j]][[i]],ind[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acmh3[j]<-cor(ram[,2], means2$Y)
  gebvmh3<-cbind(gebvmh3,ram[,2]) 
}
rm(test, x,x1,x2, ind, i, j, multitraitmh3, ram, rm)

#### Indirect Multi-trait GS MH4 #####
means2$h4<-scale(means2$h4)
gebvmh4<-1:453
acmh4<-c()
for(j in 1:length(sort)){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"M"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"h4"]<-NA
    
    multitraitmh4<-asreml(fixed=cbind(M,h4)~trait,                    
                          random=~us(trait):giv(geno),
                          rcov=~units:us(trait),
                          ginverse=list(geno=ginv), 
                          workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitmh4, pattern = 'geno')[1:453,])
    x1<-as.matrix(coef(multitraitmh4, pattern = 'geno')[454:906,])
    x2<-cbind(x, x1)
    
    x2[,1]<-x2[,1]*-0.3925930726
    x2[,2]<-x2[,2]*0.8327002
    ind<-as.matrix(apply(x2, 1, sum))
    rm[[i]]<-cbind(sort[[j]][[i]],ind[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acmh4[j]<-cor(ram[,2], means2$Y)
  gebvmh4<-cbind(gebvmh4,ram[,2]) 
}
rm(test, x,x1,x2, ind, i, j, multitraitmh4, ram, rm)


#### Indirect Multi-trait GS MA #####
means2$A<-scale(means2$A)
gebvma<-1:453
acma<-c()
for(j in 1:length(sort)){
  rm<-list()
  for(i in 1:5){
    test<-means2
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"M"]<-NA
    test[which(means2[,"geno"]%in%sort[[j]][[i]]),"A"]<-NA
    
    multitraitma<-asreml(fixed=cbind(M,A)~trait,                    
                        random=~us(trait):giv(geno),
                        rcov=~units:us(trait),
                        ginverse=list(geno=ginv), 
                        workspace=128e6, maxiter=400, data= test)
    
    x<-as.matrix(coef(multitraitma, pattern = 'geno')[1:453,])
    x1<-as.matrix(coef(multitraitma, pattern = 'geno')[454:906,])
    x2<-cbind(x, x1)
    
    x2[,1]<-x2[,1]*-0.3925930726
    x2[,2]<-x2[,2]*0.8441269
    ind<-as.matrix(apply(x2, 1, sum))
    rm[[i]]<-cbind(sort[[j]][[i]],ind[1:453,1][sort[[j]][[i]]]) 
  }
  
  ram<-rbind(rm[[1]],rm[[2]], rm[[3]], rm[[4]], rm[[5]])
  ram<-ram[order(ram[,1]),]
  acma[j]<-cor(ram[,2], means2$Y)
  gebvma<-cbind(gebvma,ram[,2]) 
}
rm(test, x,x1,x2, ind, i, j, multitraitma, ram, rm)

############## saving all GEBV files ######################
# write.csv(gebva, "gebva.csv", row.names = FALSE)
# write.csv(gebvh1, "gebvh1.csv", row.names = FALSE)
# write.csv(gebvh2, "gebvh2.csv", row.names = FALSE)
# write.csv(gebvh3, "gebvh3.csv", row.names = FALSE)
# write.csv(gebvh4, "gebvh4.csv", row.names = FALSE)
# write.csv(gebvm, "gebvm.csv", row.names = FALSE)
# write.csv(gebvma, "gebvma.csv", row.names = FALSE)
# write.csv(gebvmh3, "gebvmh3.csv", row.names = FALSE)
# write.csv(gebvmh4, "gebvmh4.csv", row.names = FALSE)
# write.csv(gebvy, "gebvy.csv", row.names = FALSE)
# write.csv(gebvya, "gebvya.csv", row.names = FALSE)
# write.csv(gebvyh1, "gebvyh1.csv", row.names = FALSE)
# write.csv(gebvyh2, "gebvyh2.csv", row.names = FALSE)
# write.csv(gebvyh3, "gebvyh3.csv", row.names = FALSE)
# write.csv(gebvyh4, "gebvyh4.csv", row.names = FALSE)
# write.csv(gebvym, "gebvym.csv", row.names = FALSE)
# write.csv(gebvymm, "gebvymm.csv", row.names = FALSE)
# write.csv(gebvyh4m, "gebvyh4m.csv", row.names = FALSE)
# write.csv(gebvyh3m, "gebvyh3m.csv", row.names = FALSE)
# write.csv(gebvyh2m, "gebvyh2m.csv", row.names = FALSE)
# write.csv(gebvyh1m, "gebvyh1m.csv", row.names = FALSE)
# write.csv(gebvyam, "gebvyam.csv", row.names = FALSE)
###############################################################

############## saving all accuracies ######################
vars<-ls()
vars<-vars[grepl("ac",vars)]
AC<-c()
for(i in 1:length(vars)){
  AC<-cbind(AC,eval(parse(text=vars[i])))
  colnames(AC)[i]<-vars[i]
}
# write.csv(AC, "accuracies.csv", row.names=F)
#### Coincidence index ####
source("CI.R")
### For each CI one must replace gebvy (standard GS) for the desired one ####
CI(x=gebvy, y=means2[,c("geno","Y")], s=0.2, top=T)
