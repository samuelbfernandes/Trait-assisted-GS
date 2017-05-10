###########################################################################
## 
## File Name: 'correlations'
##
## Authors:    SB Fernandes < samuelfernandes@agronomo.eng.br >
##             KOG Dias < kaioolimpio@hotmail.com >  
##
## Date:      May 10th, 2017
##
## Contents:  Multi-trait asreml models for obtaing genetic and
##            residual correlations
##
## input:  "means2.csv", "ginv.csv"  
##
## output: "corr.csv" 
############################################################################
#setwd("")
library(asreml)
ginv <- read.csv("ginv.csv")
attr(ginv,"rowNames")<-1:453
means2<-read.csv("means2.csv",row.names = 1)
means2<-transform(means2, geno=factor(geno))
##################################################################
#############################   Multitrait  ###############
###############################################################
multitraitm<-asreml(fixed=cbind(Y,M)~trait,                    
                    random=~corgh(trait):giv(geno),
                    rcov=~units:corgh(trait,init=c(-0.4,4.5,0.001)),
                    ginverse=list(geno=ginv), 
                    workspace=256e6, maxiter=400, data= means2)

multitraith1<-asreml(fixed=cbind(Y,h1)~trait,                    
                     random=~corgh(trait):giv(geno),
                     rcov=~units:corgh(trait),
                     ginverse=list(geno=ginv), 
                     workspace=128e6, maxiter=400, data= means2)

multitraith2<-asreml(fixed=cbind(Y,h2)~trait,                    
                     random=~corgh(trait):giv(geno),
                     rcov=~units:corgh(trait),
                     ginverse=list(geno=ginv), 
                     workspace=128e6, maxiter=500, data= means2)

multitraith3<-asreml(fixed=cbind(Y,h3)~trait,                    
                     random=~corgh(trait):giv(geno),
                     rcov=~units:corgh(trait),
                     ginverse=list(geno=ginv), 
                     workspace=128e6, maxiter=500, data= means2)

multitraith4<-asreml(fixed=cbind(Y,h4)~trait,                    
                     random=~corgh(trait):giv(geno),
                     rcov=~units:corgh(trait),
                     ginverse=list(geno=ginv), 
                     workspace=128e6, maxiter=400, data= means2)

multitraita<-asreml(fixed=cbind(Y,A)~trait,                    
                    random=~corgh(trait):giv(geno),
                    rcov=~units:corgh(trait,init=c(0.63,4,80)),
                    ginverse=list(geno=ginv), 
                    workspace=256e6, maxiter=500, data= means2)

corg<-rbind(summary(multitraitm)$varcomp["trait:giv(geno)!trait.M:!trait.Y.cor",c(2:3)],
            summary(multitraith1)$varcomp["trait:giv(geno)!trait.h1:!trait.Y.cor",c(2:3)],
            summary(multitraith2)$varcomp["trait:giv(geno)!trait.h2:!trait.Y.cor",c(2:3)],
            summary(multitraith3)$varcomp["trait:giv(geno)!trait.h3:!trait.Y.cor",c(2:3)],
            summary(multitraith4)$varcomp["trait:giv(geno)!trait.h4:!trait.Y.cor",c(2:3)],
            summary(multitraita)$varcomp["trait:giv(geno)!trait.A:!trait.Y.cor",c(2:3)])

corres<-rbind(summary(multitraitm)$varcomp["R!trait.M:!trait.Y.cor",c(2:3)],
              summary(multitraith1)$varcomp["R!trait.h1:!trait.Y.cor",c(2:3)],
              summary(multitraith2)$varcomp["R!trait.h2:!trait.Y.cor",c(2:3)],
              summary(multitraith3)$varcomp["R!trait.h3:!trait.Y.cor",c(2:3)],
              summary(multitraith4)$varcomp["R!trait.h4:!trait.Y.cor",c(2:3)],
              summary(multitraita)$varcomp["R!trait.A:!trait.Y.cor",c(2:3)])

corg<-data.frame(trait=c("M","H1","H2","H3","H4","A"),type="corg",corg)
corres<-data.frame(trait=c("M","H1","H2","H3","H4","A"),type="corres",corres)
corr<-rbind(corg,corres)
corr$sign<-ifelse(corr$component<0,"neg","pos")
corr$component<-abs(corr$component)
#write.csv(corr, "corr.csv", row.names = FALSE)
