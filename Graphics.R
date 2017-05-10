###########################################################################
## 
## File Name: 'Graphics'
##
## Author:    SB Fernandes < samuelfernandes@agronomo.eng.br >
##
## Date:      May 10th, 2017
##
## Contents:  ggplot2 codes to create all figures used in this paper
##
############## Get "accuracies.csv" file from GS_2st_stage.R ###############
################# Get "corr.csv" file from correlations.R ##################
############################################################################
library(ggplot2)
library(grid)

AC<-read.table("accuracies.txt", h=TRUE)
AC<-cbind(apply(AC, 2, mean),apply(AC, 2, sd))
AC<-data.frame(mode=rownames(AC),AC)
colnames(AC)<-c("model","mean", "sd")

#### Figure 1  #####
standard<-AC[c("acy","acm","ach1","ach2","ach3","ach4","aca"),]
standard$model<-c("Y", "M", "H1", "H2", "H3", "H4", "A")

ggplot(standard, aes(x = model, y = mean)) +
  geom_bar(colour = "black", fill = "red", stat = "identity", position = position_dodge(width = .8),size=.3, width=0.8)+
  geom_errorbar(aes(x = model, ymax=mean+sd, ymin=mean-sd) , width=.2, position=position_dodge(.8))+
  theme(panel.border = element_rect(fill=NA,color="darkred", size=0.5, linetype="solid"))+
  scale_y_continuous("Accuracy",limits=c(0,0.75), breaks=c(0, 0.25, 0.5, 0.75))  +
  scale_x_discrete(name="",limits=c("Y", "M", "H1", "H2", "H3", "H4", "A") )+
  annotate("text", family="Times",angle =90,x = standard$model, y = .15, 
           label = c("h = 0.51","h = 0.82","h = 0.54","h = 0.85","h = 0.90","h = 0.87","h = 0.94"), size = 10, parse=F)+
  theme_bw(base_size = 20) + theme( panel.border = element_rect(fill = NA, colour="black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family="Times"), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm"))+
  theme(axis.title.y = element_text(vjust=1.3,colour = "black",family="Times", size = 25),axis.text.y = element_text(  colour = "black", family="Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

ggsave(file="Fig1.eps",width=8,height=5,units=c("in"),dpi=300,family="Times")

#### Figure 2  #####
corr<-read.csv("corr.csv")
colnames(corr)<-c("trait","var", "cor","se", "sign")
pluser<- corr$cor+corr$se
minuser<-  ifelse(corr$cor-corr$se>0, corr$cor-corr$se,0)
library(scales)
ggplot(corr, aes(x = trait, y = cor, fill=as.character(interaction(corr$var,corr$sign)))) +
  geom_bar(colour="black", stat = "identity", position = position_dodge(width = .8),size=.3, width=.8)+
  geom_errorbar(aes(ymin = minuser, ymax = pluser), width = 0.2,position =position_dodge(width = .8))+
  scale_y_continuous("Correlation",limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  scale_x_discrete(name="",limits=c("M","H1","H2","H3","H4","A"), labels=c("M","H1","H2","H3","H4","A"))+
  scale_fill_manual(name="", limits=c("corg.pos", "corg.neg", "corres.pos", "corres.neg"),values = alpha(c("darkgreen","darkgreen","red","red"), c(1,.3,1,.3)),
                    labels=c("corg.pos"=expression(cor[g]*paste("","+")),"corg.neg"=expression(cor[g]*paste(""," -")),
                             "corres.pos"=expression(cor[r]*paste("","+")),"corres.neg"=expression(cor[r]*paste(""," -"))))+
  theme_bw(base_size = 20) + theme( panel.border = element_rect(fill = NA, colour="black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family="Times"), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm"))+
  theme(axis.title.y = element_text(vjust=1.3,colour = "black",family="Times", size = 25),axis.text.y = element_text(  colour = "black", family="Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+theme(legend.title=element_blank(),legend.text = element_text(family="Times", face="bold",size = 20),legend.key.width = unit(1, 'cm'))

ggsave(file="cor.pdf",width=8,height=5,units=c("in"),dpi=300,family="Times")

#### Figure 3  #####
indirect<-AC[c("acy","acm_y","ach1_y","ach2_y","ach3_y","ach4_y","aca_y", "acmh3", "acmh4","acma"),]
indirect$model<-c("Y", "M_Y", "H1_Y", "H2_Y", "H3_Y", "H4_Y", "A","MH3","MH4", "MA")
indirect$type<-as.factor(c(1,2,2,2,2,2,2,3,3,3))

ggplot(indirect, aes(x = model, y = mean, fill=type)) +
  geom_bar(colour = "black", stat = "identity", position = position_dodge(width = .8),size=.3, width=0.8)+
  geom_errorbar(aes(x = model, ymax=mean+sd, ymin=mean-sd) ,size=.3, width=.2, position = position_dodge(width = .8))+
  scale_y_continuous("Accuracy",limits=c(0,0.75), breaks=c(0, 0.25, 0.5, 0.75))  +
  scale_x_discrete(limits=c("Y", "M_Y", "H1_Y", "H2_Y", "H3_Y", "H4_Y", "A","MH3","MH4", "MA"),labels=c("Y", "M", "H1", "H2", "H3", "H4", "A","MH3", "MH4",   "MA"),name="" )+
  scale_fill_manual("",labels=c("Standard  ","Indirect  ","Multi-trait indirect"),guide = "legend",values=c("darkgreen","blue","red" ))+
  theme(panel.border = element_rect( fill=NA, size=0.5,  linetype="solid"))+
  theme_bw(base_size = 20) + theme( panel.border = element_rect(fill = NA, colour="black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family="Times"), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm"))+
  theme(axis.title.y = element_text(vjust=1.3,colour = "black",family="Times", size = 20),axis.text.y = element_text(  colour = "black", family="Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.margin.y = unit(.5, "lines"))+
  theme(strip.text.x = element_text(size=18,face="bold", hjust = 0.5, vjust = 0.5),
        strip.text.y = element_text(size=18, face="bold",  hjust = 0.5, vjust = 0.5))+
  theme(legend.title=element_text(family="Times", size = 15),legend.text = element_text(size = 20),legend.key.width = unit(0.5, 'cm'))+ theme(legend.position = "bottom")

ggsave(file="Fig3.eps",width=8,height=5,units=c("in"),dpi=300,family="Times")

#### Figure 4  #####
multitrait<-AC[c("acy","acymm","acyh1m","acyh2m","acyh3m","acyh4m","acyam", "acym","acyh1","acyh2","acyh3","acyh4","acya"),]
multitrait$model<-c("Y",   "YM" ,  "YH1",   "YH2" ,  "YH3" ,  "YH4" ,  "YA",  "YM" ,  "YH1",   "YH2" ,  "YH3" ,  "YH4" ,  "YA")
multitrait$type<-as.factor(c(1,4,4,4,4,4,4,5,5,5,5,5,5))

ggplot(multitrait, aes(x =model , y = mean, fill=type)) +
  geom_bar(colour = "black", stat = "identity", position = position_dodge(width = .8),size=.3, width=0.8)+
  geom_errorbar(aes(x = model, ymax=mean+sd, ymin=mean-sd) , width=.3, position=position_dodge(.8))+
  scale_y_continuous("Accuracy",limits=c(0,0.75), breaks=c(0, 0.25, 0.5, 0.75))  +
  scale_x_discrete(limits=c("Y",   "YM" ,  "YH1",   "YH2" ,  "YH3" ,  "YH4" ,  "YA"),name="" ,
                   labels=c("Y",   "YM" ,  "YH1",   "YH2" ,  "YH3" ,  "YH4" ,  "YA"))+
  scale_fill_manual("",labels=c("Standard  ","Multi-trait  ","Trait-assisted"),guide = "legend",values=c("darkgreen","blue","red" ))+
  theme(panel.border = element_rect( fill=NA, size=0.5,  linetype="solid"))+
  theme_bw(base_size = 20) + theme( panel.border = element_rect(fill = NA, colour="black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family="Times"), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm"))+
  theme(axis.title.y = element_text(vjust=1.3,colour = "black",family="Times", size = 20),axis.text.y = element_text(  colour = "black", family="Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.margin.y = unit(.5, "lines"))+
  theme(strip.text.x = element_text(size=18,face="bold", hjust = 0.5, vjust = 0.5),
        strip.text.y = element_text(size=18, face="bold",  hjust = 0.5, vjust = 0.5))+
  theme(legend.title=element_text(family="Times", size = 15),legend.text = element_text(size = 20),legend.key.width = unit(0.5, 'cm'))+ theme(legend.position = "bottom")

ggsave(file="Fig4.eps",width=8,height=5,units=c("in"),dpi=300,family="Times")
