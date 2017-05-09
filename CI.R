#### coincidence index (Hamblin and Zimmermann 1986) ####
#x: GEBV's from a given model
#y: reference means
#s: proportion of selection
#top: selection for top or bottom genotypes
CI<-function(x,y,s=0.2,top=T){
  ci<-c()
  for(i in 2:ncol(x)){
    x2<-as.matrix(x)
    y2<-as.matrix(y)
    size<-ceiling(nrow(x2)*s)
    x2<-x2[order(x2[,i], decreasing=top),]
    y2<-y2[order(y2[,2], decreasing=top),]
    both<-sum(x2[1:size,1]%in%y2[1:size,1])
    random<-both*s
    ci[i-1]<-(both-random)/(size-random)
  }
  return(list(ci=mean(ci), sd=sd(ci)))
}
