readResFunc<-function(fileLoc,i,burn,prmFile,thin,end=300000){
  gg<-read.csv(paste0(fileLoc,i,".csv"))[,-1]
  gg2<-read.csv(paste0(fileLoc,i+8,".csv"))[,-1]
  gg$s_Val<-as.numeric(gg$s_Val)
  gg$oDist_u<-as.numeric(gg$oDist_u)
  gg<-subset(gg,is.na(betaVer)!=T)
  gg<-subset(gg,ll>-1000)
  gg<-gg[burn:end,]
  gg2<-subset(gg2,is.na(betaVer)!=T)
  gg2<-subset(gg2,ll>-1000)
  gg2<-gg2[burn:end,]
  ggF<-Filter(var, gg)
  ggF2<-Filter(var, gg2)
  gg$chainID<-1
  gg2$chainID<-2
  gg<-rbind(gg,gg2)
  gg<-gg[seq(1, nrow(gg), by=thin), ]
}

