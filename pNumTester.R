
res<-NULL
for(i in seq(2,100,10)){
  print(i)
  particleNum<-rep(i,20)
obsPcrDens=batsBoonah
oDat=obsDataBoonah
currentParams=data.frame(paramsFunc())
obsPcrDens<-subset(obsPcrDens,pcrPos==1)
obsPcrDens<-subset(obsPcrDens,Species=="BFF")
obsPcrDens$log10Ser[is.na(obsPcrDens$log10Ser)]<-0

seroDensity=density(obsPcrDens$log10Ser)


curVal<- mcmapply(pFilt,n=particleNum, MoreArgs=list(iState,modStep,oDat,currentParams,mod=mod,obsPcrDens=obsPcrDens,seroDensity=seroDensity))
print(sd(curVal))
res<-rbind(res,sd(curVal))
}


res_365<-res
#res_365_4 down to around 2 at 52 particles
