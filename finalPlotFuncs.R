library(miscTools)
library(ggpubr)
library(binom)
library(reshape2)
library(loo)
library(cowplot)

##File locations##
if(Sys.info()['sysname']=="Darwin"){
  fileLoc="/Users/alm204/OneDrive/emma aaron/data/Boonah/hendra-virus-test-results-flying-foxes.csv"
  fileLocUrine="/Users/alm204/OneDrive/emma aaron/data/Boonah/henrda-underRoostUrine.csv"
  prmFileLoc="/Users/alm204/Documents/ModelSetups.csv"
}else{
  fileLoc="/home/aaron/hendra-virus-test-results-flying-foxes.csv"
  fileLocUrine="/home/aaron/henrda-underRoostUrine.csv"
  prmFileLoc="/home/aaron/ModelSetups.csv"

}
prmLst<-prmLstFunc(prmFileLoc)


resultsFile<-"/Users/alm204/Documents/Cambridge/results/jun22/res_"
prmFile<-read.csv("/Users/alm204/Documents/ModelSetups.csv")#read in model setups
saveLoc<-"/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/"
burn<-110000
modNums<-c(1:8)#model numbers to plot (see model setups table, prmFile)
#Plot Rt values
plotAllRt(fileLoc= resultsFile,Rt=T,saveLoc=saveLoc,burn=burn,modNums=modNums,prmFile=prmFile)
#Plot fitted simulat
plotAllRt(fileLoc= resultsFile,Rt=F,saveLoc=saveLoc,burn=burn,modNums=modNums,prmFile=prmFile)
plotEnvForc(fileLoc= resultsFile,burn=burn,modNums=modNums,saveLoc=saveLoc,prmFile=prmFile)
#run function to obtain looic vals, R0's, contributing bats etc
parmVals<-looicFunc(resultsFile=resultsFile,prmFile=prmFile,burn=burn,samples=0,prmLst=prmLst)

looicVal<-parmVals[[1]]
zetaS<-parmVals[[5]]
zetaP<-parmVals[[4]]
contBats<-parmVals[[3]]
#plotLooic(r0Val,type="R0")
#plotLooic(contBats,type="batCont")
fNm<-"/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/"

p1<-plotLooic(fileNme=paste0(fNm,"looicVals.png"),loicVals=looicVal,type="loo")
p2<-plotLooic(fileNme=paste0(fNm,"contBatsX.png"),loicVals=contBats,type="batCont")
p3<-plotLooic(fileNme=paste0(fNm,"zetaPlotsS.png"),loicVals=zetaS,type="zetaS")
p4<-plotLooic(fileNme=paste0(fNm,"zetaPlotsP.png"),loicVals=zetaP,type="zetaP")
p5<-plotLooic(fileNme=paste0(fNm,"zetaPlotsU.png"),loicVals=zetaU,type="zetaU")
p6<-plotLooic(fileNme=paste0(fNm,"r0Val.png"),loicVals=r0Val,type="R0")

gx<-ggarrange(p1,p2,p3,p4,p5,p6,common.legend = T)

ggsave(
  paste0(fNm,"loZetVals.png"),
  plot = gx,
  height = 25,
  width = 40,
  units = "cm"
)
