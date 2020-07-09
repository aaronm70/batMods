
library(magrittr)
library("dplyr")
library("odin")
library("lubridate")
library("readxl")
library("coda")
#library("mnormt")
library("parallel")
library("R.utils")
library("tidyr")
library("batMods")
library("mgcv")
library("ggplot2")
library("rmutil")
#install.packages("mcmcr")
#library("lhs")


##File locations (diff locations when running on cluster)##
if(Sys.info()['sysname']=="Darwin"){
fileLoc="data/hendra-virus-test-results-flying-foxes.csv"
fileLocUrine="data/henrda-underRoostUrine.csv"
prmFileLoc="data/ModelSetups.csv"
}else{
fileLoc="/home/aaron/hendra-virus-test-results-flying-foxes.csv"
fileLocUrine="/home/aaron/henrda-underRoostUrine.csv"
prmFileLoc="/home/aaron/ModelSetups.csv"

}
set.seed(2)

#read in bat data for boonah
obsDataBoonah<-boonahDatFunc(ret="obs",species="BFF",fileLoc=fileLoc)

dtt=4#time step - split day by
#add a running total of time steps
obsDataBoonah$NumDays<-obsDataBoonah$NumDays+yday("2013-06-19") #yday function gives the day of the year

#add a 50 year equilising period for the model to reach an equilibrium before fitting to data, each day in the year is split into .25 days (so each time step is a quater day)
obsDataBoonah$NumDays<-(obsDataBoonah$NumDays*dtt)+(50*365*dtt)
obsDataBoonah$yearlyDay<-yday(obsDataBoonah$Date)

#add under roost urine data
obsDataBoonah<-addUR(obsDataBoonah,fileLocUrine=fileLocUrine)
obsDataBoonah<-obsDataBoonah[order(obsDataBoonah$Date),]
obsDataBoonah<-obsDataBoonah[,c(1:5,10:11,14:17,6:9)] #remove some columns that are not required

#add starting date for model (50 years prior to fitting data)
obsDataBoonah<-rbind(data.frame(Date=as.Date("1993-06-19"),positives=as.numeric(0),negatives=as.numeric(0), pcrPos=as.numeric(0),meanSer=as.numeric(0)
                                ,NumDays=as.numeric(0),prev=as.numeric(0),fit=as.numeric(0),se.fit=as.numeric(0),upr=as.numeric(0),lwr=as.numeric(0),
                                PpSp=as.numeric(0),PpSn=as.numeric(0),PnSp=as.numeric(0),PnSn=as.numeric(0)), obsDataBoonah)

#aggregate data by sampling week
obsDataBoonah<-weekification(obsDataBoonah)

obsData<-obsDataBoonah

##read in parameter list
prmLst<-prmLstFunc(prmFileLoc) #list of model structures and setups with associated starting parameter values


ff2<-mcmcSampler(initParams=prmLst[[4]], #fit to one of the model structures, use mcmapply to fit to multiple at once
             sdProps=NULL,
             maxSddProps=NULL,
             niter=9900,
             particleNum=50,
             proposer = multiv.proposer,
             proposerType = "block",
             startAdapt = 1000,
             nburn = 50,
             acceptanceRate = 0.3,
             stoch=F,
             adptBurn =200,
             tell=5,
             monitoring =2,
             oDat=obsData,
             likelihoodFunc = likelihoodFuncBoonahStoch,
             priorFunc=lpriorBoonah,
             switch=500,
             switchBlock = 50000,
             juvenileInfection=F)
