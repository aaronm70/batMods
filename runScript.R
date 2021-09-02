library("magrittr")
library("dplyr")
library("odin")
library("lubridate")
library("readxl")
library("coda")
library("parallel")
library("R.utils")
library("tidyr")
library("batMods")
library("mgcv")
library("ggplot2")
library("rmutil")


##File locations (diff locations when running on cluster)##
if(Sys.info()['sysname']=="Darwin"){
fileLoc="data/hendra-virus-test-results-flying-foxes.csv"
fileLocUrine="data/henrda-underRoostUrine.csv"
prmFileLoc="data/ModelSetups.csv"
}else{
fileLoc="data\\hendra-virus-test-results-flying-foxes.csv"
fileLocUrine="data\\henrda-underRoostUrine.csv"
prmFileLoc="data\\ModelSetups.csv"

}

#library(MCMCvis)
#ff1<-(read.csv("C:\\Users\\aaron\\Documents\\batMods\\chain2\\4.csv"))
#ff2<-(read.csv("C:\\Users\\aaron\\Documents\\batMods\\chain2\\12.csv"))
#ff<-mcmc.list(as.mcmc(ff1[1:429998,],thin=10),as.mcmc(ff2[1:429998,],thin=10))
#MCMCvis::MCMCtrace(ff,iter=100000)
#

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

obsDataBoonah<-rbind(data.frame(Date=as.Date("1993-06-19"),positives=as.numeric(0),negatives=as.numeric(0), pcrPos=as.numeric(0),meanSer=as.numeric(0)
                                   ,NumDays=as.numeric(0),prev=as.numeric(0),fit=as.numeric(0),se.fit=as.numeric(0),upr=as.numeric(0),lwr=as.numeric(0),
                                   PpSp=as.numeric(0),PpSn=as.numeric(0),PnSp=as.numeric(0),PnSn=as.numeric(0)), obsDataBoonah)

#aggregate data by sampling week
obsData<-weekification(obsDataBoonah)



##read in parameter list and randomly select starting params
prmLst<-prmLstFunc(prmFileLoc) #list of model structures and setups with associated starting parameter values
#calculate maturation rates from maternal immunity etc.
#nspan <- 0.555       #duration of newborn period (i.e., duration maternal immunity) (1.800^-1)
#omegam <- 1/nspan                           #maternal immune waning rate
#jspan<-15.55/12 #15.4 approx months of juvenile lifespan
#mu <- 1/(jspan - nspan)                     #maturation rate among non-newborn juveniles
library(furrr)

#odin_package("C:\\Users\\aaron\\Documents\\batMods\\batMods")

prmLst<-list(prmLst[[2]],prmLst[[4]],prmLst[[7]],prmLst[[8]])
prmLst[[1]]$modNum<-2
prmLst[[2]]$modNum<-4
prmLst[[3]]$modNum<-7
prmLst[[4]]$modNum<-8


##fit using pMCMC
library(future.apply)
iters=1000000
plan(multisession)



gg<-mcmcSampler(initParams=prmLst[[4]], #fit to one of the model structures
             sdProps=NULL,
             maxSddProps=NULL,
             niter=iters,
             particleNum=100,
             proposer = sequential.proposer,
             proposerType = "seq",
             startAdapt = 20000,
             nburn = iters*0.1,
             acceptanceRate = 0.25,
             stoch=F,
             adptBurn =10000,
             tell=10000,
             monitoring =2,
             oDat=obsData,
             likelihoodFunc = likelihoodFuncBoonah,
             priorFunc=lpriorBoonah,
             switch=500000,
             switchBlock = 1000000000,
             juvenileInfection=F)

saveRDS(gg,"results1Mil.RDS")
