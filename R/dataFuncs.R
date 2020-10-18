library(lubridate)
library(readxl)
library(magrittr)
library(dplyr)


#'batDataFunc
#'@param ret which data type to return all or observed
#'@param species which species to return for BFF GHFF or LRFF
#'@param cutoff the cutoff value to use for serology
#'@param fileLoc location of raw data file
#'@return database
boonahDatFunc<-function(ret="all",species="BFF",cutoff=1636,fileLoc){

  #read in raw data
  batsBoonah<-read.csv(fileLoc,head=T)
  batsBoonah$Date <- as.Date(batsBoonah$Date,format="%d/%m/%Y")

  #if specific species chosen subset by species
  if(species!="all"){
    batsBoonah<-subset(batsBoonah,Species==species)
  }

  #add serological positive and negatives by cutoff value
  batsBoonah$pos <- batsBoonah$HeV.Serology >= cutoff
  batsBoonah$neg <- batsBoonah$HeV.Serology < cutoff
  #convert urine and urogential values to binomials

  ##remove bats which are missing either pcr for urine or urogential
  batsBoonah<-batsBoonah[!is.na(batsBoonah$Urine.PCR) | !is.na(batsBoonah$Urogential.PCR),]

  #convert to binary pos/neg
  batsBoonah$Urine.PCR<-as.character(batsBoonah$Urine.PCR)
  batsBoonah$Urine.PCR[batsBoonah$Urine.PCR==">40"]<-"0"
  batsBoonah$Urine.PCR[batsBoonah$Urine.PCR>"0"]<-"1"
  batsBoonah$Urogential.PCR<-as.character(batsBoonah$Urogential.PCR)
  batsBoonah$Urogential.PCR[batsBoonah$Urogential.PCR==">40"]<-"0"
  batsBoonah$Urogential.PCR[batsBoonah$Urogential.PCR>"0"]<-"1"
  batsBoonah %<>% mutate_if(is.character,as.numeric)

  #Sum PCR values to get pcr positives but convert to binom, positive yes/no
  batsBoonah$pcrPos<-rowSums(batsBoonah[,c("Urine.PCR", "Urogential.PCR")],na.rm = T)
  batsBoonah$pcrPos[batsBoonah$pcrPos>0]<-1 #in case bat is positive in both urine and urogential
  #Empirical states PCR and Sero
  batsBoonah$PpSp<-0
  batsBoonah$PpSn<-0
  batsBoonah$PnSp<-0
  batsBoonah$PnSn<-0

 #remove bats with no individual serological data
  batsBoonah<-batsBoonah[is.na(batsBoonah$HeV.Serology)==F,]


  batsBoonah$PpSp<-ifelse(batsBoonah$pos==T&batsBoonah$pcrPos==1,1,0)
  batsBoonah$PpSn<-ifelse(batsBoonah$pos==F&batsBoonah$pcrPos==1,1,0)
  batsBoonah$PnSp<-ifelse(batsBoonah$pos==T&batsBoonah$pcrPos==0,1,0)
  batsBoonah$PnSn<-ifelse(batsBoonah$pos==F&batsBoonah$pcrPos==0,1,0)

  obsDataPpSp<-aggregate(batsBoonah$PpSp~batsBoonah$Date,FUN=sum,na.rm=T)
  obsDataPpSn<-aggregate(batsBoonah$PpSn~batsBoonah$Date,FUN=sum,na.rm=T)
  obsDataPnSp<-aggregate(batsBoonah$PnSp~batsBoonah$Date,FUN=sum,na.rm=T)
  obsDataPnSn<-aggregate(batsBoonah$PnSn~batsBoonah$Date,FUN=sum,na.rm=T)


  ####
  #Aggregate serology and PCR by date to get number of positive bats on a given sampling date
  obsDataBoonah<-aggregate(batsBoonah$pos~batsBoonah$Date,FUN=sum)
  obsDataNeg<-aggregate(batsBoonah$neg~batsBoonah$Date,FUN=sum)
  obsDataPCR<-aggregate(batsBoonah$pcrPos~batsBoonah$Date,FUN=sum)
  obsDataSer<-aggregate(batsBoonah$HeV.Serology~batsBoonah$Date,FUN=mean)

  #start constructing data frame
  obsDataBoonah$neg<-obsDataNeg$`batsBoonah$neg`
  obsDataBoonah$pcrPos<-obsDataPCR$`batsBoonah$pcrPos`
  obsDataBoonah$meanSer<-obsDataSer$`batsBoonah$HeV.Serology`
  obsDataBoonah$PpSp<-obsDataPpSp$`batsBoonah$PpSp`
  obsDataBoonah$PpSn<-obsDataPpSn$`batsBoonah$PpSn`
  obsDataBoonah$PnSp<-obsDataPnSp$`batsBoonah$PnSp`
  obsDataBoonah$PnSn<-obsDataPnSn$`batsBoonah$PnSn`

  names(obsDataBoonah)<-c("Date","positives","negatives","pcrPos","meanSer","PpSp","PpSn","PnSp","PnSn")
  startdate <- as.Date(obsDataBoonah$Date[1],"%d/%m/%Y")

  #add a running value in days
  obsDataBoonah$NumDays  <- as.vector(difftime(obsDataBoonah$Date,startdate ,units="days"))
  #add serological prevalence
  obsDataBoonah$prev<-obsDataBoonah$positives/(obsDataBoonah$negatives+obsDataBoonah$positives)

  if(ret=="all") return(batsBoonah)
  if(ret=="obs") return(obsDataBoonah)

}


#set SD vals for pMCMC
#' @param prms current parameters
#' @max whether they are max SD values or starting SD values
#' @return sd values
sdFunc<-function(prms,max=F){
  if(max==F) {
    sdProps = 3*unlist(prms)

    sdProps["omega_2_Val"]<-5
    sdProps["gamma_2_Val"]<-0.5
    sdProps[c("epsilon_Val", "rho_Val" )]<-0.5
    sdProps[c( "mu_Val" , "c_Val","sigmaVer" ,"gammaVer", "betaVer" )]<-0
    sdProps[c("R0_Val","d_val"   )]<-2
    sdProps["R0_Val"]<-4
    sdProps["kappa_Val"]<-0.2
    sdProps["Phi2_val"]<-0.1
    sdProps["S2_val"]<-25
    sdProps["s_Val"]<-5
    sdProps["c_val2"]<-0
    sdProps["envOscType"]<-0
    sdProps[c("betaFX", "betaFXVal")]<-0
    sdProps[c("mj_Val", "m_Val","omega_m_Val","phi_Val","zeta_p","zeta_s","pcrProb2"  )]<-0.05
    sdProps["oDist_u"]<-1
    sdProps[c("oDist1","oDist_s")]<-0
   # sdProps["mu_Val"]<-0.05

    sdProps[which(prms==0)]<-0

  }
  else {
    sdProps = 100*unlist(prms)

    sdProps["omega_2_Val"]<-10
    sdProps["gamma_2_Val"]<-0.5
   sdProps[c("epsilon_Val", "rho_Val" )]<-0.5
    sdProps[c( "mu_Val" , "c_Val","sigmaVer" ,"gammaVer", "betaVer" )]<-0
    sdProps[c("R0_Val","d_val"   )]<-35
    sdProps["kappa_Val"]<-1
    sdProps["Phi2_val"]<-0.2
    sdProps["S2_val"]<-20
    sdProps["s_Val"]<-10
    sdProps["c_val2"]<-0
    sdProps["envOscType"]<-0
    sdProps[c("betaFX", "betaFXVal")]<-0
    sdProps[c("mj_Val", "m_Val","omega_m_Val","phi_Val","zeta_p","zeta_s","pcrProb2"  )]<-0.3
    sdProps["oDist_u"]<-3
    sdProps[c("oDist1","oDist_s")]<-0
  #  sdProps["mu_Val"]<-0.05

    sdProps[which(prms==0)]<-0

  }
  return(sdProps)
}



#Add under roost data
#'@param obsDataBoonah data cleaned by boonahDatFunc
#'@param fileLocUrine location of under roost urine data
#'@return data with under roost data added
addUR<-function(obsDataBoonah,fileLocUrine){

  #read in under roost urine data
  batsBoonahUn<-read.csv(fileLocUrine,head=T)
  batsBoonahUn$date <- as.Date(batsBoonahUn$date,format="%d/%m/%Y")

  #subset by location Boonah (22) and date to match with individual sampling period
  batsBoonahUn<-subset(batsBoonahUn,locnum==22)
  batsBoonahUn$day<-day(batsBoonahUn$date)
  batsBoonahUn$day<-paste0(month(batsBoonahUn$date),"_",batsBoonahUn$day)
  batsBoonahUn<-subset(batsBoonahUn,as.Date(date)>= as.Date("2013-06-19"))
  batsBoonahUn<-subset(batsBoonahUn,as.Date(date)<= as.Date("2014-06-04"))

  #Aggregate samples by sampling date
  batsBoonahUnSUM<-aggregate(batsBoonahUn$res~batsBoonahUn$date,FUN=sum)
  batsBoonahUn<-aggregate(batsBoonahUn$res~batsBoonahUn$date,FUN=length)
  batsBoonahUn$pos<-batsBoonahUnSUM[,2]
  batsBoonahUn$prop<-batsBoonahUnSUM[,2]/batsBoonahUn[,2]
  names(batsBoonahUn)<-c("day","total","pos","prop")

  #add day of the year
  batsBoonahUn$t<-as.numeric(yday(batsBoonahUn$day))
  batsBoonahUn$t[c(8:13)]<-batsBoonahUn$t[c(8:13)]+365

  #GAM model for predicted values
  smth<-gam(prop~s(t),data=batsBoonahUn)
  predDat<-(predict(smth,newdata = data.frame(t=c(170:555)),se.fit = T))

  predDat<-as.data.frame(predDat)
  predDat$upr <- predDat$fit + (2 * predDat$se.fit)
  predDat$lwr <- predDat$fit - (2 * predDat$se.fit)
  predDat$fit[predDat$fit<0]<-0 #remove any negative fit values

  #add dates
  predDat$day<-c(1:386)
  predDat$Date<-as.Date("2013-12-31")
  predDat$Date[c(1:211)]<- as.Date(predDat$day[c(1:211)], origin = "2013-06-18")
  predDat$Date[c(212:386)]<- as.Date(predDat$day[c(0:175)], origin = "2014-01-15")
  predDat<-predDat[order(predDat$Date),]
  predDat$day<-c(170:555)
  names(batsBoonahUn)<-c("day.y","total","pos","prop","t")

  #merge in predicted data
  predDat<-merge(predDat,batsBoonahUn,by.x="day",by.y="t",all=T)

  gUR<-ggplot(predDat,aes(x=Date,y=fit))+
          ylim(-.1,0.5)+
          geom_line(col="purple")+
          geom_point(aes(y=prop,x=Date),col="#102F47")+
          geom_ribbon(aes(ymin=lwr,ymax=upr),fill="purple",alpha=0.1)+
          ylab("Under roost prevalence")+
         # scale_colour_manual(labels=c("observed","95% CI","GAM"),values=c("#102F47","grey","#F77F00"))+
          theme_bw()+
          theme_bw(base_size = 20)

  print(gUR)
  obsDataBoonah<-merge(obsDataBoonah,predDat,by.x="Date",by.y = "Date",all=F)

  return(obsDataBoonah)
}



#'aggregate by week boonah data
#'@param data boonah data
#'@return boonah data aggregated by week
weekification<-function(obsDataBoonah){

  #add which week/period to aggregate by, checked by hand as sampling dates not consistent
  obsDataBoonah$week<-c(rep(0,1),rep(1,2),rep(2,3),rep(3,2),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,3),rep(9,4),
                        rep(10,1),rep(11,3),rep(12,3),rep(13,3),rep(14,3),rep(15,4),rep(16,2))


  tmp<-cbind(aggregate(.~week, obsDataBoonah[,c(16,2, 3,4,12,13,14,15)], sum),
             aggregate(.~week, obsDataBoonah[,c(16,5)], mean))[,-9]

  #take mid dates of weekly aggregates
  obsDataBoonahTmp<- obsDataBoonah[c(1,2,5,7,10,14,18,22,26,29,32,34,37,40,43,46,49),][,-c(2,3,4,5,12,13,14,15)]
  obsDataBoonahTmp<-merge(tmp,obsDataBoonahTmp,by.x="week",by.y="week")
  obsDataBoonahTmp<-obsDataBoonahTmp[,c(1,2,3,4,9,10:16,5:8)]
  obsDataBoonahTmp<-obsDataBoonahTmp[order(obsDataBoonahTmp$Date),]
  obsDataBoonahTmp$prev<-obsDataBoonahTmp$positives/(obsDataBoonahTmp$positives+obsDataBoonahTmp$negatives)

  #add 50 years previous

  return(obsDataBoonahTmp)
}




