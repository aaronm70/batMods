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
  batsBoonah$Urine.PCR<-as.character(batsBoonah$Urine.PCR)
  batsBoonah$Urine.PCR[batsBoonah$Urine.PCR==">40"]<-"0"
  batsBoonah$Urine.PCR[is.na(batsBoonah$Urine.PCR)]<-"0"
  batsBoonah$Urine.PCR[batsBoonah$Urine.PCR!="0"]<-"1"
  batsBoonah$Urogential.PCR<-as.character(batsBoonah$Urogential.PCR)
  batsBoonah$Urogential.PCR[batsBoonah$Urogential.PCR==">40"]<-"0"
  batsBoonah$Urogential.PCR[is.na(batsBoonah$Urogential.PCR)]<-"0"
  batsBoonah$Urogential.PCR[batsBoonah$Urogential.PCR!="0"]<-"1"
  batsBoonah %<>% mutate_if(is.character,as.numeric)

  #Sum PCR values to get pcr positives but convert to binomila, positive yes/no
  batsBoonah$pcrPos<-rowSums(batsBoonah[,c(19:20)])
  batsBoonah$pcrPos[batsBoonah$pcrPos>0]<-1 #in case bat is positive in both urine and urogential
  #Empirical states PCR and Sero
  batsBoonah$PpSp<-0
  batsBoonah$PpSn<-0
  batsBoonah$PnSp<-0
  batsBoonah$PnSn<-0

  batsBoonah$PpSp<-ifelse(batsBoonah$pos==T&batsBoonah$pcrPos==1,1,0)
  batsBoonah$PpSn<-ifelse(batsBoonah$pos==F&batsBoonah$pcrPos==1,1,0)
  batsBoonah$PnSp<-ifelse(batsBoonah$pos==T&batsBoonah$pcrPos==0,1,0)
  batsBoonah$PnSn<-ifelse(batsBoonah$pos==F&batsBoonah$pcrPos==0,1,0)

  obsDataPpSp<-aggregate(batsBoonah$PpSp~batsBoonah$Date,FUN=sum)
  obsDataPpSn<-aggregate(batsBoonah$PpSn~batsBoonah$Date,FUN=sum)
  obsDataPnSp<-aggregate(batsBoonah$PnSp~batsBoonah$Date,FUN=sum)
  obsDataPnSn<-aggregate(batsBoonah$PnSn~batsBoonah$Date,FUN=sum)


  ####
  #Aggregate serology and PCR by date to get number of positive bats on a given sampling date
  obsDataBoonah<-aggregate(batsBoonah$pos~batsBoonah$Date,FUN=sum)
  obsDataNeg<-aggregate(batsBoonah$neg~batsBoonah$Date,FUN=sum)
  obsDataPCR<-aggregate(batsBoonah$pcrPos~batsBoonah$Date,FUN=sum)
  obsDataSer<-aggregate(batsBoonah$HeV.Serology~batsBoonah$Date,FUN=mean)
  obsDataPCR<-obsDataPCR[-44,]#No serology values for this date (2014-04-29), all NA's missing data? so removed from analysis, otherwise gives unlikely 0 prevalence on this date...

  #start constructing data frame
  obsDataBoonah$neg<-obsDataNeg$`batsBoonah$neg`
  obsDataBoonah$pcrPos<-obsDataPCR$`batsBoonah$pcrPos`
  obsDataBoonah$meanSer<-obsDataSer$`batsBoonah$HeV.Serology`
  obsDataBoonah$PpSp<-obsDataPpSp$`batsBoonah$PpSp`[-44]
  obsDataBoonah$PpSn<-obsDataPpSn$`batsBoonah$PpSn`[-44]
  obsDataBoonah$PnSp<-obsDataPnSp$`batsBoonah$PnSp`[-44]
  obsDataBoonah$PnSn<-obsDataPnSn$`batsBoonah$PnSn`[-44]

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
    sdProps[2]<-0.1
    sdProps[which(prms==0)]<-0
    sdProps[c(10,15)]<-0.1
    sdProps[c(5,12,20,21,22)]<-0
    sdProps[c(19,25)]<-ifelse(prms$zeta_p!=0, 0.1,0)#pcr probs
    sdProps[c(11,16,24)]<-0.5
    sdProps[16]<-4
    sdProps[11]<-0.5
    sdProps[17]<-0.1
    sdProps[23]<-5
    sdProps[28]<-0
    sdProps[29]<-0
    sdProps[c(30,31)]<-0
    sdProps[c(6,7,8,14,19,25)]<-0.05
    sdProps[c(27)]<-0.5
    sdProps[c(26,18)]<-0

  }
  else {
    sdProps = 100*unlist(prms)
    sdProps[c(10,15)]<-1
    sdProps[2]<-1
    sdProps[3]<-0.3

    sdProps[which(prms==0)]<-0
    sdProps[c(5,12,20,21,22)]<-0
    sdProps[16]<-15
    sdProps[11]<-1
    sdProps[17]<-0.5
    sdProps[13]<-25
    sdProps[19]<-0.5
    sdProps[25]<-0.5
    sdProps[23]<-30
    sdProps[28]<-0
    sdProps[29]<-0
    sdProps[24]<-4
    sdProps[c(6,7,8,14)]<-0.1
    sdProps[c(27)]<-2
    sdProps[c(30,31)]<-0
    sdProps[c(26,18)]<-0

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
  return(obsDataBoonahTmp)
}



