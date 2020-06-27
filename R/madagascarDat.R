
#######Cambodia data
#######read in and clean data
madDatFunc<-function(fileLoc,siteLoc="all"){

  #######madagascar
  #######read in and clean data
  batsMad<-read_xls(fileLoc)
  batsMad$date <- as.Date(batsMad$date,format="%d/%m/%Y")

  Cutoff_E_dup = 402.9
  Cutoff_P_ruf = 65.55
  Cutoff_R_mad = 77.46

  batsMad<-subset(batsMad,date>=as.Date("2013-10-01"))

  batsMad<-subset(batsMad,species=="Eidolon dupreanum")
  if(siteLoc!="all")batsMad<-subset(batsMad,site==siteLoc)
  batsMad$pos<-0
  batsMad$pos[batsMad$titer>=Cutoff_E_dup]<-1
  batsMad$neg<-1
  batsMad$neg[batsMad$titer>=Cutoff_E_dup]<-0
batsMadPos<-aggregate(batsMad$pos~batsMad$date,FUN=sum)
batsMadNeg<-aggregate(batsMad$neg~batsMad$date,FUN=sum)
batsMadtmp<-merge(batsMadPos,batsMadNeg)
batsMadtmp<-batsMadtmp[-c(1:15),]#remove data from 2006 and 2007
names(batsMadtmp)<-c("date","pos","neg")

batsMadtmp<-batsMadtmp[order(batsMadtmp$date),]
batsMadtmp$month<-paste0(year(batsMadtmp$date+23),"-",month(batsMadtmp$date+23)) #aggregating from day +23 seems to give best pattern
batsMadtmp$prev<-batsMadtmp$pos/(batsMadtmp$pos+batsMadtmp$neg)
bb<-aggregate(batsMadtmp$prev~batsMadtmp$month,FUN=mean)
names(bb)<-c("Date","prev")
bb$Date<-as.Date(paste(bb$Date,"-01",sep=""))
#bb$Date<-as.Date(paste(substr(bb$Date, start = 1, stop = 4),(substr(bb$Date, start = 6, stop = 9)), 1, sep="-"), "%Y-%U-%u")

bb$pos<-aggregate(batsMadtmp$pos~batsMadtmp$month,FUN=sum)[,2]
bb$neg<-aggregate(batsMadtmp$neg~batsMadtmp$month,FUN=sum)[,2]

biCon<-binom.confint(bb$pos,(bb$neg+bb$pos),method="exact")
bb$prev<-bb$pos/(bb$neg+bb$pos)

print(ggplot(bb, aes(y=prev, x=Date))+
  geom_point()+
  geom_errorbar(aes(ymin=biCon$lower, ymax=biCon$upper), width=.2,
                position=position_dodge(0.05),alpha=0.5,colour="blue")+
  # geom_ribbon(aes(ymin=boundsPos[1,],ymax=boundsPos[2,]),alpha=0.3)+
  ylab("Sero+"))

return(bb)
}
#
#for(i in c(1:length(unique(gg$sampleID)))){
#  ss<-unique(gg$sampleID)[i]
#  h<-subset(gg,sampleID==ss)
#  if(sum(h$pos)>0&&sum(h$pos)<length(h$pos)) print(h$pos)
#
#}
#


#######Cambodia data
#######read in and clean data
camDatFunc<-function(fileLoc,siteLoc="all"){

  #######madagascar
  #######read in and clean data
  fileLoc<-"/Users/alm204/Dropbox (SPH Imperial College)/emma aaron/data/Cambodia/Capture and serology data .xlsx"
  batsCam<-read_xlsx(fileLoc)
  batsCam$SampleDate <- as.Date(batsCam$SampleDate,format="%d/%m/%Y")
  batsCam$Result[batsCam$Result=="Neg"]<-0
 batsCam$Result[batsCam$Result=="Pos"]<-1
  batsCam$Result[is.na(batsCam$Result)]<-0
  batsCam$Result[batsCam$Result=="Equivocal"]<-0
  batsCam$Result<-as.numeric(batsCam$Result)
  batsCam$Year<-year(batsCamMnth$Date)
  batsCam$Month<-month(batsCamMnth$Date)
  batsCamMnth<-aggregate(batsCam$Result~batsCam$Month+batsCam$Year,FUN=sum)
  names(batsCamMnth)<-c("Month","Year","Pos")

  batsCamMnth$Neg<-aggregate(batsCam$Result~batsCam$Month+batsCam$Year,FUN=length)[,3]
  batsCamMnth$Prev<-batsCamMnth$Pos/(batsCamMnth$Pos+batsCamMnth$Neg)

  batsCamMnth$Date<-as.Date(paste(batsCamMnth$Year,"-",batsCamMnth$Month,"-01",sep=""))
  plot(batsCamMnth$Prev~batsCamMnth$Date,xlab ="Date",ylab="Prev")

  biCon<-binom.confint(batsCamMnth$Pos,(batsCamMnth$Neg+batsCamMnth$Pos),method="exact")
  batsCamMnth$Prev<-batsCamMnth$Pos/(batsCamMnth$Neg+batsCamMnth$Pos)

  print(ggplot(batsCamMnth, aes(y=Prev, x=Date))+
    geom_point()+
    geom_errorbar(aes(ymin=biCon$lower, ymax=biCon$upper), width=.2,
                  position=position_dodge(0.05),alpha=0.5,colour="blue")+
    # geom_ribbon(aes(ymin=boundsPos[1,],ymax=boundsPos[2,]),alpha=0.3)+
    ylab("Sero+"))

}

