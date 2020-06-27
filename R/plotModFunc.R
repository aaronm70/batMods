library(reshape2)
library(parallel)
library(cowplot)


betaShapeFinder<-function(size,m,s){
  x<-rmutil::rbetabinom(n=10000,size=size,m=m,s=s)
  x = (x-min(x))/(max(x)-min(x))
  mu<-mean(x)
  vari<-var(x)
  alpha <- ((1 - mu) / vari - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#' plot model from mcmc results
#' @param mcmc mcmc results
#' @return figures
modPlotFuncOLD <- function(mcmc,assum="EI",fxd=F,osc=F) {
  #parse input data
  simALLMed<-NULL
  simTotalMed<-NULL
  simDatMed<-NULL
  simPosMed<-NULL
  simPosMed2<-NULL
  simTotalMed2<-NULL

  trajectoriesPos <- matrix(nrow=25,ncol=14)
  trajectoriesAll <- matrix(nrow=25,ncol=14)
  trajectoriesPosAll <- matrix(nrow=25,ncol=320)

  for (l in c(1:25)){

    currentTime <-0
    nextTime <- 18771
    prms<-unlist(sample_n(as.data.frame(mcmc),1))


    if(l==1){prms<-colMedians(mcmc)}
    prms<-as.data.frame(matrix(prms,ncol =length(prms),byrow = T))
    names(prms)<-names(mcmc)

    initialState <-iState(1,prms)
    prms[11]<-ifelse(prms[11]!=0,10^prms[11],0) #kappa


    pr<-as.numeric(prms)
    #initialise model
    params <-
      paramsFunc(
        gamma_1_Val = pr[1],
        # clearance rate I->S
        gamma_2_Val = pr[2],
        # clearance rate I->R
        sigma_1_Val = pr[3],
        # clearance rate E->S
        sigma_2_Val = pr[4],
        # clearance rate E->R
        mu_Val = pr[5],
        # juvenile maturation rate
        mj_Val = pr[6],
        # juvenile death rate
        m_Val = pr[7],
        # adult death rate
        omega_m_Val = pr[8],
        # maternal antibody waning rate
        omega_2_Val = pr[9],
        # immune waning rate
        epsilon_Val = 10^pr[10],
        # incubation/recurrance rate E->I
        kappa_Val = pr[11],
        # carrying capacity
        c_Val = pr[12],
        #birth pulse scalar
        s_Val = pr[13],
        #birth pulse synchronicity
        phi_Val = pr[14],
        #birth pulse timing
        rho_Val = 10^pr[15],
        R0_Val = pr[16],
        #latency I->E
        #S->E and S->I
        Phi2_val = pr[17],
        oDist_s = pr[18],
        zeta_p =pr[19],
        sigmaVer=pr[20],
        gammaVer=pr[21],
        betaVer=pr[22],
        S2_val = pr[23],
        d_val = pr[24],
        pcrProb2=pr[25],
        oDist1 = pr[26],
        oDist_u = pr[27],
        c_val2 =pr[28],
        envOscType =pr[29]
      )


    mod <-
      detSEIR(
        Sn_ini = initialState[1],
        Sj_ini = initialState[2],
        Sf_ini = initialState[3],
        Sm_ini = initialState[4],
        Ma_ini = initialState[5],
        En_ini = initialState[6],
        Ej_ini = initialState[7],
        Ef_ini = initialState[8],
        Em_ini = initialState[9],
        In_ini = initialState[10],
        Ij_ini = initialState[11],
        If_ini = initialState[12],
        Im_ini = initialState[13],
        Rn_ini = initialState[14],
        Rj_ini = initialState[15],
        Rf_ini = initialState[16],
        Rm_ini = initialState[17],
        sigma_2_Val = params$sigma_2_Val,
        gamma_2_Val = params$gamma_2_Val,
        sigma_1_Val = params$sigma_1_Val,
        gamma_1_Val = params$gamma_1_Val,
        m_Val = params$m_Val,
        rho_Val = params$rho_Val,
        R0_Val = params$R0_Val,
        phi_Val = params$phi_Val,
        s_Val = params$s_Val,
        c_Val = params$c_Val,
        omega_m_Val = params$omega_m_Val,
        omega_2_Val = params$omega_2_Val,
        kappa_Val = params$kappa_Val,
        epsilon_Val = params$epsilon_Val,
        mj_Val = params$mj_Val,
        betaVer = params$betaVer,
        gammaVer = params$gammaVer,
        sigmaVer = params$sigmaVer,
        c_val2  = params$c_val2,
        Phi2_val = params$Phi2_val,
        S2_val = params$S2_val,
        envOscType = params$envOscType,
        dt=365*4
      )
    tx = seq(currentTime, nextTime*4)


    if(fxd==T) likelihoodFunc<-likelihoodFuncBoonahFixedPr else likelihoodFuncBoonah


    #merge model output with observed data
      simDat <- as.data.frame(mod$run(tx))
      ind <- seq(1, nrow(simDat), by=4)
      simDat<- simDat[ind, ]
      simDat2<-simDat[(18453):(18771),]
      simDat3<-merge(simDat2,obsData,by.x="step",by.y="NumDays",all=T)
      simDat2<-merge(simDat2,obsData,by.x="step",by.y="NumDays")[,c(1:18)]


      obsData$total<-(obsData$positives+obsData$negatives)#(obsData$negatives+obsData$positives)
      obsData$pos<-obsData$positives#(obsData$negatives+obsData$positives)

      if(assum=="I"){
      simTotal<-rowSums(simDat2[,c(2:18)])
      simPos<-rowSums(simDat2[,c(11:14)])
      simTotal2<-rowSums(simDat3[,c(2:18)])
      simPos2<-rowSums(simDat3[,c(11:14)])
      }
      if(assum=="EI"){
        simTotal<-rowSums(simDat2[,c(2:18)])
        simPos<-rowSums(simDat2[,c(7:14)])
        simTotal2<-rowSums(simDat3[,c(2:18)])
        simPos2<-rowSums(simDat3[,c(7:14)])
      }
      if(assum=="IR"){
        simTotal<-rowSums(simDat2[,c(2:18)])
        simPos<-rowSums(simDat2[,c(11:18)])
        simTotal2<-rowSums(simDat3[,c(2:18)])
        simPos2<-rowSums(simDat3[,c(11:18)])
      }

      if(assum=="E"){
        simTotal<-rowSums(simDat2[,c(2:18)])
        simPos<-rowSums(simDat2[,c(7:10)])
        simTotal2<-rowSums(simDat3[,c(2:18)])
        simPos2<-rowSums(simDat3[,c(7:10)])
      }

      if(assum=="R"){
        simTotal<-rowSums(simDat2[,c(2:18)])
        simPos<-rowSums(simDat2[,c(14:18)])
        simTotal2<-rowSums(simDat3[,c(2:18)])
        simPos2<-rowSums(simDat3[,c(14:18)])
      }

        trajectoriesPos[l,] <- ((simPos/simTotal))


        if(l==1) {trajectoriesPosAll[l,] <- ((simPos2/simTotal2))
        simPosMed2<-cbind(simPosMed2,simPos2)
        simTotalMed2<-cbind(simTotalMed2,simTotal2)
        }

      if(l==2){simPosMed<-cbind(simPosMed,simPos)
      simTotalMed<-cbind(simTotalMed,simTotal)

      }



  }

  SDpos<-sd(obsData$pos[-1])
  SDtot<-sd(obsData$total[-1])
  simDatMedTotal<-simTotalMed
  simDatMedPos<-simPosMed
  simDatMedTotal2<-simTotalMed2
  simDatMedPos2<-simPosMed2


  boundsPos <- apply(trajectoriesPos,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))
  boundsPos2 <- apply(trajectoriesPosAll,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))


  trajectoriesPosM<-NULL
  for(v in c(1:nrow(trajectoriesPos))){
    trajectoriesPosTmp<-melt(trajectoriesPos[v,])
    trajectoriesPosTmp$variable<-obsData$Date[-1]
    trajectoriesPosTmp$id<-v
    trajectoriesPosM<-rbind(trajectoriesPosM,trajectoriesPosTmp)
  }

  ##build distribution using observed values and overdispersion parameter

  betaPrior<-mcmapply(betaShapeFinder,obsData$negatives[-1]+obsData$positives[-1], obsData$positives[-1]/(obsData$negatives[-1]+obsData$positives[-1]),params$oDist1)

  hpd <- binom.bayes(
    x = obsData$positives[-1], n = obsData$negatives[-1]+obsData$positives[-1], type = "central", conf.level = 0.9, tol = 1e-9,prior.shape1 =mean(unlist(betaPrior[1,])),prior.shape2 =mean(unlist(betaPrior[2,])) )
  print(hpd)

 # estBetaParams <- function(mu, var) {
 #   alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
 #   beta <- alpha * (1 / mu - 1)
 #   return(params = list(alpha = alpha, beta = beta))
 # }
#

#  hpd <- binom.bayes(
#    x = obsData$positives[-1], n =obsData$negatives[-1]+obsData$positives[-1], type = "central", conf.level = 0.9, tol = 1e-9,
#    prior.shape1=params$oDist1*(obsData$positives[-1]/(obsData$negatives[-1]+obsData$positives[-1])),
#    prior.shape2=params$oDist1*(1-(obsData$positives[-1]/(obsData$negatives[-1]+obsData$positives[-1]))))
#

  trajectoriesPosM$prev<-obsData$prev[-1]
  trajectoriesPosM$lower<-hpd$lower
  trajectoriesPosM$upper<-hpd$upper

  gx<-ggplot(trajectoriesPosM, aes(x = variable, y = value)) + geom_point(color = "#F77F00",alpha=0.05,size=3)+
    geom_point(aes(y=prev),col="#102F47",size=3)+
    xlab("Date")+
    ylab("Seropositive prevalence")+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
    position=position_dodge(0.05),alpha=0.01,colour="#102F47")+theme_bw(base_size = 20)




  fullData<-setNames(as.data.frame(seq(as.Date("2013-07-22"), as.Date("2014-06-05"), "days")),c("Date"))
  fullData$prev<-simDatMedPos2[-1]/simDatMedTotal2[-1]
  fullData$uppr<-boundsPos2[2,-1]
  fullData$lwr<-boundsPos2[1,-1]

  obsData2<-obsData[-1,]
  obsData2$biUp<-hpd$upper
  obsData2$bilow<-hpd$lower
  fullData<-merge(fullData,obsData2,by.x="Date",by.y="Date",all=T)

  g0<-ggplot(fullData, aes(y=prev.y, x=Date))+
    geom_point(colour="#102F47")+
    geom_errorbar(aes(ymin=bilow, ymax=biUp), width=.2,
                  position=position_dodge(0.05),alpha=1,colour="#102F47")+
    geom_line(aes(y=prev.x),se=F,color="#F77F00",alpha=0.8)+
    geom_ribbon(aes(ymin=lwr.x,ymax=uppr),alpha=0.3)+
    ylab("Sero+")+theme_bw(base_size = 20)


  ######################################################################
  #########plot individual PCR simulation against observed##############
  ######################################################################

  pcrAgg<-aggregate(obsData$pcrPos~obsData$Date,FUN=sum)
  pcrAgg$neg<-aggregate(obsData$neg~obsData$Date,FUN=sum)[,2]
  pcrAgg$inPos<-aggregate(obsData$positives~obsData$Date,FUN=sum)[,2]
  pcrAgg$fit<-aggregate(obsData$fit~obsData$Date,FUN=mean)[,2]
  names(pcrAgg)<-c("Date","pcrPos","neg","inPos","fit")
  pcrAgg$prev<-pcrAgg$pcrPos/(pcrAgg$inPos+pcrAgg$neg)
  pcrAgg<-pcrAgg[-1,]

  if(fxd==T)params$zeta_p<-0.9999
  tt<-yday(pcrAgg$Date)/365
        sDrive<-1
        if(osc==T){
          s=params$S2_val
          c=1
          phi=params$Phi2_val
          sDrive <-c* exp(-s*(cos(pi*tt- phi))^2)
        }
        simDat4<-simDat3[is.na(simDat3$Date)==F,]

        Ia<-as.vector(round(rowSums(simDat4[-1,c(11:14)])))
        Ia[is.na(Ia)]<-0

        prb<-NULL
        prb2<-NULL
        for(i in c(1:100)){
        prb<- mapply(rmutil::rbetabinom, size = Ia, m=params$zeta_p*(sDrive),MoreArgs = list(n=1,s=params$oDist1))/as.vector(round(rowSums(simDat4[-1,c(2:18)])))
        prb2<-rbind(prb,prb2)
        }

  fullData2<-as.data.frame(obsData$Date[-1])
  fullData2$resMean<-colMeans(prb2)
  names(fullData2)<-c("Date","resMean")

  biCon2<-binom.bayes(pcrAgg$pcrPos,(pcrAgg$inPos+pcrAgg$neg))

  pcrAgg$biConLwr<-biCon2$lower
  pcrAgg$biConUpr<-biCon2$upper
  fullData3<-merge(fullData2,pcrAgg,by.x="Date",by.y="Date",all=T)
  boundsPCR1<- apply(prb2,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))

  g1<-ggplot(fullData3, aes(y=prev, x=Date)) +
    geom_point(colour="#102F47")+
    ylab("PCR+ (estimated prev)")+
    geom_line(aes(y=resMean),color="#F77F00",alpha=1,se=F)+
    ylab("indv. prevalence")+
    geom_ribbon(aes(ymin=boundsPCR1[1,],ymax=boundsPCR1[2,]),alpha=0.3)+
    geom_errorbar(aes(ymin=biConLwr, ymax=biConUpr), width=.2,
                  position=position_dodge(0.05),alpha=0.5,colour="#102F47")+theme_bw(base_size = 20)

  prb2Melt<-as.data.frame(as.vector(t(prb2)))
  prb2Melt$Date<-obsData$Date[-1]
  names(prb2Melt)<-c("value","Date")
  fullDataX<-merge(prb2Melt,pcrAgg,by.x="Date",by.y="Date",all=F)

  g1All<- ggplot(fullDataX,aes(x=Date,y=prev))+
    geom_point(aes(y=value, x=Date),col="#F77F00",se=F,alpha=0.05,size=3)+
    geom_errorbar(aes(ymin=biConLwr, ymax=biConUpr), width=.2,
                  position=position_dodge(0.05),alpha=0.01,colour="#102F47")+
    geom_point(col="#102F47",alpha=0.01,size=3)+
    ylab("Individual urine RNA prevalence")+theme_bw(base_size = 20)
  #######################################################################################################################################
  #####################################################################################################################

  ######################################################################
  #########plot UR PCR simulation against observed##############
  ######################################################################

  if(fxd==T)params$pcrProb2<-0.9999
  tt<-yday(pcrAgg$Date)/365
  sDrive<-1
  if(osc==T){
    s=params$S2_val
    c=1
    phi=params$Phi2_val
    sDrive <-c* exp(-s*(cos(pi*tt- phi))^2)
  }

  Ia<-as.vector(round(rowSums(simDat4[-1,c(11:14)])))
  Ia[is.na(Ia)]<-0


  Np=46 #number of pooled samples
  x=rep(round(1+params$d_val),length(Ia)) #number of bats contributing to each pool
  Pnp= pcrAgg[,5]  #pool prevelance - prevelance in all pools
  ki=x*Np*(1-(1-Pnp)^(1/x)) #indv prevelance reflective of infection in each bat, or probability of shedding?
  #this assumes no variation in diagnostic tests with pool size

  sDrive<-1
  if(osc==T){
    tt<-(yday(pcrAgg$Date))/365
    s=params$S2_val
    c=1
    phi=params$Phi2_val
    sDrive <-c* exp(-s*(cos(pi*tt- phi))^2)
  }

  prbUR<-((Ia*params$pcrProb2*(sDrive))/rowSums(simDat2[,c(2:18)])) #each simulated bat has a probability of shedding prb = prevelance based on random shedding events
  Npx=x*Np #number bats contributing to pooled samples
  prbUR<-ifelse(prbUR>=1,0.99,prbUR)
  prbUR<-ifelse(prbUR<=0,0.00001,prbUR)


  prbx<-NULL
  prb2x<-NULL
  toteski<-NULL
  totesN<-NULL
  for (o in c(1:100)){
    prbx<-rmutil::rbetabinom(1,size=Npx,m=prbUR,s=params$oDist1)/Npx
  prb2x<-rbind(prbx,prb2x)

  toteski<-rbind(ki,toteski)
  totesN<-rbind(Npx,totesN)
  }


  res2<-data.frame(t(prb2x))
  res2$Date<-pcrAgg$Date
  resXS2<-as.data.frame(t(res2[,-ncol(res2)]))
  boundsPCR2<- apply(resXS2,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))


  res2X<-data.frame(colMeans(prb2x))
  res2X$totesKi<-colMeans(toteski)
  res2X$totesN<-colMeans(totesN)
  res2X$Date<-obsData$Date[-1]


  biCon3 <- binom.bayes(
    x = res2X$totesKi, n = res2X$totesN, type = "central", conf.level = 0.95, tol = 1e-9)
  print(biCon3)

  g2<- ggplot(res2X,aes(x=Date,y=totesKi/totesN))+
    geom_errorbar(aes(ymin=biCon3$lower, ymax=biCon3$upper), width=.2,
                  position=position_dodge(0.05),alpha=0.5,colour="#102F47")+
    geom_point(col="#102F47",alpha=0.7)+
    ylab("Under-roost urine prevalence")+
    geom_smooth(aes(y=res2X[,1], x=Date),col="#F77F00",se=F,alpha=0.9)+
    geom_ribbon(aes(ymin=boundsPCR2[1,],ymax=boundsPCR2[2,]),alpha=0.3)+
    theme_bw(base_size = 18)+theme_bw(base_size = 20)

  resXS2all<-melt(t(prb2x))

  resXS2all$Date<-pcrAgg$Date
  resXS2all$totes<-res2X$totesKi/res2X$totesN
  resXS2all$low<-biCon3$lower
  resXS2all$up<-biCon3$upper

  g2All<- ggplot(resXS2all,aes(x=Date,y=totes))+
    geom_point(aes(y=value, x=Date),col="#F77F00",se=F,alpha=0.05,size=3)+
    geom_errorbar(aes(ymin=low, ymax=up), width=.2,
                  position=position_dodge(0.05),alpha=0.01,colour="#102F47")+
    geom_point(col="#102F47",alpha=0.01,size=3)+
    ylab("under roost prevalence")+theme_bw(base_size = 20)
  #######################################################################################################################################
  #####################################################################################################################

  print(ggarrange(g0,g1,g2))
  print(ggarrange(gx,g1All,g2All,ncol=1))
}

#modPlotFunc(gg,assum="I",osc=T,fxd=T)


plotDistsFunc<-function(aicvl,llT=F){
  gf<-NULL
 indx<-aicvl[order(aicvl$V1),2]
for(i in c(indx[1:10])){

  print(i)
  gg<-read.csv(paste0("/Users/alm204/Documents/Cambridge/results/feb_10/res_",i,".csv"))[,-1]
  gg<-subset(gg,is.na(gamma_1_Val)!=T)
  gg<-gg[,-(which(colSums(gg)==0))]
  gg<-if (mean(gg$zeta_p)==0.5480735) gg[ , -which(names(gg) %in% c("zeta_p","pcrProb2","sigmaVer", "gammaVer", "betaVer"))] else  gg[ , -which(names(gg) %in% c("sigmaVer", "gammaVer", "betaVer"))]
  if(prmFile$lFunc[i]!="Osc"&&prmFile$lFunc[i]!="prFxOsc") gg<- gg[ , -which(names(gg) %in% c("c_val2","d_val","Phi2_val"))]
  #modPlotFunc(gg)
  aicval<--2*(gg$ll) + 2*ncol(gg)
  if(llT==T) aicval<-gg$ll
  gf<-cbind(gf,aicval)
}
gf<-as.data.frame(gf)
names(gf)<-(indx[1:ncol(gf)])
dataG<- melt(gf)
print(ggplot(dataG,aes(x=value, fill=variable)) + geom_density(alpha=0.3,adjust=2)+
  guides(fill=guide_legend(title="Model Number"))+
  xlab("AIC")+
  scale_fill_viridis_d(option="C",direction=-1)+
  theme_bw(base_size = 20))
}




plotDistsFunc<-function(aicvl,parmNme="R0_Val"){
  gf<-NULL
  indx<-aicvl[order(aicvl$V1),2]
  for(i in c(indx[1:10])){

    print(i)
    gg<-read.csv(paste0("/Users/alm204/Documents/Cambridge/results/feb_10/res_",i,".csv"))[,-1]
    gg<-subset(gg,is.na(gamma_1_Val)!=T)
    gf<-cbind(gf,gg[ ,which(names(gg)==parmNme)])
  }
  gf<-as.data.frame(gf)
  biCon<-as.data.frame(HPDinterval(as.mcmc(gf)))
  biCon$med<-colMedians(gf)
  biCon$indx<-c(1:ncol(gf))
  print(ggplot(biCon,aes(y=med,x=indx))+
          geom_errorbar(aes(ymin=biCon$lower, ymax=biCon$upper), width=.2,
                        position=position_dodge(0.05),alpha=0.5,colour="blue")+
          geom_point(col="black",alpha=0.7)+
          xlab("Model rank")+
          ylab(parmNme)+
          scale_fill_viridis_d(option="C",direction=-1)+
          theme_bw(base_size = 20))
}


##calc DIC
calc.dic <- function(gg,prmFil) {
  lik<-gg$ll
  gg<-subset(gg,is.na(gamma_1_Val)!=T)
  gg<-gg[,-(which(colSums(gg)<=0))]
  gg<-if(prmFil$lFunc=="prFx"&&prmFil$lFunc=="prFxOsc") gg[ , -which(names(gg) %in% c("zeta_p","pcrProb2","sigmaVer", "gammaVer", "betaVer"))] else  gg[ , -which(names(gg) %in% c("sigmaVer", "gammaVer", "betaVer"))]
  if(prmFil$lFunc!="Osc"&&prmFil$lFunc!="prFxOsc") gg<- gg[ , -which(names(gg) %in% c("c_val2","d_val","Phi2_val"))]
  #modPlotFunc(gg)

  D.bar <- -2*mean(lik)
  theta.bar <- colMeans(gg)
  theta.bar<-as.data.frame(matrix(theta.bar,ncol =length(theta.bar),byrow = T))
  names(theta.bar)<-names(gg)
  D.hat <- -2*detModFunc(theta.bar,obsData,assump=prmFil$assump,likelihoodFuncBoonah)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}


##calc LOO
calc.loo <- function(data=obsData,gg=gg,assump=assump,nDraws,likelihoodFunc,birthType) {
  library(doParallel)

  gg1<-gg[gg$chainID==1,]
  gg2<-gg[gg$chainID==2,]

  post<-gg1#[sample(nrow(gg1), nDraws), ]
  post2<-gg2#[sample(nrow(gg2), nDraws), ]

# llMat<-pFiltMat(data=data,gg=post,assump=assump,likelihoodFunc=likelihoodFuncBoonahStoch,birthType=birthType,nDraws=nDraws)
# llMat2<-pFiltMat(data=data,gg=post2,assump=assump,likelihoodFunc=likelihoodFuncBoonahStoch,birthType=birthType,nDraws=nDraws)
#

  post$mu_Val<-as.numeric(as.character(post$mu_Val))
  post2$mu_Val<-as.numeric(as.character(post2$mu_Val))

  llMat<-detModFuncLOO(data,post,assump,likelihoodFunc,birthType)
  llMat2<-detModFuncLOO(data,post2,assump,likelihoodFunc,birthType)

  nlength<-ifelse(nrow(llMat)<nrow(llMat2),nrow(llMat),nrow(llMat2))


  llMat10<-rbind(llMat,llMat2)[1:(2*nlength),]

  rel_n_eff <- relative_eff(exp(llMat10),chain_id = c(rep(1,nlength),rep(2,nlength)) )

  #llMat<-detModFuncLOO_ind(data[2,],post,assump,likelihoodFunc)
  loo_3<-loo(llMat10, r_eff = rel_n_eff, cores = 11)
  #loo_3 <- loo(detModFuncLOO_ind, data = obsData[-1,], draws = post, r_eff = NA,likelihoodFunc=likelihoodFunc,assump=assump,birthType=birthType,cores=10)
#  return(loo(llMat, r_eff = rel_n_eff, cores = 11))
return(list(loo_3,llMat,rel_n_eff))
}

##calc LOO
calc.looStoch <- function(data=obsData,gg=gg,assump=assump,nDraws,likelihoodFunc,birthType) {

  gg1<-gg[gg$chainID==1,]
  gg2<-gg[gg$chainID==2,]

    llMat<-pFiltMat(data=data,gg=gg1,assump=assump,likelihoodFunc=likelihoodFuncBoonahStoch,birthType=birthType,nDraws=nDraws)
    llMat2<-pFiltMat(data=data,gg=gg2,assump=assump,likelihoodFunc=likelihoodFuncBoonahStoch,birthType=birthType,nDraws=nDraws)

     llMat<-rbind(gg1[,33:48],gg2[,33:48])


    rel_n_eff <- relative_eff(exp(llMat),chain_id=c(rep(1,nrow(gg1)),rep(2,nrow(gg2))))

  loo3<- loo(llMat, r_eff = rel_n_eff, cores = 11)

  return(list(loo3,llMat,rel_n_eff))
}



##plot environmental force
plotEnvForc<-function(){

  for (i in c(4,6,8)){
    prmFile<-read.csv("/Users/alm204/Documents/ModelSetups.csv")

    print(i)
    gg<-read.csv(paste0("/Users/alm204/Documents/Cambridge/results/jun3/res_",i,".csv"))[,-1]
    gg<-subset(gg,is.na(betaVer)!=T)
    if (prmFile$lFunc[i]=="prFxOsc") gg$zeta_p<-1
    if (prmFile$lFunc[i]=="prFxOsc") gg$pcrProb2<-1
    gg<-subset(gg,ll>-500)
    gg<-gg[94998:nrow(gg),]

    s=median(gg$S2_val)
    phi=median(gg$Phi2_val)
    c=median(gg$c_val2)
    x<-seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")
    x<-yday(x)/365

    if(i == 4) {


      medVal<-10^median(gg$epsilon_Val)
    ciValsEp<- 10^HPDinterval(as.mcmc(gg$epsilon_Val))

    sDriveEpOrig <-c* exp(-s*(cos(pi*x- phi))^2)
    sDriveEp<-medVal*sDriveEpOrig
    sDriveEp<-as.data.frame(sDriveEp)
    sDriveEp$date<-seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

    G4<- ggplot(sDriveEp,aes(y=sDriveEp,x=date))+
      geom_line(colour="purple",size=1.5)+
      geom_rect(data=head(sDriveEp),aes(xmin = as.Date("2013-06-01"), xmax = as.Date("2013-08-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
      geom_rect(data=head(sDriveEp),aes(xmin = as.Date("2014-05-01"), xmax = as.Date("2014-06-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
      theme_bw(base_size = 20)+
      scale_x_date(date_breaks = "2 months" , date_labels = "%b")+
      ylab(expression(epsilon_t~(year^-1)))+
      geom_ribbon(aes(ymin=sDriveEpOrig*ciValsEp[1],ymax=sDriveEpOrig*ciValsEp[2]),fill="purple",alpha=0.3)

    }
    if(i == 6){
      medValGam<-10^median(gg$gamma_2_Val)

      ciValsGam<- 10^HPDinterval(as.mcmc(gg$gamma_2_Val))

      sDriveGamOrig <-c* exp(-s*(cos(pi*x- phi))^2)
      sDriveGam<-medValGam*sDriveGamOrig
      sDriveGam<-as.data.frame(sDriveGam)
      sDriveGam$date<-seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G6<- ggplot(sDriveGam,aes(y=sDriveGam,x=date))+
        geom_line(colour="purple",size=1.5)+
        geom_rect(data=head(sDriveGam),aes(xmin = as.Date("2013-06-01"), xmax = as.Date("2013-08-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
        geom_rect(data=head(sDriveGam),aes(xmin = as.Date("2014-05-01"), xmax = as.Date("2014-06-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
        theme_bw(base_size = 20)+
        scale_x_date(date_breaks = "2 months" , date_labels = "%b")+
        ylab(expression(gamma_t~(year^-1)))+
        geom_ribbon(aes(ymin=sDriveGamOrig*ciValsGam[1],ymax=sDriveGamOrig*ciValsGam[2]),fill="purple",alpha=0.3)

    }


    if(i == 8) {
      medValOm<-median(gg$omega_2_Val)

      ciValsOmeg<- HPDinterval(as.mcmc(gg$omega_2_Val))

      sDriveOmegOrig <-c* exp(-s*(cos(pi*x- phi))^2)
      sDriveOmeg<-medValOm*sDriveOmegOrig
      sDriveOmeg<-as.data.frame(sDriveOmeg)
      sDriveOmeg$date<-seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G8<- ggplot(sDriveOmeg,aes(y=sDriveOmeg,x=date))+
        geom_line(colour="purple",size=1.5)+
        geom_rect(data=head(sDriveOmeg),aes(xmin = as.Date("2013-06-01"), xmax = as.Date("2013-08-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
        geom_rect(data=head(sDriveOmeg),aes(xmin = as.Date("2014-05-01"), xmax = as.Date("2014-06-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
        theme_bw(base_size = 20)+
        scale_x_date(date_breaks = "2 months" , date_labels = "%b")+
        ylab(expression(omega_t~(year^-1)))+
        geom_ribbon(aes(ymin=sDriveOmegOrig*ciValsOmeg[1],ymax=sDriveOmegOrig*ciValsOmeg[2]),fill="purple",alpha=0.3)

    }


  }

  weath<-data.frame(month=as.Date(c("2013-07-01","2013-08-01","2013-09-01","2013-10-01","2013-11-01","2013-12-01","2014-01-01","2014-02-01","2014-03-01","2014-04-01","2014-05-01","2014-06-01")),
                    val=as.numeric(c("21","22","25","27","29","30","31","30","28","26","23","21")),Mn=c(1:12))

  gW<- ggplot(weath,aes(x=month,y=val))+
    ylab("Average Temp. Celsius")+
    xlab("Date")+
    geom_rect(data=head(weath),aes(xmin = as.Date("2013-06-01"), xmax = as.Date("2013-08-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
    geom_rect(data=head(weath),aes(xmin = as.Date("2014-05-01"), xmax = as.Date("2014-06-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
    theme_bw(base_size = 20)+
    ylim(0,40)+
    geom_bar(stat = 'identity',fill="#102F47",width=15,alpha=0.8)+
    ggtitle("Boonah Temperature")+
    scale_x_date(date_breaks = "2 months" , date_labels = "%b")


  ggRt<- ggarrange(G4+ggtitle("SILI (Maternal Immunity)"),G6+ggtitle("SIR"),G8+ggtitle("SIRS"),gW)
  print(ggRt)
  ggsave("/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/envForceVals.png",plot = ggRt,width = 18, height = 13)

}

##plot environmental force
#getR0<-function(params,rev=F){
#  s=params$S2_val
#  phi=params$Phi2_val
#  c=params$c_val2
#  x<-seq(as.Date("2013-07-22"), as.Date("2014-07-21"), "days")
#  x<-yday(x)/365
#  sDrive <-c* exp(-s*(cos(pi*x- phi))^2)
#  epsilon<-10^params$epsilon_Val*sDrive
#  R0_Val<-params$R0_Val
#  sigma_2_Val<-params$sigma_2_Val
#  m_Val<-params$m_Val
#  N<-params$kappa_Val
#  rho_Val<-params$rho_Val
#  gamma_2_Val<-params$gamma_2_Val
#
#  beta<-R0_Val*((epsilon+sigma_2_Val+m_Val)*(gamma_2_Val+m_Val+rho_Val)
#          -epsilon*rho_Val)/(N*(epsilon+sigma_2_Val+m_Val))
#
#  ((beta*N)*(epsilon+m_Val))/((epsilon+m_Val)*(gamma_2_Val+m_Val+rho_Val)-epsilon*rho_Val)
#
#
#  sDrive<-(params$R0_Val*sDrive)
#  sDrive<-as.data.frame(sDrive)
#  sDrive$date<-seq(as.Date("2013-07-22"), as.Date("2014-07-21"), "days")
#  if(rev==T ) plot(sDrive,type="l",xlab="time",ylab="sDrive")
#  else ggplot(sDrive,aes(y=sDrive,x=date))+
#    geom_line(colour="#F77F00",size=1.5)+
#    ylab(expression(nu))+
#    theme_bw(base_size = 20)
#
#}
#


#loicVals<-data.frame(Looic=rnorm(27,200,2))
#loicVals$se<-round(rnorm(27,5,2))
#loicVals$modType<-c(rep("SILI",9),rep("SIR",9),rep("SIRS",9))
#loicVals$Structure<-rep(c("std","fxd","osc"),9)
#library(cowplot)
plotLooic<-function(fileNme,loicVals,loicVals2=NULL,type="loo"){
  loicVals$Structure[loicVals$Structure==1]<-"EF"
  loicVals$Structure[loicVals$Structure==2]<-"No EF"
  loicVals$Structure[loicVals$Structure==3]<-"DC"
  loicVals$modType[loicVals$modType==3]<-"SIRS"
  loicVals$modType[loicVals$modType==2]<-"SIR"
  loicVals$modType[loicVals$modType==1]<-"SILI"

  loicVals2$Structure[loicVals2$Structure==1]<-"EF"
  loicVals2$Structure[loicVals2$Structure==2]<-"No EF"
  loicVals2$Structure[loicVals2$Structure==3]<-"DC"
  loicVals2$modType[loicVals2$modType==3]<-"SIRS"
  loicVals2$modType[loicVals2$modType==2]<-"SIR"
  loicVals2$modType[loicVals2$modType==1]<-"SILI"

 if(type=="loo") {loicVals$seUp<-loicVals$Looic+loicVals$se
  loicVals$seLow<-loicVals$Looic-loicVals$se}


  loicVals$location<-c(4,1,5,2,7,8,10,11)
  loicVals2$location<-c(4,1,5,2,7,8,10,11)

  if(type=="loo") yName<- "Looic"
  if(type=="R0") yName<-expression("R"[0])
  if(type=="batCont") yName<-"Bats contributing to pooled samples (1+)"
  if(type=="loo") yLim<- c(200,350)
  if(type=="R0") yLim<-c(0,50)
  if(type=="batCont") yLim<-c(0,25)
  if(type=="zeta") yLim<-c(0,1)
  if(type=="zeta") yName<- expression(zeta[s])
  if(type=="zeta") yName2<- expression(zeta[p])


   g1<-ggplot(loicVals,aes(y=Looic,x=location))+
     ylim(yLim)+
          geom_errorbar(aes(ymin=seLow, ymax=seUp), width=0.2,
                        position=position_dodge(0.05),alpha=1,colour="#102F47")+
          geom_point(aes(colour=Structure),size=5,alpha=1)+
          xlab("")+
     ylab(yName)+
          scale_fill_viridis_d(option="C",direction=-1)+
      scale_x_continuous(breaks=c(1.5,4.5,7.5,10.5),
                       labels=c("SILI (Maternal immunity)", "SILI", "SIR","SIRS"))+  theme_bw(base_size = 20)+
     panel_border( size = 0.5,color = "black")+
     scale_color_viridis_d(option="D")+
     geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="lightskyblue")+
   geom_rect(aes(xmin = 3, xmax = 6, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="darkorange")+
   geom_rect(aes(xmin = 6, xmax = 9, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="lightskyblue")+
   geom_rect(aes(xmin = 9, xmax = 12, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="darkorange")

     ggsave(fileNme,plot=g1,height=25,width=40,units="cm")


if(type=="zeta"){
g2<-ggplot(loicVals2,aes(y=Looic,x=location))+
  ylim(yLim)+
  geom_errorbar(aes(ymin=seLow, ymax=seUp), width=0.2,
                position=position_dodge(0.05),alpha=1,colour="#102F47")+
  geom_point(aes(colour=Structure),size=5,alpha=1)+
  xlab("")+
  ylab(yName2)+
  scale_fill_viridis_d(option="C",direction=-1)+
  scale_x_continuous(breaks=c(1.5,4.5,7.5,10.5),
                     labels=c("SILI (Maternal immunity)", "SILI", "SIR","SIRS"))+  theme_bw(base_size = 20)+
  panel_border( size = 0.5,color = "black")+
  scale_color_viridis_d(option="D")+
  geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="lightskyblue")+
  geom_rect(aes(xmin = 3, xmax = 6, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="darkorange")+
  geom_rect(aes(xmin = 6, xmax = 9, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="lightskyblue")+
  geom_rect(aes(xmin = 9, xmax = 12, ymin = -Inf, ymax = Inf), alpha = 0.01,fill="darkorange")
print(ggarrange(g1,g2,common.legend = T))
ggsave(fileNme,plot=ggarrange(g1,g2,common.legend = T),height=25,width=40,units="cm")
 return(ggarrange(g1,g2,common.legend = T))
}
}

#plotLooic(r0Val,type="R0")
#plotLooic(contBats,type="batCont")
#
#plotLooic(looicVal,type="loo")
#fNm<-"/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/"
#plotLooic(fileNme=paste0(fNm,"zetaPlots.png"),loicVals=zetaU,loicVals2=zetaP,type="zeta")
#plotLooic(fileNme=paste0(fNm,"looicVals.png"),loicVals=looicVal,type="loo")
#plotLooic(fileNme=paste0(fNm,"contBatsX.png"),loicVals=contBats,type="batCont")




