
#' plot model from mcmc results
#' @param mcmc mcmc results
#' @return figures
modPlotFuncContinuous <- function(mcmc,assum="EI",fxd=F,iters=100,pName,birthType,Rt=F) {
  #parse input data
  simALLMed<-NULL
  simTotalMed<-NULL
  simDatMed<-NULL
  simPosMed<-NULL
  simPosMed2<-NULL
  simTotalMed2<-NULL

  trajectoriesPos <- matrix(nrow=iters,ncol=366)
  trajectoriesPOP <- matrix(nrow=iters,ncol=366)
  r0<-matrix(nrow=iters,ncol=366)

  trajectoriesPcrPos <- matrix(nrow=iters,ncol=366)

  trajectoriesAll <- matrix(nrow=iters,ncol=16)
  trajectoriesPosAll <- matrix(nrow=iters,ncol=320)

  currentTime <-0
  nextTime <- 18771

  prms<-colMedians(mcmc)
  prms<-as.data.frame(matrix(prms,ncol =length(prms),byrow = T))
  names(prms)<-names(mcmc)

  prms$betaFXVal<- betaMean(prms)

  initialState <-iState(1,prms)
  prms[11]<-ifelse(prms[11]!=0,10^prms[11],0) #kappa
  prms[10]<-ifelse(prms[10]!=0,10^prms[10],0) #kappa
  prms[15]<-ifelse(prms[15]!=0,10^prms[15],0) #kappa
  prms[2]<-ifelse(prms[2]!=0,10^prms[2],0) #kappa


  pr<-as.numeric(prms)
  #initialise model
  params <-
    paramsFunc(
      gamma_1_Val = pr[1],
      # clearance rate I->S
      gamma_2_Val = pr[2],
      # clearance rate I->R
      zeta_s = pr[3],
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
      epsilon_Val = pr[10],
      # incubation/recurrance rate E->I
      kappa_Val = pr[11],
      # carrying capacity
      c_Val = pr[12],
      #birth pulse scalar
      s_Val = pr[13],
      #birth pulse synchronicity
      phi_Val = pr[14],
      #birth pulse timing
      rho_Val = pr[15],
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
      envOscType =pr[29],
      betaFX=prms$betaFX,
      betaFXVal=prms$betaFXVal
    )


  mod <-
    stochSEIR(
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
      zeta_s = params$zeta_s,
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
      Phi2_val = params$Phi2_val,
      S2_val  = params$S2_val,
      c_val2  = params$c_val2,
      envOscType = params$envOscType,
      dt=365*4,
      birthType = birthType,
      betaFX=prms$betaFX,
      betaFXVal=prms$betaFXVal
    )
  tx = seq(currentTime, nextTime*4)



  #merge model output with observed data
  for (l in c(1:iters)){
    simDat <- as.data.frame(mod$run(tx))
    ind <- seq(1, nrow(simDat), by=4)
    simDat<- simDat[ind, ]
    simDat2<-simDat[(18406):(18770),]

    r0[l,]<-simDat3$R0_out


    simDat3<-merge(simDat2,obsData,by.x="step",by.y="NumDays",all=T)
    simDat2<-merge(simDat2,obsData,by.x="step",by.y="NumDays")[,c(1:18)]


    obsData$total<-(obsData$positives+obsData$negatives)#(obsData$negatives+obsData$positives)
    obsData$pos<-obsData$positives#(obsData$negatives+obsData$positives)

    I<-rowSums(simDat3[,c(11:14)])
    E<-rowSums(simDat3[,c(7:10)])
    R<-rowSums(simDat3[,c(6,15:18)])
    S<-rowSums(simDat3[,c(2:5)])

    PpSp<-I*params$zeta_p*params$zeta_s
    PnSp<-(I+E+R)*params$zeta_s*(1-params$zeta_p)
    PpSn<-I*params$zeta_p*(1-params$zeta_s)
    PnSn<-I*(1-params$zeta_p)*(1-params$zeta_s) +S

    simPos<-PpSp+PnSp
    simPCRPos<-PpSp+PpSn
    simTotal<-PpSp+PnSp+PpSn+PnSn

    trajectoriesPcrPos[l,]<-(simPCRPos/simTotal)
    trajectoriesPos[l,] <- ((simPos/simTotal))
    trajectoriesPOP[l,]<-simTotal
  }



  SDpos<-sd(obsData$pos[-1])
  SDtot<-sd(obsData$total[-1])


  r0M<-melt(t(r0))
  r0M$Date<-seq(as.Date("2013-06-03"), as.Date("2014-06-03"), "days")



  gRT<-ggplot(r0M, aes(x = Date, y = value,group=Var2)) +
    geom_line(col="#F77F00",size=.5,alpha=.2)+
    xlab("Date")+
    ylab(expression(R[t]))+
    theme_bw(base_size = 20)+
    geom_rect(data=head(r0M),aes(xmin = as.Date("2013-06-04"), xmax = as.Date("2013-08-01"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")+
    geom_rect(data=head(r0M),aes(xmin = as.Date("2014-05-01"), xmax = as.Date("2014-06-03"), ymin = -Inf, ymax = Inf), alpha = 0.05,fill="lightskyblue")


  print(gRT)
  if(Rt==T) return(gRT)

  if(Rt!=T){

    betaPrior<-mcmapply(betaShapeFinder,obsData$negatives[-1]+obsData$positives[-1], obsData$positives[-1]/(obsData$negatives[-1]+obsData$positives[-1]),10^params$oDist1)

    hpd <- binom.bayes(
      x = obsData$positives[-1], n = obsData$negatives[-1]+obsData$positives[-1], type = "central", conf.level = 0.9, tol = 1e-9,prior.shape1 =mean(unlist(betaPrior[1,])),prior.shape2 =mean(unlist(betaPrior[2,])) )
    print(hpd)

    obsDataX<-obsData[-1,]
    # obsDataX$simPop<-colMeans(trajectoriesPOP)
    # obsDataX<-merge(obsDataX,s,all=T)

    obsDataX$lower<-hpd$lower
    obsDataX$upper<-hpd$upper
    biCon2<-binom.bayes(obsDataX$pcrPos,(obsDataX$positives+obsDataX$negatives))

    obsDataX$PCRLwr<-biCon2$lower
    obsDataX$PCRUpr<-biCon2$upper

    trajectoriesPosM<-melt(t(trajectoriesPos))
    trajectoriesPosM$population<-melt(t(trajectoriesPOP))[,3]
    trajectoriesPosM$pcrPrevSim<-melt(t(trajectoriesPcrPos))[,3]

    trajectoriesPosM$Date<-seq(ymd('"2013-06-18"'),ymd('2014-06-18'), by = '1 day')
    trajectoriesPosM<-merge(trajectoriesPosM,obsDataX,by.x="Date",by.y="Date",all=T)

    trajectoriesPosM$pcrPrev<-trajectoriesPosM$pcrPos/(trajectoriesPosM$positives+trajectoriesPosM$negatives)






    gx<-ggplot(trajectoriesPosM, aes(x = Date, y = value,group=Var2)) +
      geom_line(color = "#F77F00",alpha=0.1,size=1)+
      geom_point(aes(y=prev),col="#102F47",size=3,alpha=0.01)+
      xlab("Date")+
      ylab("Seropositive prevalence")+
      geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                    position=position_dodge(0.05),alpha=0.01,colour="#102F47")+theme_bw(base_size = 20)

    g1All<-ggplot(trajectoriesPosM, aes(x = Date, y = pcrPrevSim,group=Var2)) +
      geom_line(color = "#F77F00",alpha=0.1,size=1)+
      geom_point(aes(y=pcrPrev),col="#102F47",size=3,alpha=0.01)+
      xlab("Date")+
      ylab("Indv. PCR prevalence")+
      geom_errorbar(aes(ymin=PCRLwr, ymax=PCRUpr), width=.2,
                    position=position_dodge(0.05),alpha=0.01,colour="#102F47")+theme_bw(base_size = 20)


    #######################################################################################################################################
    #####################################################################################################################

    ######################################################################
    #########plot UR PCR simulation against observed##############
    ######################################################################

    if(fxd==T)params$pcrProb2<-0.9999
    prbx<-NULL
    prb2x<-NULL
    toteski<-NULL
    totesN<-NULL

    for (o in c(1:iters)){

      simDat <- as.data.frame(mod$run(tx))
      ind <- seq(1, nrow(simDat), by=4)
      simDat<- simDat[ind, ]
      simDat2<-simDat[(18313):(18771),]
      simDat3<-merge(simDat2,obsData,by.x="step",by.y="NumDays",all=T)
      simDat2<-merge(simDat2,obsData,by.x="step",by.y="NumDays")[,c(1:18)]


      Ia<-as.vector(round(rowSums(simDat2[,c(11:14)])))
      Ia[is.na(Ia)]<-0

      Np=46 #number of pooled samples
      x=rep(round(1+params$d_val),length(Ia)) #number of bats contributing to each pool
      Pnp= obsData[-1,9]  #pool prevelance - prevelance in all pools
      ki=x*Np*(1-(1-Pnp)^(1/x)) #number of infected bats contributing to the pooled samples
      #this assumes no variation in diagnostic tests with pool size
      prbUR<-((Ia*params$pcrProb2)/rowSums(simDat2[,c(2:18)])) #each simulated bat has a probability of shedding prb = prevelance based on random shedding events
      Npx=x*Np #number bats contributing to pooled samples
      prbUR<-ifelse(prbUR>=1,0.99,prbUR)
      prbUR<-ifelse(prbUR<=0,0.00001,prbUR)


      prbx<-rmutil::rbetabinom(1,size=Npx,m=prbUR,s=10^params$oDist1)/Npx
      prb2x<-rbind(prbx,prb2x)

      toteski<-rbind(ki,toteski)
      totesN<-rbind(Npx,totesN)
    }

    res2<-data.frame(t(prb2x))
    res2$Date<-obsData$Date[-1]
    resXS2<-as.data.frame(t(res2[,-ncol(res2)]))
    boundsPCR2<- apply(resXS2,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))


    res2X<-data.frame(colMeans(prb2x))
    res2X$totesKi<-colMeans(toteski)
    res2X$totesN<-colMeans(totesN)
    res2X$Date<-obsData$Date[-1]


    biCon3 <- binom.bayes(
      x = res2X$totesKi, n = res2X$totesN, type = "central", conf.level = 0.95, tol = 1e-9)
    print(biCon3)


    resXS2all<-melt(t(prb2x))

    resXS2all$Date<-obsData$Date[-1]
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
    pt<-ggarrange(gx,g1All,g2All,ncol=1)
    print(pt)
    ggsave(paste0('/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/',pName,"_plot.png"),plot = pt,width = 18, height = 13)
  }

}

#modPlotFuncX(gg,assum="R",fxd=F,iters=100)
plotAll<-function(){
  for (i in c(7,15,19)){
    prmFile<-read.csv("/Users/alm204/Documents/ModelSetups.csv")

    print(i)
    gg<-read.csv(paste0("/Users/alm204/Documents/Cambridge/results/may18/res_",i,".csv"))[,-1]
    gg<-subset(gg,is.na(betaVer)!=T)
    if (prmFile$lFunc[i]=="prFxOsc") gg$zeta_p<-1
    if (prmFile$lFunc[i]=="prFxOsc") gg$pcrProb2<-1
    gg<-subset(gg,ll>-500)
    gg<-gg[245000:nrow(gg),]

    fxdY<-if(prmFile$lFunc[i]=="prFxOsc") T else F
    birthType<-if(prmFile$birthType[i]=="noImmune") 1 else 0

    modPlotFuncX(gg,assum=prmFile$assump[i],fxd=fxdY,iters=100,pName="SILI",birthType = birthType)
  }
}


