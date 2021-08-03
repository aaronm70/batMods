


#' model step function - runs model in steps using Odin, returning the model state at a subsequent time point
#'@param particles number of particles in run
#'@param mod  redundant - needs removing
#'@param startTime current time in model
#'@param endTime time of next datapoint
#'@param birthtype maternal immunity on off
#'@return simulated data for time step
modStepJI <- function(particles,mod,birthType,startTime,endTime) {
  #parse input data

  initialState <-particles#state of bats
  currentTime <-startTime #time in model to run from
  nextTime <- endTime#time in model to run to
  birthType<-birthType

  #initialise model
  particles<-as.numeric(particles)
  N = sum(particles)

  #assign("params", params, env = .GlobalEnv)
  mod$set_user(
    Sn_ini = particles[1],
    Sj_ini = particles[2],
    Sf_ini = particles[3],
    Sm_ini = particles[4],
    Ma_ini = particles[5],
    En_ini = particles[6],
    Ej_ini = particles[7],
    Ef_ini = particles[8],
    Em_ini = particles[9],
    In_ini = particles[10],
    Ij_ini = particles[11],
    If_ini = particles[12],
    Im_ini = particles[13],
    Rn_ini = particles[14],
    Rj_ini = particles[15],
    Rf_ini = particles[16],
    Rm_ini = particles[17],
    timeOscIni = currentTime/(365*4)
  )


  tx = seq(currentTime, nextTime)
  simDat <- as.data.frame(mod$run(tx))
  return(simDat)
}



##function for running deterministic model
#'@param params model parameters
#'@param obsData  observed data
#'@param assump redundant - needs removing
#'@param likelihoodFunc likelihood function to use
#'@param birthtype maternal immunity on off
detModFunc<-function(params,obsData,assump,likelihoodFunc,birthType){
  currentTime=obsData$NumDays[1]
  nextTime=obsData$NumDays[nrow(obsData)]


#Some parameter values log10 is used to keep them on a similar scale to other parameters, makes fitting easier
  params$kappa_Val<-10^params$kappa_Val

  initialState<- as.data.frame(iState(prms=params,n=2))[1,]
  initialState<-as.numeric(initialState)

  if(params$gammaVer==0){
    params$rho_Val<-ifelse(params$rho_Val>-1,10^params$rho_Val,-1)
    params$epsilon_Val<-ifelse(params$epsilon_Val>-1,10^params$epsilon_Val,-1)
  }
  if(params$gammaVer>0){
    params$gamma_2_Val<-ifelse(params$gamma_2_Val>-1,10^params$gamma_2_Val,-1)
  }


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
      birthType=birthType,
      betaFX=params$betaFX,
      betaFXVal=params$betaFXVal,
      mu_Val=params$mu_Val

    )
  tx = seq(currentTime, nextTime)
  simDat <- as.data.frame(mod$run(tx))

  ll=NULL
  for (i in c(1:(nrow(obsData)-1))){
    detMod2<-simDat[obsData$NumDays[i+1],-1]
    ll=rbind(ll,mapply(likelihoodFunc,
                 S=sum(round(as.vector(detMod2[c(1:4)]))),
                 E=sum(round(as.vector(detMod2[c(6:9)]))),
                 I=sum(round(as.vector(detMod2[c(10:13)]))),
                 R=sum(round(as.vector(detMod2[c(14:17,5)]))),
                 MoreArgs = list(
                   parms = params,
                   N = obsData$negatives[i + 1] + obsData$positives[i + 1],
                   k_PpSp = obsData$PpSp[i + 1],
                   k_PpSn = obsData$PpSn[i + 1],
                   k_PnSp = obsData$PnSp[i + 1],
                   k_PnSn = obsData$PnSn[i + 1],
                   PCR = obsData$pcrPos[i + 1],
                   urPredDat = tryCatch({obsData$fit[i + 1]}, error =function(ex){NULL}),
                   time = obsData$Date[i + 1],
                   assump = assump
                 )))
  }
return(ll)

}


####function for running deterministic model and getting output for LOO and WAIC
#'@param data_i observed data
#'@param draws  posterior draws
#'@param assump redundant - needs removing
#'@param likelihoodFunc likelihood function to use
#'@param birthtype maternal immunity on off
detModFuncLOO<-function(data_i,draws,assump="I",likelihoodFunc,birthType){

  llMat<-NULL
  for(j in c(1:nrow(draws))) {

  params<-draws[j,]

  currentTime=data_i$NumDays[1]
  nextTime=data_i$NumDays[nrow(data_i)]
  initialState<- as.data.frame(iState(prms=params,n=2))[1,]
  initialState<-as.numeric(initialState)
  assign("params", params, env = .GlobalEnv)

  params$kappa_Val<-params$kappa_Val
  if(params$gammaVer==0){
    params$rho_Val<-ifelse(params$rho_Val>-1,params$rho_Val,-1)
    params$epsilon_Val<-ifelse(params$epsilon_Val>-1,params$epsilon_Val,-1)
  }
  if(params$gammaVer>0){
    params$gamma_2_Val<-ifelse(params$gamma_2_Val>-1,params$gamma_2_Val,-1)
  }
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
      birthType=birthType,
      betaFX=params$betaFX,
      betaFXVal=params$betaFXVal
    )
  tx = seq(currentTime, nextTime)
  simDat <- as.data.frame(mod$run(tx))


    ll=NULL
    for (i in c(1:(nrow(data_i)-1))){
      detMod2<-simDat[data_i$NumDays[i+1],-1]
      ll=rbind(ll,mapply(likelihoodFunc,
                         S=sum((as.vector(detMod2[c(1:4)]))),
                         E=sum((as.vector(detMod2[c(6:9)]))),
                         I=sum((as.vector(detMod2[c(10:13)]))),
                         R=sum((as.vector(detMod2[c(14:17,5)]))),
                         MoreArgs = list(
                           parms = params,
                           N = obsData$negatives[i + 1] + obsData$positives[i + 1],
                           k_PpSp = obsData$PpSp[i + 1],
                           k_PpSn = obsData$PpSn[i + 1],
                           k_PnSp = obsData$PnSp[i + 1],
                           k_PnSn = obsData$PnSn[i + 1],
                           PCR = obsData$pcrPos[i + 1],
                           urPredDat = tryCatch({obsData$fit[i + 1]}, error =function(ex){NULL}),
                           time = obsData$Date[i + 1],
                           assump = assump
                         )))
    }

    llMat<-rbind(llMat,as.vector(ll))

  }

    return(as.matrix(llMat))

}




##function for running deterministic model in shiny - old and needs updating
detModFuncShiny<-function(params,times,assump=NULL,likelihoodFunc=NULL){
  currentTime=times[1]
  nextTime=times[2]
  params$kappa_Val<-log10(params$kappa_Val)
  initialState<- as.data.frame(iState(prms=params,n=2))[1,]
  params$kappa_Val<-params$kappa_Val

  initialState<-as.numeric(initialState)
  assign("params", params, env = .GlobalEnv)

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
      #   In_ini = initialState[10],
      #  Ij_ini = initialState[11],
      If_ini = initialState[10],
      Im_ini = initialState[11],
      Rn_ini = initialState[12],
      Rj_ini = initialState[13],
      Rf_ini = initialState[14],
      Rm_ini = initialState[15],
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
      betaFX=params$betaFX,
      betaFXVal=params$betaFXVal
    )
  tx = seq(currentTime, nextTime)
  simDat <- as.data.frame(mod$run(tx))

 # ll=0
 # for (i in c(1:(nrow(obsData)-1))){
 #   detMod2<-simDat[obsData$NumDays[i+1],-1]
 #   ll=ll+mapply(likelihoodFunc,
 #                S=sum(round(as.vector(detMod2[c(1:4)]))),
 #                E=sum(round(as.vector(detMod2[c(6:9)]))),
 #                I=sum(round(as.vector(detMod2[c(10:13)]))),
 #                R=sum(round(as.vector(detMod2[c(14:17,5)]))),
 #                MoreArgs = list(
 #                  parms = params,
 #                  N = obsData[i + 1, 3] + obsData[i + 1, 2],
 #                  k = obsData[i + 1, 2],
 #                  PCR = obsData[i + 1, 4],
 #                  urPredDat = tryCatch({obsData[i + 1, 8]}, error =function(ex){NULL}),
 #                  time = obsData$Date[i + 1],
 #                  assump = assump
 #                ))
 # }
  simDat<-simDat[c((nrow(simDat)-365*4):nrow(simDat)),]
  simDat2<-as.data.frame(rowSums(simDat[,2:5]))
  names(simDat2)<-"S_f"
  simDat2$L_f<-rowSums(simDat[,7:10])
  simDat2$I_f<-rowSums(simDat[,11:14])
  simDat2$R_f<-rowSums(simDat[,15:18,6])
  simDat2$S<-rowSums(simDat[,2:5])/rowSums(simDat[,2:18])
  simDat2$L<-rowSums(simDat[,7:10])/rowSums(simDat[,2:18])
  simDat2$I<-rowSums(simDat[,11:14])/rowSums(simDat[,2:18])
  simDat2$R<-rowSums(simDat[,15:18,6])/rowSums(simDat[,2:18])
  simDat2$time<-c(1:nrow(simDat2))
  simDat2$step<-simDat$step

  return(simDat2)

}




####function for running pFilt and obtaining pointwise log likelihood values for old runs where this was not explicitly implemented
#'@param data observed data
#'@param gg  mcmc results
#'@param assump redundant - needs removing
#'@param likelihoodFunc likelihood function to use
#'@param birthtype maternal immunity on off
#'@param nDraws if drawing from posterior
pFiltMat<-function(data,gg,assump="I",likelihoodFunc=NULL,birthType=0,nDraws){

  draws<-gg[sample(nrow(gg), nDraws), ]
  draws$mu_Val<-as.numeric(as.character(draws$mu_Val))




  pDraws<-NULL
  cat("Percent complete = ")
  for(i in c(1:nrow(draws))){
    percentage<-(100/nDraws)*i
if(percentage %% 5 == 0)    {cat(paste(percentage,"% "))}

    draws_i=draws[i,]

    pDrawsTmp <- -1000
    attempt <- 1
    while( mean(pDrawsTmp)<=-250 && attempt <= 1000 ) {
      attempt <- attempt + 1
      if(attempt>=10) draws_i = gg[sample(nrow(gg), 1), ] #if pFilt fails too many times, try a new sample
      try(
        pDrawsTmp<-pFilt(likelihoodFunc=likelihoodFuncBoonahStoch,n=25,iState=iState,stepFun=modStepJI,obsData=data,
                                             prms=draws_i,assump=assump,mat=T,birthType=birthType)#+lpriorBoonah(draws[i,])
        )
    }

    pDraws<- rbind(pDraws,as.vector(pDrawsTmp))

  }

  return(pDraws)

}
