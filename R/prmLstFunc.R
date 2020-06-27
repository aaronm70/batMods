prmLstFunc<-function(fileLoc){
  prmFile<-read.csv(fileLoc)
  prmLst<-NULL
  for(i in c(1:nrow(prmFile))){
    prmTmp<-paramsFunc(
      omega_2_Val = prmFile$Omega[i],
      gamma_1_Val = ifelse(prmFile$GammaVer[i]==1,runif(1,0.1,2),0),
      gamma_2_Val = ifelse(prmFile$GammaVer[i]==2,runif(1,0.1,1),0),
      gammaVer = ifelse(prmFile$GammaVer[i]==2,2,0),
      zeta_s = ifelse(prmFile$SigmaVer[i]==1,runif(1,0.1,1),0),
      sigma_2_Val = ifelse(prmFile$SigmaVer[i]==2,1,0),
      epsilon_Val = ifelse(prmFile$Epsilon[i]>=1,runif(1,0.1,1), 0),
      rho_Val = ifelse(prmFile$Rho[i]>=1,runif(1,0.1,1), 0),
      S2_val=ifelse(prmFile$EnvOsc[i]==T,1,0),
      zeta_p=ifelse(prmFile$zeta[i]==T,runif(1,0.1,1),0),
      pcrProb2=ifelse(prmFile$zeta[i]==T,runif(1,0.1,1),0),
      Phi2_val = ifelse(prmFile$EnvOsc[i]==T,runif(1,0.1,1),0),
      c_val2 = ifelse(prmFile$EnvOsc[i]==T,1,0),
      envOscType=prmFile$envOscType[i],
      betaVer = prmFile$BetaVer[i],
      betaFX=prmFile$betaFX[i]
      )

    prmTmp<-as.list(prmTmp)

    if(prmFile$lFunc[i]=="std") prmTmp$lFunc<-likelihoodFuncBoonah
    if(prmFile$lFunc[i]=="prFx") prmTmp$lFunc<-likelihoodFuncBoonahFixedPr
    if(prmFile$lFunc[i]=="prFxOsc") prmTmp$lFunc<-likelihoodFuncBoonahFixedPr
    if(prmFile$lFunc[i]=="Osc") prmTmp$lFunc<-likelihoodFuncBoonah



    prmTmp$modNum<-prmFile$ModelNum[i]
    prmTmp$assump<-as.character(prmFile$assump[i])
    prmTmp$birthType<-ifelse(as.character(prmFile$birthType[i])=="noImmune",1,0)

    prmLst[[i]] <- prmTmp

  }

return(prmLst)
}
