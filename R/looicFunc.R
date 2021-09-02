# looic function, fairly messy and needs tidying slightly
#' @param resultFile location of mcmc results
#' @param prmFile file for model setups
#' @param burn any extra chain burning
#' @param samples if sampling from posterior for pointwise likelihoods
#' @param prmLst more parameter values
#' @return looic and associated data for all models
  looicFunc<-function(resultsFile,prmFile,burn,samples,prmLst,thin=10){
    contBats<-NULL
    r0Val<-NULL
    pval<-NULL
    pvalUp<-NULL
    pvalLw<-NULL
    looicVal<-NULL
    zetaP<-NULL
    zetaU<-NULL
    zetaS<-NULL

for(i in c(2,4,8)){
  print(i)
  gg<-readResFunc(fileLoc=resultsFile,i=i,burn=burn,prmFile=prmFile,thin=10,end=300000)
median(gg$ll)

i<-if(i == 24) 8 else i


  birthType<-ifelse(prmFile$birthType[i]=="noImmune",1,0)

    looTmp<-calc.loo(gg)

  print(looTmp[[1]])
  write.csv(c(looTmp$estimates,looTmp$diagnostics$pareto_k),paste0("C:\\Users\\aaron\\Documents\\batMods\\test\\",i,".csv"))
  nam <- paste("Mod", i, sep = "")
  assign(nam, looTmp)

  looicVal<-as.data.frame(rbind(looicVal, cbind(looTmp[[1]]$estimates[[3]],looTmp[[1]]$estimates[[6]],prmFile$ModelType[i],prmFile$lFunc[i])))
  contBats<-as.data.frame(rbind(contBats,cbind(median(gg$d_val),HPDinterval(as.mcmc(gg$d_val),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaP<-as.data.frame(rbind(zetaP,cbind(median(gg$zeta_p),HPDinterval(as.mcmc(gg$zeta_p),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaS<-as.data.frame(rbind(zetaS,cbind(median(gg$zeta_s),HPDinterval(as.mcmc(gg$zeta_s),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaU<-as.data.frame(rbind(zetaU,cbind(median(gg$pcrProb2),HPDinterval(as.mcmc(gg$pcrProb2),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  r0Val<-as.data.frame(rbind(r0Val,cbind(median(gg$R0_Val),HPDinterval(as.mcmc(gg$R0_Val),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))

}

    names(looicVal)<-c("Looic","se","modType","Structure")
    looicVal
    names(r0Val)<-c("Looic","seLow","seUp","modType","Structure")
    r0Val$mod<-c(1:nrow(r0Val))
    names(contBats)<-c("Looic","seLow","seUp","modType","Structure")
    contBats$mod<-c(1:nrow(contBats))

    names(zetaP)<-c("Looic","seLow","seUp","modType","Structure")
    zetaP$mod<-c(1:nrow(zetaP))

    names(zetaU)<-c("Looic","seLow","seUp","modType","Structure")
    zetaU$mod<-c(1:nrow(zetaU))

    names(zetaS)<-c("Looic","seLow","seUp","modType","Structure")
    zetaS$mod<-c(1:nrow(zetaS))

    looData<-list(looicVal,r0Val,contBats,zetaP,zetaU,zetaS)
    saveRDS(looData, file = "C:\\Users\\aaron\\Documents\\batMods\\looDat.rds")

    return(looData)
  }




  ##calc LOO
  calc.loo <- function(gg) {
    gg1<-gg[gg$chainID==1,]
    gg2<-gg[gg$chainID==2,]

    nlength<-ifelse(nrow(gg1)<nrow(gg2),nrow(gg1),nrow(gg2))

    gg1<-gg1[1:nlength,]
    gg2<-gg2[1:nlength,]


    llMat<-rbind(gg1,gg2)[1:(2*nlength),]
    post<-as.matrix(llMat)
    rel_n_eff <- relative_eff(exp(post[,c(33:48)]),chain_id =llMat$chainID )
    loo_res<-loo(post[,c(33:48)], r_eff = rel_n_eff, cores = 11, moment_match = TRUE)

    return(list(loo_res,llMat,rel_n_eff))
  }


  #trace plots of pMCMC and gelman rubin test
  plotMCMC<-function(gg){
    #convert logged values back to full vals
   # gg$epsilon_Val<-10^gg$epsilon_Val
    #gg$rho_Val<-10^gg$rho_Val
    #gg$gamma_2_Val<-10^gg$gamma_2_Val
    gg1<-subset(gg,chainID==1)
    gg2<-subset(gg,chainID==2)
    gg1<-gg1[1:nrow(gg2),]

    ggF1<-Filter(var, gg1)
    ggF2<-Filter(var, gg2)

    mlst<-mcmc.list(as.mcmc(ggF1[,1:19]),as.mcmc(ggF2[,1:19]))
    gelman.diag(mlst,autoburnin = F,transform=F)


    color_scheme_set("viridisC")
    st=1
    sf=nrow(ggF1)
    mcmc_trace(list(as.matrix(ggF1[st:sf,1:20]),as.matrix(ggF2[st:sf,1:20])))
  }
#6,20,22
