# looic function, fairly messy and needs tidying slightly
#' @param resultFile location of mcmc results
#' @param prmFile file for model setups
#' @param burn any extra chain burning
#' @param samples if sampling from posterior for pointwise likelihoods
#' @param prmLst more parameter values
#' @return looic and associated data for all models
  looicFunc<-function(resultsFile,prmFile,burn,samples,prmLst){
    contBats<-NULL
    r0Val<-NULL
    pval<-NULL
    pvalUp<-NULL
    pvalLw<-NULL
    looicVal<-NULL
    zetaP<-NULL
    zetaU<-NULL
    zetaS<-NULL

for(i in c(1:8)){
  print(i)
  gg<-readResFunc(fileLoc=resultsFile,i=i,burn=burn,prmFile=prmFile,thin=10)
median(gg$ll)
  birthType<-ifelse(prmFile$birthType[i]=="noImmune",1,0)

  #SIRS: Immune waning at the height of rainy season when food is abundant is unlikely
  #also SIRS the later the seasonal force wanes immunity, the faster infection needs to occur (i.e. really high R0), less problem with SILI as straight back to infectious and not susceptible
  #SIR: Increased recovery rate during the winter season to reduce infectious numbers during the summer...but this would unlikely, recovery drops during summer and susceptibles occur as maternal immunity wanes
  #switch from latent to infectious at height of winter season most probable

    looTmp<-calc.loo(data=obsData,gg=gg,assump=prmFile$assump[i],likelihoodFunc = prmLst[[i]]$lFunc,nDraws=samples,birthType = birthType,threshold=NULL)


  print(looTmp[[1]])
  write.csv(c(looTmp$estimates,looTmp$diagnostics$pareto_k),paste0("/Users/alm204/Documents/test",i,".csv"))
  nam <- paste("Mod", i, sep = "")
  assign(nam, looTmp)

  looicVal<-as.data.frame(rbind(looicVal, cbind(looTmp[[1]]$estimates[[3]],looTmp[[1]]$estimates[[6]],prmFile$ModelType[i],prmFile$lFunc[i])))
  contBats<-as.data.frame(rbind(contBats,cbind(median(gg$d_val),HPDinterval(as.mcmc(gg$d_val),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaP<-as.data.frame(rbind(zetaP,cbind(median(gg$zeta_p),HPDinterval(as.mcmc(gg$zeta_p),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaS<-as.data.frame(rbind(zetaS,cbind(median(gg$zeta_s),HPDinterval(as.mcmc(gg$zeta_s),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  zetaU<-as.data.frame(rbind(zetaU,cbind(median(gg$pcrProb2),HPDinterval(as.mcmc(gg$pcrProb2),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))
  r0Val<-as.data.frame(rbind(r0Val,cbind(median(gg$R0_Val),HPDinterval(as.mcmc(gg$R0_Val),prob=0.89),prmFile$ModelType[i],prmFile$lFunc[i])))

  gg$epsilon_Val<-10^gg$epsilon_Val
  gg$rho_Val<-10^gg$rho_Val
  gg$gamma_2_Val<-10^gg$gamma_2_Val
  gg1<-subset(gg,chainID==1)
  gg2<-subset(gg,chainID==2)
  gg1<-gg1[1:nrow(gg2),]

  ggF1<-Filter(var, gg1)
  ggF2<-Filter(var, gg2)
  #ggF3<-Filter(var, gg3)
  mlst<-mcmc.list(as.mcmc(ggF1[,1:19]),as.mcmc(ggF2[,1:19]))
  gelman.diag(mlst,autoburnin = F,transform=F)

 # ggF1<- ggF1[seq(1, nrow(ggF1), by=10), ]
 # ggF2<- ggF2[seq(1, nrow(ggF2), by=10), ]

color_scheme_set("viridisC")
st=1
sf=nrow(ggF1)
mcmc_trace(list(as.matrix(ggF1[st:sf,1:20]),as.matrix(ggF2[st:sf,1:20])))

}

    log_lik_list <- list(Mod4[[2]], Mod6[[2]], Mod8[[2]])
    r_eff_list<-list(Mod4[[3]], Mod6[[3]], Mod8[[3]])
    loo_list <- lapply(1:length(log_lik_list), function(j) {
      loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
    })

   modWeights<- loo_model_weights(loo_list, method = "pseudobma")

   wts2 <- loo_model_weights(
     loo_list,
     method = "stacking",
     optim_control = list(reltol=1e-10)
   )

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
    saveRDS(looData, file = "/Users/alm204/Documents/looDat.rds")

    return(looData)
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

    post<-as.matrix(post)

    rel_n_eff <- relative_eff(exp(post[,c(33:48)]),chain_id =rep(1,nrow(post)) )

    #llMat<-detModFuncLOO_ind(data[2,],post,assump,likelihoodFunc)
    loo_3<-loo(post[,c(33:48)], r_eff = rel_n_eff, cores = 11, moment_match = TRUE)
    #loo_3 <- loo(detModFuncLOO_ind, data = obsData[-1,], draws = post, r_eff = NA,likelihoodFunc=likelihoodFunc,assump=assump,birthType=birthType,cores=10)
    #  return(loo(llMat, r_eff = rel_n_eff, cores = 11))
    return(list(loo_3,llMat,rel_n_eff))
  }


