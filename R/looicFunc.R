
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

  gg$epsilon_Val<-gg$epsilon_Val
  gg$rho_Val<-gg$rho_Val
  gg$gamma_2_Val<-gg$gamma_2_Val
  gg1<-subset(gg,chainID==1)
  gg2<-subset(gg,chainID==2)
  gg2<-gg2[1:nrow(gg1),]

  ggF1<-Filter(var, gg1)
  ggF2<-Filter(var, gg2)
  #ggF3<-Filter(var, gg3)
  mlst<-mcmc.list(as.mcmc(ggF1[,1:19]),as.mcmc(ggF2[,1:19]))
  gelman.diag(mlst,autoburnin = F,transform=F)

 # ggF1<- ggF1[seq(1, nrow(ggF1), by=10), ]
 # ggF2<- ggF2[seq(1, nrow(ggF2), by=10), ]

color_scheme_set("viridisC")
mcmc_trace(list(as.matrix(ggF1[,1:19]),as.matrix(ggF2[,1:19])))

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
  calc.loo <- function(data=obsData,gg=gg,assump=assump,nDraws,likelihoodFunc,birthType,threshold) {

    gg1<-gg[gg$chainID==1,]
    gg2<-gg[gg$chainID==2,]

    if(nDraws>0){
      post<-gg1[sample(nrow(gg1), nDraws), ]
      post2<-gg2[sample(nrow(gg2), nDraws), ]
    }
    else {
      post<-gg1
      post2<-gg2
    }

    post$mu_Val<-as.numeric(as.character(post$mu_Val))
    post2$mu_Val<-as.numeric(as.character(post2$mu_Val))

    llMat<-detModFuncLOO(data,post,assump,likelihoodFunc,birthType)
    llMat2<-detModFuncLOO(data,post2,assump,likelihoodFunc,birthType)

    nlength<-ifelse(nrow(llMat)<nrow(llMat2),nrow(llMat),nrow(llMat2))


    llMat10<-rbind(llMat,llMat2)[1:(2*nlength),]

    rel_n_eff <- relative_eff(exp(llMat10),chain_id = c(rep(1,nlength),rep(2,nlength)) )

    #llMat<-detModFuncLOO_ind(data[2,],post,assump,likelihoodFunc)
    loo_3<-loo(llMat10, r_eff = rel_n_eff, cores = 11, k_threshold = 0.7,save_psis = T)
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
