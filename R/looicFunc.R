
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


    looTmp<-calc.loo(data=obsData,gg=gg,assump=prmFile$assump[i],likelihoodFunc = prmLst[[i]]$lFunc,nDraws=samples,birthType = birthType)


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
  #  ggFm<-list(as.mcmc(ggF[250000:nrow(ggF2),]),as.mcmc(ggF2[250000:nrow(ggF2),]))
  #gelman.diag(ggFm)

 # ggF1<- ggF1[seq(1, nrow(ggF1), by=10), ]
 # ggF2<- ggF2[seq(1, nrow(ggF2), by=10), ]

 # color_scheme_set("viridisC")
  #mcmc_trace(list(as.matrix(ggF1[,1:19]),as.matrix(ggF2[,1:19])))

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


    return(list(looicVal,r0Val,contBats,zetaP,zetaU,zetaS))
  }



 # gg$epsilon_Val<-10^gg$epsilon_Val
 # gg$rho_Val<-10^gg$rho_Val
 # gg$gamma_2_Val<-10^gg$gamma_2_Val
#
 # gg2$epsilon_Val<-10^gg2$epsilon_Val
 # gg2$rho_Val<-10^gg2$rho_Val
 # gg2$gamma_2_Val<-10^gg2$gamma_2_Val

 # ind <- seq(1, nrow(gg), by=5)
 # gg<- gg[ind, ]
 # gg2<- gg2[ind, ]

 # gg3<-subset(gg3,is.na(betaVer)!=T)
 # if (prmFile$lFunc[i]=="prFxOsc") gg3$zeta_p<-1
 # if (prmFile$lFunc[i]=="prFxOsc") gg3$pcrProb2<-1
 # gg3<-subset(gg3,ll>-1000)
 # gg3<-gg3[150000:200000,]
 # gg4<-subset(gg4,is.na(betaVer)!=T)
 # if (prmFile$lFunc[i]=="prFxOsc") gg4$zeta_p<-1
 # if (prmFile$lFunc[i]=="prFxOsc") gg4$pcrProb2<-1
 # gg4<-subset(gg4,ll>-1000)
 # gg4<-gg4[150000:200000,]
#
 # gg3$epsilon_Val<-10^gg3$epsilon_Val
 # gg3$rho_Val<-10^gg3$rho_Val
 # gg3$gamma_2_Val<-10^gg3$gamma_2_Val
#
 # gg4$epsilon_Val<-10^gg4$epsilon_Val
 # gg4$rho_Val<-10^gg4$rho_Val
 # gg4$gamma_2_Val<-10^gg4$gamma_2_Val
#
 # ind <- seq(1, nrow(gg3), by=5)
 # gg3<- gg3[ind, ]
 # gg4<- gg4[ind, ]
#

#
# ggsave(filename=paste0("/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/mcmc_",i,".png"),plot=f,width = 18, height = 13)

  #ll<-gg$ll
  #gg<-gg[,-30]
  #gg$betaFX<-0
  #gg$betaFXVal<-0
  #gg$ll<-ll


  #gg<-gg[,-(which(colSums(gg)==0))]
 #gg<-if (prmFile$lFunc[i]=="prFxOsc"&&prmFile$lFunc[i]=="prFx") gg[ , -which(names(gg) %in% c("zeta_p","pcrProb2","sigmaVer", "gammaVer", "betaVer"))] else  gg[ , -which(names(gg) %in% c("sigmaVer", "gammaVer", "betaVer"))]
 #gg<- if(prmFile$lFunc[i]!="Osc"&&prmFile$lFunc[i]!="prFxOsc")  gg[ , -which(names(gg) %in% c("c_val2","d_val","Phi2_val"))] else  gg[ , -which(names(gg) %in% c("c_val2"))]
  #modPlotFunc(gg)
#print(paste(round(median(1e-20+gg$R0_Val),3),",",round(median(1e-20+gg$gamma_1_Val),3),",",round(median(1e-20+gg$gamma_2_Val),3),",",
#            round(median(1e-20+gg$zeta_s),3),",",round(median(1e-20+gg$sigma_2_Val),3),",",round(median(1e-20+gg$epsilon_Val),3),",",
#            round(median(1e-20+gg$omega_2_Val),3),",", round(median(1e-20+gg$rho_Val),3),",", round(median(1e-20+10^gg$kappa_Val),3)))
 #gg<- gg[, colSums(gg != 0) > 0]

 #
#ciV<- HPDinterval(as.mcmc(gg$pcrProb2))
# pval<-rbind(pval,median(gg$pcrProb2))
# pvalUp<-rbind(pvalUp,ciV[2])
# pvalLw<-rbind(pvalLw,ciV[1])
#
# gg2<-colMedians(gg)
# gg2<-as.list(gg2)
#
# gg2$gamma_1_Val<-ifelse(is.null(gg2$gamma_1_Val)==T,0,gg2$gamma_1_Val)
# gg2$gamma_2_Val<-ifelse(is.null(gg2$gamma_2_Val)==T,0,gg2$gamma_2_Val)
# gg2$zeta_s<-ifelse(is.null(gg2$zeta_s)==T,0,gg2$zeta_s)
# gg2$sigma_2_Val<-ifelse(is.null(gg2$sigma_2_Val)==T,0,gg2$sigma_2_Val)
# gg2$omega_2_Val<-ifelse(is.null(gg2$omega_2_Val)==T,0,gg2$omega_2_Val)
# gg2$epsilon_Val<-ifelse(is.null(gg2$epsilon_Val)==T,0,gg2$epsilon_Val)
#
#beta_Val=((gg2$R0_Val*((gg2$m_Val^2)+(gg2$rho_Val*gg2$m_Val)+(gg2$sigma_2_Val*gg2$m_Val)+(gg2$m_Val*gg2$gamma_1_Val)+
#                     (gg2$m_Val*gg2$gamma_2_Val)+(gg2$m_Val*gg2$sigma_2_Val)+(gg2$m_Val*gg2$epsilon_Val)+
#                     (gg2$rho_Val*gg2$sigma_2_Val)+(gg2$sigma_2_Val*gg2$gamma_1_Val)+(gg2$sigma_2_Val*gg2$gamma_2_Val)+
#                     (gg2$rho_Val*gg2$zeta_s)+(gg2$gamma_1_Val*gg2$zeta_s)+(gg2$gamma_2_Val*gg2$zeta_s)+
##                     (gg2$gamma_1_Val*gg2$epsilon_Val)+(gg2$gamma_2_Val*gg2$epsilon_Val))/(10^gg2$kappa_Val*gg2$epsilon_Val)))
##print(beta_Val)
##
##
##beta_Val2<-gg2$R0_Val*((gg2$epsilon_Val+(gg2$zeta_s+gg2$sigma_2_Val)+gg2$m_Val)*((gg2$gamma_1_Val+gg2$gamma_2_Val)+gg2$m_Val+gg2$rho_Val)
##                   -gg2$epsilon_Val*gg2$rho_Val)/(10^gg2$kappa_Val*gg2$epsilon_Val+(gg2$zeta_s+gg2$sigma_2_Val)+gg2$m_Val)
##print(beta_Val2)
#


#
#aicval<-as.data.frame(aicval)
#aicval$modNum<-c(1:22)
#aicval[order(aicval$V1),]
#
#llVal<-as.data.frame(llVal)
#llVal$modNum<-c(1:22)
#llVal[order(llVal$V1),]
##aicval$mod<-c(1:nrow(aicval))
#
##pval<-as.data.frame(pval)
##pval$aic<-aicval
##pval$Up<-pvalUp
##pval$Lw<-pvalLw
##pval<-pval[order(pval$aic),]
##pval$rank<-c(1:nrow(pval))
##
##ggplot(pval,aes(x=rank,y=V1))+
##  geom_errorbar(aes(ymin=Lw, ymax=Up,col="blue"), width=.2,
##                position=position_dodge(0.05),alpha=0.5)+
##  geom_point(aes(col="black"),alpha=0.7)+
##  ylab("Probability of shedding UR")+
##  scale_color_manual(labels = c("R0", "95% CI"), values = c("black", "blue"))+
##  theme_bw(base_size = 20 )+
##  theme(legend.title = element_blank())
##
##
##
##
##
##1
##7
##11
##10
##9
##16
##5
##14
##8
##15
##x <- data.frame(v1=rnorm(100),v2=rnorm(100,1,1),v3=rnorm(100,0,2))
##library(ggplot2);library(reshape2)
##data<- melt(x)
##ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
#plot(density(gg$ll,adjust=2))
#for(i in c(7,11,10,9,16,5,14,8,15)){
#
#  print(i)
#  gg<-read.csv(paste0("/Users/alm204/Documents/Cambridge/res",i,".csv"))[,-1]
#  gg<- gg[, colSums(gg != 0) > 0]
#
#  aicval<--2*(gg$ll) + 2*ncol(gg)
#  gf<-cbind(gf,aicval)
#}
#
##names(gf)<-c("1","2","3","4","5","6","7","8","9","10")
##dataG<- melt(gf)
##ggplot(dataG,aes(x=value, fill=variable)) + geom_density(alpha=0.3,adjust=2)+
##  guides(fill=guide_legend(title="Model Ranking"))+
##  xlab("AIC")+
##  scale_fill_viridis_d(option="C",direction=-1)+
##  theme_bw(base_size = 20)
##
#
##ggplot(gg, aes(x=zeta_p, y=rho_Val) ) +
##  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
##  theme_bw()+
##  xlab("zeta_i")+
##  ylab("rho")
##
