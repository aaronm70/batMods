lpriorBoonah <- function(parms) with(parms, {
  gamma2_prior<-ifelse(parms$gamma_2_Val!=0,dunif(parms$gamma_2_Val,min=-1 ,max=2.56,log=T),dunif(0.5,min=-1,max=2.56,log=T))
  zeta_s_prior<-dunif(parms$zeta_s,min=0 ,max=1,log=T)
  omega2_prior<-ifelse(parms$omega_2_Val!=0,dunif(parms$omega_2,min=0,max=365,log=T),dunif(1,min=0,max=365,log=T))
  epsilon_prior<-ifelse(parms$epsilon_Val!=0, dunif(parms$epsilon,min=-1,max=2.56,log=T),dunif(0.5,min=-1,max=2.56,log=T))
  kappa_prior<-dnorm(parms$kappa,mean=3.60206,sd=0.15,log=T)
  R0_prior<-dunif(parms$R0_Val,min=1,max=150,log=T)
  rho_prior<-ifelse(parms$rho_Val!=0,dunif(parms$rho_Val,min=-1,max=2.56,log=T),dunif(0.5,min=-1,max=2.56,log=T))
  s_prior<-dnorm(parms$s_Val,mean=130,sd=15,log=T)
  omega_m_prior<-dnorm(parms$omega_m_Val, mean = 1.431, sd = 0.087,log=T)
  mj_Val_prior<-dnorm(parms$mj_Val,mean=0.5,sd=0.07,log=T)
  m_Val_prior<-dnorm(parms$m_Val,mean=0.135,sd=0.03,log=T)
  d_val_prior<-dunif(parms$d_val,min=0,max=30,log=T)
  S2_val_prior<-ifelse(parms$c_val2==1, dunif(parms$S2_val,min=1,max=50,log=T),dunif(1,min=1,max=50,log=T))
  zeta_p_prior<-dunif(parms$zeta_p,min=0,max=1,log=T)
  pcrProb_prior2<-dunif(parms$pcrProb2,min=0,max=1,log=T)
  phi_prior<-dnorm(parms$phi_Val,mean=7.18,sd=0.1,log=T)
  #mu_Val_prior<-dnorm(parms$mu_Val ,mean=1.37,sd=0.02,log=T)


  # sigmaProbPrior<-dunif(parms$sigmaVerProbVal,min=0.01,max=0.99,log=T)
  oDist_uPrior<-dunif(parms$oDist_u,min=0,max=10,log=T)

  Phi2_valPrior<-dunif(parms$Phi2_val,min=0 ,max=3.5,log=T)

  priorSum<-(S2_val_prior+Phi2_valPrior+omega_m_prior+s_prior+phi_prior+m_Val_prior+mj_Val_prior+rho_prior+R0_prior+kappa_prior+
               omega2_prior+gamma2_prior+zeta_s_prior+pcrProb_prior2+
               epsilon_prior+oDist_uPrior+d_val_prior+zeta_p_prior)
  return(priorSum)
})






createPriors_sitka<-function(parms){

  nm<-c("gamma_2_Val","zeta_s","omega_2_Val","epsilon_Val","kappa","R0_Val","rho_Val","s_Val","omega_m_Val","mj_Val","m_Val",
        "d_val","c_val2","zeta_p","pcrProb2","phi_Val","oDist_u","Phi2_val")


  gamma2_prior<-ifelse(parms$gamma_2_Val!=0,c(-1 ,2.56),c(0.5,-1,max=2.56,log=T))
  zeta_s_prior<-dunif(parms$zeta_s,min=0 ,max=1,log=T)
  omega2_prior<-ifelse(parms$omega_2_Val!=0,dunif(parms$omega_2,min=0,max=365,log=T),dunif(1,min=0,max=365,log=T))
  epsilon_prior<-ifelse(parms$epsilon_Val!=0, dunif(parms$epsilon,min=-1,max=2.56,log=T),dunif(0.5,min=-1,max=2.56,log=T))
  kappa_prior<-dnorm(parms$kappa,mean=3.60206,sd=0.15,log=T)
  R0_prior<-dunif(parms$R0_Val,min=1,max=150,log=T)
  rho_prior<-ifelse(parms$rho_Val!=0,dunif(parms$rho_Val,min=-1,max=2.56,log=T),dunif(0.5,min=-1,max=2.56,log=T))
  s_prior<-dnorm(parms$s_Val,mean=130,sd=15,log=T)
  omega_m_prior<-dnorm(parms$omega_m_Val, mean = 1.431, sd = 0.087,log=T)
  mj_Val_prior<-dnorm(parms$mj_Val,mean=0.5,sd=0.07,log=T)
  m_Val_prior<-dnorm(parms$m_Val,mean=0.135,sd=0.03,log=T)
  d_val_prior<-dunif(parms$d_val,min=0,max=30,log=T)
  S2_val_prior<-ifelse(parms$c_val2==1, dunif(parms$S2_val,min=1,max=50,log=T),dunif(1,min=1,max=50,log=T))
  zeta_p_prior<-dunif(parms$zeta_p,min=0,max=1,log=T)
  pcrProb_prior2<-dunif(parms$pcrProb2,min=0,max=1,log=T)
  phi_prior<-dnorm(parms$phi_Val,mean=7.18,sd=0.1,log=T)
  #mu_Val_prior<-dnorm(parms$mu_Val ,mean=1.37,sd=0.02,log=T)

  # sigmaProbPrior<-dunif(parms$sigmaVerProbVal,min=0.01,max=0.99,log=T)
  oDist_uPrior<-dunif(parms$oDist_u,min=0,max=10,log=T)

  Phi2_valPrior<-dunif(parms$Phi2_val,min=0 ,max=3.5,log=T)



  ##Need to check what priors we are using!
  pMaxima <- as.vector(unlist(parms[nm])*(1+(f.increase)))
  pMinima <- as.vector(unlist(parms[nm])*(1-(f.decrease)))
  pValues <- as.vector(unlist(parms[nm]))

  pMaxima[1:11] <- f.increase[1:11]
  pMinima[1:11] <- f.decrease[1:11]

  pMaxima[[30]]<-0.5
  pMinima[[31]]<-0.01

  pMaxima[[48]]<-0.65
  pMinima[[48]]<-0.4

  pMaxima[[49]]<-7
  pMinima[[49]]<-3

  pMaxima[[34]]<-0.3
  pMinima[[34]]<-0.01

  sdVals<-(pMaxima-pMinima)*0.8
  #ES1 and 2 on very different scales so set manually
  sdVals[7]<-15
  sdVals[8]<-15


  sdVals[[48]]<-0.1
  sdVals[[49]]<-1.3

  sdVals[[50]]<-20
  sdVals[[51]]<-0.05
  pMaxima[[50]]<--5
  pMinima[[50]]<--95

  pMaxima[[51]]<-0.9
  pMinima[[51]]<-0.7


  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)

  return(priorVals)

}


