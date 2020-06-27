lpriorBoonah <- function(parms) with(parms, {
  gamma2_prior<-dunif(parms$gamma_2_Val,min=-1 ,max=2.5,log=T)
  zeta_s_prior<-dunif(parms$zeta_s,min=0 ,max=1,log=T)
  omega2_prior<-dunif(parms$omega_2,min=0,max=100,log=T)
  epsilon_prior<-dunif(parms$epsilon,min=-1,max=2.5,log=T)
  kappa_prior<-dnorm(parms$kappa,mean=3.6,sd=0.1,log=T)
  R0_prior<-dunif(parms$R0_Val,min=1,max=50,log=T)
  rho_prior<-dunif(parms$rho_Val,min=-1,max=2.5,log=T)
  s_prior<-dnorm(parms$s_Val,mean=130,sd=10,log=T)
  omega_m_prior<-dnorm(parms$omega_m_Val, mean = 0.8, sd = 0.03,log=T)
  mj_Val_prior<-dnorm(parms$mj_Val,mean=0.5,sd=0.01,log=T)
  m_Val_prior<-dnorm(parms$m_Val,mean=0.186,sd=0.02,log=T)
  d_val_prior<-dunif(parms$d_val,min=0,max=20,log=T)
  S2_val_prior<-dunif(parms$S2_val,min=0,max=200,log=T)
  zeta_p_prior<-dunif(parms$zeta_p,min=0,max=1,log=T)
  pcrProb_prior2<-dunif(parms$pcrProb2,min=0,max=1,log=T)
  phi_prior<-dnorm(parms$phi_Val,mean=7.18,sd=0.01,log=T)
  c_val2Prior<-dunif(parms$c_val2,min=0,max=1,log=T)
  # sigmaProbPrior<-dunif(parms$sigmaVerProbVal,min=0.01,max=0.99,log=T)
  oDist_uPrior<-dunif(parms$oDist_u,min=0,max=8,log=T)

  Phi2_valPrior<-dunif(parms$Phi2_val,min=0 ,max=100,log=T)

  priorSum<-(S2_val_prior+Phi2_valPrior+omega_m_prior+c_val2Prior+s_prior+phi_prior+m_Val_prior+mj_Val_prior+rho_prior+R0_prior+kappa_prior+
               omega2_prior+gamma2_prior+zeta_s_prior+pcrProb_prior2+
               epsilon_prior+oDist_uPrior+d_val_prior+zeta_p_prior)
  return(priorSum)
})
