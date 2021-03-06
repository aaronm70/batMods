#'initialState
#'@param n number of particles
#'@return initial conditions of each particle

iState<-function(n,prms){
  b <- (prms$omega_m_Val+prms$mj_Val)*((prms$mj_Val+prms$mu_Val)/prms$omega_m_Val)*prms$m_Val/prms$mu_Val
  N<-prms$kappa_Val
  Na=(round(N/(1+prms$m_Val/prms$mu_Val+b/(prms$omega_m_Val+prms$mj_Val))))
  Nj <- round(Na*prms$m_Val/prms$mu_Val)
  Nn <- N - Na - Nj

  conds<-cbind(  Sn_ini = Nn,
                 Sj_ini = Nj,
                 Sf_ini = 0.5*Na,
                 Sm_ini = 0.5*Na,
                 Ma_ini = 0,
                 En_ini = 0,
                 Ej_ini = 0,
                 Ef_ini = 0,
                 Em_ini = 0,
                 In_ini = 0,
                 Ij_ini = 0,
                 If_ini = 0.5*(Na/100*5),
                 Im_ini = 0.5*(Na/100*5),
                 Rn_ini = 0,
                 Rj_ini = 0,
                 Rf_ini = 0,
                 Rm_ini = 0
                 )
  conds <- conds[rep(seq_len(nrow(conds)), each = n),]
  return(conds)
}

scalarFunc<-function(currentParams,m,mu,b,omegam,mj,s,phi){
jspan=15.55/12#juvenile stage span of ~15.55 months
  omegam <- currentParams$omega_m_Val                           #maternal immune waning rate
  mu <- 1/((jspan) - (1/omegam))
  mj_Val<-currentParams$mj_Val
  m_Val<-currentParams$m_Val
  m <- (omegam+mj_Val)*((mj_Val+mu)/omegam)*m_Val/mu

  #demographic parameter values
  #scaling factor (maximum birth rate)
  bessI <- function(z){
    bessInner <- function(x){
      return(exp(z*cos(x)))
    }
    return(integrate(bessInner,lower=0,upper=pi)$value/pi)
  }


  #caclulate initial value of birth rate constant (to optimize)
  c0 <- m/bessI(s/2)*exp(s/2)*sqrt(s/pi)
  #optimize birth rate constant to minimize adult population growth
  grow <- function(k){
    inner <- function(t){
      return(k*exp(-s*(cos(pi*t-phi))^2)*exp(-jspan*mj))
    }
    return((integrate(inner,lower=0,upper=1)$value - m)^2)
  }
  c<-optim(c0,grow,method="Brent",lower=0,upper=10^10)$par[1]
  return(c)
}
