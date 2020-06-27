#   __ _ _        _ _ _                     _     __                  _   _
#  / /(_) | _____| (_) |__   ___   ___   __| |   / _|_   _ _ __   ___| |_(_) ___  _ __  ___
# / / | | |/ / _ \ | | '_ \ / _ \ / _ \ / _` |  | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#/ /__| |   <  __/ | | | | | (_) | (_) | (_| |  |  _| |_| | | | | (__| |_| | (_) | | | \__ \
#\____/_|_|\_\___|_|_|_| |_|\___/ \___/ \__,_|  |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

#' Beta binomial likelihood function alt formulation
#' @param k Observed data point
#' @param n Simulated data point
#' @param p fraction of population
#' @param w overdisperion parameter
#' @return log likelihood
betaBinom <- function(k, n, p, w) {

  if(is.na(n)) return(-10000) else {

  if (n >= k) {
    if (n <= 0) return (-10000)
    else {
      a <- p * ((1 / w) - 1)
      b <- (1 - p) * ((1 / w) - 1)
      res=lbeta(k + a, n - k + b) - lbeta(a, b) + lchoose(n, k)
      return(res)
    }
  }
  else return (-10000)
  }
}


#' likelihood function
#' @param parms parameter values
#' @param N observed total population
#' @param k observed sero positive
#' @param SN total susceptable
#' @param I total infected
#' @param Ia infected adult males and adult females
#' @param E total exposed/latent
#' @param R total recovered
#' @param PCR observed PCR positive counts - out of N
#' @param time week of the year 1-52 for seasonal information
#' @param urPredDat predicted values for under roost from GAM
#' @return log-likelihood value
#'
#'

likelihoodFuncBoonah<-function(parms,N,k_PpSp,k_PnSp,k_PpSn,k_PnSn,S,I,Ia,E,R,PCR,time,assump="EIR",urPredDat){

  PpSp<-I*parms$zeta_p*parms$zeta_s
  PnSp<-(I+E+R)*parms$zeta_s*(1-parms$zeta_p)
  PpSn<-I*parms$zeta_p*(1-parms$zeta_s)
  PnSn<-(I*(1-parms$zeta_p)*(1-parms$zeta_s)) +S

z<-c(PpSp,PnSp,PpSn,PnSn)
k<-c(k_PpSp,k_PnSp,k_PpSn,k_PnSn)

out<-tryCatch({dmultinom(x=k, prob=z/(S+I+E+R), size=N,log=T)}, error =function(ex){-100000})#


x=1+parms$d_val #number of bats contributing to each pool (d), minimum of 1

p=(I*parms$pcrProb2)/(S+I+E+R)# simulated prevalence?

Pt = 1 -(1 - p)^x #prob of positive in pools or bats?

N_u=46

X_u<-(urPredDat*N_u)

   #betabinomial from rmutil package
  out<-out + tryCatch({ rmutil::dbetabinom(y=round(X_u),m=Pt,size=round(N_u),s=parms$oDist_u,log=T)}, error =function(ex){-100000})

  out<- tryCatch({ if(N>1e5) out-1e10 else out }, error =function(ex){-1e10})

  return(out)
}


#' likelihood function
#' @param parms parameter values
#' @param N observed total population
#' @param k observed sero positive
#' @param SN total susceptable
#' @param I total infected
#' @param Ia infected adult males and adult females
#' @param E total exposed/latent
#' @param R total recovered
#' @param PCR observed PCR positive counts - out of N
#' @param time week of the year 1-52 for seasonal information
#' @param urPredDat predicted values for under roost from GAM
#' @return log-likelihood value
#'
#'

likelihoodFuncBoonahStoch<-function(parms,N,k_PpSp,k_PnSp,k_PpSn,k_PnSn,S,I,Ia,E,R,PCR,time,assump="EIR",urPredDat){

  PpSp<-I*parms$zeta_p*parms$zeta_s
  PnSp<-(I+E+R)*parms$zeta_s*(1-parms$zeta_p)
  PpSn<-I*parms$zeta_p*(1-parms$zeta_s)
  PnSn<-(I*(1-parms$zeta_p)*(1-parms$zeta_s)) +S

  z<-c(PpSp,PnSp,PpSn,PnSn)
  k<-c(k_PpSp,k_PnSp,k_PpSn,k_PnSn)

  out<-tryCatch({dmultinom(x=k, prob=z/(S+I+E+R), size=N,log=T)}, error =function(ex){-100000})#

  x=1+parms$d_val #number of bats contributing to each pool (d), minimum of 1

  p=(I*parms$pcrProb2)/(S+I+E+R)# simulated prevalence?

  Pt = 1 -(1 - p)^x #prob of positive in pools or bats?

  N_u=46

  X_u<-(urPredDat*N_u)

  #betabinomial from rmutil package
  out<-out + tryCatch({ rmutil::dbetabinom(y=round(X_u),m=Pt,size=round(N_u),s=parms$oDist_u,log=T)}, error =function(ex){-100000})

  out<- tryCatch({ if(N>1e5) out-1e10 else out }, error =function(ex){-1e10})

  return(out)
}

