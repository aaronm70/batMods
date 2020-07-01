
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    Metropolis-Hastings MCMC Sampler                                                            #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################


#add pcr data
#check serology cutoff
#stratify by species
#add seasonality to carrying capacity
#spatial temporal aspects of hendra virus hume field et al. - shows above roost and underroost
#population meta data in the paper
#population = 25,000?
#check birth pulse timing
#little reds different birth pulse timings
#blacks and greys similar but different shedding
#THINK CAREFULLY ABOUT LIKELIHOOD EG dbinom(50, prob=point.prev, size=100,log=T) POINT.PREV=SIMULATED PREVELANCE




#' Sequential proposal function
#' @param current current parameters
#' @param sdTune s.d. tuning from tuner function
#' @param prmNum which parameter to update
#' @param sdProps
#' @param tune
#' @return proposal
sequential.proposer <- function(current, prmNum, sdProps) {
  proposal <- current
  propVal <-
    proposal[prmNum] + (rnorm(1, mean = 0, sd = sdProps[prmNum]))
  propVal<-if(propVal<=-0.69897 && prmNum==10) -0.69897 else propVal
  propVal<-if(propVal<=-0.69897 && prmNum==15) -0.69897 else propVal
  propVal<-if(propVal<=-0.69897 && prmNum==2) -0.69897 else propVal

  proposal[prmNum] <- propVal
  proposal[proposal<=-0.69897]<--0.69897
  return(proposal)
}



#' function for adaptive blocked proposals based on var-covar matrix
#' @param current current parameters
#' @param sdTune tuned s.d.
#' @return proposal parameters
multiv.proposer <-  function(current, sdProps) {
    proposal <-
      current[c(which(sdProps !=0))] + (rmnorm(1, mean = 0, varcov = covar))*sdProps[c(which(sdProps !=0))]

           current[c(which(sdProps !=0))]<-proposal
proposal<-current
  proposal
}


#' proposal sd tuning function
#' @param current s.d.
#' @param target acceptance ratio
#' @param current acceptance ratio
#' @param maximum proposal s.d.
#' @return tuned s.d.
tuner <- function(curSd, acptR, curAcptR, maxSddProps) {
    curAcptR[curAcptR==1] <- 0.99
    curAcptR[curAcptR==0] <- 0.01
  curSd = (curSd * qnorm(acptR / 2)) / qnorm(curAcptR / 2)
  curSd[c(which(curSd > maxSddProps))] <-maxSddProps[c(which(curSd > maxSddProps))]
  curSd[c(which(maxSddProps ==0))] <-0

  return(curSd)
}


#'sequentialTuner
#' @param current s.d.
#' @param target acceptance ratio
#' @param current acceptance ratio
#' @param maximum proposal s.d.
#' @param i which parameter to work on
#' @return tuned s.d.
tunerSeq <- function(curSd, acptR, curAcptR, maxSddProps, i) {
 if (curAcptR[i] == 1)
    curAcptR[i] <- 0.99
  if (curAcptR[i] == 0)
    curAcptR[i] <- 0.01
  curSd[i] = (curSd[i] * qnorm(acptR[i] / 2)) / qnorm(curAcptR[i] / 2)
  curSd[i][curSd[i] > maxSddProps[i]] <- maxSddProps[i]

  curSd[c(which(curSd > maxSddProps))] <-maxSddProps[c(which(curSd > maxSddProps))]
  curSd[c(which(maxSddProps ==0))] <-0

  return(curSd[i])
}

###hyperfunc
#hypFunc<-function(currentParams,nRws,maxSddProps){
#
#  X <- randomLHS(nRws, length(which(maxSddProps>0)))
#  for(i in c(1:ncol(X))){
#    st<-maxSddProps[which(maxSddProps>0)][i]
#    X[,i]<-if(st>1) X[,i]*st else X[,i]
#  }
#  LS<-matrix(0, ncol=length(currentParams),nrow=nrow(X))
#  LS[,which(maxSddProps>0)]<-X
#  LS[,which(maxSddProps==0)]<-matrix(unlist(data.frame(currentParams)[,which(maxSddProps==0)]),ncol=length(which(maxSddProps==0)),nrow=nRws,byrow=T)
#  return(data.frame(LS))
#}

##calculates the mean beta from oscilating variables
betaMean<-function(params){
  s=params$S2_val
  phi=params$Phi2_val
  c=params$c_val2
  x<-seq(as.Date("2013-07-22"), as.Date("2014-07-21"), "days")
  x<-yday(x)/365
  sDrive <-c* exp(-s*(cos(pi*x- phi))^2)

  epsilon_Val<-params$epsilon_Val
  rho_Val<-params$rho_Val

  R0_Val<-params$R0_Val
  gamma_2_Val<-params$gamma_2_Val

  if(params$gammaVer==0){
    epsilon_Val<-params$epsilon_Val*sDrive
    rho_Val<-params$rho_Val
  }
  if(params$gammaVer>0){
    gamma_2_Val<-params$gamma_2_Val*sDrive
  }

 #
 # r0<-((beta*params$kappa_Val)*(epsilon_Val+params$m_Val))/
 #   ((epsilon_Val+params$m_Val)*(gamma_2_Val+params$m_Val+rho_Val)-epsilon_Val*rho_Val)
 #
  return(mean(params$R0_Val*((epsilon_Val+params$sigma_2_Val+params$m_Val)*
                                       (params$gamma_2_Val+params$m_Val+rho_Val)
                                     -epsilon_Val*rho_Val)/(params$kappa_Val*(epsilon_Val+params$sigma_2_Val+params$m_Val))))
}


#' pMCMCsampler
#' @param initial parameter guess
#' @param randInit T then randomly sample initial parameters instead of initParams value
#' @param proposer proposal function, multivariate block (adaptiveMCMC must = T) or sequential can have adaptive or not adaptive tuning
#' @param sdProps standard deviation for proposal distributions - in adaptive this is the starting sd
#' @param maxSddProps maximum values for sd proposals if using adaptive MCMC
#' @param niter number of iterations to run the MCMC for
#' @param particles number of particles for particle filter
#' @param nburn number of mcmc iterations to burn
#' @param monitoring 0 = no monitoring, > 0 prints more progress information
#' @param adaptiveMCMC T/F if true uses tuning for s.d. of parmater proposal distributions based on acceptance ratios
#' @param proposerType "seq" or "block", seq for sequential proposing, block for blocked proposing, blocks based on var-covar matrix - block proposing only available when using adaptiveMCMC
#' @param startAdapt starting iteration for adapting
#' @param adaptBurn burn n number of iterations for defining var-covar matrix in block proposing
#' @param acceptanceRate acceptance rates, for adaptive sequential use string, one rate for each param, for block use single value
#' @param tell print monitoring information every x number of iterations
#' @param cluster T/F if true, using dide cluster, if false running locally
#' @return mcmc

mcmcSampler <- function(initParams,
                        randInit = F,
                        proposer = sequential.proposer,
                        sdProps,
                        maxSddProps,
                        niter = 100,
                        particleNum = 100,
                        nburn = 100,
                        monitoring = 2,
                        adaptiveMCMC = T,
                        proposerType = 'seq',
                        startAdapt = 1000,
                        adptBurn = 200,
                        acceptanceRate = 0.25,
                        tell = 5,
                        cluster = F,
                        oDat,
                        stoch=T,
                        priorFunc=NA,
                        switch=2500,
                        switchBlock=2500,
                        likelihoodFunc=likelihoodFunc1,
                        juvenileInfection=F) {


if(is.null(initParams$lFunc)!=T )likelihoodFunc<-initParams$lFunc else print("No likelihood function in param list, using default from function parameters")
  assump=initParams$assump
  birthType=initParams$birthType
  modNum<-initParams$modNum
  initParams<-initParams[c(1:31)] #remove extra functions from initial parameter list
  sdProps<-sdFunc(initParams,F)
  maxSddProps<-sdFunc(initParams,T)

  stepFun = modStepJI#ifelse(juvenileInfection==T,modStep,modStepJI)
  covar<-matrix(0,length(which(maxSddProps!=0)),length(which(maxSddProps!=0)))
  diag(covar)<-0.1
  assign('covar', covar, envir = .GlobalEnv)
  ##glmm for predicting shedding
  oDat$Day<-day(oDat$Date)
  mod<-NULL#glmer(pcrPos ~ log10Ser+ (1|Age)+(1|Day), obsPcrDens, family="binomial")
  ###


  aratio <- rep(0, length(sdProps))#starting acceptance ratio
  acceptRseq <- rep(0, length(sdProps))
  iterR <- rep(0, length(sdProps))
  sdp <- sdProps

  prmNum <- 1 # which paramter is currently being fitted
  if (randInit==T) initParams <- initRand(initParams)
  currentParams <- as.data.frame(initParams)
  nfitted <- length(currentParams[maxSddProps!=0]) ## number fitted parameters
  iter <- 2 ## mcmc iteration started at 1 so already on 2
  accept <- 0 ## initialize proportion of iterations accepted
  acceptR <- 0 #number of accepts for aratio

  #calc mean beta if R0 is fixed
  currentParams$betaFXVal<-if(currentParams$betaFX==1) betaMean(currentParams) else 0
  ##calc scaling factor
  currentParams$c_Val<- scalarFunc(currentParams,currentParams$m_Val,currentParams$mu_Val,b,currentParams$omega_m_Val,currentParams$mj_Val,currentParams$s_Val,phi = 0)



    ## Calculate log likelihood for first value
  if(stoch==T){
    curValX<- pFilt(likelihoodFunc=likelihoodFuncBoonahStoch,n=particleNum,iState,stepFun=stepFun,oDat,
                   currentParams,assump=assump,birthType = birthType,mat=T)
  }else {
    curValX <- as.vector(detModFunc(params=currentParams,obsData=oDat,assump=assump,likelihoodFunc=likelihoodFunc,birthType=birthType))
  }

  #curvalX is pointwise output
  curVal<-sum(curValX)+priorFunc(currentParams)

  curVal[is.na(curVal)]<--10000000
  ## Initialize matrix to store MCMC chain
  out <- matrix(NA, nr = niter, nc = length(currentParams) + 1+(nrow(obsData)-1))
  out[1,] <- c(as.numeric(currentParams), ll = curVal+priorFunc(currentParams),curValX) ## add first value
  colnames(out) <- c(names(currentParams), 'll',c(paste0("pW",(1:(nrow(obsData)-1))))) ## name columns
    originalCovar <-
    get('covar', envir = .GlobalEnv)## Store original covariance matrix

  if (adaptiveMCMC == T & proposerType == 'seq') acceptanceRate = rep(acceptanceRate,length(currentParams))
  if (adaptiveMCMC == T & proposerType == 'block') acceptanceRate = rep(acceptanceRate,length(currentParams))


  while (iter <= niter) {



if (iter >= switch)  stoch= T else stoch =F
if (iter == switch)  curVal= -1000000 else curVal =curVal
if (iter >= switchBlock)  proposerType= "block" else proposerType ="seq"
if (iter >= switchBlock)  proposer= multiv.proposer else proposer = sequential.proposer

    ##var covar matrix update - currently every 50 iterations
    if (adaptiveMCMC == T &
        proposerType == 'block' & iter > startAdapt & iter %% 50 == 0) {
      ##modulur division of 50, update covar every 50 iterations
      adptBurn <- min((startAdapt - 50), adptBurn)
      assign('out', out, envir = .GlobalEnv)
      adaptedCovar <-
        (2.38 ^ 2 / nfitted) * cov(log(out[adptBurn:(iter - 1), which(maxSddProps!=0)]))
      adaptedCovar[is.na(adaptedCovar)]<-0

      assign('adaptedCovar', adaptedCovar, envir = .GlobalEnv)
      adaptedCovar <-
        adaptedCovar * .95 + originalCovar * .05 ## 95% adapted & 5% original
      rownames(adaptedCovar) <-
        colnames(adaptedCovar) <-names(currentParams[maxSddProps!=0])
      assign('covar', adaptedCovar, envir = .GlobalEnv)

    }


    if (adaptiveMCMC == T & proposerType == 'block') {
      sdp[prmNum] <-
        tunerSeq(sdp, acceptanceRate, aratio, maxSddProps, prmNum)
      proposal <-
        proposer(currentParams, sdProps = sdp)
    }

    if (adaptiveMCMC == T & proposerType == 'seq') {
      if (iter >= startAdapt)
        sdp[prmNum] <-
          tunerSeq(sdp, acceptanceRate, aratio, maxSddProps, prmNum)
      proposal <-
        proposer(currentParams, prmNum, sdp)

    }

    if (adaptiveMCMC == F & proposerType == 'seq') {
      proposal <-
        proposer(currentParams, prmNum, sdp)
    }

proposal$c_Val<- scalarFunc(proposal,proposal$m_Val,proposal$mu_Val,b,proposal$omega_m_Val,proposal$mj_Val,proposal$s_Val,phi = 0)
proposal$betaFXVal<-if(proposal$betaFX==1) betaMean(proposal) else 0


    if(proposerType=="block" || maxSddProps[prmNum]!=0){
      propValX <-
          if(stoch==T){
            tryCatch({pFilt(likelihoodFunc=likelihoodFuncBoonahStoch,n=particleNum,iState,stepFun=stepFun,oDat,
                            proposal,assump=assump,birthType = birthType,mat=T)}, error =function(ex){-10000000})
          }else {
            detModFunc(params=proposal,obsData=oDat,assump=assump,likelihoodFunc=likelihoodFunc,birthType = birthType)
            }

      propVal<-sum(propValX)+priorFunc(proposal)
      lmh <-
        propVal - curVal ## likelihood ratio = log likelihood difference

      lmh = if(is.na(lmh)==T) -1e10 else lmh
        ## if it's not NA then do acception/rejection algorithm
        if ((lmh >= 0) | (runif(1, 0, 1) <= exp(lmh))) {
          currentParams <- proposal
          if (adaptiveMCMC == T & proposerType == 'block') {
            acceptRseq <- acceptRseq + 1
          } else{
            acceptRseq[prmNum] <- acceptRseq[prmNum] + 1
          }
          if (iter > nburn)
            accept <- accept + 1 ## track acceptance after burn-in
          curVal <- propVal
          curValX<-propValX
        }


      out[iter, ] <- c(as.numeric(currentParams), ll = curVal,as.vector(curValX))
      iter <- iter + 1

    if ((monitoring > 1 && iter %% tell == 0)){
      print(paste0("current likelihood = ",curVal))
      print(paste0("proposed likelihood = ",propVal))
      print(paste0("Assumption = ",assump))
      print(paste0("ModNum = ",modNum))

      print(sdp)
      print('sd Vals#############')
      print(aratio)
      print('aratio#############')
      print(paste0("prm num = ",prmNum))
      print('Current Params#############')
      print(currentParams)
      print('Proposed Params#############')
      print(proposal)
      print(paste0("prm num = ",prmNum))
      print(paste0("iteration  = ",iter," out of ",niter))

    }


    if (adaptiveMCMC == T & proposerType == 'block') {
      aratio <-
        acceptRseq/ (iter)
      aratio<-ifelse(aratio>1,1,aratio)
    } else {
      iterR[prmNum] <- iterR[prmNum] + 1
      if(iter>=startAdapt){
        aratio[prmNum] <-
          acceptRseq[prmNum] / (iterR[prmNum])#acceptance ratio change for specific parameter number
      }

    }

    }
    prmNum <- prmNum + 1#progress parameter number
    if (prmNum > length(sdProps))
      prmNum <-1#if parameter number reaches end of parameters, switch back to start


    if(iter %% 10000 == 0) write.csv(as.mcmc(out[1:nrow(out) > (nburn + 1),]),paste0("/home/aaron/res_",modNum,".csv"))
    }


  colnames(out) <- c(names(currentParams), 'll',c(paste0("pW",(1:(nrow(obsData)-1)))))
  results <- as.mcmc(out[1:nrow(out) > (nburn + 1),])
  return(list(
    initParams = initParams
    ,
    aratio = aratio
    ,
    results = results
  ))
}

