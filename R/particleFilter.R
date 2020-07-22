


#  @@@@@@@    @@@@@@   @@@@@@@   @@@@@@@  @@@   @@@@@@@  @@@       @@@@@@@@      @@@@@@@@  @@@  @@@       @@@@@@@  @@@@@@@@  @@@@@@@
#  @@@@@@@@  @@@@@@@@  @@@@@@@@  @@@@@@@  @@@  @@@@@@@@  @@@       @@@@@@@@      @@@@@@@@  @@@  @@@       @@@@@@@  @@@@@@@@  @@@@@@@@
#  @@!  @@@  @@!  @@@  @@!  @@@    @@!    @@!  !@@       @@!       @@!           @@!       @@!  @@!         @@!    @@!       @@!  @@@
#  !@!  @!@  !@!  @!@  !@!  @!@    !@!    !@!  !@!       !@!       !@!           !@!       !@!  !@!         !@!    !@!       !@!  @!@
#  @!@@!@!   @!@!@!@!  @!@!!@!     @!!    !!@  !@!       @!!       @!!!:!        @!!!:!    !!@  @!!         @!!    @!!!:!    @!@!!@!
#  !!@!!!    !!!@!!!!  !!@!@!      !!!    !!!  !!!       !!!       !!!!!:        !!!!!:    !!!  !!!         !!!    !!!!!:    !!@!@!
#  !!:       !!:  !!!  !!: :!!     !!:    !!:  :!!       !!:       !!:           !!:       !!:  !!:         !!:    !!:       !!: :!!
#  :!:       :!:  !:!  :!:  !:!    :!:    :!:  :!:       :!:       :!:           :!:       :!:  :!:         :!:    :!:       :!:  !:!
#  ::       ::   :::  ::   :::     ::     ::   ::: :::   :: ::::   :: ::::       ::        ::   :: ::::     ::     :: ::::   ::   :::
#  :         :   : :   :   : :     :      :    :: :: :   : :: : :  : :: ::       :         :    : :: : :     :     : :: ::   :    : :


#' Particle filter function - runs model in steps between each observed data point
#' before using weighted re-sampling of model state at each step to start next
#' @param n number of particles
#' @param iState function for calculating initial state of particles
#' @param stepFun function for taking parmeters and time point, running model in steps and outputting model state
#' @param likeFunc likelihood function - obsolete
#' @param obsData observed data to fit model to
#' @param prms model parameters
#' @param resM True/False whether to output likelihood values or results of simulation for plotting figures
#' @param cluster T/F if true, using the cluster if false, run locally
#' @return if resM = F log likelihood value, if resM = T mean results of simulation
pFilt <-
  function (n,
            iState,
            stepFun,
            obsData,
            prms,
            resM = F,
            cluster = F,
            assump,
            likelihoodFunc,
            full = F,
            mat = F,
            birthType)
  {
    times = obsData$NumDays
    dayObs = obsData$Date
    prms$kappa_Val <- 10 ^ prms$kappa_Val

    particles = round(as.data.frame(iState(prms = prms, n = n))) #initial state

    ll = 0# only need start value if fitting initial conditions?
    llMat <- NULL
    rMeans <- NULL
    particlesFull <- NULL

    if(prms$gammaVer==0){
      prms$rho_Val<-ifelse(prms$rho_Val>-1,10^prms$rho_Val,-1)
      prms$epsilon_Val<-ifelse(prms$epsilon_Val>-1,10^prms$epsilon_Val,-1)
    }
    if(prms$gammaVer>0){
      prms$gamma_2_Val<-ifelse(prms$gamma_2_Val>-1,10^prms$gamma_2_Val,-1)
    }

    #calculate mu from juvenile lifespan and mat immune waning
    params$mu_Val<- 1/((15.55/12) - (1/params$omega_m_Val))


    mod <-
      stochSEIR(
        sigma_2_Val = prms$sigma_2_Val,
        gamma_2_Val = prms$gamma_2_Val,
        zeta_s = prms$zeta_s,
        gamma_1_Val = prms$gamma_1_Val,
        m_Val = prms$m_Val,
        rho_Val = prms$rho_Val,
        R0_Val = prms$R0_Val,
        phi_Val = prms$phi_Val,
        s_Val = prms$s_Val,
        c_Val = prms$c_Val,
        omega_m_Val = prms$omega_m_Val,
        omega_2_Val = prms$omega_2_Val,
        kappa_Val = prms$kappa_Val,
        epsilon_Val = prms$epsilon_Val,
        mj_Val = prms$mj_Val,
        betaVer = prms$betaVer,
        gammaVer = prms$gammaVer,
        sigmaVer = prms$sigmaVer,
        Phi2_val = prms$Phi2_val,
        S2_val  = prms$S2_val,
        c_val2  = prms$c_val2,
        envOscType = prms$envOscType,
        dt=365*4,
        birthType=birthType,
        betaFX=prms$betaFX,
        betaFXVal=prms$betaFXVal
      )

    for (i in 1:length(times[-length(times)])) {
      #Start/end time for model run
      startTime <- times[i]
      endTime <- times[i + 1]
      #use mapply with low particle numbers as overheads not worth it - currently using parralel for multi models
      if (n < 50) {
        particlesTemp = mapply(
          stepFun,
          particles = as.list(as.data.frame(t(particles))),
          MoreArgs = list(
            mod = mod,
            birthType = birthType,
            startTime = startTime,
            endTime = endTime
          ),
          SIMPLIFY = F
        )
      } else {
        particlesTemp = mcmapply(
          stepFun,
          particles = as.list(as.data.frame(t(particles))),
          MoreArgs = list(
            mod = mod,
            birthType = birthType,
            startTime = startTime,
            endTime = endTime
          ),
          SIMPLIFY = F,
          mc.cores = detectCores() - 1
        )
      }

      #format results of model runs to get DF of final states of each particle
      pTemp <- lapply(particlesTemp, tail, 1)
      pTemp <- as.data.frame(matrix(unlist(pTemp),nrow=length(pTemp),byrow=TRUE))
      pTemp <- pTemp[, -1]

      #weight particles using likelihood function based on final states of particles
      weights <- mapply(
        likelihoodFunc,
        S = round(rowSums(as.vector(pTemp[c(1:4)]))),
        E = round(rowSums(pTemp[c(6:9)])),
        I = round(rowSums(pTemp[c(10:13)])),
        Ia = round(rowSums(pTemp[c(10, 13)])),
        R = round(rowSums(pTemp[c(14:17, 5)])),
        MoreArgs = list(
          parms = prms,
          N = obsData$negatives[i + 1] + obsData$positives[i + 1],
          k_PpSp = obsData$PpSp[i + 1],
          k_PpSn = obsData$PpSn[i + 1],
          k_PnSp = obsData$PnSp[i + 1],
          k_PnSn = obsData$PnSn[i + 1],
          PCR = obsData$pcrPos[i + 1],
          urPredDat = tryCatch({
            obsData$fit[i + 1]
          }, error = function(ex) {
            NULL
          }),
          time = dayObs[i + 1],
          assump = assump
        )
      )
      if (mat == T) llMat = rbind(llMat, mean(weights))#return matrix of pointwise likelihoods for loo

       ll = ll + mean(weights)#add mean likelihood value

      #normalise weights
      mxWeight <- max(weights)
      weights <- weights - mxWeight
      swP = sum(exp(weights))
      weights = exp(weights) / swP
      weights[is.na(weights)] <- 1e-50

      rows = sample(1:n, n, replace = TRUE, prob = weights)
      particles = pTemp[rows,]
      #add mean to full particle results
      if (full == T) {
        mostCom = as.numeric(names(sort(table(rows), decreasing = TRUE))[1])
        if (i == 1)
          particlesFull <- particlesTemp[rows]
        if (i > 1)
          for (c in 1:length(particlesTemp)) {
            particlesFull[[c]] <-
              rbind(particlesFull[[c]], particlesTemp[[rows[c]]])
          }
      }
    }

    if (mat == T) {
      return(as.vector(llMat))
    } else {
      if (full == T)
        return(particlesFull[[mostCom]])
      else
        return (ll)
    }


  }
