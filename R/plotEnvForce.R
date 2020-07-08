



##plot environmental force
plotEnvForc <- function(resultsFile,burn,modNums,saveLoc,prmFile) {
  for (i in modNums) {
    prmFile <- read.csv("/Users/alm204/Documents/ModelSetups.csv")

    gg<-readResFunc(fileLoc=resultsFile,i=i,burn=burn,prmFile=prmFile,thin=10)
   params<- colMedians(gg)
    params<-as.data.frame(matrix(params,ncol =length(params),byrow = T))
    names(params)<-names(gg)


    s = median(gg$S2_val)
    phi = median(gg$Phi2_val)
    c = median(gg$c_val2)
    x <- seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")
    x <- yday(x) / 365


    if (i == 3) {
      medVal <- 10 ^ median(gg$epsilon_Val)
      ciValsEp <- 10 ^ HPDinterval(as.mcmc(gg$epsilon_Val),prob=.89)

      sDriveEpOrig <- c * exp(-s * (cos(pi * x - 1)) ^ 2)
      sDriveEp <- medVal * sDriveEpOrig
      sDriveEp <- as.data.frame(sDriveEp)
      sDriveEp$date <-
        seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G3 <- ggplot(sDriveEp, aes(y = sDriveEp, x = date)) +
        geom_line(colour = "purple", size = 1.5) +
        theme_bw(base_size = 20) +
        scale_x_date(date_breaks = "2 months" , date_labels = "%b") +
        ylab(expression(epsilon[t] ~ (year ^ -1))) +
        geom_ribbon(
          aes(ymin = sDriveEpOrig * ciValsEp[1], ymax = sDriveEpOrig * ciValsEp[2]),
          fill = "purple",
          alpha = 0.3
        )

    }

    if (i == 4) {
      medVal <- 10 ^ median(gg$epsilon_Val)
      ciValsEp2 <- 10 ^ HPDinterval(as.mcmc(gg$epsilon_Val),prob=.89)

      sDriveEp2Orig <- c * exp(-s * (cos(pi * x - phi)) ^ 2)
      sDriveEp2 <- medVal * sDriveEp2Orig
      sDriveEp2 <- as.data.frame(sDriveEp2)
      sDriveEp2$date <-
        seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G4 <- ggplot(sDriveEp2, aes(y = sDriveEp2, x = date)) +
        geom_line(colour = "purple", size = 1.5) +
       theme_bw(base_size = 20) +
        scale_x_date(date_breaks = "2 months" , date_labels = "%b") +
        ylab(expression(epsilon[t] ~ (year ^ -1))) +
        geom_ribbon(
          aes(ymin = sDriveEp2Orig * ciValsEp2[1], ymax = sDriveEp2Orig * ciValsEp2[2]),
          fill = "purple",
          alpha = 0.3
        )

    }
    if (i == 6) {
      medValGam <- 10 ^ median(gg$gamma_2_Val)

      ciValsGam <- 10 ^ HPDinterval(as.mcmc(gg$gamma_2_Val),prob=.89)

      sDriveGamOrig <- c * exp(-s * (cos(pi * x - phi)) ^ 2)
      sDriveGam <- medValGam * sDriveGamOrig
      sDriveGam <- as.data.frame(sDriveGam)
      sDriveGam$date <-
        seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G6 <- ggplot(sDriveGam, aes(y = sDriveGam, x = date)) +
        geom_line(colour = "purple", size = 1.5) +
        theme_bw(base_size = 20) +
        scale_x_date(date_breaks = "2 months" , date_labels = "%b") +
        ylab(expression(gamma[t] ~ (year ^ -1))) +
        geom_ribbon(
          aes(
            ymin = sDriveGamOrig * ciValsGam[1],
            ymax = sDriveGamOrig * ciValsGam[2]
          ),
          fill = "purple",
          alpha = 0.3
        )

    }


    if (i == 8) {
      medValOm <- median(gg$omega_2_Val)

      ciValsOmeg <- HPDinterval(as.mcmc(gg$omega_2_Val),prob=.89)

      sDriveOmegOrig <- c * exp(-s * (cos(pi * x - phi)) ^ 2)
      sDriveOmeg <- medValOm * sDriveOmegOrig
      sDriveOmeg <- as.data.frame(sDriveOmeg)
      sDriveOmeg$date <-
        seq(as.Date("2013-06-01"), as.Date("2014-06-01"), "days")

      G8 <- ggplot(sDriveOmeg, aes(y = sDriveOmeg, x = date)) +
        geom_line(colour = "purple", size = 1.5) +
        theme_bw(base_size = 20) +
        scale_x_date(date_breaks = "2 months" , date_labels = "%b") +
        ylab(expression(omega[t] ~ (year ^ -1))) +
        geom_ribbon(
          aes(
            ymin = sDriveOmegOrig * ciValsOmeg[1],
            ymax = sDriveOmegOrig * ciValsOmeg[2]
          ),
          fill = "purple",
          alpha = 0.3
        )



    }


  }
  ggRt <-
    ggarrange(G4 + ggtitle("SILI (Maternal Immunity)"),
              G3 + ggtitle("SILI"),
              G6 + ggtitle("SIR"),
              G8 + ggtitle("SIRS"),
              common.legend = T)
  print(ggRt)
  ggsave(
   paste0(saveLoc,"envForceVals.png"),
    plot = ggRt,
    width = 18,
    height = 13
  )

}
