


plotLooic <- function(fileNme,
                      loicVals,
                      loicVals2 = NULL,
                      type = "loo") {
  loicVals$Structure[loicVals$Structure == 1] <- "EF"
  loicVals$Structure[loicVals$Structure == 2] <- "No EF"
  loicVals$Structure[loicVals$Structure == 3] <- "DC"
  loicVals$modType[loicVals$modType == 3] <- "SIRS"
  loicVals$modType[loicVals$modType == 1] <- "SILI"

 # loicVals2$Structure[loicVals2$Structure == 1] <- "EF"
 # loicVals2$Structure[loicVals2$Structure == 2] <- "No EF"
 # loicVals2$Structure[loicVals2$Structure == 3] <- "DC"
 # loicVals2$modType[loicVals2$modType == 3] <- "SIRS"
 # loicVals2$modType[loicVals2$modType == 2] <- "SIR"
 # loicVals2$modType[loicVals2$modType == 1] <- "SILI"

  loicVals[,1]<-as.numeric(loicVals[,1])
  loicVals[,2]<-as.numeric(loicVals[,2])
  if (type != "loo") loicVals[,3]<-as.numeric(loicVals[,3])

  if (type == "loo") {
    loicVals$seUp <- loicVals$Looic + loicVals$se
    loicVals$seLow <- loicVals$Looic - loicVals$se
  }


  loicVals$location <- c(4, 2, 6)
  #loicVals2$location <- c(4, 1, 5, 2, 7, 8, 10, 11)

  if (type == "loo")
    yName <- "Looic"
  if (type == "R0")
    yName <- expression("R"[0])
  if (type == "batCont")
    yName <- expression("d")
  if (type == "loo")
    yLim <- c(200, 350)
  if (type == "R0")
    yLim <- c(0, 50)
  if (type == "batCont")
    yLim <- c(0, 25)
  if (type == "zetaP"||type == "zetaS"||type == "zetaU")
    yLim <- c(0, 1)
  if (type == "zetaS")
    yName <- expression(zeta[s])
  if (type == "zetaP")
    yName <- expression(zeta[p])
  if (type == "zetaU")
    yName <- expression(zeta[u])

  g1 <- ggplot(loicVals, aes(y = Looic, x = location)) +
    ylim(yLim) +
    geom_errorbar(
      aes(ymin = seLow, ymax = seUp),
      width = 0.2,
      position = position_dodge(0.05),
      alpha = 1,
      colour = "#102F47"
    ) +
    geom_point(aes(colour = Structure), size = 5, alpha = 1) +
    xlab("") +
    ylab(yName) +
    scale_fill_viridis_d(option = "C", direction = -1) +
    scale_x_continuous(
      breaks = c(2, 4, 6),
      labels = c("SILI","SILI", "SIRS")
    ) +  theme_bw(base_size = 20) +
    panel_border(size = 0.5, color = "black") +
    scale_color_viridis_d(option = "D") +
    geom_rect(aes(
      xmin = 0,
      xmax = 3,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = 0.01,
    fill = "lightskyblue") +
    geom_rect(aes(
      xmin = 3,
      xmax = 6,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = 0.01,
    fill = "darkorange") +
    geom_rect(aes(
      xmin = 6,
      xmax = 9,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = 0.01,
    fill = "lightskyblue") +
    geom_rect(aes(
      xmin = 9,
      xmax = 12,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = 0.01,
    fill = "darkorange")

  ggsave(
    fileNme,
    plot = g1,
    height = 25,
    width = 40,
    units = "cm"
  )
  print(g1)
return(g1)

  }


