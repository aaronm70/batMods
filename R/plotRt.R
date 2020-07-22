

plotAllRt <- function(fileLoc,Rt = F,saveLoc,prmFile,burn,modNums) {
  for (i in modNums) {

    print(i)

    gg<-readResFunc(fileLoc=fileLoc,i=i,burn=burn,prmFile=prmFile)

    fxdY <- if (prmFile$lFunc[i] == "prFxOsc") T  else  F
    birthType <- if (prmFile$birthType[i] == "noImmune") 1 else 0
    g1 <-
      modPlotFuncX(
        gg,
        assum = prmFile$assump[i],
        fxd = fxdY,
        iters = 10,
        pName = i,
        birthType = birthType,
        Rt = Rt,
        saveLoc=saveLoc
      )


    assign(paste("G", i, sep = ""), g1)

  }

if(Rt==F) stop

  ggRt <-
    ggarrange(G4 + ggtitle("SILI (Mat. Immunity) (EF)"),
              G3 + ggtitle("SILI (EF)"),
              G6 + ggtitle("SIR (EF)"),
              G8 + ggtitle("SIRS (EF)"),
              G2 + ggtitle("SILI (Mat. Immunity)"),
              G1 + ggtitle("SILI"),
              G5 + ggtitle("SIR"),
              G7 + ggtitle("SIRS"),
              common.legend = T,legend="right"
              )
  ggsave(
    paste0(saveLoc,"RtAll_plot.png"),
    plot = ggRt,
    width = 18,
    height = 13
  )


}
