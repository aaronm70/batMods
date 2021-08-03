

plotAllRt <- function(resultsFile,Rt = F,saveLoc,prmFile,burn,modNums,thin=10) {
  for (i in modNums) {

    print(i)



    gg<-readResFunc(fileLoc=resultsFile,i=i,burn=burn,prmFile=prmFile,thin=thin)

    i<-if(i == 24) 8 else i

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
              G2 + ggtitle("SILI (Mat. Immunity)"),
              G1 + ggtitle("SILI"),
              G6 + ggtitle("SIR (EF)"),
              G5 + ggtitle("SIR"),

              G8 + ggtitle("SIRS (EF)"),

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
