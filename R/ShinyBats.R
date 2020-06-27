#
#library("shiny")
#library("deSolve")
#library("cowplot")
#library("ggplot2")
#library("tidyverse")
#library("ggrepel")
#library("shinydashboard")
#
#
##betaPrior<-mcmapply(betaShapeFinder,obsData$negatives[-1]+obsData$positives[-1], obsData$positives[-1]/(obsData$negatives[-1]+obsData$positives[-1]),params$oDist1)
##
##hpd <- binom.bayes(
##  x = obsData$positives[-1], n = obsData$negatives[-1]+obsData$positives[-1], type = "central", conf.level = 0.9, tol = 1e-9,prior.shape1 =mean(unlist(betaPrior[1,])),prior.shape2 =mean(unlist(betaPrior[2,])) )
##print(hpd)
##
##obsDataX<-obsData[-1,]
##obsDataX$upper<-hpd$upper
##obsDataX$lower<-hpd$lower
#
#
## Define UI
##
#
#ui <- dashboardPage(
#  dashboardHeader(disable = F,title="ShinyBats"),
#  dashboardSidebar(
#    selectInput("BetaVer", h3("BetaVer"),
#                choices = list(1, 2), selected = 1
#                ),
#    selectInput("envOscType", h3("envOscType"),
#                choices = list(0, 1), selected = 0
#    ),
#    sliderInput("epsilon",
#                "Epsilon:",
#                min = 0, max = 3, value = 1.8,step=0.1
#    ),
#    sliderInput("phi",
#                "phi:",
#                min = 0, max = 50, value = 7.8,step=0.1
#    ),
#    sliderInput("s",
#                "s:",
#                min = 0, max = 200, value = 130,step=0.5
#    ),
#    sliderInput("c",
#                "c:",
#                min = 0, max = 50, value = 15.9,step=0.1
#    ),
#    sliderInput("pV",
#                "pV:",
#                min = 0, max = 50, value = 7.8,step=0.1
#    ),
#    sliderInput("sV",
#                "sV:",
#                min = 0, max = 200, value = 130,step=0.5
#    ),
#    sliderInput("cV",
#                "cV:",
#                min = 0, max = 50, value = 15.9,step=0.1
#    ),
#    sliderInput("omega",
#                "omega:",
#                min = 0, max = 10, value = 0.8,step=0.1
#    ),
#    sliderInput("gamma",
#                "gamma:",
#                min = 0, max = 3, value = 0,step=0.1
#    ),
#    sliderInput("rho",
#                "Rho:",
#                min = 0, max = 3, value = 1,step=0.1
#    ),
#    sliderInput("R0",
#                "R0:",
#                min = 0.5, max = 50, value = 10.1,step=0.1
#    ),
#    sliderInput("pcrProb",
#                "zeta_1:",
#                min = 0.01, max = 1, value = 0.1
#    ),
#    sliderInput("m_Val",
#                "m_Val:",
#                min = 0.01, max = 1, value = 0.1
#    ),
#    sliderInput("kappa",
#                "kappa:",
#                min = 1, max = 5, value = 3.59,step=0.1
#    ),
#    sliderInput("S_osc",
#                     "S_osc:",
#                     min = 0.01, max = 20, value = 0.2,step=0.1
#    ),
#    sliderInput("Osc_Phi",
#                "Osc_Phi:",
#                min = 0.01, max = 10, value = 1.8,step=0.1
#    ),
#    sliderInput("oDist1",
#                "Overdispersion_indv:",
#                min = 0.01, max = 100, value = 73
#    )
#
#  ),
#  dashboardBody(
#    tags$head(tags$style(HTML('
#                              /* body */
#                              .content-wrapper, .right-side {
#                              background-color: #fffff8;
#                              }
#                              '))),
#
#    #    mainPanel(
#    titlePanel("ShinyBats", windowTitle = "ShinyBats"),
#    fluidRow(plotOutput("distPlot")),
#    br(),
#    br(),
#    fluidRow(plotOutput("IndvPlot")),
#    br()
#  )
#)
#
##
## Define server
##
#server <- function(input, output) {
#  # Create reactive input
#  dataInput <- reactive({
#    init       <-
#      c(
#        S = 1000,
#        I = 10,
#        R = 1,
#        E = 1
#      )
#    ## beta: infection parameter; gamma: recovery parameter
#
#    params <-
#      paramsFunc(
#        epsilon_Val = 10^input$epsilon,
#        # incubation/recurrance rate E->I
#        kappa_Val = ifelse(10^input$kappa>1,10^input$kappa,0),
#        rho_Val = ifelse(10^input$rho>1,10^input$rho,0),
#        gamma_2_Val  = ifelse(10^input$gamma>1,10^input$gamma,0),
#        R0_Val = input$R0,
#        pcrProb =input$pcrProb,
#        pcrProb2=input$pcrProb2,
#        m_Val=input$m_Val,
#        omega_2_Val = input$omega,
#        betaVer = as.numeric(input$BetaVer),
#        phi_Val = input$phi,
#        c_Val = input$c,
#        s_Val = input$s,
#        envOscType=as.numeric(input$envOscType),
#        lambda_val  = input$sV,
#        P_val = input$pV,
#        c_val2  = input$cV
#        )
#
#
#    ## Time frame
#    times  <- c(0, 75080)
#
#    ## Solve using ode (General Solver for Ordinary Differential Equations)
#    out <-detModFuncShiny(params,times,assump=NULL,likelihoodFunc=NULL)
#    #    out
#   as.data.frame(out)
#  })
#
#  output$distPlot <- renderPlot({
#    datF<- dataInput()[seq(1, nrow(dataInput()), 4),-c(1:4)]
#    datF$time<-datF$time
#    datF$time<-datF$time/4
#    datF2<-merge(datF,obsDataX,by.x="step",by.y="NumDays",all=T)
#    datF2<-datF2[,c(1:7,12,17,18)]
#
#    out <-
#      datF %>%
#      gather(key, value, -time) %>%
#      mutate(
#        id = row_number(),
#        key2 = recode(
#          key,
#          S = "Susceptible (S)",
#          I = "Infected (I)",
#          R = "Recovered (R)",
#          L = "Latent (L)"
#        ),
#        keyleft = recode(
#          key,
#          S = "Susceptible (S)",
#          I = "",
#          R = "",
#          L = "Latent (L)"
#        ),
#        keyright = recode(
#          key,
#          S = "",
#          I = "Infected (I)",
#          R = "Recovered (R)",
#          L = ""
#        )
#      )
#
#    out1<-merge(out,datF2,by.x="time",by.y="time")
#    out1<-out1[-c(1:245),]
#    out1$step<-(out1$step-(50*(365*4)))/4
#
#    out1$Date<- as.Date( out1$step, origin = "2012-12-31")
#
#    ggplot(data = out1,
#           aes(
#             x = Date,
#             y = value,
#             group = key2,
#             col = key2,
#             label = key2,
#             data_id = id
#           )) +
#      ylab("Proportion of population sero+") + xlab("Time (days)") +
#      geom_line(size = 1.5,alpha=0.9) +
#      geom_point(y=out1$prev,col="red")+
#      geom_errorbar(aes(ymin=out1$lower, ymax=out1$upper), width=.2,
#                    position=position_dodge(0.05),alpha=0.5,colour="blue")+
#      geom_text_repel(
#        data = subset(out1, time == max(time)),
#        aes(label = keyright),
#        size = 6,
#        segment.size  = 0.2,
#        segment.color = "grey50",
#        nudge_x = 0,
#        hjust = 1,
#        direction = "y"
#      ) +
#      geom_text_repel(
#        data = subset(out1, time == min(time)),
#        aes(label = keyleft),
#        size = 6,
#        segment.size  = 0.2,
#        segment.color = "grey50",
#        nudge_x = 0,
#        hjust = 0,
#        direction = "y"
#      ) +
#      theme(legend.position = "none") +
#      scale_color_viridis_d()+
#      scale_y_continuous(labels = scales::percent, limits = c(0, 1))
#  })
#
#
#
#  output$IndvPlot<- renderPlot({
#    out <-
#    dataInput()[,c(1:4,10)]
#
#  tt<-yday(obsData$Date[-1])/365
#  sDrive<-0
#
#    s=input$S_osc
#    c=1
#    phi=input$Osc_Phi
#    sDrive <-c* exp(-s*(cos(pi*tt- phi))^2)
#
#    simDat3<-merge(out,obsData,by.x="step",by.y="NumDays",all=T)
#
#  simDat4<-simDat3[is.na(simDat3$Date)==F,]
#
#  Ia<-as.vector(round((simDat4[-1,c(4)])))
#  Ia[is.na(Ia)]<-0
#
#  prb<-NULL
#  prb2<-NULL
#  for(i in c(1:50)){
#    prb<- mapply(rmutil::rbetabinom, size = Ia, m=input$pcrProb*(1-sDrive),MoreArgs = list(n=1,s=input$oDist1))/as.vector(round(rowSums(simDat4[-1,c(2:5)])))
#    prb2<-rbind(prb,prb2)
#  }
#
#  fullData2<-as.data.frame(obsData$Date[-1])
#  fullData2$resMean<-colMeans(prb2)
#  names(fullData2)<-c("Date","resMean")
#
#  biCon2<-binom.bayes(obsData$pcrPos,(obsData$positives+obsData$neg))
#
#  obsData$biConLwr<-biCon2$lower
#  obsData$biConUpr<-biCon2$upper
#  fullData3<-merge(fullData2,obsData,by.x="Date",by.y="Date",all=T)
#  boundsPCR1<- apply(prb2,2,function(x) quantile(x, c(0.05,0.95),na.rm=T))
#
#
#  prb2Melt<-as.data.frame(as.vector(t(prb2)))
#  prb2Melt$Date<-obsData$Date[-1]
#  names(prb2Melt)<-c("value","Date")
#  fullDataX<-merge(prb2Melt,obsData,by.x="Date",by.y="Date",all=F)
#
#ggplot(fullDataX,aes(x=Date,y=pcrPos/(positives+negatives)))+
#    geom_point(aes(y=value, x=Date),col="red",se=F,alpha=0.1)+
#    geom_errorbar(aes(ymin=biConLwr, ymax=biConUpr), width=.2,
#                  position=position_dodge(0.05),alpha=0.01,colour="blue")+
#    geom_point(col="black",alpha=0.01)+
#    ylab("indv. prevalence")})
#
#
#}
#
## Run the application
#shinyApp(ui = ui, server = server)
#
