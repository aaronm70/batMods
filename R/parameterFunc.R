
#                                                   __                         ___                         __
#                                                  /\ \__                    /'___\                       /\ \__  __
# _____      __     _ __    __      ___ ___      __\ \ ,_\    __   _ __     /\ \__/  __  __    ___     ___\ \ ,_\/\_\    ___     ___     ____
#/\ '__`\  /'__`\  /\`'__\/'__`\  /' __` __`\  /'__`\ \ \/  /'__`\/\`'__\   \ \ ,__\/\ \/\ \ /' _ `\  /'___\ \ \/\/\ \  / __`\ /' _ `\  /',__\
#\ \ \L\ \/\ \L\.\_\ \ \//\ \L\.\_/\ \/\ \/\ \/\  __/\ \ \_/\  __/\ \ \/     \ \ \_/\ \ \_\ \/\ \/\ \/\ \__/\ \ \_\ \ \/\ \L\ \/\ \/\ \/\__, `\
# \ \ ,__/\ \__/.\_\\ \_\\ \__/.\_\ \_\ \_\ \_\ \____\\ \__\ \____\\ \_\      \ \_\  \ \____/\ \_\ \_\ \____\\ \__\\ \_\ \____/\ \_\ \_\/\____/
#  \ \ \/  \/__/\/_/ \/_/ \/__/\/_/\/_/\/_/\/_/\/____/ \/__/\/____/ \/_/       \/_/   \/___/  \/_/\/_/\/____/ \/__/ \/_/\/___/  \/_/\/_/\/___/
#   \ \_\
#    \/_/
#

#annual juvenile mortality rates of 0.57 for females and 0.47 for males,
#given that females reach adult size at 14.8 months of age and males at 16.3 months of age
#Vardon, Michael J., and Christopher R. Tidemann. "The black flying-fox (Pteropus alecto) in north Australia:
#juvenile mortality and longevity." Australian Journal of Zoology 48.1 (2000): 91-97.

#'parameter function
#'@param x any parameter
#'@return paramerter values
paramsFunc <-
  function(#parameters
    gamma_1_Val = 0,
    # clearance rate I->S
    gamma_2_Val =        0,
    # clearance rate I->R
    zeta_s =   rnorm(1,0.2,0.1),
    # clearance rate E->S
    sigma_2_Val =     0,
    # clearance rate E->R
    mu_Val = 1.37,# 2.27 for straw col bats,
    # juvenile maturation rate
    mj_Val =   rnorm(1,0.5,0.1), #0.796 for straw coloureds, split to male and female?
    # juvenile death rate
    m_Val =  rnorm(1,0.186,0.05) ,#cant get a good value in the lit so putting loose prior on it
    # adult death rate
    omega_m_Val = rnorm(1,0.8,0.1),
    # maternal antibody waning rate
    omega_2_Val =    0,
    # immune waning rate
    epsilon_Val =    0 ,
    # incubation/recurrance rate E->I
    kappa_Val =      rnorm(1,mean=3.60206,sd=0.1)  ,
    # carrying capacity
    c_Val = 2.113776  , # p.alecto 0.4? 1.53 for ghana,
    #birth pulse scalar
    s_Val = rnorm(1,130,10)  , #14.3 for straw col,
    #birth pulse synchronicity
    phi_Val = rnorm(1,7.18,0.1), #7.18 for p.alecto mid october
    #birth pulse timing
    rho_Val =  0,
    #latency I->E
    R0_Val = runif(1,1,25),

    Phi2_val= 1,
    #proportion of population sampled
    oDist_s = 2  ,
    #overdispersion parameter
    zeta_p =   runif(1,0.1,1),
    pcrProb2= runif(1,0.1,1),
    sigmaVer=ifelse(zeta_s>0,1,2),
    gammaVer=ifelse(gamma_1_Val>0,1,2),
    betaVer=2,
    S2_val =runif(1,1,25) ,
    d_val =   runif(1,1,25),
    oDist1=2,
    oDist_u=runif(1,0,3),
    c_val2=1,
    envOscType=0,
    betaFX=0,
    betaFXVal=0

    )
list(
  gamma_1_Val = gamma_1_Val,
  gamma_2_Val = gamma_2_Val,
  zeta_s = zeta_s,
  sigma_2_Val = sigma_2_Val,
  mu_Val = mu_Val,
  mj_Val = mj_Val,
  m_Val = m_Val,
  omega_m_Val = omega_m_Val,
  omega_2_Val = omega_2_Val,
  epsilon_Val = epsilon_Val,
  kappa_Val = kappa_Val,
  c_Val = c_Val,
  s_Val = s_Val,
  phi_Val = phi_Val,
  rho_Val = rho_Val,
  R0_Val=  R0_Val,
  Phi2_val = Phi2_val,
  oDist_s = oDist_s,
  zeta_p =zeta_p,
    sigmaVer=sigmaVer,
    gammaVer=gammaVer,
  betaVer=betaVer,
  S2_val=S2_val,
  d_val=d_val,
  pcrProb2=pcrProb2,
  oDist1=oDist1,
  oDist_u=oDist_u,
  c_val2=c_val2,
  envOscType= envOscType,
  betaFX=betaFX,
  betaFXVal=betaFXVal

)

#bbJ<-subset(batsBoonah,Age=="J")
#bbJ$day<-yday(bbJ$Date)
#bbJ<-aggregate(bbJ$Age~bbJ$day,FUN=length)
#plot(bbJ)
#res<-NULL
#for (i in 1:365){
#  b <-c * exp(-70*(cos(3.14*i - 4.5))^2)
#  res<-rbind(res,b)
#}
#lines(res,type="l")



