
## Initial states
initial(Sn) <- Sn_ini
initial(Sj) <- Sj_ini
initial(Sm) <- Sm_ini
initial(Sf) <- Sf_ini
initial(Ma) <- Ma_ini

initial(En) <- En_ini
initial(Ej) <- Ej_ini
initial(Em) <- Em_ini
initial(Ef) <- Ef_ini

initial(In) <- In_ini
initial(Ij) <- Ij_ini
initial(Im) <- Im_ini
initial(If) <- If_ini

initial(Rn) <- Rn_ini
initial(Rj) <- Rj_ini
initial(Rm) <- Rm_ini
initial(Rf) <- Rf_ini
initial(timeOsc)<-timeOscIni

#parameters
Sn_ini <- user(20) # susceptibles
Sj_ini <- user(20) # susceptibles
Sf_ini <- user(20) # susceptibles
Sm_ini <- user(20) # susceptibles
Ma_ini<-user(0)
En_ini <- user(5) # infected
Ej_ini <- user(5) # infected
Ef_ini <- user(5) # infected
Em_ini <- user(5) # infected

In_ini<- user(0)
Ij_ini<- user(0)
Im_ini<- user(0)
If_ini<- user(0)

Rn_ini<- user(0)
Rj_ini<- user(0)
Rm_ini<- user(0)
Rf_ini<- user(0)
timeOscIni<-user(0)

gamma_1_Val <- user(0.8) # clearance rate I->S
gamma_2_Val <- user(0) # clearance rate I->R
zeta_s <-user(0.8) # clearance rate E->S
sigma_2_Val <-user(0) # clearance rate E->R
mu_Val <- user(0.44) # juvenile maturation rate
mj_Val <- user(0.796) # juvenile death rate
m_Val <- user(0.186) # adult death rate
omega_m_Val <- user(0.4) # maternal antibody waning rate
omega_2_Val <-user(0.4) # immune waning rate
epsilon_Val <- user(0.9) # incubation/recurrance rate E->I
kappa_Val <- user(1000)# carrying capacity
c_Val <- user(1.53) #birth pulse scalar
s_Val <- user(14.3) #birth pulse synchronicity
phi_Val <- user(4.5) #birth pulse timing
rho_Val <- user(0.9) #latency I->E
R0_Val <-user(0.5)
betaVer<-user(0)
gammaVer<-user(0)
sigmaVer<-user(0)
envOscType<-user(0)
initial(R0_out)<-0
birthType<-user(0)
betaFX<-user(0)
betaFXVal<-user(0)
##sDrive component
S2_val<-user(0)
c_val2<-user(0)
Phi2_val<-user(0)
sV<-S2_val
cV<-c_val2
phiV<-Phi2_val
sDrive <- if(envOscType==1) cV* exp(-sV*(cos(3.141593*timeOsc- phiV))^2) else 1

epsilon_ValS<-if (epsilon_Val>0) epsilon_Val*sDrive else epsilon_Val
omega_2<-if (omega_2_Val>0) (omega_2_Val*sDrive)/dt else omega_2_Val/dt
gamma_2_ValS<-if(betaVer==1)(gamma_2_Val*sDrive) else gamma_2_Val

###

dt<-user(365)
## Total population size, and number of infected bats
N <- Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If + Rn +
  Rj + Rm + Rf + Ma

# beta_Val=((R0_Val*((m_Val^2)+(rho_Val*m_Val)+(sigma_2_Val*m_Val)+(m_Val*gamma_1_Val)+
#                      (m_Val*gamma_2_Val)+(m_Val*sigma_2_Val)+(m_Val*epsilon_ValS)+
#                      (rho_Val*sigma_2_Val)+(sigma_2_Val*gamma_1_Val)+(sigma_2_Val*gamma_2_Val)+
#                      (rho_Val*zeta_s)+(gamma_1_Val*zeta_s)+(gamma_2_Val*zeta_s)+
#                      (gamma_1_Val*epsilon_ValS)+(gamma_2_Val*epsilon_ValS))/(N*epsilon_ValS)))
#
beta_2x<-if(betaFX==0) R0_Val*((epsilon_ValS+sigma_2_Val+m_Val)*(gamma_2_ValS+m_Val+rho_Val)
                               -epsilon_ValS*rho_Val)/(N*(epsilon_ValS+sigma_2_Val+m_Val)) else betaFXVal


#parameters
gamma_1 <- if(gammaVer==1 )gamma_1_Val/dt else 0
gamma_2 <- if(gammaVer==2 )gamma_2_ValS/dt else 0
sigma_1 <- 0#if(sigmaVer==1 )zeta_s/dt else 0
sigma_2 <- if(sigmaVer==2 )sigma_2_Val/dt else 0
beta_1<- 0#if(beta_1x>1) 1 else beta_1x
beta_2 <-beta_2x/dt

mu <- mu_Val/dt
mj <- mj_Val/dt
m <-  m_Val/dt
omega_m <-omega_m_Val/dt
epsilon <-epsilon_ValS/dt
rho <- rho_Val/dt

kappa <- kappa_Val
c <- c_Val
s <- s_Val
phi <- phi_Val

update(R0_out)<- beta_2*((Sn+Sj+Sf+Sm)*(In+Ij+If+Im)) #((beta_2x*N)*(epsilon_ValS+m_Val))/
  #((epsilon_ValS+m_Val)*(gamma_2_ValS+m_Val+rho_Val)-epsilon_ValS*rho_Val)

update(timeOsc)<-timeOsc+1/dt

b <-c * exp(-s*(cos(3.141593*timeOsc - phi))^2)

## transitions between neonatal compartments:
Sn_births <- if(birthType==1) (b*(Sf + If+Ef) )/dt else (b*(Sf + If) )/dt
Ma_births <- if(birthType==1) (b*(Rf))/dt else (b*(Rf+Ef))/dt


SNB<-(Sn_births)
MAB<-(Ma_births)

update(Sn) <- if ((Sn + SNB - B_snSumBIN + n_En_Sn + n_Rn_Sn+n_In_Sn)<=0) 0 else (Sn + SNB - B_snSumBIN + n_En_Sn + n_Rn_Sn+n_In_Sn)
update(En) <- if ((En - B_enSumBIN + n_Sn_En+n_In_En)<=0) 0 else (En - B_enSumBIN + n_Sn_En+n_In_En)
update(In) <- if (( In - B_inSumBIN+ n_Sn_In + n_En_In)<=0) 0 else (In - B_inSumBIN+ n_Sn_In + n_En_In)
update(Rn) <- if ((Rn - B_RnSumBIN+ n_En_Rn+n_In_Rn)<=0) 0 else (Rn - B_RnSumBIN+ n_En_Rn+n_In_Rn)
update(Ma) <- if ((Ma + MAB - B_maSumBIN)<=0) 0 else (Ma + MAB - B_maSumBIN)

##Transition probabilities
P_death<-mj * (N / kappa)#prob juvenile death
P_deathA<-m * (N / kappa)#prob adult death
P_toE<-(beta_1 * (If + Im + Ij + In))#transmission to E
P_toI<-(beta_2 * (If + Im + Ij + In))#transmission to I
P_matWane<-omega_m
P_E_S<-sigma_1
P_E_R<-sigma_2
P_I_S<-gamma_1
P_I_R<-gamma_2
P_E_I<-  epsilon #Incubation/recurance
P_I_E<-rho
P_J_A<-mu#juvenile maturation rate
P_ImWane<-omega_2

#B_sn
B_snSum <- P_matWane + P_death + P_toE + P_toI
B_snPr <- 1-exp(-(P_matWane + P_death + P_toE + P_toI))
B_snSumX<-if (B_snPr>1) 1 else B_snPr
B_snSumBIN<-rbinom(Sn,B_snSumX)
B_sn[] <- rmultinom(B_snSumBIN, p)
p[1] <-  P_matWane/B_snSum
p[2] <- P_death/B_snSum
p[3] <- P_toE/B_snSum
p[4] <-  P_toI/B_snSum

dim(p) <- 4
dim(B_sn) <- 4
n_Sn_Sj <- B_sn[1]

n_Sn_En <- 0
n_Sn_In <- B_sn[3] + B_sn[4]

#n_Sn_En <- if (betaVer == 1)
#  B_sn[3] + B_sn[4] else  0
#n_Sn_In <- if (betaVer == 2)
#  B_sn[3] + B_sn[4] else 0



#B_en
B_en[] <- rmultinom(B_enSumBIN, p1)
B_enSum <- (P_matWane +P_death + P_E_S + P_E_R + P_E_I)
B_enPr <- 1-exp(-(P_matWane +P_death + P_E_S + P_E_R + P_E_I))

B_enSumX<-if (B_enPr>1) 1 else B_enPr
B_enSumBIN<-rbinom(En,B_enSumX)

p1[1] <-   P_matWane/B_enSum
p1[2] <-  P_death/B_enSum
p1[3] <- P_E_S/B_enSum
p1[4] <- P_E_R/B_enSum
p1[5] <- P_E_I/B_enSum

dim(p1) <- 5
dim(B_en) <- 5
n_En_Ej <- B_en[1]
n_En_In <- B_en[5]
n_En_Sn <- if (sigmaVer == 1)
  B_en[3] + B_en[4] else 0
n_En_Rn <- if (sigmaVer == 2)
  B_en[3] + B_en[4] else 0



#B_in
B_in[] <- rmultinom(B_inSumBIN, p2)
B_inSum <- P_matWane +P_death +  P_I_S + P_I_R + P_I_E
B_inPr<-1-exp(-(P_matWane +P_death +  P_I_S + P_I_R + P_I_E))
B_inSumX<- if (B_inPr>1) 1 else B_inPr
B_inSumBIN<-rbinom(In,B_inSumX)
p2[1] <-  P_matWane/B_inSum
p2[2] <-  P_death/B_inSum
p2[3] <-  P_I_S/B_inSum
p2[4] <-  P_I_R/B_inSum
p2[5] <-  P_I_E/B_inSum
dim(p2) <- 5
dim(B_in) <- 5
n_In_Sn <- if (gammaVer == 1)
  B_in[3] + B_in[4] else 0
n_In_Rn <- if (gammaVer == 2)
  B_in[3] + B_in[4] else 0
n_In_En <- B_in[5]
n_In_Ij <- B_in[1]

#B_rn
B_rn[] <- rmultinom(B_RnSumBIN, p3)
B_rnSum <- P_matWane +  P_death + P_ImWane
B_rnPr<-1-exp(-(P_matWane +  P_death + P_ImWane))
B_rnSumX <- if (B_rnPr>1) 1 else B_rnPr
B_RnSumBIN<-rbinom(Rn,B_rnSumX)

p3[1] <-  P_matWane/B_rnSum
p3[2] <-  P_death/B_rnSum
p3[3] <- P_ImWane/B_rnSum
dim(p3) <- 4
dim(B_rn) <- 4
n_Rn_Sn <- B_rn[3]
n_Rn_Rj <- B_rn[1]

#B_ma
#MaL<-if(Ma>0) Ma else 1
B_ma[] <- rmultinom(B_maSumBIN, p4)
B_maSum <- P_matWane + P_death
B_maPr <- 1-exp(-(P_matWane + P_death))
B_maSumX<-if(B_maPr>1) 1 else B_maPr
B_maSumBIN<-rbinom(Ma,B_maSumX)
p4[1] <- P_matWane/B_maSum
p4[2] <- P_death/B_maSum
dim(p4) <- 2
dim(B_ma) <- 2
n_Ma_Sj <- B_ma[1]


## transitions between juvenile compartments:
update(Sj) <- if ((Sj - B_sjSumBIN + n_Sn_Sj + n_Ma_Sj + n_Ej_Sj  + n_Rj_Sj+ n_Ij_Sj)<=0) 0 else (Sj - B_sjSumBIN + n_Sn_Sj + n_Ma_Sj + n_Ej_Sj  + n_Rj_Sj+ n_Ij_Sj)
update(Ej) <- if ((Ej -  B_ejSumBIN+ n_En_Ej + n_Sj_Ej + n_Ij_Ej)<=0) 0 else (Ej -  B_ejSumBIN+ n_En_Ej + n_Sj_Ej+n_Ij_Ej)
update(Ij) <-if (( Ij - B_ijSumBIN + n_In_Ij + n_Sj_Ij + n_Ej_Ij)<=0) 0 else (Ij - B_ijSumBIN + n_In_Ij + n_Sj_Ij + n_Ej_Ij)
update(Rj) <- if ((Rj - B_rjSumBIN + n_Rn_Rj + n_Ej_Rj + n_Ij_Rj)<=0) 0 else (Rj - B_rjSumBIN + n_Rn_Rj + n_Ej_Rj+n_Ij_Rj)



#B_sj
B_sj[] <- rmultinom(B_sjSumBIN, pj)
B_sjSum <-P_J_A +P_death + P_toE + P_toI
B_sjPr <-1-exp(-(P_J_A +P_death + P_toE + P_toI))

B_sjSumX<- if (B_sjPr>1) 1 else B_sjPr
B_sjSumBIN<-rbinom(Sj,B_sjSumX)
pj[1] <- P_J_A/B_sjSum
pj[2] <- P_death/B_sjSum
pj[3] <-  P_toE/B_sjSum
pj[4] <-  P_toI/B_sjSum
dim(pj) <- 4
dim(B_sj) <- 4
n_Sj_SmSf <- B_sj[1]

n_Sj_Ej <- 0
n_Sj_Ij <-B_sj[3] + B_sj[4]

#
# n_Sj_Ej <- if (betaVer == 1)
#   B_sj[3] + B_sj[4] else 0
# n_Sj_Ij <- if (betaVer == 2)
#   B_sj[3] + B_sj[4] else  0

#B_ej
B_ejSum <-P_J_A +P_death + P_E_S + P_E_R + P_E_I
B_ejPr <-1-exp(-(P_J_A +P_death + P_E_S + P_E_R + P_E_I))

B_ej[] <- rmultinom(B_ejSumBIN, pj1)
B_ejSumX <- if (B_ejPr>1) 1 else B_ejPr
B_ejSumBIN<-rbinom(Ej,B_ejSumX)
pj1[1] <- P_J_A/B_ejSum
pj1[2] <-  P_death/B_ejSum
pj1[3] <- P_E_S/B_ejSum
pj1[4] <- P_E_R/B_ejSum
pj1[5] <-  P_E_I/B_ejSum
dim(pj1) <- 5
dim(B_ej) <- 5
n_Ej_EmEf <- B_ej[1]
n_Ej_Ij <- B_ej[5]
n_Ej_Sj <- if (sigmaVer == 1)
  B_ej[3] + B_ej[4] else 0
n_Ej_Rj <- if (sigmaVer == 2)
  B_ej[3] + B_ej[4] else 0


#B_ij
B_ijSum <- P_J_A + P_death + P_I_S + P_I_R + P_I_E
B_ijPr <- 1-exp(-(P_J_A + P_death + P_I_S + P_I_R + P_I_E))

B_ij[] <- rmultinom(B_ijSumBIN, pj2)
B_ijSumX <- if (B_ijPr>1) 1 else B_ijPr
B_ijSumBIN<-rbinom(Ij,B_ijSumX)

pj2[1] <- P_J_A/B_ijSum
pj2[2] <- P_death/B_ijSum
pj2[3] <- P_I_S/B_ijSum
pj2[4] <- P_I_R/B_ijSum
pj2[5] <- P_I_E/B_ijSum
dim(pj2) <- 5
dim(B_ij) <- 5
n_Ij_ImIf <- B_ij[1]
n_Ij_Ej <- B_ij[5]
n_Ij_Sj <- if (gammaVer == 1)
  B_ij[3] + B_ij[4] else 0
n_Ij_Rj <- if (gammaVer == 2)
  B_ij[3] + B_ij[4] else 0

#B_rj
B_rj[] <- rmultinom(B_rjSumBIN, pj3)
B_rjSumX<-if(B_rjPr>1) 1 else B_rjPr
B_rjSumBIN<-rbinom(Rj,B_rjSumX)
B_rjSum <- P_J_A + P_death + P_ImWane
B_rjPr <- 1-exp(-(P_J_A + P_death + P_ImWane))

pj3[1] <-P_J_A/B_rjSum
pj3[2] <-P_death/B_rjSum
pj3[3] <- P_ImWane/B_rjSum
dim(pj3) <- 3
dim(B_rj) <- 3
n_Rj_RmRf <- B_rj[1]
n_Rj_Sj <- B_rj[3]


## transitions between adult male compartments:
update(Sm) <- if((Sm - B_smSumBIN + (n_Sj_SmSf/2) + n_Em_Sm + n_Im_Sm + n_Rm_Sm)<=0) 0 else (Sm - B_smSumBIN + (n_Sj_SmSf/2) + n_Em_Sm + n_Im_Sm + n_Rm_Sm)
update(Em) <- if((Em - B_emSumBIN + (n_Ej_EmEf/2) + n_Sm_Em + n_Im_Em)<=0) 0 else (Em - B_emSumBIN + (n_Ej_EmEf/2) + n_Sm_Em + n_Im_Em)
update(Im) <- if((Im - B_imSumBIN + n_Sm_Im + n_Em_Im+(0.5*n_Ij_ImIf))<=0) 0 else (Im - B_imSumBIN + n_Sm_Im + n_Em_Im+(0.5*n_Ij_ImIf))
update(Rm) <- if((Rm - B_rmSumBIN+ (n_Rj_RmRf/2) + n_Em_Rm + n_Im_Rm)<=0) 0 else (Rm - B_rmSumBIN+ (n_Rj_RmRf/2) + n_Em_Rm + n_Im_Rm)

#B_sm
B_sm[] <- rmultinom(B_smSumBIN, pm)
B_smSumX<- if (B_smPr>1) 1 else B_smPr
B_smSumBIN<-rbinom((Sm),B_smSumX)
B_smSum <-  P_deathA + P_toE + P_toI
B_smPr <-  1-exp(-(P_deathA + P_toE + P_toI))

pm[1] <- P_deathA/B_smSum
pm[2] <-  P_toE/B_smSum
pm[3] <- P_toI/B_smSum
dim(pm) <- 3
dim(B_sm) <- 3

n_Sm_Em <- 0
n_Sm_Im <- B_sm[2] + B_sm[3]

#
#n_Sm_Em <- if (betaVer == 1)
#  B_sm[2] + B_sm[3] else 0
#n_Sm_Im <- if (betaVer == 2)
#  B_sm[2] + B_sm[3] else 0

#B_em
B_em[] <- rmultinom(B_emSumBIN, pm1)
B_emSumX<- if (B_emPr>1) 1 else B_emPr
B_emSumBIN<-rbinom((Em),B_emSumX)
B_emSum <- P_deathA + P_E_S + P_E_R + P_E_I
B_emPr <- 1-exp(-(P_deathA + P_E_S + P_E_R + P_E_I))

pm1[1] <- P_deathA/B_emSum
pm1[2] <-  P_E_S/B_emSum
pm1[3] <- P_E_R/B_emSum
pm1[4] <- P_E_I/B_emSum
dim(pm1) <- 4
dim(B_em) <- 4
n_Em_Im <- B_em[4]
n_Em_Sm <- if (sigmaVer == 1)
  B_em[2] + B_em[3] else 0
n_Em_Rm <- if (sigmaVer == 2)
  B_em[2] + B_em[3] else 0

#B_im
B_im[] <- rmultinom(B_imSumBIN, pm2)
B_imSumX<- if (B_imSPr>1) 1 else B_imSPr
B_imSumBIN<-rbinom((Im),B_imSumX)
B_imSum <- P_deathA + P_I_S + P_I_R + P_I_E
B_imSPr <- 1-exp(-(P_deathA + P_I_S + P_I_R + P_I_E))

pm2[1] <- P_deathA/B_imSum
pm2[2] <- P_I_S/B_imSum
pm2[3] <-  P_I_R/B_imSum
pm2[4] <-  P_I_E/B_imSum
dim(pm2) <- 4
dim(B_im) <- 4
n_Im_Em <- B_im[4]
n_Im_Sm <- if (gammaVer == 1)
  B_im[2] + B_im[3] else 0
n_Im_Rm <- if (gammaVer == 2)
  B_im[2] + B_im[3] else 0

#B_rm
B_rm[] <- rmultinom(B_rmSumBIN, pm3)
B_rmSumX<- if (B_rmPr>1) 1 else B_rmPr
B_rmSumBIN<-rbinom((Rm),B_rmSumX)
B_rmSum <-P_deathA + P_ImWane
B_rmPr <-1-exp(-(P_deathA + P_ImWane))

pm3[1] <-  P_deathA/B_rmSum
pm3[2] <- P_ImWane/B_rmSum
dim(pm3) <- 2
dim(B_rm) <- 2
n_Rm_Sm <- B_rm[2]


## transitions between adult male compartments:
update(Sf) <- if((Sf - B_sfSumBIN+  (n_Sj_SmSf/2) + n_Ef_Sf + n_If_Sf + n_Rf_Sf)<=0) 0 else (Sf - B_sfSumBIN +  (n_Sj_SmSf/2) + n_Ef_Sf + n_If_Sf + n_Rf_Sf)
update(Ef) <- if((Ef - B_efSumBIN + (n_Ej_EmEf/2) + n_Sf_Ef + n_If_Ef)<=0) 0 else (Ef - B_efSumBIN + (n_Ej_EmEf/2) + n_Sf_Ef + n_If_Ef)
update(If) <- if((If - B_ifSumBIN + n_Sf_If + n_Ef_If+(0.5*n_Ij_ImIf))<=0) 0 else (If - B_ifSumBIN + n_Sf_If + n_Ef_If+(0.5*n_Ij_ImIf))
update(Rf) <- if((Rf - B_rfSumBIN + (n_Rj_RmRf/2) + n_Ef_Rf + n_If_Rf)<=0) 0 else (Rf - B_rfSumBIN + (n_Rj_RmRf/2) + n_Ef_Rf + n_If_Rf)

#B_s
B_sf[] <- rmultinom(B_sfSumBIN, pf)
B_sfSumX <- if (B_sfPr>1) 1 else B_sfPr
B_sfSumBIN<-rbinom((Sf),B_sfSumX)
B_sfSum <-P_deathA + P_toE + P_toI
B_sfPr <-1-exp(-(P_deathA + P_toE + P_toI))

pf[1] <- P_deathA/B_sfSum
pf[2] <-   P_toE/B_sfSum
pf[3] <-P_toI/B_sfSum
dim(pf) <- 3
dim(B_sf) <- 3
n_Sf_Ef <- 0
n_Sf_If <-  B_sf[2] + B_sf[3]

#
#n_Sf_Ef <- if (betaVer == 1)
#  B_sf[2] + B_sf[3] else 0
#n_Sf_If <- if (betaVer == 2)
#  B_sf[2] + B_sf[3] else 0

#B_ef
B_ef[] <- rmultinom(B_efSumBIN, pf1)
B_efSumX<- if (B_efPr>1) 1 else B_efPr
B_efSumBIN<-rbinom((Ef),B_efSumX)
B_efSum <- P_deathA + P_E_S + P_E_R + P_E_I
B_efPr <- 1-exp(-(P_deathA + P_E_S + P_E_R + P_E_I))

pf1[1] <- P_deathA/B_efSum
pf1[2] <- P_E_S/B_efSum
pf1[3] <- P_E_R/B_efSum
pf1[4] <- P_E_I/B_efSum
dim(pf1) <- 4
dim(B_ef) <- 4
n_Ef_If <- B_ef[4]
n_Ef_Sf <- if (sigmaVer == 1)
  B_ef[2] + B_ef[3] else 0
n_Ef_Rf <- if (sigmaVer == 2)
  B_ef[2] + B_ef[3] else 0

#B_if
B_if[] <- rmultinom(B_ifSumBIN, pf2)
B_ifSumX<- if (B_ifPr>1) 1 else B_ifPr
B_ifSumBIN<-rbinom((If),B_ifSumX)
B_ifSum <- P_deathA + P_I_S + P_I_R + P_I_E
B_ifPr <- 1-exp(-(P_deathA + P_I_S + P_I_R + P_I_E))

pf2[1] <-  P_deathA/B_ifSum
pf2[2] <- P_I_S/B_ifSum
pf2[3] <-  P_I_R/B_ifSum
pf2[4] <-  P_I_E/B_ifSum
dim(pf2) <- 4
dim(B_if) <- 4
n_If_Ef <- B_if[4]
n_If_Sf <- if (gammaVer == 1)
  B_if[2] + B_if[3] else 0
n_If_Rf <- if (gammaVer == 2)
  B_if[2] + B_if[3] else 0

#B_rf
B_rf[] <- rmultinom(B_rfSumBIN, pf3)
B_rfSumX<- if (B_rfPr>1 ) 1 else B_rfPr
B_rfSumBIN<-rbinom((Rf),B_rfSumX)
B_rfSum <- P_deathA + P_ImWane
B_rfPr <- 1-exp(-(P_deathA + P_ImWane))
pf3[1] <-P_deathA/B_rfSum
pf3[2] <-  P_ImWane/B_rfSum
dim(pf3) <- 2
dim(B_rf) <- 2
n_Rf_Sf <- B_rf[2]
