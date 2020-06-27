
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
initial(R0_out)<-0

#parameter input
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
birthType<-user(0)
betaFX<-user(0)
betaFXVal<-user(0)

##sDrive component
envOscType<-user(0)
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

dt<-user(365)
## Total population size, and number of infected bats
N <- (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If + Rn +
        Rj + Rm + Rf + Ma)




#beta calculations
#beta_Val=((R0_Val*((m_Val^2)+(rho_Val*m_Val)+(sigma_2_Val*m_Val)+(m_Val*gamma_1_Val)+
#                     (m_Val*gamma_2_Val)+(m_Val*sigma_2_Val)+(m_Val*epsilon_ValS)+
#                     (rho_Val*sigma_2_Val)+(sigma_2_Val*gamma_1_Val)+(sigma_2_Val*gamma_2_Val)+
#                     (rho_Val*zeta_s)+(gamma_1_Val*zeta_s)+(gamma_2_Val*zeta_s)+
#                     (gamma_1_Val*epsilon_ValS)+(gamma_2_Val*epsilon_ValS))/(N*epsilon_ValS)))



beta_2x<-if(betaFX==0) R0_Val*((epsilon_ValS+sigma_2_Val+m_Val)*(gamma_2_ValS+m_Val+rho_Val)
                               -epsilon_ValS*rho_Val)/(N*(epsilon_ValS+sigma_2_Val+m_Val)) else betaFXVal



# ((epsilon_ValS+m_Val)*(gamma_2_Val+m_Val+rho_Val)-epsilon_ValS*rho_Val)
#beta_Val=R0_Val
#beta_Val2=R0_Val

#parameters
gamma_1 <- if(gammaVer==1 )gamma_1_Val/dt else 0
gamma_2 <- if(gammaVer==2 )gamma_2_ValS/dt else 0
sigma_1 <-0 # if(sigmaVer==1 )zeta_s/dt else 0
sigma_2 <- if(sigmaVer==2 )sigma_2_Val/dt else 0
beta_1<- 0#if(beta_1x>1) 1 else beta_1x
beta_2 <- beta_2x/dt

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


update(timeOsc)<-timeOsc+1/dt

b <-c * exp(-s*(cos(3.141593*timeOsc - phi))^2)

Sbirths<- if(birthType==1) (b*(Sf + If+Ef) )/dt else (b*(Sf + If) )/dt
MaBirths<- if(birthType==1) (b*(Rf))/dt else (b*(Rf+Ef))/dt


update(R0_out)<-R0_Val*((Sn+Sj+Sf+Sm)/N)#((beta_2x*N)*(epsilon_ValS+m_Val))/
  #((epsilon_ValS+m_Val)*(gamma_2_Val+m_Val+rho_Val)-epsilon_ValS*rho_Val)


betSN<-omega_m+mj*(N/kappa) + ((beta_1+beta_2)*(If+Im+Ij+In))
betSj<-mu+mj*(N/kappa)+(beta_1+beta_2)*(If+Im+In+Ij)
betSmf<-m*(N/kappa)+(beta_1+beta_2)*(If+Im+In+Ij)
betE<- beta_1*(If+Im+In+Ij)
betI<-beta_2*(If+Im+In+Ij)


## transitions between neonatal compartments:
update(Sn) <-if(Sn+(Sbirths)- (betSN*Sn)+sigma_1*En+gamma_1*In+omega_2*Rn <=0) 0 else Sn+(Sbirths)- (betSN*Sn)+sigma_1*En+gamma_1*In+omega_2*Rn
update(En) <-if(En - ((omega_m+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*En)+betE*Sn+rho*In <=0) 0 else En - ((omega_m+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*En)+betE*Sn+rho*In
update(In) <- if(In -((omega_m+mj*(N/kappa)+gamma_1+gamma_2+rho)*In)+(betI*Sn)+epsilon*En <=0) 0 else In -((omega_m+mj*(N/kappa)+gamma_1+gamma_2+rho)*In)+(betI*Sn)+epsilon*En
update(Rn) <- if(Rn -(omega_m+mj*(N/kappa)+omega_2)*Rn+sigma_2*En+gamma_2*In <=0) 0 else Rn -(omega_m+mj*(N/kappa)+omega_2)*Rn+sigma_2*En+gamma_2*In
update(Ma) <- if(Ma+ (MaBirths) - (omega_m+mj*(N/kappa))*Ma <=0) 0 else Ma+ (MaBirths) - (omega_m+mj*(N/kappa))*Ma

## transitions between juvenile compartments:
update(Sj) <-if(Sj+omega_m*(Sn+Ma) -betSj*Sj + sigma_1*Ej+gamma_1*Ij+omega_2*Rj <=0) 0 else Sj+omega_m*(Sn+Ma) -betSj*Sj + sigma_1*Ej+gamma_1*Ij+omega_2*Rj
update(Ej) <-if(Ej+omega_m*En-(mu+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*Ej+betE*Sj+rho*Ij <=0) 0 else Ej+omega_m*En-(mu+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*Ej+betE*Sj+rho*Ij
update(Ij) <-if(Ij+omega_m*In-(mu+mj*(N/kappa)+gamma_1+gamma_2+rho)*Ij+betI*Sj+epsilon*Ej <=0) 0 else Ij+omega_m*In-(mu+mj*(N/kappa)+gamma_1+gamma_2+rho)*Ij+betI*Sj+epsilon*Ej
update(Rj) <-if(Rj+omega_m*Rn-(mu+mj*(N/kappa)+omega_2)*Rj+sigma_2*Ej+gamma_2*Ij <=0) 0 else Rj+omega_m*Rn-(mu+mj*(N/kappa)+omega_2)*Rj+sigma_2*Ej+gamma_2*Ij

## transitions between adult male compartments:
update(Sm) <-if(Sm+ mu*(Sj/2)-betSmf*Sm+sigma_1*Em+gamma_1*Im+omega_2*Rm <=0) 0 else Sm+ mu*(Sj/2)-betSmf*Sm+sigma_1*Em+gamma_1*Im+omega_2*Rm
update(Em) <-if(Em+ mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Em+betE*Sm+rho*Im <=0) 0 else Em+ mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Em+betE*Sm+rho*Im
update(Im) <-if(Im+ mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*Im+(betI*Sm)+epsilon*Em <=0) 0 else Im+ mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*Im+(betI*Sm)+epsilon*Em
update(Rm) <-if(Rm+ mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rm+sigma_2*Em+gamma_2*Im <=0) 0 else Rm+ mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rm+sigma_2*Em+gamma_2*Im

## transitions between adult female compartments:
update(Sf) <-if (Sf+  mu*(Sj/2)-betSmf*Sf+sigma_1*Ef+gamma_1*If+omega_2*Rf <=0) 0 else Sf+  mu*(Sj/2)-betSmf*Sf+sigma_1*Ef+gamma_1*If+omega_2*Rf
update(Ef) <-if (Ef+  mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Ef+(betE*Sf)+rho*If <=0) 0 else Ef+  mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Ef+(betE*Sf)+rho*If
update(If) <-if (If+  mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*If+(betI*Sf)+epsilon*Ef <=0) 0 else If+  mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*If+(betI*Sf)+epsilon*Ef
update(Rf) <-if (Rf+  mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rf+sigma_2*Ef+gamma_2*If <=0) 0 else Rf+  mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rf+sigma_2*Ef+gamma_2*If
