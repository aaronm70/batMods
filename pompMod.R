skel <- Csnippet("
double b = b0*(1+b1*cos(2*M_PI*t/365));
DS = S - Beta*S*I/N+ b*N - d*S ;
DI = I + Beta*S*I/N - gamma*I- d*I;
DR = R + gamma*I- d*R;
")





skel <- Csnippet("
N <- (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If + Rn +
        Rj + Rm + Rf + Ma);

b <-(c * exp(-s*(cos(3.141593*timeOsc - phi))^2));

DSn=Sn+(b*(Rf+Ef))- (omega_m+mj*(N/kappa) + (beta*(If+Im+Ij+In))*Sn)+sigma_1*En+gamma_1*In+omega_;
DEn=En - ((omega_m+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*En)+rho*In;
DIn=In -((omega_m+mj*(N/kappa)+gamma_1+gamma_2+rho)*In)+((beta*(If+Im+In+Ij))*Sn)+epsilon*En;
DRn=Rn -(omega_m+mj*(N/kappa)+omega_2)*Rn+sigma_2*En+gamma_2*In;
DMa=Ma+ (b**(Sf+If)) - (omega_m+mj*(N/kappa))*Ma;

DSj=Sj+omega_m*(Sn+Ma) -mu+mj*(N/kappa)+(beta*(If+Im+In+Ij))*Sj + sigma_1*Ej+gamma_1*Ij+omega_2*R;
DEj=Ej+omega_m*En-(mu+mj*(N/kappa)+sigma_1+sigma_2+epsilon)*Ej+rho*Ij;
DIj=Ij+omega_m*In-(mu+mj*(N/kappa)+gamma_1+gamma_2+rho)*Ij+((beta*(If+Im+In+Ij))*Sj)+epsilon*Ej;
DRj=Rj+omega_m*Rn-(mu+mj*(N/kappa)+omega_2)*Rj+sigma_2*Ej+gamma_2*Ij;

DSm=Sm+ mu*(Sj/2)-m*(N/kappa)+(beta*(If+Im+In+Ij))*Sm+sigma_1*Em+gamma_1*Im+omega_2*Rm;
DEm=Em+ mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Em+rho*Im;
DIm=Im+ mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*Im+((beta*(If+Im+In+Ij))*Sm)+epsilon*Em;
DRm=Rm+ mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rm+sigma_2*Em+gamma_2*Im;

DSf=Sf+  mu*(Sj/2)-m*(N/kappa)+(beta*(If+Im+In+Ij))*Sf+sigma_1*Ef+gamma_1*If+omega_2*Rf;
DEf=Ef+  mu*(Ej/2) -(m*(N/kappa)+sigma_1+sigma_2+epsilon)*Ef+rho*If;
DIf=If+  mu*(Ij/2)-(m*(N/kappa)+gamma_1+gamma_2+rho)*If+((beta*(If+Im+In+Ij))*Sf)+epsilon*Ef;
DRf=Rf+  mu*(Rj/2)-(m*(N/kappa)+omega_2)*Rf+sigma_2*Ef+gamma_2*If;

")




N <- (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If + Rn +
        Rj + Rm + Rf + Ma)

##Transition probabilities

step.fun <- Csnippet("
P_death = 1-exp(-(mj * (N / kappa)));
P_deathA = 1-exp(-(m * (N / kappa)));
P_toI = 1-exp(-(beta * (If + Im + Ij + In)));
P_matWane = 1-exp(-omega_m);
P_I_R = 1-exp(-gamma_2);
P_E_I = 1-exp(-epsilon);
P_I_E = 1-exp(-rho);
P_J_A = 1-exp(-mu);
P_ImWane = 1-exp(-omega_2);


double Snrates[3];
double Sn_DN[3];
Snrates[0]=P_death;
Snrates[1]=P_toI;
Snrates[2]=P_matWane;

reulermultinom(1, Sn, &Snrates[0], dt, &Sn_DN[0]);
reulermultinom(1, Sn, &Snrates[1], dt, &Sn_DN[1]);
reulermultinom(1, Sn, &Snrates[2], dt, &Sn_DN[2]);
d_Sn += - Sn_DN[0] - Sn_DN[1] -Sn_DN[2];


double SJrates[3];
double SJ_DN[3];
SJrates[0]=P_death;
SJrates[1]=P_toI;
SJrates[2]=P_J_A;

reulermultinom(1, SJ, &SJrates[0], dt, &SJ_DN[0]);
reulermultinom(1, SJ, &SJrates[1], dt, &SJ_DN[1]);
reulermultinom(1, SJ, &SJrates[2], dt, &SJ_DN[2]);

d_Sj += Sn_DN[2] - SJ_DN[0] - SJ_DN[1] -SJ_DN[2];


")




stochStep <- Csnippet("
double b = b0*(1+b1*cos(2*M_PI*t/365));
double db_rates[2];
double db_DN[2];

p_SI = 1 - exp(-Beta*I/N);
p_IR = 1 - exp(-gamma);
p_b = 1 - exp(-b);
p_d = 1 - exp(-d);

db_rates[0] = p_b;
db_rates[1] = p_d;

reulermultinom(1, N, &db_rates[0], dt, &db_DN[0]);
reulermultinom(1, I, &db_rates[1], dt, &db_DN[1]);


  n_SI = rbinom(S,p_SI);
  n_IR = rbinom(I,p_IR);
  b_N  = db_DN[0];
  d_I  = db_DN[1];
  d_R  = rbinom(R, p_d);
  d_S  = rbinom(S, p_d);
  S=S-n_SI+b_N -d_S;
I = I + n_SI - n_IR- d_I;
R = R + n_IR- d_R;
")


stochStep <- Csnippet("
double b = b0*(1+b1*cos(2*M_PI*t/365));

p_SI = 1 - exp(-Beta*I/N);
p_IR = 1 - exp(-gamma);
p_b = 1 - exp(-b);
p_d = 1 - exp(-d);
  n_SI = rbinom(S,p_SI);
  n_IR = rbinom(I,p_IR);
  b_N  = rbinom(N, p_b);
  d_I  = rbinom(I, p_d);
  d_R  = rbinom(R, p_d);
  d_S  = rbinom(S, p_d);
  S=S-n_SI+b_N -d_S;
I = I + n_SI - n_IR- d_I;
R = R + n_IR- d_R;
")

# initial values of states
init1 <- Csnippet(" S = N-5;
I = 5; R = 0; ")

params1 <- c(Beta=1, b0=0.1, b1=0.4, gamma = 1/13, N= 763, d = 0.1)


closed_sir <- pomp(data=data.frame(time=1:700,data=NA), times="time",
                   t0=0,
                   rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                   skeleton = map(skel),
                   rinit=init1,
                   statenames=c("S","I","R", "p_SI", "p_IR", "p_b", "p_d", "n_SI", "n_IR","b_N", "d_I","d_R","d_S"),
                   paramnames=c("Beta","b0","b1","gamma","N", "d"))


sim <- simulate(closed_sir,params=params1,format = "data.frame", nsim = 100)
traj <- trajectory(closed_sir,params=params1,format="data.frame")

birth_plot <- ggplot(NULL)+ geom_line(data=sim,aes(x=time, y = I+S+R,))+ylim(0,NA)+ labs(title="Seasonal birth pulse")+xlab("Days")+ylab("Size of population")+ theme(plot.title = element_text(hjust = 0.5))
plot(birth_plot)
