Readme
================
Aaron M.
27 June, 2020

# Overview

batMods fits a number of discrete time stochastic models of varying
structures with and without seasonal forces, to observed bat virus data
(currently from boonah australia (Field et al. 2015)), using particle
MCMC based methods. The goal is to identify which dynamical model best
represents the observed viral samples from wild populations and gain
further insight into between-host viral dynamics in bats. Model
comparison is conducted used an approximate leave one out cross
validation algorithm, incorporating Pareto smoothed importance sampling
(Vehtari et al. 2019). This algorithm uses pointwise likelihood values
to compute the log pointwise predictive density and its Monte Carlo
standard error, the effective number of parameters, Pareto k diagnostic
values (which can help assess if a model is well specified) and an
information criterion “looic” (lower values suggest a better model fit).
(Vehtari, Gelman, and Gabry 2017a, 2017b; Vehtari et al. 2015).

<br>

<div class="figure">

<img src="/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/adultMod-Paper.png" alt="Figure 1: Model structures for SILI, SIR and SIRS type models, each model is built on top of an age structured bat population model and transitions occur between variable states as probablistic draws from binoial distributions, see suppplementary materials for full details" width="100%" height="75%" />

<p class="caption">

Figure 1: Model structures for SILI, SIR and SIRS type models, each
model is built on top of an age structured bat population model and
transitions occur between variable states as probablistic draws from
binoial distributions, see suppplementary materials for full details

</p>

</div>

<br>

  - Currently batMods fits three primary models structures with and
    without maternal immmunity and seasonal forces (figure 1) to
    multiple data-types, including serology and PCR data.

  - The analysis can be run from the runscript.R file.

  - The metropilis hastings and particle filter algorithms run in R,
    whilst the model itself runs in C code, which is implemented via the
    Odin package.

  - The model and fitting methods are described below.

  - This is a work in progress as part of a paper on bat virus dynamics,
    as such should not be seen as a final analysis <br>

# Model description

<br> **Discrete time model equations** <br>

The age structured viral dynamics model, with all potential transitions
is stated in the equations below, with parameters corresponding to table
1 (currently in main text of paper), transitions between states occurs
in discrete time \((t\to t+1)\), in each model discrete time setps
\((t)\) are set to 0.25 days. 

\(S_{N}(t+1)=S_{N}(t) + b(t)(S_F(t)+I_F(t))-(\omega_m+m_j\frac{N}{\kappa}+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t)))S_N(t)+\omega R_N(t)\)

<br>
\(S_{J}(t+1)=S_{J}(t)+\omega_m(S_N(t)+Ma)-(\mu+m_j\frac{N}{\kappa}+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t)))S_J(t)+\omega R_J(t)\)

<br>
\(S_{M}(t+1)=S_{M}+(t)\mu \frac{S_{J}}{2}-(m \frac{N}{\kappa}+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t))) S_{M}(t) I_M(t)+\omega R_M(t)\)

<br>
\(S_{F}(t+1)=S_{F}(t)+\mu \frac{S_{J}}{2}-(m \frac{N}{\kappa}+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t))) S_{F}(t)+\omega R_F(t)\)

<br>
\(L_{N}(t+1)= L_{N}(t)-( \omega_m+m_j\frac{N}{\kappa}+\epsilon)L_N+\rho I_N\)

<br>
\(L_{J}(t+1)= L_{J}(t) + \omega_mL_N(t) - (\mu +m_j\frac{N}{\kappa}+\epsilon)L_J(t)+\rho L_J(t)\)

<br>
\(L_{M}(t+1)= L_{M}(t) - (\mu +m_j\frac{N}{\kappa}+\epsilon)L_M(t)+\rho L_M(t)\)

<br>
\(L_{F}(t+1)= L_{F}(t) - (\mu +m_j\frac{N}{\kappa}+\epsilon)L_F(t)+\rho L_F(t)\)

<br>
\(I_{N}(t+1)= I_{J}(t) - (\omega_m+m_j\frac{N}{\kappa}+\gamma+\rho)I_J+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t)))S_N +\epsilon L_N\)

<br>
\(I_{J}(t+1)= I_{J}(t) - (\mu+m_j\frac{N}{\kappa}+\gamma+\rho)I_J+\beta(I_J(t)+I_N(t)+I_{F}(t)+I_{M}(t)))S_J +\epsilon L_J\)

<br>
\(I_{F}(t+1)=I_{F}(t)-(m \frac{N}{\kappa}\gamma+\rho) I_{F}(t)+\beta (I_{F}(t)+I_{M}(t))S_{F}(t)+\epsilon L_{F}(t)\)

<br>
\(I_{M}(t+1)=I_{M}(t)-(m \frac{N}{\kappa}+\gamma+\rho) I_{M}(t)+\beta (I_{M}(t)+I_{F}(t))S_{M}(t)+\epsilon L_{M}(t)\)

<br> \(R_N(t+1) = R_N(t)-(\omega_m +m_j\frac{N}{\kappa}+\omega)R_N(t)\)

<br>
\(R_J(t+1) = R_J(t)+\omega_m R_N(t) -(\mu+m_j\frac{N}{\kappa}+\omega)R_J(t)\)

<br>
\(R_M(t+1) = R_M(t)+\mu \frac{R_J(t)}{2} -(m\frac{N}{\kappa}+\omega)R_M(t)+\gamma I_M(t)\)

<br>
\(R_F(t+1) = R_F(t)+\mu \frac{R_J(t)}{2} -(m\frac{N}{\kappa}+\omega)R_F(t)+\gamma I_F(t)\)

\(M_a(t+1)=M_a(t) +b(R_F(t)+L_F(t))-(\omega_m+m_j\frac{N}{\kappa})M_a(t)\)  

Births occur in seasonal pulses, as previously described (Peel et al.
2014), with timing, amplitude and seasonality derived from three
parameters in a modified Gaussian function.

\(\beta\) is derived from a rearrangement of equation (2) in the main
text, here \(N\) is considered an upper limit to the population, thus
\(N\) is equal to the carrying capacity \(\kappa\):  

**Stochastic transitions between states**

Movement between states at each discrete time-step (set as 0.25 days)
for individuals is a stochastic process, using a combination of
binomial, and multinomial distributions where multiple transition
possibilities occur.

e.g. for bats exiting the variable state \(R\), transitions will occur
with the following steps:

  - The total number of individuals transitioning between states is
    determined by first summing all of the rates of exiting a current
    state, e.g. \((R \to exit)=\omega+\frac{N}{\kappa}\).

  - We then convert the sum of rates to a probability for each possible
    transition occurring during a single time-step, e.g. for a rate
    \(x\), the equivalent probability of transition in one time-step is
    \(p = 1-e^{-x}\), thus \(Pr(R \to exit) = 1-e^{-(R \to exit)}\)

  - For each individual, a random draw from a binomial distribution is
    then performed at each time-step to obtain the number of individuals
    exiting \(R\) e.g. \(R_e = bin(R,Pr(R \to exit))\).

  - The relative probabilities are then computed for transition to each
    of the potential states <br>  <br>and <br> 

  - A multinomial distribution is used with each probability to
    determine the destination of each individual
    \(M (R_e,Pr(R\to S),Pr(R\to death))\) 

# Model fitting

The models are considered state-space (hidden Markov) models, with
observations \(k\) and hidden state variables \(Z\) at time \(t\) and we
seek to identify an unbiased estimate of the likelihood value, (also
known as a marginal likelihood) \(pr(k|\theta)\), by essentially
“marginalising out” (\(Z\)). Here we do this using a particle markov
chain monte carlo (pMCMC) algorithm. Which combines a particle filter as
an unbiased estimator of the marginal likelihood value, and an MCMC
acceptence/rejection algorithm (Andrieu, Doucet, and Holenstein 2010).

# pMCMC generic example

Particle filters act as an unbiased estimator of marginal likelihood by
applying a form of importance re-sampling, to generate an approximate
sample from, and make inferences about an unobserved Markov process
(Smith 2013). They are becoming popular in disciplines where state-space
models are common such as ecology and epidemiology (Kantas et al. 2015;
Peters, Hosack, and Hayes 2010; Knape and De Valpine 2012; Fasiolo et
al. 2016; Sheinson, Niemi, and Meiring 2014).

In brief, the marginal likelihood \(pr(k|\theta)\) for the observed data
\(k\), the hidden state \(Z\) and the parameter vector \(\theta\) is
considered the joint posterior probability
\(pr(k|Z,\theta)*p(Z|\theta)\). Using Monte-Carlo approximation for
\(p(k|\theta)\), a particle filter with \(j\) particles which have the
possible trajectories/evolution of \(Z_j\), the marginal likelihood can
be considered as
\(pr(k | \theta) \approx \sum_{J} pr\left(k | Z_{j}, \theta\right) * pr\left(Z_{j} | \theta\right)\).

To run the particle filter algorithm requires the following steps:

1.  Initialise the particles with equal weights.  

2.  For each particle j at time t, simulate the initial conditions at
    the first observed data point.  

3.  Calculate a probability weighting for each particle based on the
    results of the simulation and the observed data.  

4.  The marginal likelihood for this data point can be considered the
    average of each particle weighting.  

5.  Normalise the weightings \(\frac{w_{jt}}{(\sum_Jw_{jt} )}\),
    resample with replacement each particle based on their weighting and
    simulate forward to the next observed data point.  

Repeat steps 3 to 5 for all observed data, an estimate for the total
marginal likelihood value can be considered the product of the average
particle weights at each observed timepoint.  

This marginal likelihood value is then used in a traditional MCMC
algorithm, in this instance Metropolis-Hastings, where rejection or
acceptance of a parameter proposal is based on the current and previous
marginal likelihood values.

**Model fitting with pMCMC to Boonah bat data**

To fit the model with pMCMC to our observed data, we developed a
likelihood function which incorporated the three unique longitudinal
data-types; Individual serology for Hendra virus antibodies, individual
PCR of urine for Hendra virus RNA, and under roost PCR of urine for
Hendra virus RNA.

  - *Firstly, we consider the observed data:*
    
    If the sampling time points are \(i_1...i_n\), the number of
    individual samples through time are \(N_{i...n}\) of which
    \(k^p_{i...n}\) are PCR positive and \(k^s_{i...n}\) are
    seropositive, the number of under-roost urine samples is
    \(N^u_{i...n}\) of which \(k^u_{i...n}\) are PCR positive.

  - *Secondly, the hidden state:*
    
    If the system state of the model at each time point is
    \(Z_{i...n}\), the transition probability density function for \(Z\)
    from the stochastic model, conditioned on the parameters \(\theta\),
    is \(pr(Z_i|Z_{i-1},\theta)\).

Below we define how each data type contributes to the likelihood
function in detail: 

***Individual serological and PCR of urine samples***

In the observed data, individual bats are found to be in all four
empirical states of serological and PCR positivity:

  - PCR positive and seropositive (\(k^{P^+S^+}\))
  - PCR negative and seropositive (\(k^{P^-S^+}\))
  - PCR positive and seronegative (\(k^{P^+S^-}\))
  - PCR negative and seronegative (\(k^{P^-S^-}\))

Logically, a relationship between the PCR and serology states is likely
and so we consider a joint conditional probability function rather than
treating the data as independent. Transitory states and test detection
failure are also considered by fitting the coefficients \(\zeta_s\) and
\(\zeta_p\) which are bounded between 0 and 1. Thus, the simulated
implementation of each observed empirical state are as follows:

  - \(k^{P^+S^+}\): Bats are assumed to be in the I state of the model
    and are equivalent to the number of I state bats proportionate to
    the value of the coefficients; \(z^{P^+S^+}=I\zeta_p\zeta_s\)

  - \(k^{P^-S^+}\): Bats are assumed to be in either the L or R states
    of the models, or in the I state of the model;
    \(z^{P^-S^+}=(I+L+R)\zeta_s(1-\zeta_p)\)

  - \(k^{P^+S^-}\): Bats are considered to be in the I state of the
    model; \(z^{P^+S^-}=I(1-\zeta_s)\zeta_p\)

  - \(k^{P^-S^-}\): Bats are assumed to be in the S state of the models
    or in the I state; \(z^{P^-S^-}=I(1-\zeta_s)(1-\zeta_p)+S\)

If \(k\) is a vector of the observed counts of bats in each empirical
state, \(k=(k^{P^+S^+},k^{P^-S^+},k^{P^+S^-},k^{P^-S^-})\), \(z\) is a
vector of the simulated prevalences of each state
\(z=(z^{P^+S^+},z^{P^-S^+},z^{P^+S^-},z^{P^-S^-})\) and \(N\) is the
total observed bats, we can then propose a joint probability function
based on a multinomial distribution, \(pr(k^p,k^s|Z,N)\).

***Pooled under roost samples***

If the probability of at least one bat within a pooled urine sample of
\(x\) individuals being PCR positive for virus RNA can be considered
\(P_t=1-(1-p)^x\), where \(p\) is the probability of each individual bat
contributing to a pool being positive (Chiang and Reeves 1962). The
result of PCR testing on a pooled sample (\(T_i\)) will be one of two
states, with the following probabilities: 

<!-- Nu(t) samples, of whichXu(t)were positive -->

<!--   #prob of Nu positives, given Xu samples and state of model Z -->

<!--   #Pt is the probability of a pooled sample being positive -->

<!--   #p is the probability of a bat being positive, x is the number of bats contributing -->

<!--   #p is sampled from Z, as the prevalence of infectious (I/N) -->

<!--   #x is a poison dist variable -->

Therefore, assuming that any one under-roost sample comes from \(x\)
bats with a prevalence of \(p\), considering the model state \(Z(t)\),
we aim to identify the likelihood of the contribution of positive
under-roost samples at time (\(t\)) given (from the predicted GAM
values) the number of samples \(N^u(t)\), of which \(k^u(t)\) were
positive, \(pr(k^u|N^u,Z)\). As the predicted values from the GAM are
for prevalence, the total number of bats (\(N^u(t)\)) was considered the
mean number of samples collected on each sampling occasion.

Due to the nature and challenges of field-sampling, under-roost data is
likely to contain a degree of overdispersion, as it is highly exposed to
external factors which are difficult to account for in study design.
Therefore, the probability is unlikely to follow a standard binomial
distribution, where the variance is defined by the mean. To account for
this, we used a probability function based around a beta-binomial
distribution, which allows for an additional variance parameter to be
fitted to account for any overdispersion in the observations. A
beta-binomial distribution, is a binomial distribution such
\(pr(k^u|N^u,Z)\) would be:

However, \(p_t\) is not constant and is generated from a beta
distribution, which takes two shape parameters \(\alpha\) and \(\beta\),
\(beta(\alpha,\beta)\) thus;

Here, \(O_u\) accounts for overdispersion in the data, and the
prevalence of infectious bats \(p\) in \(p_t\) is derived from the model
state \(Z(t)\) as \(\frac{I(t)\zeta_u}{N(t)}\), where \(\zeta_u\) is a
fitted coefficient accounting for detection failure.

As there is no precise measure of exact bat numbers contributing to a
pooled sample, \(x\) is a fitted parameter with a prior based on expert
knowledge and estimates from the field. For simplicity, only one value
of \(x\) is derived for all pools on a single sampling occasion,
assuming that each pooled sample on average has the same number of bats
contributing. We also assume that there is no effect of pooling on the
diagnostic tests, and any that does occur should be primarily accounted
for by the \(\zeta_u\) coefficient.

***full likelihood function***

Considering the above, and the sampling time points (\(i\)), we can
define the full likelihood function as:

<br> Table 1: Informative priors are based on normal distributions and
are shown with 95% credible intervals in brackets, uninformative priors
use a uniform distribution and are shown as a minimum and maximum value

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Parameter

</th>

<th style="text-align:left;">

Prior

</th>

<th style="text-align:left;">

Description

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

\(S\)

</td>

<td style="text-align:left;">

  - 
    
    </td>
    
    <td style="text-align:left;">
    
    Number of susceptible bats
    
    </td>
    
    </tr>
    
    <tr>
    
    <td style="text-align:left;">
    
    \(I\)
    
    </td>
    
    <td style="text-align:left;">
    
      - 
        
        </td>
        
        <td style="text-align:left;">
        
        Number of infectious bats
        
        </td>
        
        </tr>
        
        <tr>
        
        <td style="text-align:left;">
        
        \(R\)
        
        </td>
        
        <td style="text-align:left;">
        
          - 
            
            </td>
            
            <td style="text-align:left;">
            
            Number of recovered bats
            
            </td>
            
            </tr>
            
            <tr>
            
            <td style="text-align:left;">
            
            \(L\)
            
            </td>
            
            <td style="text-align:left;">
            
              - 
                
                </td>
                
                <td style="text-align:left;">
                
                Number of latently infected bats
                
                </td>
                
                </tr>
                
                <tr>
                
                <td style="text-align:left;">
                
                \(N\)
                
                </td>
                
                <td style="text-align:left;">
                
                  - 
                    
                    </td>
                    
                    <td style="text-align:left;">
                    
                    Total population of bats
                    
                    </td>
                    
                    </tr>
                    
                    <tr>
                    
                    <td style="text-align:left;">
                    
                    \(\theta\)
                    
                    </td>
                    
                    <td style="text-align:left;">
                    
                      - 
                        
                        </td>
                        
                        <td style="text-align:left;">
                        
                        Parameter vector
                        
                        </td>
                        
                        </tr>
                        
                        <tr>
                        
                        <td style="text-align:left;">
                        
                        \(Z\)
                        
                        </td>
                        
                        <td style="text-align:left;">
                        
                          - 
                            
                            </td>
                            
                            <td style="text-align:left;">
                            
                            System state
                            
                            </td>
                            
                            </tr>
                            
                            <tr>
                            
                            <td style="text-align:left;">
                            
                            \(z\)
                            
                            </td>
                            
                            <td style="text-align:left;">
                            
                              - 
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                Vector of simulated individual bat
                                states
                                
                                </td>
                                
                                </tr>
                                
                                <tr>
                                
                                <td style="text-align:left;">
                                
                                \(\zeta_s\)
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                0-1
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                Individual serology coefficient
                                
                                </td>
                                
                                </tr>
                                
                                <tr>
                                
                                <td style="text-align:left;">
                                
                                \(\zeta_p\)
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                0-1
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                Individual PCR coefficient
                                
                                </td>
                                
                                </tr>
                                
                                <tr>
                                
                                <td style="text-align:left;">
                                
                                \(\zeta_u\)
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                0-1
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                Under roost PCR coefficient
                                
                                </td>
                                
                                </tr>
                                
                                <tr>
                                
                                <td style="text-align:left;">
                                
                                \(k^s\)
                                
                                </td>
                                
                                <td style="text-align:left;">
                                
                                  - 
                                    
                                    </td>
                                    
                                    <td style="text-align:left;">
                                    
                                    Observed individual seropositive
                                    bats
                                    
                                    </td>
                                    
                                    </tr>
                                    
                                    <tr>
                                    
                                    <td style="text-align:left;">
                                    
                                    \(k^p\)
                                    
                                    </td>
                                    
                                    <td style="text-align:left;">
                                    
                                      - 
                                        
                                        </td>
                                        
                                        <td style="text-align:left;">
                                        
                                        Observed individual PCR positive
                                        bats
                                        
                                        </td>
                                        
                                        </tr>
                                        
                                        <tr>
                                        
                                        <td style="text-align:left;">
                                        
                                        \(k^u\)
                                        
                                        </td>
                                        
                                        <td style="text-align:left;">
                                        
                                          - 
                                            
                                            </td>
                                            
                                            <td style="text-align:left;">
                                            
                                            Predicted no. positive bats
                                            contributing to pooled under
                                            roost samples
                                            
                                            </td>
                                            
                                            </tr>
                                            
                                            <tr>
                                            
                                            <td style="text-align:left;">
                                            
                                            \(k\)
                                            
                                            </td>
                                            
                                            <td style="text-align:left;">
                                            
                                              - 
                                                
                                                </td>
                                                
                                                <td style="text-align:left;">
                                                
                                                Vector of observed
                                                individual bat states
                                                
                                                </td>
                                                
                                                </tr>
                                                
                                                <tr>
                                                
                                                <td style="text-align:left;">
                                                
                                                \(N^u\)
                                                
                                                </td>
                                                
                                                <td style="text-align:left;">
                                                
                                                  - 
                                                    
                                                    </td>
                                                    
                                                    <td style="text-align:left;">
                                                    
                                                    Number of pooled
                                                    under-roost samples
                                                    (mean)
                                                    
                                                    </td>
                                                    
                                                    </tr>
                                                    
                                                    <tr>
                                                    
                                                    <td style="text-align:left;">
                                                    
                                                    \(p\)
                                                    
                                                    </td>
                                                    
                                                    <td style="text-align:left;">
                                                    
                                                      - 
                                                        
                                                        </td>
                                                        
                                                        <td style="text-align:left;">
                                                        
                                                        Probability of
                                                        bat being
                                                        positive
                                                        (infection
                                                        prevalence)
                                                        
                                                        </td>
                                                        
                                                        </tr>
                                                        
                                                        <tr>
                                                        
                                                        <td style="text-align:left;">
                                                        
                                                        \(P_t\)
                                                        
                                                        </td>
                                                        
                                                        <td style="text-align:left;">
                                                        
                                                          - 
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            Probability
                                                            of pooled
                                                            sample being
                                                            positive
                                                            
                                                            </td>
                                                            
                                                            </tr>
                                                            
                                                            <tr>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            \(O_u\)
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            \>0
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            Overdispersion
                                                            parameter
                                                            for under
                                                            roost data
                                                            
                                                            </td>
                                                            
                                                            </tr>
                                                            
                                                            <tr>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            \(x\)
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            0-10
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            Number of
                                                            bats
                                                            contributing
                                                            to each
                                                            pooled under
                                                            roost sample
                                                            
                                                            </td>
                                                            
                                                            </tr>
                                                            
                                                            <tr>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                            \(T_i\)
                                                            
                                                            </td>
                                                            
                                                            <td style="text-align:left;">
                                                            
                                                              - 
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Pooled
                                                                under
                                                                roost
                                                                sample
                                                                test
                                                                state
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(\beta\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \> 0
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Infection
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(\gamma\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \> 1 day
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Recovery
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(\rho\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \> 1 day
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Latency
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(\epsilon\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \> 1 day
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Recurrence
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(m\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                0.186
                                                                (0.146-0.225)
                                                                year\(^{-1}\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Adult
                                                                mortality
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(m_j\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                0.500
                                                                (0.480-0.520)
                                                                year\(^{-1}\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                Juvenile
                                                                mortality
                                                                rate
                                                                
                                                                </td>
                                                                
                                                                </tr>
                                                                
                                                                <tr>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                \(b\)
                                                                
                                                                </td>
                                                                
                                                                <td style="text-align:left;">
                                                                
                                                                  - 
                                                                    
                                                                    </td>
                                                                    
                                                                    <td style="text-align:left;">
                                                                    
                                                                    Birth
                                                                    rate
                                                                    
                                                                    </td>
                                                                    
                                                                    </tr>
                                                                    
                                                                    <tr>
                                                                    
                                                                    <td style="text-align:left;">
                                                                    
                                                                    \(c\)
                                                                    
                                                                    </td>
                                                                    
                                                                    <td style="text-align:left;">
                                                                    
                                                                      - 
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Birth
                                                                        pulse
                                                                        scalar
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(s\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        130
                                                                        (111-150)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Birth
                                                                        pulse
                                                                        synchronicity
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(\phi\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        7.180
                                                                        (6.787-7.571)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Birth
                                                                        pulse
                                                                        timing
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(\mu\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        1.37
                                                                        year\(^{-1}\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Maturation
                                                                        rate
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(\omega_m\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        0.800
                                                                        (0.741-0.859)
                                                                        year\(^{-1}\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Maternal
                                                                        antibody
                                                                        waning
                                                                        rate
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(\omega\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \>
                                                                        1
                                                                        day
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Antibody
                                                                        waning
                                                                        rate
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(R_0\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \>
                                                                        1
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(R_{0}\)
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(d\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        0-25
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Pooled
                                                                        sample
                                                                        contributing
                                                                        bats
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(c_ u\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        1
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        Env.
                                                                        force
                                                                        scalar
                                                                        
                                                                        </td>
                                                                        
                                                                        </tr>
                                                                        
                                                                        <tr>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                        \(s_ u\)
                                                                        
                                                                        </td>
                                                                        
                                                                        <td style="text-align:left;">
                                                                        
                                                                          - 
                                                                            
                                                                            </td>
                                                                            
                                                                            <td style="text-align:left;">
                                                                            
                                                                            Env.
                                                                            force
                                                                            synchronicity
                                                                            
                                                                            </td>
                                                                            
                                                                            </tr>
                                                                            
                                                                            <tr>
                                                                            
                                                                            <td style="text-align:left;">
                                                                            
                                                                            \(\phi_ u\)
                                                                            
                                                                            </td>
                                                                            
                                                                            <td style="text-align:left;">
                                                                            
                                                                              - 
                                                                                
                                                                                </td>
                                                                                
                                                                                <td style="text-align:left;">
                                                                                
                                                                                Env.
                                                                                force
                                                                                timing
                                                                                
                                                                                </td>
                                                                                
                                                                                </tr>
                                                                                
                                                                                <tr>
                                                                                
                                                                                <td style="text-align:left;">
                                                                                
                                                                                \(\kappa\)
                                                                                
                                                                                </td>
                                                                                
                                                                                <td style="text-align:left;">
                                                                                
                                                                                4000
                                                                                (2570-6300)
                                                                                
                                                                                </td>
                                                                                
                                                                                <td style="text-align:left;">
                                                                                
                                                                                Environmental
                                                                                carrying
                                                                                capacity
                                                                                
                                                                                </td>
                                                                                
                                                                                </tr>
                                                                                
                                                                                </tbody>
                                                                                
                                                                                </table>

<br>

**Initial conditions**

Starting parameters are estimated by running a Metropolis-Hastings MCMC
algorithm with a deterministic version of the model for 50,000
iterations with a 10,000 iteration burn-in. The median parameter set
from each posterior distribution was then used in the pMCMC.

Initial conditions for the bat population are derived from the carrying
capacity parameter (\(\kappa\)). A population equal to the carrying
capacity is first assumed and this is then split into demographic
compartments as follows:

Infectious and exposed individuals are added as 5% of the population and
the model is run for 50 years to reach an equilibrium before fitting to
the first data point.

# Example output

<br> 

<div class="figure">

<img src="/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/2_plot.png" alt="Figure 1: 100 repeat runs of the maternal immunity SILI model without seasonal forcing (model 2 in table 1), parameterised using the median values from the posterior estimates of each parameter. Dark blue points show the observed data, purple points show the model outputs for each run of the model at the corresponding observed time period. The top figure shows the simulated and observed serological prevalence, the middle figure shows the simulated and observed virus RNA prevalence from individual urine samples, and the bottom figure shows the simulated and predicted (using observed data and contributing bats parameter $d$) virus RNA prevalence from the under-roost samples. Confidence intervals on observed data are calculated as binomial confidence intervals, with a beta prior on the binomial distribution; the shape of the beta prior for individual samples is uninformitive and for under-roost data is derived from the corresponding fitted overdispersion parameters for each data-type." width="100%" height="100%" />

<p class="caption">

Figure 1: 100 repeat runs of the maternal immunity SILI model without
seasonal forcing (model 2 in table 1), parameterised using the median
values from the posterior estimates of each parameter. Dark blue points
show the observed data, purple points show the model outputs for each
run of the model at the corresponding observed time period. The top
figure shows the simulated and observed serological prevalence, the
middle figure shows the simulated and observed virus RNA prevalence from
individual urine samples, and the bottom figure shows the simulated and
predicted (using observed data and contributing bats parameter \(d\))
virus RNA prevalence from the under-roost samples. Confidence intervals
on observed data are calculated as binomial confidence intervals, with a
beta prior on the binomial distribution; the shape of the beta prior for
individual samples is uninformitive and for under-roost data is derived
from the corresponding fitted overdispersion parameters for each
data-type.

</p>

</div>

<br><br>

<div class="figure">

<img src="/Users/alm204/OneDrive/Cambridge/Projects/model_comparisons/figures/4_plot.png" alt="Figure 2: 100 repeat runs of the maternal immunity SILI model with seasonal forcing (model 4 in table 1), parameterised using the median values from the posterior estimates of each parameter. Dark blue points show the observed data, purple points show the model outputs for each run of the model at the corresponding observed time period. The top figure shows the simulated and observed serological prevalence, the middle figure shows the simulated and observed virus RNA prevalence from individual urine samples, and the bottom figure shows the simulated and predicted (using observed data and contributing bats parameter $d$) virus RNA prevalence from the under-roost samples. Confidence intervals on observed data are calculated as binomial confidence intervals, with a beta prior on the binomial distribution; the shape of the beta prior for individual samples is uninformitive and for under-roost data is derived from the corresponding fitted overdispersion parameters for each data-type." width="100%" height="100%" />

<p class="caption">

Figure 2: 100 repeat runs of the maternal immunity SILI model with
seasonal forcing (model 4 in table 1), parameterised using the median
values from the posterior estimates of each parameter. Dark blue points
show the observed data, purple points show the model outputs for each
run of the model at the corresponding observed time period. The top
figure shows the simulated and observed serological prevalence, the
middle figure shows the simulated and observed virus RNA prevalence from
individual urine samples, and the bottom figure shows the simulated and
predicted (using observed data and contributing bats parameter \(d\))
virus RNA prevalence from the under-roost samples. Confidence intervals
on observed data are calculated as binomial confidence intervals, with a
beta prior on the binomial distribution; the shape of the beta prior for
individual samples is uninformitive and for under-roost data is derived
from the corresponding fitted overdispersion parameters for each
data-type.

</p>

</div>

<br> 

# References

<div id="refs" class="references hanging-indent">

<div id="ref-andrieu2010particle">

Andrieu, Christophe, Arnaud Doucet, and Roman Holenstein. 2010.
“Particle Markov Chain Monte Carlo Methods.” *Journal of the Royal
Statistical Society: Series B (Statistical Methodology)* 72 (3):
269–342.

</div>

<div id="ref-chiang1962statistical">

Chiang, Chin-long, and William C et al. Reeves. 1962. “Statistical
Estimation of Virus Infection Rates in Mosquito Vector Populations.”
*American Journal of Hygiene* 75 (3).

</div>

<div id="ref-fasiolo2016comparison">

Fasiolo, Matteo, Natalya Pya, Simon N Wood, and others. 2016. “A
Comparison of Inferential Methods for Highly Nonlinear State Space
Models in Ecology and Epidemiology.” *Statistical Science* 31 (1):
96–118.

</div>

<div id="ref-field2015spatiotemporal">

Field, Hume, David Jordan, Daniel Edson, Stephen Morris, Debra Melville,
Kerryn Parry-Jones, Alice Broos, et al. 2015. “Spatiotemporal Aspects of
Hendra Virus Infection in Pteropid Bats (Flying-Foxes) in Eastern
Australia.” *PloS One* 10 (12): e0144055.

</div>

<div id="ref-kantas2015particle">

Kantas, Nikolas, Arnaud Doucet, Sumeetpal S Singh, Jan Maciejowski,
Nicolas Chopin, and others. 2015. “On Particle Methods for Parameter
Estimation in State-Space Models.” *Statistical Science* 30 (3): 328–51.

</div>

<div id="ref-knape2012fitting">

Knape, Jonas, and Perry De Valpine. 2012. “Fitting Complex Population
Models by Combining Particle Filters with Markov Chain Monte Carlo.”
*Ecology* 93 (2): 256–63.

</div>

<div id="ref-peel2014effect">

Peel, Alison J, JRC Pulliam, AD Luis, RK Plowright, TJ O’Shea, DTS
Hayman, JLN Wood, CT Webb, and O Restif. 2014. “The Effect of Seasonal
Birth Pulses on Pathogen Persistence in Wild Mammal Populations.”
*Proceedings of the Royal Society B: Biological Sciences* 281 (1786):
20132962.

</div>

<div id="ref-peters2010ecological">

Peters, Gareth W, Geoff R Hosack, and Keith R Hayes. 2010. “Ecological
Non-Linear State Space Model Selection via Adaptive Particle Markov
Chain Monte Carlo (Adpmcmc).” *arXiv Preprint arXiv:1005.2238*.

</div>

<div id="ref-sheinson2014comparison">

Sheinson, Daniel M, Jarad Niemi, and Wendy Meiring. 2014. “Comparison of
the Performance of Particle Filter Algorithms Applied to Tracking of a
Disease Epidemic.” *Mathematical Biosciences* 255: 21–32.

</div>

<div id="ref-smith2013sequential">

Smith, Adrian. 2013. *Sequential Monte Carlo Methods in Practice*.
Springer Science & Business Media.

</div>

<div id="ref-VehtariLooPackage">

Vehtari, Aki, Jonah Gabry, Mans Magnusson, Yuling Yao, and Andrew
Gelman. 2019. “Loo: Efficient Leave-One-Out Cross-Validation and Waic
for Bayesian Models.” <https://mc-stan.org/loo>.

</div>

<div id="ref-vehtari2017practical">

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017a. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and Waic.”
*Statistics and Computing* 27 (5): 1413–32.

</div>

<div id="ref-AkiLoo">

———. 2017b. “Practical Bayesian Model Evaluation Using Leave-One-Out
Cross-Validation and Waic.” *Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

</div>

<div id="ref-vehtari2015pareto">

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2015. “Pareto Smoothed Importance Sampling.” *arXiv Preprint
arXiv:1507.02646*.

</div>

</div>
