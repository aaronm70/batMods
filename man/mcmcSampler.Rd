% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metHast.R
\name{mcmcSampler}
\alias{mcmcSampler}
\title{pMCMCsampler}
\usage{
mcmcSampler(
  initParams,
  randInit = F,
  proposer = sequential.proposer,
  sdProps,
  maxSddProps,
  niter = 100,
  particleNum = 100,
  nburn = 100,
  monitoring = 2,
  adaptiveMCMC = T,
  proposerType = "seq",
  startAdapt = 1000,
  adptBurn = 200,
  acceptanceRate = 0.25,
  tell = 5,
  cluster = F,
  oDat,
  stoch = T,
  priorFunc = NA,
  switch = 2500,
  switchBlock = 2500,
  likelihoodFunc = likelihoodFunc1,
  juvenileInfection = F
)
}
\arguments{
\item{randInit}{T then randomly sample initial parameters instead of initParams value}

\item{proposer}{proposal function, multivariate block (adaptiveMCMC must = T) or sequential can have adaptive or not adaptive tuning}

\item{sdProps}{standard deviation for proposal distributions - in adaptive this is the starting sd}

\item{maxSddProps}{maximum values for sd proposals if using adaptive MCMC}

\item{niter}{number of iterations to run the MCMC for}

\item{nburn}{number of mcmc iterations to burn}

\item{monitoring}{0 = no monitoring, > 0 prints more progress information}

\item{adaptiveMCMC}{T/F if true uses tuning for s.d. of parmater proposal distributions based on acceptance ratios}

\item{proposerType}{"seq" or "block", seq for sequential proposing, block for blocked proposing, blocks based on var-covar matrix - block proposing only available when using adaptiveMCMC}

\item{startAdapt}{starting iteration for adapting}

\item{acceptanceRate}{acceptance rates, for adaptive sequential use string, one rate for each param, for block use single value}

\item{tell}{print monitoring information every x number of iterations}

\item{cluster}{T/F if true, using dide cluster, if false running locally}

\item{initial}{parameter guess}

\item{particles}{number of particles for particle filter}

\item{adaptBurn}{burn n number of iterations for defining var-covar matrix in block proposing}
}
\value{
mcmc
}
\description{
pMCMCsampler
}
