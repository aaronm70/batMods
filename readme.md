Readme
================
Aaron M.
28 June, 2020

batMods fits a number of discrete time stochastic models of varying
structures with and without seasonal forces, to observed bat virus data
(currently from boonah australia (Field et al. 2015)), using particle
MCMC based methods. The goal is to identify which dynamical model best
represents the observed viral samples from wild populations and gain
further insight into between-host viral dynamics in bats. Model
comparison is conducted using an approximate leave one out cross
validation algorithm, incorporating Pareto smoothed importance sampling
(Vehtari et al. 2019). This algorithm uses pointwise likelihood values
to compute the log pointwise predictive density and its Monte Carlo
standard error, the effective number of parameters, Pareto k diagnostic
values (which can help assess if a model is well specified) and an
information criterion “looic” (lower values suggest a better model fit).
(Vehtari, Gelman, and Gabry 2017a, 2017b; Vehtari et al. 2015).

<br>

<div class="figure">

<img src="https://github.com/aaronm70/batMods/blob/master/figures/adultMod-Paper.png?raw=true" alt="Figure 1: Model structures for SILI, SIR and SIRS type models, each model is built on top of an age structured bat population model and transitions occur between variable states as probablistic draws from binoial distributions, see suppplementary materials for full details" width="100%" height="75%" />

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

  - The model and fitting methods are described in modelMethods.pdf
    <https://github.com/aaronm70/batMods/blob/master/modelMethods.pdf>

  - This is a work in progress as part of a paper on bat virus dynamics,
    as such should not be seen as a final analysis <br>

## References

<div id="refs" class="references hanging-indent">

<div id="ref-field2015spatiotemporal">

Field, Hume, David Jordan, Daniel Edson, Stephen Morris, Debra Melville,
Kerryn Parry-Jones, Alice Broos, et al. 2015. “Spatiotemporal Aspects of
Hendra Virus Infection in Pteropid Bats (Flying-Foxes) in Eastern
Australia.” *PloS One* 10 (12): e0144055.

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

Vehtari, A., Gelman, A. and Gabry, J. 2017b. “Practical Bayesian Model Evaluation Using Leave-One-Out
Cross-Validation and Waic.” *Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

</div>

<div id="ref-vehtari2015pareto">

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2015. “Pareto Smoothed Importance Sampling.” *arXiv Preprint
arXiv:1507.02646*.

</div>

</div>
