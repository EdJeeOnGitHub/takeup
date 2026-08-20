# Bayesian Data-to-Parameter Sensitivity

To clarify which sources of variation inform each structural parameter, we
report a Bayesian data-to-parameter sensitivity analysis. This analysis is
analogous in spirit to the moment-sensitivity calculations commonly reported
for minimum-distance estimators, but it is adapted to our Bayesian estimator,
which is fit directly to microdata rather than to a predetermined vector of
empirical moments.

Let \(\theta\) denote the vector of structural parameters. We partition the
likelihood into four components:

\[
L(\theta)
=
L_T(\theta)
L_B(\theta)
L_W(\theta)
L_D(\theta),
\]

where \(L_T\) is the likelihood contribution of the experimental take-up data,
\(L_B\) is the contribution of the beliefs data, \(L_W\) is the contribution
of the willingness-to-pay data, and \(L_D\) is the contribution of the
observed distance distribution. The posterior is therefore

\[
p(\theta\mid y)
\propto
p(\theta)
L_T(\theta)
L_B(\theta)
L_W(\theta)
L_D(\theta).
\]

To measure the influence of each data source, introduce a likelihood weight
\(\alpha_b\) for block \(b\):

\[
p_{\boldsymbol{\alpha}}(\theta\mid y)
\propto
p(\theta)
\prod_{b\in\{T,B,W,D\}}
L_b(\theta)^{\alpha_b}.
\]

The estimated model corresponds to \(\alpha_b=1\) for every block. Increasing
\(\alpha_b\) places marginally more weight on data source \(b\), while
decreasing it places less weight on that source.

For structural parameter \(\theta_k\), the local sensitivity of its posterior
mean to the weight on data block \(b\) is

\[
\left.
\frac{\partial
E_{\boldsymbol{\alpha}}[\theta_k\mid y]}
{\partial\alpha_b}
\right|_{\boldsymbol{\alpha}=\boldsymbol{1}}
=
\operatorname{Cov}_{post}
\left(
\theta_k,\log L_b(\theta)
\right).
\]

Thus, the sensitivity can be calculated directly from the posterior draws and
the likelihood contribution of each data block. Intuitively, a positive value
means that posterior draws assigning a relatively high likelihood to data
block \(b\) also tend to have a relatively high value of parameter
\(\theta_k\). Marginally increasing the weight on that data block would
therefore raise the posterior mean of the parameter. A negative value has the
converse interpretation.

Because the structural parameters have different units, we report
sensitivities normalized by the posterior standard deviation of each
parameter:

\[
S_{kb}
=
\frac{
\operatorname{Cov}_{post}
\left(
\theta_k,\log L_b(\theta)
\right)
}{
\operatorname{sd}_{post}(\theta_k)
}.
\]

An entry \(S_{kb}\) gives the approximate change in the posterior mean of
parameter \(k\), measured in posterior-standard-deviation units, resulting
from a one-unit marginal increase in the likelihood weight on data source
\(b\).

The sensitivity table reports these quantities for every structural
parameter. Its columns correspond to the take-up, beliefs,
willingness-to-pay, and distance likelihoods. This produces a quantitative
map from the different sources of evidence to the parameters reported in the
structural-parameter table. Importantly, the calculation allows a data source
to inform parameters outside the submodel in which it appears directly. For
example, the beliefs data directly determine predicted visibility, but the
visibility parameters also enter reputational utility and hence predicted
take-up. The joint posterior therefore allows both beliefs and take-up data to
inform these parameters.

The sensitivity matrix is a local diagnostic. It describes small changes in
likelihood weights around the estimated posterior and need not accurately
approximate complete removal of a data source when the posterior is highly
nonlinear. Moreover, a large sensitivity does not by itself imply weak
identification: it indicates that changing the relative weight on a data
source moves the posterior mean. We therefore interpret the matrix as a
transparent description of how the joint model combines its different
sources of evidence, rather than as a test of global identification.

### Posterior-uncertainty sensitivity

For dispersion parameters such as \(\sigma_u\), and for nonlinear derived
quantities, a likelihood block may affect posterior precision without moving
the posterior mean. We therefore also report the local response of the
posterior standard deviation:

\[
\left.
\frac{\partial\log \operatorname{sd}_{\boldsymbol{\alpha}}(\theta_k\mid y)}
{\partial\alpha_b}
\right|_{\boldsymbol{\alpha}=\boldsymbol{1}}
=
\frac{
\operatorname{Cov}_{post}
\left(
[\theta_k-E_{post}(\theta_k)]^2,\log L_b(\theta)
\right)
}{
2\operatorname{Var}_{post}(\theta_k)
}.
\]

A negative entry means that marginally increasing the weight on data block
\(b\) contracts the posterior standard deviation of the row quantity. A
positive entry means that it expands locally. As with the posterior-mean
sensitivity, this is a derivative at \(\alpha_b=1\), not the result of
completely adding or removing a data source.

### Derived structural quantities

The same covariance identities apply to any scalar function of a posterior
draw. We calculate both posterior-mean and posterior-uncertainty sensitivity
for two paper-facing structural objects:

- the social multiplier, using the paper's sign and scale convention,
  \[
  M_z(d)
  =
  -\frac{
  E_c[\texttt{cluster\_social\_multiplier}_{cz}(d)]
  }{\delta};
  \]
- the net social-image return,
  \[
  SI_z(d)
  =
  E_c\!\left[
  \lambda p_{\mathrm{Observed},cz}(d)\Delta(w^*_{cz}(d))
  \right],
  \]
  stored by Stan as the cluster-average of `cluster_rep_return`.

We evaluate these quantities for Control, Ink, Calendar, and Bracelet at
0.5, 1.5, and 2.5 km, spanning the main experimental distance range.

## Proposed table

**Title:** Bayesian Data-to-Parameter Sensitivity

Use the same parameter rows and group headings as the existing structural
parameter table. The four columns are:

1. Take-up
2. Beliefs
3. WTP
4. Distance

### Table note

*Notes:* Each cell reports the local sensitivity of the posterior mean of the
row parameter to the likelihood weight on the column data source, normalized
by the posterior standard deviation of the parameter. For parameter
\(\theta_k\) and likelihood block \(b\), the reported quantity is
\(\operatorname{Cov}_{post}(\theta_k,\log L_b(\theta))/
\operatorname{sd}_{post}(\theta_k)\). Positive entries indicate that
increasing the relative weight on the data source would raise the posterior
mean of the parameter; negative entries indicate that it would lower it. The
four likelihood blocks contain the experimental take-up data, beliefs data,
willingness-to-pay data, and observed distance distribution, respectively.
The calculation is local to the estimated posterior.

## Implementation

The calculation is implemented in two stages. The first stage evaluates the
four likelihood-block totals at every saved posterior draw using standalone
Stan generated quantities. It does **not** rerun MCMC:

```sh
Rscript scripts/appendix/generate_bayesian_sensitivity_quantities.R 104 \
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
```

The second stage joins the likelihood totals to the corresponding posterior
draws, calculates the posterior covariances, divides each parameter row by its
posterior standard deviation, and writes CSV, LaTeX, and Markdown tables:

```sh
Rscript scripts/appendix/bayesian_sensitivity.R 104 \
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
  --derived
```

By default, generated quantities and tables are written to
`temp-data/bayesian-sensitivity/`. The generated table files begin with
`bayesian-parameter-sensitivity-fit104-` for reported parameters and
`bayesian-derived-sensitivity-fit104-` for the multiplier and social-image
returns. Files ending in `-uncertainty` report the derivative of the log
posterior standard deviation; files without that suffix report standardized
posterior-mean sensitivity.

The relevant implementation files are:

- `stan_models/bayesian_sensitivity_generated_quantities.stan`, which
  reproduces the take-up, beliefs, WTP, and distance likelihood contributions;
- `stan_models/takeup_struct_sensitivity.stan`, the lightweight standalone
  generated-quantities model;
- `scripts/appendix/generate_bayesian_sensitivity_quantities.R`, which reconstructs the
  original Stan data and evaluates the likelihood blocks at the saved draws;
- `scripts/appendix/bayesian_sensitivity.R`, which calculates and formats the sensitivity
  matrix.

The scripts process the original fit CSVs one chain at a time because these
files are large. Draws remain in their original within-chain order, ensuring
that each parameter draw is paired with the likelihood totals evaluated at
that draw.

For the main fit, \(\alpha_b=1\). A marginal change
\(\mathrm{d}\alpha_b\) changes the posterior mean by approximately
\(S_{kb}\mathrm{d}\alpha_b\) posterior standard deviations. Thus, a full
one-unit increase in \(\alpha_b\) would double the log-likelihood contribution
of block \(b\). The reported derivative is local at \(\alpha_b=1\); it should
not be interpreted as an exact prediction for such a large change.
