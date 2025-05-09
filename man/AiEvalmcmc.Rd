% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{AiEvalmcmc}
\alias{AiEvalmcmc}
\title{Gibbs sampler for the main analysis}
\usage{
AiEvalmcmc(
  data,
  rho = 0,
  Sigma0.beta.inv = NULL,
  Sigma0.alpha.inv = NULL,
  sigma0 = NULL,
  beta = NULL,
  alpha = NULL,
  theta = NULL,
  delta = NULL,
  n.mcmc = 5 * 10,
  verbose = FALSE,
  out.length = 10,
  beta.zx.off = FALSE,
  theta.z.off = FALSE
)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.}

\item{rho}{A sensitivity parameter. The default is  \code{0} which implies the unconfoundedness assumption (Assumption 4).}

\item{Sigma0.beta.inv}{Inverse of the prior covariance matrix of beta. The default is a diagonal matrix with  \code{0.01} diagonal entries.}

\item{Sigma0.alpha.inv}{Inverse of the prior covariance matrix of alpha. The default is a diagonal matrix with  \code{0.01} diagonal entries.}

\item{sigma0}{Prior variance of the cutoff points (theta and delta)}

\item{beta}{Initial value for beta.}

\item{alpha}{Initial value for alpha.}

\item{theta}{Initial value for theta.}

\item{delta}{Initial value for delta.}

\item{n.mcmc}{The total number of MCMC iterations. The default is \code{50}.}

\item{verbose}{A logical argument specified to print the progress on the screen. The default is \code{FALSE}.}

\item{out.length}{An integer to specify the progress on the screen. If \code{verbose = TRUE}, every \code{out.length}-th iteration is printed on the screen. The default is \code{10}.}

\item{beta.zx.off}{A logical argument specified to exclude the interaction terms (Z by X) from the model. The default is \code{FALSE}.}

\item{theta.z.off}{A logical argument specified to set same cutoffs theta for treatment and control group. The default is \code{FALSE}.}
}
\value{
An object of class \code{mcmc} containing the posterior samples.
}
\description{
See Appendix S5 for more details.
}
\examples{
data(synth)
sample_mcmc <- AiEvalmcmc(data = synth, n.mcmc = 2)

}
