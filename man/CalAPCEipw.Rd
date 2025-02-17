% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequentist.R
\name{CalAPCEipw}
\alias{CalAPCEipw}
\title{Compute APCE using frequentist analysis}
\usage{
CalAPCEipw(data)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.}
}
\value{
An object of class \code{list} with the following elements:
\item{P.D1}{An array with dimension (k+1) by (k+2) for quantity P(D(1)=d| R=r), dimension 1 is (k+1) values of D from 0 to k, dimension 2 is (k+2) values of R from 0 to k+1.}
\item{P.D0}{An array with dimension (k+1) by (k+2) for quantity P(D(0)=d| R=r).}
\item{APCE}{An array with dimension (k+1) by (k+2) for quantity P(D(1)=d| R=r)-P(D(0)=d| R=r).}
\item{P.R}{An array with dimension (k+2) for quantity P(R=r) for r from 0 to (k+1).}
\item{alpha}{An array with estimated alpha.}
\item{delta}{An array with estimated delta.}
}
\description{
Estimate propensity score and use Hajek estimator to compute APCE. See S7 for more details.
}
\examples{
data(synth)
freq_apce <- CalAPCEipw(synth)

}
