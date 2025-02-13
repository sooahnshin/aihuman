% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_itt.R
\name{CalDIMsubgroup}
\alias{CalDIMsubgroup}
\title{Calculate diff-in-means estimates}
\usage{
CalDIMsubgroup(
  data,
  subgroup,
  name.group = c("Overall", "Female", "Male", "Non-white\\nMale", "White\\nMale")
)
}
\arguments{
\item{data}{A \code{data.frame} of which columns includes a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y).}

\item{subgroup}{A list of numeric vectors for the index of each of the five subgroups.}

\item{name.group}{A character vector including the labels of five subgroups.}
}
\value{
A \code{data.frame} of diff-in-means estimates for each value of D and Y for each subgroup.
}
\description{
Calculate average causal effect based on diff-in-means estimator.
}
\examples{
data(synth)
subgroup_synth <- list(
  1:nrow(synth), which(synth$Sex == 0), which(synth$Sex == 1),
  which(synth$Sex == 1 & synth$White == 0), which(synth$Sex == 1 & synth$White == 1)
)
CalDIMsubgroup(synth, subgroup = subgroup_synth)

}
