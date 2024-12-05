#' Gibbs sampler for the main analysis
#'
#' See Appendix S5 for more details.
#'
#' @param data A \code{data.frame} or \code{matrix} of which columns consists of pre-treatment covariates, a binary treatment (Z), an ordinal decision (D), and an outcome variable (Y). The column names of the latter three should be specified as "Z", "D", and "Y" respectively.
#' @param rho A sensitivity parameter. The default is  \code{0} which implies the unconfoundedness assumption (Assumption 4).
#' @param Sigma0.beta.inv Inverse of the prior covariance matrix of beta. The default is a diagonal matrix with  \code{0.01} diagonal entries.
#' @param Sigma0.alpha.inv Inverse of the prior covariance matrix of alpha. The default is a diagonal matrix with  \code{0.01} diagonal entries.
#' @param sigma0 Prior variance of the cutoff points (theta and delta)
#' @param beta Initial value for beta.
#' @param alpha Initial value for alpha.
#' @param theta Initial value for theta.
#' @param delta Initial value for delta.
#' @param n.mcmc The total number of MCMC iterations. The default is \code{50}.
#' @param verbose A logical argument specified to print the progress on the screen. The default is \code{FALSE}.
#' @param out.length An integer to specify the progress on the screen. If \code{verbose = TRUE}, every \code{out.length}-th iteration is printed on the screen. The default is \code{10}.
#' @param beta.zx.off A logical argument specified to exclude the interaction terms (Z by X) from the model. The default is \code{FALSE}.
#' @param theta.z.off A logical argument specified to set same cutoffs theta for treatment and control group. The default is \code{FALSE}.
#'
#' @return An object of class \code{mcmc} containing the posterior samples.
#'
#' @importFrom coda mcmc
#' @importFrom stats model.matrix
#'
#' @examples
#' data(synth)
#' sample_mcmc <- AiEvalmcmc(data = synth, n.mcmc = 2)
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
AiEvalmcmc <- function(data,
                       rho = 0,
                       Sigma0.beta.inv = NULL, Sigma0.alpha.inv = NULL, sigma0 = NULL,
                       beta = NULL, alpha = NULL, theta = NULL, delta = NULL,
                       n.mcmc = 5 * 10, verbose = FALSE, out.length = 10,
                       beta.zx.off = FALSE, theta.z.off = FALSE) {
  if (sum(is.na(data)) > 0) {
    stop("NA in the data")
  }

  Z <- data$Z
  D <- data$D
  Y <- data$Y

  if (length(unique(Z)) != 2) {
    stop("Non-binary treatment")
  }
  if (length(unique(Y)) != 2) {
    stop("Non-binary outcome")
  }

  X <- as.matrix(subset(data, select = -c(Z, D, Y)))
  zX <- as.matrix(subset(data, select = -c(Z, D, Y)))
  lZX <- ncol(zX)

  D.factor <- model.matrix(~ 0 + as.factor(D))

  if (beta.zx.off) {
    ZX <- cbind(Z, X)
  } else {
    ZX <- cbind(Z, X, Z * zX)
  }

  p <- dim(X)[2]
  k <- length(unique(D)) - 1

  C <- matrix(0, (k + 1), (k + 1))
  C[lower.tri(C)] <- 1
  C <- t(C) %*% C # for delta
  out <- .AiEvalProbitOrdinal(Z, D, Y, X, D.factor, ZX, rho, Sigma0.beta.inv, Sigma0.alpha.inv, sigma0, beta, alpha, theta, delta, n.mcmc, verbose, out.length, beta.zx.off, theta.z.off, lZX, C)

  #### combine the posterior samples in an mcmc object
  MCMC <- mcmc(out)


  ### name each variable
  colnames(MCMC) <- rep("name", dim(MCMC)[2])
  colnames(MCMC)[1] <- "betaZ"
  if (beta.zx.off & (!theta.z.off)) {
    colnames(MCMC)[2 * p + 2 * k + 2 + k] <- paste("delta", k + 1, sep = "")
    for (l in 1:p) {
      colnames(MCMC)[1 + l] <- paste("beta", l, sep = "")
      colnames(MCMC)[p + k + 1 + l + k] <- paste("alpha", l, sep = "")
    }

    for (r in 1:k) {
      colnames(MCMC)[p + 1 + r] <- paste("theta1", r, sep = "")
      colnames(MCMC)[p + 1 + r + k] <- paste("theta0", r, sep = "")
      colnames(MCMC)[2 * p + 1 + r + 2 * k] <- paste("delta", r, sep = "")
    }
  } else if ((!beta.zx.off) & (!theta.z.off)) {
    colnames(MCMC)[2 * p + 2 * k + 2 + k + lZX] <- paste("delta", k + 1, sep = "")
    for (l in 1:p) {
      colnames(MCMC)[1 + l] <- paste("beta", l, sep = "")
      if (l < (1 + lZX)) {
        colnames(MCMC)[1 + l + p] <- paste("Zbeta", l, sep = "")
      }
      colnames(MCMC)[p + k + 1 + l + k + lZX] <- paste("alpha", l, sep = "")
    }

    for (r in 1:k) {
      colnames(MCMC)[p + 1 + r + lZX] <- paste("theta1", r, sep = "")
      colnames(MCMC)[p + 1 + r + k + lZX] <- paste("theta0", r, sep = "")
      colnames(MCMC)[2 * p + 1 + r + 2 * k + lZX] <- paste("delta", r, sep = "")
    }
  } else if (beta.zx.off & theta.z.off) {
    colnames(MCMC)[2 * p + 2 * k + 2] <- paste("delta", k + 1, sep = "")
    for (l in 1:p) {
      colnames(MCMC)[1 + l] <- paste("beta", l, sep = "")
      colnames(MCMC)[p + k + 1 + l] <- paste("alpha", l, sep = "")
    }

    for (r in 1:k) {
      colnames(MCMC)[p + 1 + r] <- paste("theta", r, sep = "")
      colnames(MCMC)[2 * p + 1 + r + k] <- paste("delta", r, sep = "")
    }
  } else {
    colnames(MCMC)[2 * p + 2 * k + 2 + lZX] <- paste("delta", k + 1, sep = "")
    for (l in 1:p) {
      colnames(MCMC)[1 + l] <- paste("beta", l, sep = "")
      if (l < (1 + lZX)) {
        colnames(MCMC)[1 + l + p] <- paste("Zbeta", l, sep = "")
      }
      colnames(MCMC)[p + k + 1 + l + lZX] <- paste("alpha", l, sep = "")
    }

    for (r in 1:k) {
      colnames(MCMC)[p + 1 + r + lZX] <- paste("theta", r, sep = "")
      colnames(MCMC)[2 * p + 1 + r + k + lZX] <- paste("delta", r, sep = "")
    }
  }
  return(MCMC)
}
