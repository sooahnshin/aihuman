#' Crossfitting for nuisance functions
#'
#' Implement crossfitting with boosting methods and get predicted values for outcome/decision regression or propensity score models
#'
#' @param data A \code{data.frame} or \code{matrix}  to fit on.
#' @param include_for_fit Boolean vector for whether or not a unit should be included in fitting (e.g. treated/control).
#' @param form Formula for outcome regression/propensity score models.
#' @param ... Additional arguments to be passed to \code{gbm} function.
#'
#' @aliases gbm
#'
#' @return A vector of predicted values
#'
#' @importFrom gbm gbm
#' @importFrom gbm gbm.perf
#' @importFrom stats predict
#'
#' @useDynLib aihuman, .registration=TRUE
#'
crossfit <- function(data, include_for_fit, form, ...) {
  n <- nrow(data)
  # split into 3 folds
  folds <- split(sample(n), cut(1:n, 3))

  pred <- rep(NA, n)
  for (fold in folds) {
    # training set
    idx <- setdiff(1:n, fold)
    idx <- idx[include_for_fit[idx]]
    train_data <- data[idx, ]
    # fit model using boosting method
    gb <- gbm(form, data = train_data, ...)
    best_iter <- gbm.perf(gb, plot.it = FALSE)
    # compute predictions
    pred[fold] <- predict(gb, data[fold, ], n.trees = best_iter, type = "response")
  }
  return(pred)
}

#' Fit outcome/decision and propensity score models
#'
#' Fit (1) the decision model \eqn{m^{D}(z, X_i) := \Pr(D = 1 \mid Z = z, X = X_i)} and
#' (2) the outcome model \eqn{m^{Y}(z, X_i) := \Pr(Y = 1 \mid D = 0, Z = z, X = X_i)}
#' for each treatment group \eqn{z \in \{0,1\}} and (3) the propensity score model
#' \eqn{e(1, X_i) := \Pr(Z = 1 \mid X = X_i)}.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param V Pretreatment covariates for nuisance functions. A vector, a matrix, or a data frame.
#' @param d_form A formula for decision model where the dependent variable is \code{D}.
#' @param y_form A formula for outcome model where the dependent variable is \code{Y}.
#' @param ps_form A formula for propensity score model.
#' @param distribution A distribution argument used in \code{gbm} function. Default is \code{"bernoulli"}.
#' @param n.trees Integer specifying the total number of trees to fit used in \code{gbm} function.
#' @param shrinkage A shrinkage parameter used in \code{gbm} function.
#' @param interaction.depth Integer specifying the maximum depth of each tree used in \code{gbm} function.
#' @param ... Additional arguments to be passed to \code{gbm} function called in \code{crossfit}
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{z_models}}{A \code{data.frame} with the following columns:
#'     \describe{
#'       \item{\code{idx}}{Index of observation.}
#'       \item{\code{d_pred}}{Predicted probability of decision.}
#'       \item{\code{y_pred}}{Predicted probability of outcome.}
#'       \item{\code{Z}}{Treatment group.}
#'     }
#'   }
#'   \item{\code{pscore}}{A vector of predicted propensity scores.}
#' }
#' 
#' @importFrom magrittr %>%
#' @importFrom tidyselect starts_with
#' @importFrom dplyr select bind_rows mutate if_else
#'
#' @examples
#' compute_nuisance_functions(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   V = NCAdata[, c("Sex", "White", "Age")],
#'   shrinkage = 0.01,
#'   n.trees = 1000
#' )
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
compute_nuisance_functions <- function(Y, D, Z, V,
                                       d_form = D ~ .,
                                       y_form = Y ~ .,
                                       ps_form = Z ~ .,
                                       distribution = "bernoulli",
                                       n.trees = 1000,
                                       shrinkage = 0.01,
                                       interaction.depth = 1, ...) {
  if (!is.vector(Y) || !is.vector(D) || !is.vector(Z)) {
    stop("Y, D, and Z must be vectors")
  }
  if (is.vector(V)) {
    if (!isTRUE(all.equal(length(V), length(Y), length(D), length(Z)))) {
      stop("V must have the same length as Y, D, and Z")
    }
  } else if (is.matrix(V) || is.data.frame(V)) {
    if (!isTRUE(all.equal(nrow(V), length(Y), length(D), length(Z)))) {
      stop("V must have the same number of rows as Y, D, and Z")
    }
  } else {
    stop("V must be either a vector, a matrix, or a data frame")
  }
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary")
  }
  if (!all(D %in% c(0, 1))) {
    stop("D must be binary")
  }
  if (!all(Z %in% c(0, 1))) {
    stop("Z must be binary")
  }
  # fit different outcome/decision models for each Z
  dat <- data.frame(Y = Y, D = D, Z = Z, V = V)
  zs <- unique(dat$Z)
  models <- lapply(zs, function(z) {
    d_preds <- crossfit(
      dat %>%
        select(starts_with("V"), D),
      (Z == z),
      d_form,
      distribution = distribution,
      n.trees = n.trees, shrinkage = shrinkage,
      interaction.depth = interaction.depth, ...
    )
    y_preds <- crossfit(
      dat %>%
        select(starts_with("V"), Y),
      (Z == z & D == 0),
      y_form,
      distribution = distribution,
      n.trees = n.trees, shrinkage = shrinkage,
      interaction.depth = interaction.depth, ...
    )

    return(data.frame(
      idx = 1:length(Y),
      d_pred = d_preds,
      y_pred = y_preds,
      Z = z
    ))
  }) %>%
    bind_rows()
  pscores <- crossfit(dat %>% select(starts_with("V"), Z),
    rep(TRUE, length(Z)),
    ps_form,
    distribution = distribution,
    n.trees = n.trees, shrinkage = shrinkage,
    interaction.depth = interaction.depth, ...
  )

  res <- list(z_models = models, pscore = pscores)
  class(res) <- "nuisance_functions"

  return(res)
}
#' Fit outcome/decision and propensity score models conditioning on the AI recommendation
#'
#' Fit (1) the decision model \eqn{m^{D}(z, a, X_i) := \Pr(D = 1 \mid Z = z, A = a, X = X_i)} and
#' (2) the outcome model \eqn{m^{Y}(z, a, X_i) := \Pr(Y = 1 \mid D = 0, Z = z, A = a, X = X_i)}
#' for each treatment group \eqn{z \in \{0,1\}} and AI recommendation \eqn{a \in \{0,1\}},
#' and (3) the propensity score model \eqn{e(1, X_i) := \Pr(Z = 1 \mid X = X_i)}.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param V A \code{matrix} of pretreatment covariates for nuisance functions.
#' @param d_form A formula for decision model where the dependent variable is \code{D}.
#' @param y_form A formula for outcome model where the dependent variable is \code{Y}.
#' @param ps_form A formula for propensity score model.
#' @param distribution A distribution argument used in \code{gbm} function. Default is \code{"bernoulli"}.
#' @param n.trees Integer specifying the total number of trees to fit used in \code{gbm} function.
#' @param shrinkage A shrinkage parameter used in \code{gbm} function.
#' @param interaction.depth Integer specifying the maximum depth of each tree used in \code{gbm} function.
#' @param ... Additional arguments to be passed to \code{gbm} function called in \code{crossfit}
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{z_models}}{A \code{data.frame} with the following columns:
#'     \describe{
#'       \item{\code{idx}}{Index of observation.}
#'       \item{\code{d_pred}}{Predicted probability of decision.}
#'       \item{\code{y_pred}}{Predicted probability of outcome.}
#'       \item{\code{Z}}{Treatment group.}
#'       \item{\code{A}}{AI recommendation.}
#'     }
#'   }
#'   \item{\code{pscore}}{A vector of predicted propensity scores.}
#' }
#' 
#' @importFrom magrittr %>%
#' @importFrom tidyselect starts_with
#' @importFrom dplyr select bind_rows mutate if_else
#'
#' @examples
#' compute_nuisance_functions_ai(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   V = NCAdata[, c("Sex", "White", "Age")],
#'   shrinkage = 0.01,
#'   n.trees = 1000
#' )
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
compute_nuisance_functions_ai <- function(Y, D, Z, A, V,
                                          d_form = D ~ .,
                                          y_form = Y ~ .,
                                          ps_form = Z ~ .,
                                          distribution = "bernoulli",
                                          n.trees = 1000, shrinkage = 0.01, interaction.depth = 1, ...) {
  if (!is.vector(Y) || !is.vector(D) || !is.vector(Z) || !is.vector(A)) {
    stop("Y, D, Z, and A must be vectors")
  }
  if (is.vector(V)) {
    if (!isTRUE(all.equal(length(V), length(Y), length(D), length(Z), length(A)))) {
      stop("V must have the same length as Y, D, Z, and A")
    }
  } else if (is.matrix(V) || is.data.frame(V)) {
    if (!isTRUE(all.equal(nrow(V), length(Y), length(D), length(Z), length(A)))) {
      stop("V must have the same number of rows as Y, D, Z, and A")
    }
  } else {
    stop("V must be either a vector, a matrix, or a data frame")
  }
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary")
  }
  if (!all(D %in% c(0, 1))) {
    stop("D must be binary")
  }
  if (!all(Z %in% c(0, 1))) {
    stop("Z must be binary")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary")
  }
  # fit different outcome/decision models for each Z
  dat <- data.frame(Y = Y, D = D, Z = Z, A = A, V = V)
  zs <- unique(dat$Z)
  as <- unique(dat$A)
  params <- expand.grid(zs, as)
  models <- purrr::map(1:nrow(params), function(i) {
    z <- params$Var1[i]
    a <- params$Var2[i]
    d_preds <- crossfit(
      dat %>%
        select(starts_with("V"), D),
      (Z == z & A == a),
      d_form,
      distribution = distribution,
      n.trees = n.trees, shrinkage = shrinkage,
      interaction.depth = interaction.depth, ...
    )
    y_preds <- crossfit(
      dat %>%
        select(starts_with("V"), Y),
      (Z == z & A == a & D == 0),
      y_form,
      distribution = distribution,
      n.trees = n.trees, shrinkage = shrinkage,
      interaction.depth = interaction.depth, ...
    )

    return(data.frame(
      idx = 1:length(Y),
      d_pred = d_preds,
      y_pred = y_preds,
      Z = z,
      A = a
    ))
  }) %>%
    bind_rows()
  pscores <- crossfit(dat %>% select(starts_with("V"), Z),
    rep(TRUE, length(Z)),
    ps_form,
    distribution = distribution,
    n.trees = n.trees, shrinkage = shrinkage,
    interaction.depth = interaction.depth, ...
  )
  
  res <- list(z_models = models, pscore = pscores)
  class(res) <- "nuisance_functions"
  return(res)
}
