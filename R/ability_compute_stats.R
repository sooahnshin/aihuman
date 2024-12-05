#' Compute Risk (Human+AI v. Human)
#'
#' Compute the difference in risk between human+AI and human decision makers using AIPW estimators.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}
#' @param true.pscore A vector of true propensity scores (numeric), if available. Optional.
#' @param X Pretreatment covariate used for subgroup analysis (vector). Must be the same length as Y, D, Z, and A if provided. Default is NULL.
#' @param l01 Ratio of the loss between false positives and false negatives
#'
#' @return A tibble the following columns:
#' \itemize{
#'  \item \code{Z_focal}: The focal treatment indicator. `1` indicates the treatment group.
#'  \item \code{Z_compare}: The comparison treatment indicator. `0` indicates the control group.
#'  \item \code{X}: Pretreatment covariate (if provided).
#'  \item \code{loss_diff}: The difference in loss between human+AI and human decision
#'  \item \code{loss_diff_se}: The standard error of the difference in loss
#'  \item \code{fn_diff}: The difference in false negatives between human+AI and human decision
#'  \item \code{fn_diff_se}: The standard error of the difference in false negatives
#'  \item \code{fp_diff}: The difference in false positives between human+AI and human decision
#'  \item \code{fp_diff_se}: The standard error of the difference in false positives
#'  }
#'
#'
#' @examples
#' compute_stats_aipw(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   nuis_funcs = nuis_func,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   X = NULL,
#'   l01 = 1
#' )
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols group_by summarise select mutate n
#' @importFrom tidyr pivot_wider
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
compute_stats_aipw <- function(Y, D, Z, nuis_funcs, true.pscore = NULL, X = NULL, l01 = 1) {
  if (!is.vector(Y) || !is.vector(D) || !is.vector(Z)) {
    stop("Y, D, and Z must be vectors")
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
  flag <- is.null(X)

  # Check lengths using all.equal
  if (flag) {
    if (!isTRUE(all.equal(length(Y), length(D), length(Z)))) {
      stop("Y, D, and Z must have the same length")
    }
  } else if (is.vector(X)) {
    if (!isTRUE(all.equal(length(X), length(Y), length(D), length(Z)))) {
      stop("X must have the same length as Y, D, and Z")
    }
  } else {
    stop("X must be either NULL or a vector")
  }

  if (flag) {
    X <- rep(1, length(Y))
  }
  
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
  }
  
  dat <- data.frame(Y = Y, D = D, Z = Z, X = X, pscore = nuis_funcs$pscore)

  preds <- nuis_funcs$z_models |>
    pivot_wider(names_from = c(Z), values_from = c(-Z, -idx), names_sep = "") |>
    select(-idx)

  dat <- dat %>%
    bind_cols(preds)

  props <- dat %>%
    mutate(
      # FN Pr(Y(0) = 1, D(1) = 0) = Pr(Y=1,D=0|Z=1)
      p10.1i = (1 - d_pred1) * y_pred1 + Z * (1 - D) * (Y - y_pred1) / pscore - y_pred1 * Z * (D - d_pred1) / pscore,
      # FN Pr(Y(0) = 1, D(0) = 0) = Pr(Y=1,D=0|Z=0)
      p10.0i = (1 - d_pred0) * y_pred0 + (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - y_pred0 * (1 - Z) * (D - d_pred0) / (1 - pscore),
      # Pr(Y=0,D=0|Z=1)
      p00.1i = (1 - d_pred1) * (1 - y_pred1) - Z * (1 - D) * (Y - y_pred1) / pscore - (1 - y_pred1) * Z * (D - d_pred1) / pscore,
      # Pr(Y=0,D=0|Z=0)
      p00.0i = (1 - d_pred0) * (1 - y_pred0) - (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - (1 - y_pred0) * (1 - Z) * (D - d_pred0) / (1 - pscore),
      # Pr(W|Z=1) where W = Y(1-D) - ell_01 (1-Y)(1-D)
      p.1i = (1 - d_pred1) * ((1 + l01) * y_pred1 - l01) + (1 + l01) * Z * (1 - D) * (Y - y_pred1) / pscore - ((1 + l01) * y_pred1 - l01) * Z * (D - d_pred1) / pscore,
      # Pr(W|Z=0)
      p.0i = (1 - d_pred0) * ((1 + l01) * y_pred0 - l01) + (1 + l01) * (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - ((1 + l01) * y_pred0 - l01) * (1 - Z) * (D - d_pred0) / (1 - pscore)
    ) %>%
    group_by(X) %>%
    summarise(
      n = n(),
      p1 = mean(p.1i),
      p0 = mean(p.0i),
      p1_p0_se = sqrt(sum(((p.1i - p.0i) - (p1 - p0))^2) / n^2),
      p10.1 = mean(p10.1i),
      p10.0 = mean(p10.0i),
      p10.1_p10.0_se = sqrt(sum(((p10.1i - p10.0i) - (p10.1 - p10.0))^2) / n^2),
      p00.1 = mean(p00.1i),
      p00.0 = mean(p00.0i),
      p00.0_p00.1_se = sqrt(sum(((p00.0i - p00.1i) - (p00.0 - p00.1))^2) / n^2),
      .groups = "drop"
    )
  out <- props %>%
    mutate(
      loss_diff = p1 - p0,
      loss_diff_se = p1_p0_se,
      fn_diff = p10.1 - p10.0,
      fn_diff_se = p10.1_p10.0_se,
      fp_diff = p00.0 - p00.1, # Pr(Y(0)=0,D(1)=1) - Pr(Y(0)=0,D(0)=1) = Pr(Y(0)=0,D(0)=0) - Pr(Y(0)=0,D(1)=0) = Pr(Y=0,D=0|Z=0) - Pr(Y=0,D=0|Z=1)
      fp_diff_se = p00.0_p00.1_se
    ) %>%
    mutate(Z_focal = 1, Z_compare = 0) %>%
    select(Z_focal, Z_compare, X, contains("diff"))

  if (flag) {
    return(out %>% select(-X))
  } else {
    return(out)
  }
}
#' Compute Risk (Human+AI v. Human)
#'
#' Compute the difference in risk between human+AI and human decision makers using difference-in-means estimators.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param X Pretreatment covariate used for subgroup analysis (vector). Must be the same length as Y, D, Z, and A if provided. Default is NULL.
#' @param l01 Ratio of the loss between false positives and false negatives
#'
#' @return A tibble the following columns:
#' \itemize{
#'  \item \code{Z_focal}: The focal treatment indicator. `1` indicates the treatment group.
#'  \item \code{Z_compare}: The comparison treatment indicator. `0` indicates the control group.
#'  \item \code{X}: Pretreatment covariate (if provided).
#'  \item \code{loss_diff}: The difference in loss between human+AI and human decision
#'  \item \code{loss_diff_se}: The standard error of the difference in loss
#'  \item \code{fn_diff}: The difference in false negatives between human+AI and human decision
#'  \item \code{fn_diff_se}: The standard error of the difference in false negatives
#'  \item \code{fp_diff}: The difference in false positives between human+AI and human decision
#'  \item \code{fp_diff_se}: The standard error of the difference in false positives
#'  }
#'
#' @examples
#' compute_stats(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   X = NULL,
#'   l01 = 1
#' )
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
compute_stats <- function(Y, D, Z, X = NULL, l01 = 1) {
  if (!is.vector(Y) || !is.vector(D) || !is.vector(Z)) {
    stop("Y, D, and Z must be vectors")
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
  flag <- is.null(X)
  # Check lengths using all.equal
  if (flag) {
    if (!isTRUE(all.equal(length(Y), length(D), length(Z)))) {
      stop("Y, D, and Z must have the same length")
    }
  } else if (is.vector(X)) {
    if (!isTRUE(all.equal(length(X), length(Y), length(D), length(Z)))) {
      stop("X must have the same length as Y, D, and Z")
    }
  } else {
    stop("X must be either NULL or a vector")
  }

  if (flag) {
    X <- rep(1, length(Y))
  }
  props <- data.frame(Y = Y, D = D, Z = Z, X = X, l01 = l01) %>%
    mutate(
      W0 = ifelse(D == 0 & Y == 0, 1, 0),
      W1 = ifelse(D == 0 & Y == 1, 1, 0),
      W = W1 - l01 * W0
    ) %>%
    group_by(Z, X) %>%
    summarise(
      Wbar = mean(W),
      n = n(),
      sigma_W = sd(W),
      p00 = mean(W0),
      p10 = mean(W1),
      p00_se = sd(W0),
      p10_se = sd(W1),
      .groups = "drop"
    )

  out <- cross_join(props, props) %>%
    filter(X.y == X.x) %>%
    mutate(
      loss_diff = Wbar.x - Wbar.y,
      loss_diff_se = sqrt(sigma_W.x^2 / n.x + sigma_W.y^2 / n.y),
      fn_diff = p10.x - p10.y,
      fn_diff_se = sqrt(p10_se.x^2 / n.x + p10_se.y^2 / n.y),
      fp_diff = -(p00.x - p00.y),
      fp_diff_se = sqrt(p00_se.x^2 / n.x + p00_se.y^2 / n.y)
    ) %>%
    filter(Z.x > Z.y) %>%
    select(Z.x, Z.y, X.x, contains("diff")) %>%
    rename(Z_focal = Z.x, Z_compare = Z.y, X = X.x)

  if (flag) {
    return(out %>% select(-X))
  } else {
    return(out)
  }
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Z_compare", "Z_focal", "d_pred0", "d_pred1", "idx", "p.0i", "p.1i", "p0", "p00.0",
    "p00.0_p00.1_se", "p00.0i", "p00.1", "p00.1i", "p1", "p10.0", "p10.0i", "p10.1",
    "p10.1_p10.0_se", "p10.1i", "p1_p0_se", "pscore", "y_pred0", "y_pred1",
    "W0", "W1", "p00.x", "p00.y", "p00_se.x", "p00_se.y", "p10.x", "p10.y", "p10_se.x", "p10_se.y"
  ))
}
