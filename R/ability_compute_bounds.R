#' Compute Risk (AI v. Human)
#'
#' Compute the difference in risk between AI and human decision makers using AIPW estimators.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param A An observed AI recommendation (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param X Pretreatment covariate used for subgroup analysis (vector). Must be the same length as Y, D, Z, and A if provided. Default is NULL.
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}
#' @param nuis_funcs_ai output from \code{\link{compute_nuisance_functions_ai}}
#' @param true.pscore A vector of true propensity scores (numeric), if available. Optional.
#' @param l01 Ratio of the loss between false positives and false negatives
#'
#' @return A tibble the following columns:
#' \itemize{
#'  \item \code{Z_focal}: The focal treatment indicator. `1` indicates the treatment group.
#'  \item \code{Z_compare}: The comparison treatment indicator. `0` indicates the control group.
#'  \item \code{X}: Pretreatment covariate (if provided).
#'  \item \code{fn_diff_lb}: The lower bound of difference in false negatives
#'  \item \code{fn_diff_ub}: The upper bound of difference in false negatives
#'  \item \code{fp_diff_lb}: The lower bound of difference in false positives
#'  \item \code{fp_diff_ub}: The upper bound of difference in false positives
#'  \item \code{loss_diff_lb}: The lower bound of difference in loss
#'  \item \code{loss_diff_ub}: The upper bound of difference in loss
#'  \item \code{fn_diff_lb_se}: The standard error of the difference in false negatives
#'  \item \code{fn_diff_ub_se}: The standard error of the difference in false negatives
#'  \item \code{fp_diff_lb_se}: The standard error of the difference in false positives
#'  \item \code{fp_diff_ub_se}: The standard error of the difference in false positives
#'  \item \code{loss_diff_lb_se}: The standard error of the difference in loss
#'  \item \code{loss_diff_ub_se}: The standard error of the difference in loss
#'  }
#'
#'
#'
#' @examples
#' compute_bounds_aipw(
#'   Y = NCAdata$Y,
#'   A = PSAdata$DMF,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   nuis_funcs = nuis_func,
#'   nuis_funcs_ai = nuis_func_ai,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   X = NULL,
#'   l01 = 1
#' )
#' @importFrom tidyr separate
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
#'
compute_bounds_aipw <- function(Y, A, D, Z, X = NULL, nuis_funcs, nuis_funcs_ai, true.pscore = NULL, l01 = 1) {
  if (!is.vector(Y) || !is.vector(A) || !is.vector(D) || !is.vector(Z)) {
    stop("Y, A, D, and Z must be vectors")
  }
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary")
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
    if (!isTRUE(all.equal(length(Y), length(A), length(D), length(Z)))) {
      stop("Y, D, and Z must have the same length")
    }
  } else if (is.vector(X)) {
    if (!isTRUE(all.equal(length(X), length(Y), length(A), length(D), length(Z)))) {
      stop("X must have the same length as Y, A, D, and Z")
    }
  } else {
    stop("X must be either NULL or a vector")
  }
  
  if (flag) {
    X <- rep(1, length(Y))
  }
  
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
    nuis_funcs_ai$pscore <- true.pscore
  }
  
  # fit different outcome/decision models for each Z
  data <- data.frame(Y = Y, D = D, Z = Z, A = A, X = X, pscore = nuis_funcs$pscore)
  
  preds1 <- nuis_funcs$z_models |>
    pivot_wider(names_from = c(Z), values_from = c(-Z, -idx), names_sep = "") |>
    select(-idx)
  
  preds2 <- nuis_funcs_ai$z_models |>
    pivot_wider(names_from = c(Z, A), values_from = c(-Z, -A, -idx), names_sep = "") |>
    select(-idx)
  
  data <- data %>%
    bind_cols(preds1, preds2)
  # d_pred.z, y_pred.z: Pr(D=1|Z=z,X=x) and Pr(Y=1|D=0,Z=z,X=x)
  # d_pred.za, y_pred.za: Pr(D=1|Z=z,A=a,X=x) and Pr(Y=1|D=0,Z=z,A=a,X=x)
  
  if(ncol(preds1)==4 & ncol(preds2)==8){
    props <- data %>%
      mutate(
        # gL1.i = 1 if z'=0 i.e. Pr(Y=1,D=0|A=0,Z=0,X=x) > Pr(Y=1,D=0|A=0,Z=1,X=x)
        gL1.i = 1 * ((1 - d_pred00) * y_pred00 >= (1 - d_pred10) * y_pred10),
        # gU1.i = 1 if z'=0 i.e. Pr(Y=0,D=0|A=0,Z=0,X=x) > Pr(Y=0,D=0|A=0,Z=1,X=x)
        gU1.i = 1 * ((1 - d_pred00) * (1 - y_pred00) >= (1 - d_pred10) * (1 - y_pred10)),
        # gL0.i = 1 if z'=1 i.e. Pr(Y=1,D=0|A=0,Z=1,X=x) > Pr(Y=1,D=0|A=0,Z=0,X=x)
        gL0.i = 1 * ((1 - d_pred10) * y_pred10 >= (1 - d_pred00) * y_pred00),
        # gU0.i = 1 if z'=1 i.e. Pr(Y=0,D=0|A=0,Z=1,X=x) > Pr(Y=0,D=0|A=0,Z=0,X=x)
        gU0.i = 1 * ((1 - d_pred10) * (1 - y_pred10) >= (1 - d_pred00) * (1 - y_pred00))
      ) %>%
      mutate(
        # phi_z(x): Pr(Y=1,D=0|Z=z,X=x)
        p1.i = (1 - d_pred1) * y_pred1 + Z * (1 - D) * (Y - y_pred1) / pscore - y_pred1 * Z * (D - d_pred1) / pscore,
        p0.i = (1 - d_pred0) * y_pred0 + (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - y_pred0 * (1 - Z) * (D - d_pred0) / (1 - pscore),
        # phi_z1(x): Pr(Y=1,D=0,A=0|Z=z,X=x)
        p11.i = (1 - A) * (1 - d_pred10) * y_pred10 + Z * (1 - A) * (1 - D) * (Y - y_pred10) / pscore - y_pred10 * Z * (1 - A) * (D - d_pred10) / pscore,
        p01.i = (1 - A) * (1 - d_pred00) * y_pred00 + (1 - Z) * (1 - A) * (1 - D) * (Y - y_pred00) / (1 - pscore) - y_pred00 * (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore),
        # phi_z0(x): Pr(Y=0,D=0,A=0|Z=z,X=x)
        p10.i = (1 - A) * (1 - d_pred10) * (1 - y_pred10) - Z * (1 - A) * (1 - D) * (Y - y_pred10) / pscore - (1 - y_pred10) * Z * (1 - A) * (D - d_pred10) / pscore,
        p00.i = (1 - A) * (1 - d_pred00) * (1 - y_pred00) - (1 - Z) * (1 - A) * (1 - D) * (Y - y_pred00) / (1 - pscore) - (1 - y_pred00) * (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore),
        # phi_z1^D(x): Pr(D=0,A=1|Z=z,X=x)
        p11.D.i = A * (1 - d_pred11) - Z * A * (D - d_pred11) / pscore,
        p01.D.i = A * (1 - d_pred01) - (1 - Z) * A * (D - d_pred01) / (1 - pscore),
        # phi_z0^D(x): Pr(D=0,A=0|Z=z,X=x)
        p10.D.i = (1 - A) * d_pred10 + Z * (1 - A) * (D - d_pred10) / pscore,
        p00.D.i = (1 - A) * d_pred00 + (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore)
      ) %>%
      mutate(
        p11_gL1.i = p11.i * gL1.i,
        p11_gL0.i = p11.i * gL0.i,
        p01_gL1.i = p01.i * gL1.i,
        p01_gL0.i = p01.i * gL0.i,
        p10_gL1.i = p10.i * gL1.i,
        p10_gL0.i = p10.i * gL0.i,
        p00_gL1.i = p00.i * gL1.i,
        p00_gL0.i = p00.i * gL0.i,
        p11_gU1.i = p11.i * gU1.i,
        p11_gU0.i = p11.i * gU0.i,
        p01_gU1.i = p01.i * gU1.i,
        p01_gU0.i = p01.i * gU0.i,
        p10_gU1.i = p10.i * gU1.i,
        p10_gU0.i = p10.i * gU0.i,
        p00_gU1.i = p00.i * gU1.i,
        p00_gU0.i = p00.i * gU0.i
      ) %>%
      mutate(
        fn_diff_lb.1i = p11.i - p1.i + p01_gL1.i - p11_gL1.i,
        fn_diff_ub.1i = p11.i - p1.i + p10.D.i - p00_gU1.i + p10_gU1.i,
        fn_diff_lb.0i = p01.i - p0.i + p11_gL0.i - p01_gL0.i,
        fn_diff_ub.0i = p01.i - p0.i + p00.D.i - p10_gU0.i + p00_gU0.i,
        fp_diff_lb.1i = p11.i - p1.i + p11.D.i - p10.D.i + p01_gL1.i - p11_gL1.i,
        fp_diff_ub.1i = p11.i - p1.i + p11.D.i - p00_gU1.i + p10_gU1.i,
        fp_diff_lb.0i = p01.i - p0.i + p01.D.i - p00.D.i + p11_gL0.i - p01_gL0.i,
        fp_diff_ub.0i = p01.i - p0.i + p01.D.i - p10_gU0.i + p00_gU0.i,
        loss_diff_lb.1i = fn_diff_lb.1i + l01 * fp_diff_lb.1i,
        loss_diff_ub.1i = fn_diff_ub.1i + l01 * fp_diff_ub.1i,
        loss_diff_lb.0i = fn_diff_lb.0i + l01 * fp_diff_lb.0i,
        loss_diff_ub.0i = fn_diff_ub.0i + l01 * fp_diff_ub.0i
      )
    out <- props %>%
      group_by(X) %>%
      summarise(
        n = n(),
        fn_diff_lb.1 = mean(fn_diff_lb.1i),
        fn_diff_ub.1 = mean(fn_diff_ub.1i),
        fn_diff_lb.0 = mean(fn_diff_lb.0i),
        fn_diff_ub.0 = mean(fn_diff_ub.0i),
        fp_diff_lb.1 = mean(fp_diff_lb.1i),
        fp_diff_ub.1 = mean(fp_diff_ub.1i),
        fp_diff_lb.0 = mean(fp_diff_lb.0i),
        fp_diff_ub.0 = mean(fp_diff_ub.0i),
        loss_diff_lb.1 = mean(loss_diff_lb.1i),
        loss_diff_ub.1 = mean(loss_diff_ub.1i),
        loss_diff_lb.0 = mean(loss_diff_lb.0i),
        loss_diff_ub.0 = mean(loss_diff_ub.0i),
        fn_diff_lb_se.1 = sqrt(sum((fn_diff_lb.1i - fn_diff_lb.1)^2) / n^2),
        fn_diff_ub_se.1 = sqrt(sum((fn_diff_ub.1i - fn_diff_ub.1)^2) / n^2),
        fn_diff_lb_se.0 = sqrt(sum((fn_diff_lb.0i - fn_diff_lb.0)^2) / n^2),
        fn_diff_ub_se.0 = sqrt(sum((fn_diff_ub.0i - fn_diff_ub.0)^2) / n^2),
        fp_diff_lb_se.1 = sqrt(sum((fp_diff_lb.1i - fp_diff_lb.1)^2) / n^2),
        fp_diff_ub_se.1 = sqrt(sum((fp_diff_ub.1i - fp_diff_ub.1)^2) / n^2),
        fp_diff_lb_se.0 = sqrt(sum((fp_diff_lb.0i - fp_diff_lb.0)^2) / n^2),
        fp_diff_ub_se.0 = sqrt(sum((fp_diff_ub.0i - fp_diff_ub.0)^2) / n^2),
        loss_diff_lb_se.1 = sqrt(sum((loss_diff_lb.1i - loss_diff_lb.1)^2) / n^2),
        loss_diff_ub_se.1 = sqrt(sum((loss_diff_ub.1i - loss_diff_ub.1)^2) / n^2),
        loss_diff_lb_se.0 = sqrt(sum((loss_diff_lb.0i - loss_diff_lb.0)^2) / n^2),
        loss_diff_ub_se.0 = sqrt(sum((loss_diff_ub.0i - loss_diff_ub.0)^2) / n^2),
        .groups = "drop"
      ) %>%
      select(-n) %>%
      pivot_longer(cols = -X, names_to = "metric", values_to = "value") %>%
      separate(metric, into = c("metric", "Z_compare"), sep = "\\.") %>%
      pivot_wider(names_from = "metric", values_from = "value") %>%
      mutate(Z_focal = "AI") %>%
      relocate(Z_focal, Z_compare, X, contains("diff"))
  } else if(ncol(preds1)==2 & ncol(preds2)==4 & all(unique(Z)==0)) {
    props <- data %>%
      mutate(
        # phi_z(x): Pr(Y=1,D=0|Z=z,X=x)
        p0.i = (1 - d_pred0) * y_pred0 + (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - y_pred0 * (1 - Z) * (D - d_pred0) / (1 - pscore),
        # phi_z1(x): Pr(Y=1,D=0,A=0|Z=z,X=x)
        p01.i = (1 - A) * (1 - d_pred00) * y_pred00 + (1 - Z) * (1 - A) * (1 - D) * (Y - y_pred00) / (1 - pscore) - y_pred00 * (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore),
        # phi_z0(x): Pr(Y=0,D=0,A=0|Z=z,X=x)
        p00.i = (1 - A) * (1 - d_pred00) * (1 - y_pred00) - (1 - Z) * (1 - A) * (1 - D) * (Y - y_pred00) / (1 - pscore) - (1 - y_pred00) * (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore),
        # phi_z1^D(x): Pr(D=0,A=1|Z=z,X=x)
        p01.D.i = A * (1 - d_pred01) - (1 - Z) * A * (D - d_pred01) / (1 - pscore),
        # phi_z0^D(x): Pr(D=0,A=0|Z=z,X=x)
        p00.D.i = (1 - A) * d_pred00 + (1 - Z) * (1 - A) * (D - d_pred00) / (1 - pscore)
      ) %>%
      mutate(
        fn_diff_lb.0i = p01.i - p0.i,
        fn_diff_ub.0i = p01.i - p0.i + p00.D.i,
        fp_diff_lb.0i = p01.i - p0.i + p01.D.i - p00.D.i,
        fp_diff_ub.0i = p01.i - p0.i + p01.D.i,
        loss_diff_lb.0i = fn_diff_lb.0i + l01 * fp_diff_lb.0i,
        loss_diff_ub.0i = fn_diff_ub.0i + l01 * fp_diff_ub.0i
      )
    out <- props %>%
      group_by(X) %>%
      summarise(
        n = n(),
        fn_diff_lb.0 = mean(fn_diff_lb.0i),
        fn_diff_ub.0 = mean(fn_diff_ub.0i),
        fp_diff_lb.0 = mean(fp_diff_lb.0i),
        fp_diff_ub.0 = mean(fp_diff_ub.0i),
        loss_diff_lb.0 = mean(loss_diff_lb.0i),
        loss_diff_ub.0 = mean(loss_diff_ub.0i),
        fn_diff_lb_se.0 = sqrt(sum((fn_diff_lb.0i - fn_diff_lb.0)^2) / n^2),
        fn_diff_ub_se.0 = sqrt(sum((fn_diff_ub.0i - fn_diff_ub.0)^2) / n^2),
        fp_diff_lb_se.0 = sqrt(sum((fp_diff_lb.0i - fp_diff_lb.0)^2) / n^2),
        fp_diff_ub_se.0 = sqrt(sum((fp_diff_ub.0i - fp_diff_ub.0)^2) / n^2),
        loss_diff_lb_se.0 = sqrt(sum((loss_diff_lb.0i - loss_diff_lb.0)^2) / n^2),
        loss_diff_ub_se.0 = sqrt(sum((loss_diff_ub.0i - loss_diff_ub.0)^2) / n^2),
        .groups = "drop"
      ) %>%
      select(-n) %>%
      pivot_longer(cols = -X, names_to = "metric", values_to = "value") %>%
      separate(metric, into = c("metric", "Z_compare"), sep = "\\.") %>%
      pivot_wider(names_from = "metric", values_from = "value") %>%
      mutate(Z_focal = "AI") %>%
      relocate(Z_focal, Z_compare, X, contains("diff"))
  } else if(ncol(preds1)==2 & ncol(preds2)==4 & all(unique(Z)==1)) {
    props <- data %>%
      mutate(
        # phi_z(x): Pr(Y=1,D=0|Z=z,X=x)
        p1.i = (1 - d_pred1) * y_pred1 + Z * (1 - D) * (Y - y_pred1) / pscore - y_pred1 * Z * (D - d_pred1) / pscore,
        # phi_z1(x): Pr(Y=1,D=0,A=0|Z=z,X=x)
        p11.i = (1 - A) * (1 - d_pred10) * y_pred10 + Z * (1 - A) * (1 - D) * (Y - y_pred10) / pscore - y_pred10 * Z * (1 - A) * (D - d_pred10) / pscore,
        # phi_z0(x): Pr(Y=0,D=0,A=0|Z=z,X=x)
        p10.i = (1 - A) * (1 - d_pred10) * (1 - y_pred10) - Z * (1 - A) * (1 - D) * (Y - y_pred10) / pscore - (1 - y_pred10) * Z * (1 - A) * (D - d_pred10) / pscore,
        # phi_z1^D(x): Pr(D=0,A=1|Z=z,X=x)
        p11.D.i = A * (1 - d_pred11) - Z * A * (D - d_pred11) / pscore,
        # phi_z0^D(x): Pr(D=0,A=0|Z=z,X=x)
        p10.D.i = (1 - A) * d_pred10 + Z * (1 - A) * (D - d_pred10) / pscore,
      ) %>%
      mutate(
        fn_diff_lb.1i = p11.i - p1.i,
        fn_diff_ub.1i = p11.i - p1.i + p10.D.i,
        fp_diff_lb.1i = p11.i - p1.i + p11.D.i - p10.D.i,
        fp_diff_ub.1i = p11.i - p1.i + p11.D.i,
        loss_diff_lb.1i = fn_diff_lb.1i + l01 * fp_diff_lb.1i,
        loss_diff_ub.1i = fn_diff_ub.1i + l01 * fp_diff_ub.1i,
      )
    out <- props %>%
      group_by(X) %>%
      summarise(
        n = n(),
        fn_diff_lb.1 = mean(fn_diff_lb.1i),
        fn_diff_ub.1 = mean(fn_diff_ub.1i),
        fp_diff_lb.1 = mean(fp_diff_lb.1i),
        fp_diff_ub.1 = mean(fp_diff_ub.1i),
        loss_diff_lb.1 = mean(loss_diff_lb.1i),
        loss_diff_ub.1 = mean(loss_diff_ub.1i),
        fn_diff_lb_se.1 = sqrt(sum((fn_diff_lb.1i - fn_diff_lb.1)^2) / n^2),
        fn_diff_ub_se.1 = sqrt(sum((fn_diff_ub.1i - fn_diff_ub.1)^2) / n^2),
        fp_diff_lb_se.1 = sqrt(sum((fp_diff_lb.1i - fp_diff_lb.1)^2) / n^2),
        fp_diff_ub_se.1 = sqrt(sum((fp_diff_ub.1i - fp_diff_ub.1)^2) / n^2),
        loss_diff_lb_se.1 = sqrt(sum((loss_diff_lb.1i - loss_diff_lb.1)^2) / n^2),
        loss_diff_ub_se.1 = sqrt(sum((loss_diff_ub.1i - loss_diff_ub.1)^2) / n^2),
        .groups = "drop"
      ) %>%
      select(-n) %>%
      pivot_longer(cols = -X, names_to = "metric", values_to = "value") %>%
      separate(metric, into = c("metric", "Z_compare"), sep = "\\.") %>%
      pivot_wider(names_from = "metric", values_from = "value") %>%
      mutate(Z_focal = "AI") %>%
      relocate(Z_focal, Z_compare, X, contains("diff"))
  } else {
    stop("The number of columns in the nuisance functions output is not as expected.")
  }
  
  if (flag) {
    return(out %>% select(-X))
  } else {
    return(out)
  }
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "d_pred00", "d_pred01", "d_pred10", "d_pred11", 
    "fn_diff_lb.0", "fn_diff_lb.0i", "fn_diff_lb.1", "fn_diff_lb.1i", 
    "fn_diff_ub.0", "fn_diff_ub.0i", "fn_diff_ub.1", "fn_diff_ub.1i", 
    "fp_diff_lb.0", "fp_diff_lb.0i", "fp_diff_lb.1", "fp_diff_lb.1i", 
    "fp_diff_ub.0", "fp_diff_ub.0i", "fp_diff_ub.1", "fp_diff_ub.1i", 
    "gL0.i", "gL1.i", "gU0.i", "gU1.i", 
    "loss_diff_lb.0", "loss_diff_lb.0i", "loss_diff_lb.1", "loss_diff_lb.1i", 
    "loss_diff_ub.0", "loss_diff_ub.0i", "loss_diff_ub.1", "loss_diff_ub.1i", 
    "p0.i", "p00.D.i", "p00.i", "p00_gU0.i", 
    "p00_gU1.i", "p01.D.i", "p01.i", "p01_gL0.i", 
    "p01_gL1.i", "p1.i", "p10.D.i", "p10.i", 
    "p10_gU0.i", "p10_gU1.i", "p11.D.i", "p11.i", 
    "p11_gL0.i", "p11_gL1.i", "y_pred00", "y_pred10"
  ))
}
