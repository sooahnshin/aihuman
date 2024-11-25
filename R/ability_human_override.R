#' Compute Risk (Human+AI v. Human) for a Subgroup Defined by AI Recommendation
#'
#' Compute the difference in risk between human+AI and human decision makers, for a subgroup \eqn{\{A_i = a\}}, using AIPW estimators.
#' This can be used for computing how the decision maker overrides the AI recommendation.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param a A specific AI recommendation value to create the subset (numeric: 0 or 1).
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}. If NULL, the function will compute the nuisance functions using the provided data. Note that \code{V} must be provided if \code{nuis_funcs} is NULL.
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
#'  \item \code{tn_fn_diff}: The difference in true negatives and false negatives between human+AI and human decision
#'  \item \code{tn_fn_diff_se}: The standard error of the difference in true negatives and false negatives
#'  \item \code{tp_diff}: The difference in true positives between human+AI and human decision
#'  \item \code{tp_diff_se}: The standard error of the difference in true positives
#'  \item \code{tn_diff}: The difference in true negatives between human+AI and human decision
#'  \item \code{tn_diff_se}: The standard error of the difference in true negatives
#'  \item \code{fn_diff}: The difference in false negatives between human+AI and human decision
#'  \item \code{fn_diff_se}: The standard error of the difference in false negatives
#'  \item \code{fp_diff}: The difference in false positives between human+AI and human decision
#'  \item \code{fp_diff_se}: The standard error of the difference in false positives
#'  }
#'
#' @examples
#' compute_stats_subgroup(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   a = 1,
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
compute_stats_subgroup <- function(Y, D, Z, A, a = 1, nuis_funcs, true.pscore = NULL, X = NULL, l01 = 1) {
  if (!is.vector(Y) || !is.vector(D) || !is.vector(Z) || !is.vector(A)) {
    stop("Y, D, Z, and A must be vectors")
  }
  # Check that Y, D, Z, and A are binary
  if (!all(Y %in% c(0, 1))) stop("Y must be binary")
  if (!all(D %in% c(0, 1))) stop("D must be binary")
  if (!all(Z %in% c(0, 1))) stop("Z must be binary")
  if (!all(A %in% c(0, 1))) stop("A must be binary")
  if (!all(a %in% c(0, 1))) stop("a must be binary")
  
  # Check lengths using all.equal
  if (is.null(X)) {
    if (!isTRUE(all.equal(length(Y), length(D), length(Z), length(A)))) {
      stop("Y, D, Z, and A must have the same length")
    }
  } else if (is.vector(X)) {
    if (!isTRUE(all.equal(length(X), length(Y), length(D), length(Z), length(A)))) {
      stop("X must have the same length as Y, D, Z, and A")
    }
  } else {
    stop("X must be either NULL or a vector")
  }
  
  # Use a placeholder if X is NULL
  flag <- is.null(X)
  if (flag) {
    X <- rep(1, length(Y))
  }
  
  index <- (A == a)
  
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
  }
  
  Y <- Y[index]
  D <- D[index]
  Z <- Z[index]
  A <- A[index]
  X <- X[index]
  nuis_funcs$pscore <- nuis_funcs$pscore[index]
  nuis_funcs$z_models <- nuis_funcs$z_models |> filter(idx %in% index)
  
  dat <- data.frame(Y = Y, D = D, Z = Z, X = X, pscore = nuis_funcs$pscore)
  
  preds <- nuis_funcs$z_models |>
    pivot_wider(names_from = c(Z), values_from = c(-Z, -idx), names_sep = "") |>
    select(-idx)
  
  dat <- dat %>%
    bind_cols(preds)
  
  props <- dat %>%
    mutate(
      p10.1i = (1 - d_pred1) * y_pred1 + Z * (1 - D) * (Y - y_pred1) / pscore - y_pred1 * Z * (D - d_pred1) / pscore,
      p10.0i = (1 - d_pred0) * y_pred0 + (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - y_pred0 * (1 - Z) * (D - d_pred0) / (1 - pscore),
      # the equation below is not exactly the prediction of p00.1, but p00.1 - p00.0 = - TP.1 + TP.0
      p00.1i = (1 - d_pred1) * (y_pred1 - l01) + Z * (1 - D) * (Y - y_pred1) / pscore - (y_pred1 - l01) * Z * (D - d_pred1) / pscore,
      p00.0i = (1 - d_pred0) * (y_pred0 - l01) + (1 - Z) * (1 - D) * (Y - y_pred0) / (1 - pscore) - (y_pred0 - l01) * (1 - Z) * (D - d_pred0) / (1 - pscore),
      p.1i = (1 - d_pred1) * ((1 + l01) * y_pred1 - l01) + (1 + l01) * Z * (1 - D) * (Y - y_pred1) / pscore - ((1 + l01) * y_pred1 - l01) * Z * (D - d_pred1) / pscore,
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
      p00.1_p00.0_se = sqrt(sum(((p00.1i - p00.0i) - (p00.1 - p00.0))^2) / n^2),
      .groups = "drop"
    )
  out <- props %>%
    mutate(
      loss_diff = p1 - p0,
      loss_diff_se = p1_p0_se,
      tn_fn_diff = -p1 + p0,
      tn_fn_diff_se = p1_p0_se,
      tp_diff = -p10.1 + p10.0,
      tp_diff_se = p10.1_p10.0_se,
      tn_diff = -p00.1 + p00.0,
      tn_diff_se = p00.1_p00.0_se,
      fn_diff = p10.1 - p10.0,
      fn_diff_se = p10.1_p10.0_se,
      fp_diff = p00.1 - p00.0,
      fp_diff_se = p00.1_p00.0_se
    ) %>%
    mutate(Z_focal = 1, Z_compare = 0) %>%
    select(Z_focal, Z_compare, X, contains("diff"))
  
  if (flag) {
    return(out %>% select(-X))
  } else {
    return(out)
  }
}
#' Visualize Difference in Risk (Human+AI v. Human) for a Subgroup Defined by AI Recommendation
#'
#' Visualize the the difference in risk between human+AI and human decision makers using AIPW estimators, for a subgroup defined by AI recommendation.
#' Generate a plot based on the overall agreement and subgroup-specific agreement.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param A An AI recommendation (binary: numeric vector of 0 or 1).
#' @param a A specific AI recommendation value to create the subset (numeric: 0 or 1).
#' @param V A matrix of pretreatment covariates (numeric matrix). Optional.
#' @param l01 Ratio of the loss between false positives and false negatives. Default 1.
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}. If NULL, the function will compute the nuisance functions using the provided data. Note that \code{V} must be provided if \code{nuis_funcs} is NULL.
#' @param true.pscore A vector of true propensity scores (numeric), if available. Optional.
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -1.
#' @param p.ub An upper bound for the y-axis (numeric). Default 1.
#' @param y.lab A label for the y-axis (character). Default "Impact of PSA".
#' @param p.label A vector of two labels for the annotations (character). Default c("Human correct", "PSA correct").
#' @param label A label for the plot (character). Default "TNP - FNP".
#' @param metrics A vector of metrics to include in the plot (character). Default c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion").
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_diff_subgroup(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   a = 1,
#'   l01 = 1,
#'   nuis_funcs = nuis_func,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male"),
#'   p.title = NULL, p.lb = -0.5, p.ub = 0.5,
#'   label = "TNP - FNP",
#'   metrics = c("True Negative Proportion (TNP)", "False Negative Proportion (FNP)", "TNP - FNP")
#' )
#'
#' @useDynLib aihuman, .registration=TRUE
#' @export
plot_diff_subgroup <- function(Y, D, Z, A, a = 1, V = NULL, l01 = l01,
                               nuis_funcs = NULL, true.pscore = NULL,
                               subgroup1, subgroup2,
                               label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                               x.order = NULL,
                               p.title = NULL, p.lb = -1, p.ub = 1, y.lab = "Impact of PSA", p.label = c("Human correct", "PSA correct"),
                               label = "TNP - FNP", metrics = c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion")) {
  if (is.null(nuis_funcs)) {
    nuis_funcs <- compute_nuisance_functions(Y = Y, D = D, Z = Z, V = V)
  }
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
  }
  stats_subset <- compute_stats_subgroup(Y = Y, D = D, Z = Z, A = A, a = a, nuis_funcs = nuis_funcs, X = NULL, l01 = l01)
  stats_subset_race <- compute_stats_subgroup(Y = Y, D = D, Z = Z, A = A, a = a, nuis_funcs = nuis_funcs, X = subgroup1, l01 = l01)
  stats_subset_gender <- compute_stats_subgroup(Y = Y, D = D, Z = Z, A = A, a = a, nuis_funcs = nuis_funcs, X = subgroup2, l01 = l01)
  
  df_subset <- stats_subset %>%
    mutate(cov = "Overall", X = "Overall") %>%
    bind_rows(stats_subset_race %>% mutate(cov = label.subgroup1)) %>%
    bind_rows(stats_subset_gender %>% mutate(cov = label.subgroup2))
  
  df_subset <- df_subset |>
    pivot_longer(cols = -c("X", "cov", "Z_focal", "Z_compare")) |>
    mutate(
      metric = case_when(
        str_detect(name, "tn_fn_diff") ~ label,
        str_detect(name, "tp_diff") ~ "True Positive Proportion (TPP)",
        str_detect(name, "tn_diff") ~ "True Negative Proportion (TNP)",
        str_detect(name, "fn_diff") ~ "False Negative Proportion (FNP)",
        str_detect(name, "fp_diff") ~ "False Positive Proportion (FPP)",
        TRUE ~ ""
      ),
      se = ifelse(str_detect(name, "_se"), "se", "est")
    ) |>
    select(-name) |>
    pivot_wider(values_from = "value", names_from = "se")
  
  vis_diff_subgroup(df_subset,
                    label.subgroup1 = label.subgroup1, label.subgroup2 = label.subgroup2, x.order = x.order,
                    p.title = p.title, p.lb = p.lb, p.ub = p.ub, y.lab = y.lab, p.label = p.label, metrics = metrics
  )
}

#' Visualize Risk (Human+AI v. Human; internal)
#'
#' Internal function to visualize the difference in risk between human+AI and human decision makers, for a subgroup defined by AI recommendation.
#'
#' @param df A data frame generated by \code{compute_stats_aipw} and reshaped.
#' @param label.subgroup1 A label for subgroup1 (character).
#' @param label.subgroup2 A label for subgroup2 (character).
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -0.2.
#' @param p.ub An upper bound for the y-axis (numeric). Default 0.2.
#' @param y.lab A label for the y-axis (character). Default "Impact of PSA".
#' @param p.label A vector of two labels for the annotations (character). Default c("PSA harms", "PSA helps").
#' @param metrics A vector of metrics to include in the plot (character). Default c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion").
#'
#' @return A ggplot object.
#'
#' @importFrom forcats fct_relevel
#'
#'
#' @keywords internal
#' @useDynLib aihuman, .registration = TRUE
#'
vis_diff_subgroup <- function(df,
                              label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                              x.order = NULL,
                              p.title = NULL, p.lb = -0.2, p.ub = 0.2,
                              y.lab = "Impact of PSA", p.label = c("Human correct", "PSA correct"),
                              metrics = c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion")) {
  annotations <- data.frame(
    metric = c(metrics[3], metrics[3]),
    x = c(1 - 0.2, 1 - 0.2),
    y.start = c(p.ub - 0.3 * p.ub, p.lb + 0.3 * p.ub),
    y.end = c(p.ub, p.lb),
    label = p.label
  )
  annotations_bottom <- annotations %>%
    filter(label == p.label[1]) # Assuming this is your bottom annotation
  
  annotations_top <- annotations %>%
    filter(label == p.label[2]) # Assuming this is your top annotation
  
  df %>%
    filter(
      metric %in% metrics,
      cov %in% c("Overall", label.subgroup1, label.subgroup2)
    ) %>%
    mutate(
      cov = factor(cov, levels = c("Overall", label.subgroup1, label.subgroup2)),
      X = factor(X, levels = x.order)
    ) %>%
    ggplot(aes(x = X, color = cov, shape = cov)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_pointrange(aes(y = est, ymin = est - qnorm(0.975) * se, ymax = est + qnorm(0.975) * se),
                    position = position_dodge(.5), lineend = "round", linewidth = 1.2, size = 1.2
    ) +
    facet_grid(. ~ fct_relevel(metric, metrics), space = "free_x", scales = "free") +
    ylab(y.lab) +
    scale_color_brewer("", type = "qual", palette = "Set1") +
    ggtitle(p.title) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_y_continuous(limits = c(p.lb, p.ub)) +
    geom_text(
      data = annotations_bottom, aes(x = x, y = y.end, label = label),
      inherit.aes = FALSE, vjust = -1, hjust = 1.2, angle = 90
    ) +
    geom_text(
      data = annotations_top, aes(x = x, y = y.end, label = label),
      inherit.aes = FALSE, vjust = -1, hjust = -0.2, angle = 90
    ) +
    geom_segment(
      data = annotations, aes(x = x, xend = x, y = y.start, yend = y.end),
      arrow = arrow(length = unit(0.1, "inches")), inherit.aes = FALSE
    )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "p00.1_p00.0_se"
  ))
}