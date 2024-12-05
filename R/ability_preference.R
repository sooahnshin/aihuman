#' Visualize Preference
#'
#' Compute the difference in risk between AI and human decision makers using AIPW estimators
#' over a set of loss ratios, and then visualize when we prefer human over AI decision makers.
#' Generate a plot based on the overall and subgroup-specific results.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param V A matrix of pretreatment covariates (numeric matrix). Optional.
#' @param A An observed AI recommendation (binary: numeric vector of 0 or 1).
#' @param z_compare A compare treatment indicator (numeric). Default 0.
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}. If NULL, the function will compute the nuisance functions using the provided data. Note that \code{V} must be provided if \code{nuis_funcs} is NULL.
#' @param nuis_funcs_ai output from \code{\link{compute_nuisance_functions_ai}}
#' @param true.pscore A vector of true propensity scores (numeric), if available. Optional.
#' @param l01_seq A candidate list of ratio of the loss between false positives and false negatives. Default \code{10^seq(-2, 2, length.out = 100)}.
#' @param alpha A significance level (numeric). Default 0.05.
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param p.title A title for the plot (character). Default NULL.
#' @param legend.position Position of the legend (character).
#' @param p.label A vector of three labels for the annotations (character). Default c("AI-alone preferred", "Human-alone preferred", "Ambiguous").
#'
#' @importFrom purrr map_dfr
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_preference(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   z_compare = 0,
#'   nuis_funcs = nuis_func,
#'   nuis_funcs_ai = nuis_func_ai,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   l01_seq = 10^seq(-2, 2, length.out = 10),
#'   alpha = 0.05,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male"),
#'   p.title = NULL, legend.position = "none",
#'   p.label = c("AI-alone preferred", "Human-alone preferred", "Ambiguous")
#' )
#' 
#' @useDynLib aihuman, .registration = TRUE
#' @export

plot_preference <- function(Y, D, Z, V = NULL, A,
                            z_compare = 0,
                            true.pscore = NULL, nuis_funcs = NULL, nuis_funcs_ai = NULL,
                            l01_seq = 10^seq(-2, 2, length.out = 100),
                            alpha = 0.05,
                            subgroup1, subgroup2,
                            label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                            x.order = NULL,
                            p.title = NULL, legend.position = "none",
                            p.label = c("AI-alone preferred", "Human-alone preferred", "Ambiguous")) {
  if (is.null(nuis_funcs)) {
    nuis_funcs <- compute_nuisance_functions(Y = Y, D = D, Z = Z, V = V)
  }
  if (is.null(nuis_funcs_ai)) {
    nuis_funcs_ai <- compute_nuisance_functions_ai(Y = Y, D = D, Z = Z, A = A, V = V)
  }
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
    nuis_funcs_ai$pscore <- true.pscore
  }
  ## Use future_map_dfr to parallelize the computation if needed
  diff_df <- map_dfr(l01_seq, function(l01) {
    diff_bounds <- compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, l01 = l01, X = NULL)
    diff_bounds_sub1 <- compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, X = subgroup1, l01 = l01)
    diff_bounds_sub2 <- compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, X = subgroup2, l01 = l01)

    diff_bounds |>
      filter(Z_compare == z_compare) |>
      mutate(tstat_ub = loss_diff_ub / loss_diff_ub_se, tstat_lb = loss_diff_lb / loss_diff_lb_se, l01 = l01, cov = "Overall", X = "Overall") |>
      relocate(Z_focal, Z_compare, X) |>
      bind_rows(diff_bounds_sub1 |> filter(Z_compare == z_compare) |> mutate(tstat_ub = loss_diff_ub / loss_diff_ub_se, tstat_lb = loss_diff_lb / loss_diff_lb_se, l01 = l01, cov = label.subgroup1)) |>
      bind_rows(diff_bounds_sub2 |> filter(Z_compare == z_compare) |> mutate(tstat_ub = loss_diff_ub / loss_diff_ub_se, tstat_lb = loss_diff_lb / loss_diff_lb_se, l01 = l01, cov = label.subgroup2))
  })

  vis_df <- diff_df |>
    mutate(pref = factor(
      case_when(
        tstat_ub <= qnorm(alpha) ~ "Prefer AI",
        tstat_lb >= qnorm(1 - alpha) ~ "Prefer Human",
        TRUE ~ "Ambiguous"
      ),
      levels = c("Prefer AI", "Prefer Human", "Ambiguous")
    ))

  p_human <- vis_df |>
    filter(Z_compare == z_compare) |>
    mutate(X = factor(X, levels = x.order)) |>
    vis_preference(p.title = p.title, legend.position = legend.position, p.label = p.label)

  return(p_human)
}

#' Visualize Preference (internal)
#'
#' Internal function to visualize preference for Human over AI decision makers.
#'
#'
#' @param df Data frame with columns `X`, `l01`, and `pref`.
#' @param p.title Plot title (character).
#' @param legend.position Position of the legend (character).
#' @param p.label Labels for the lines (character).
#'
#' @return A ggplot object.
#'
#'
#' @keywords internal
#' @useDynLib aihuman, .registration = TRUE
vis_preference <- function(df,
                           p.title = NULL,
                           legend.position = "none",
                           p.label = c("AI-alone preferred", "Human-alone preferred", "Ambiguous")) {
  ggplot(df, aes(x = X)) +
    geom_hline(yintercept = 1, lty = 2, color = "grey50") +
    geom_line(aes(y = l01, color = pref, lty = pref), linewidth = 2) +
    xlab(NULL) +
    ggtitle(p.title) +
    # scale_y_continuous(expression(frac("Loss of False Positive", "Loss of False Negative") ~ ~ ("\u2113"["01"])),
    scale_y_continuous(expression(frac("Loss of False Positive", "Loss of False Negative") ~ ~ ("l"["01"])),
      trans = "log10",
      breaks = c(1 / 100, 1 / 50, 1 / 25, 1 / 10, 1 / 5, 1 / 2, 1, 2, 5, 10, 25, 50, 100),
      labels = c(
        paste(c(100, 50, 25, 10, 5, 2), "x FN"), "1-1",
        paste(rev(c(100, 50, 25, 10, 5, 2)), "x FP")
      )
    ) +
    scale_color_manual("",
      values = c("#c51b7d", "#4d9221", "grey30"),
      labels = p.label, drop = FALSE
    ) +
    scale_linetype_manual("", values = c(1, 1, 3), labels = p.label, drop = FALSE) +
    theme(legend.position = legend.position)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "l01", "loss_diff_lb", "loss_diff_lb_se",
    "loss_diff_ub", "loss_diff_ub_se", "pref"
  ))
}