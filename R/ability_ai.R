#' Visualize Difference in Risk (AI v. Human)
#'
#' Visualize the difference in risk between AI and human decision makers using AIPW estimators.
#' Generate a plot based on the overall and subgroup-specific results.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param V A matrix of pretreatment covariates (numeric matrix). Optional.
#' @param z_compare A compare treatment indicator (numeric). Default 0.
#' @param A An observed AI recommendation (binary: numeric vector of 0 or 1).
#' @param l01 Ratio of the loss between false positives and false negatives. Default 1.
#' @param nuis_funcs output from \code{\link{compute_nuisance_functions}}. If NULL, the function will compute the nuisance functions using the provided data. Note that \code{V} must be provided if \code{nuis_funcs} is NULL.
#' @param nuis_funcs_ai output from \code{\link{compute_nuisance_functions_ai}}
#' @param true.pscore A vector of true propensity scores (numeric), if available. Optional.
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param zero.line A logical indicating whether to include a zero line. Default TRUE.
#' @param arrows A logical indicating whether to include arrows. Default TRUE.
#' @param y.min A lower bound for the y-axis (numeric). Default -Inf.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -0.2.
#' @param p.ub An upper bound for the y-axis (numeric). Default 0.2.
#' @param y.lab A label for the y-axis (character). Default "PSA versus Human".
#' @param p.label A vector of two labels for the annotations (character). Default c("PSA harms", "PSA helps").
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_diff_ai_aipw(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   A = PSAdata$DMF,
#'   z_compare = 0,
#'   nuis_funcs = nuis_func,
#'   nuis_funcs_ai = nuis_func_ai,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   l01 = 1,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male"),
#'   zero.line = TRUE, arrows = TRUE, y.min = -Inf,
#'   p.title = NULL, p.lb = -0.3, p.ub = 0.3
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
plot_diff_ai_aipw <- function(Y, D, Z, V = NULL, A, 
                              z_compare = 0, l01 = 1,
                              nuis_funcs = NULL, nuis_funcs_ai = NULL,
                              true.pscore = NULL,
                              subgroup1, subgroup2,
                              label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                              x.order = NULL,
                              zero.line = TRUE, arrows = TRUE, y.min = -Inf,
                              p.title = NULL, p.lb = -1, p.ub = 1, y.lab = "PSA versus Human", p.label = c("PSA worse", "PSA better")) {
  if(!z_compare %in% c(0, 1)) {
    stop("z_compare must be 0 or 1")
  }
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
  diff_bounds <- compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, l01 = l01, X = NULL) |>
    mutate(cov = "Overall", X = "Overall") |>
    relocate(Z_focal, Z_compare, X) |>
    bind_rows(compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, X = subgroup1, l01 = l01) |> mutate(cov = label.subgroup1)) |>
    bind_rows(compute_bounds_aipw(Y = Y, D = D, Z = Z, A = A, nuis_funcs = nuis_funcs, nuis_funcs_ai = nuis_funcs_ai, X = subgroup2, l01 = l01) |> mutate(cov = label.subgroup2))

  diff_ai_human <- diff_bounds |>
    filter(Z_compare == z_compare)

  diff_ai_human_vis <- diff_ai_human |>
    pivot_longer(cols = -c("Z_focal", "Z_compare", "X", "cov")) |>
    mutate(
      metric = case_when(
        str_detect(name, "loss") ~ "Misclassification Rate",
        str_detect(name, "fn") ~ "False Negative Proportion",
        str_detect(name, "fp") ~ "False Positive Proportion",
        TRUE ~ ""
      ),
      se = case_when(
        str_detect(name, "_se") ~ "se",
        TRUE ~ "est"
      ),
      ublb = case_when(
        str_detect(name, "ub") ~ "ub",
        str_detect(name, "lb") ~ "lb",
        TRUE ~ NA
      )
    ) |>
    select(-name) |>
    pivot_wider(values_from = "value", names_from = c("se", "ublb"))

  vis_diff_ai(diff_ai_human_vis,
    label.subgroup1 = label.subgroup1,
    label.subgroup2 = label.subgroup2,
    x.order = x.order,
    zero.line = zero.line, arrows = arrows, y.min = y.min,
    p.title = p.title, p.lb = p.lb, p.ub = p.ub, y.lab = y.lab, p.label = p.label
  )
}


#' Visualize Risk (AI v. Human; internal)
#'
#' Internal function to visualize the difference in risk between AI and human decision makers.
#'
#' @param df A data frame generated by \code{compute_stats_aipw} and reshaped.
#' @param label.subgroup1 A label for subgroup1 (character).
#' @param label.subgroup2 A label for subgroup2 (character).
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param zero.line A logical indicating whether to include a zero line. Default TRUE.
#' @param arrows A logical indicating whether to include arrows. Default TRUE.
#' @param y.min A lower bound for the y-axis (numeric). Default -Inf.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -0.2.
#' @param p.ub An upper bound for the y-axis (numeric). Default 0.2.
#' @param y.lab A label for the y-axis (character). Default "PSA versus Human".
#' @param p.label A vector of two labels for the annotations (character). Default c("PSA harms", "PSA helps").
#'
#' @return A ggplot object.
#'
#'
#' @keywords internal
#' @useDynLib aihuman, .registration = TRUE
#'
vis_diff_ai <- function(df,
                        label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                        x.order = NULL,
                        zero.line = TRUE, arrows = TRUE, y.min = -Inf,
                        p.title = NULL, p.lb = -1, p.ub = 1, y.lab = "PSA versus Human", p.label = c("PSA worse", "PSA better")) {
  p <- df %>%
    filter(
      metric %in% c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion"),
      cov %in% c("Overall", label.subgroup1, label.subgroup2)
    ) %>%
    mutate(
      cov = factor(cov, levels = c("Overall", label.subgroup1, label.subgroup2)),
      X = factor(X, levels = x.order)
    ) %>%
    ggplot(aes(x = X, color = cov, shape = cov)) +
    geom_linerange(aes(ymin = pmax(est_lb - qnorm(1 - 0.05) * se_lb, y.min), ymax = est_ub + qnorm(1 - 0.05) * se_ub),
      position = position_dodge(.5), lineend = "round", linewidth = 1.2
    ) +
    geom_linerange(aes(ymin = pmax(est_lb, y.min), ymax = est_ub),
      lineend = "round", linewidth = 2.5,
      position = position_dodge(.5)
    )
  if (zero.line) {
    p <- p + geom_hline(yintercept = 0, lty = 2)
  }
  if ("est_NA" %in% colnames(df)) {
    p <- p + geom_pointrange(aes(y = est_NA, ymin = pmax(est_NA - 2 * se_NA, y.min), ymax = est_NA + 2 * se_NA), position = position_dodge(.5), lineend = "round")
  }

  p <- p + facet_grid(. ~ fct_relevel(metric, "Misclassification Rate"), space = "free_x", scales = "free") +
    ylab(y.lab) +
    scale_color_brewer("", type = "qual", palette = "Set1") +
    ggtitle(p.title) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_y_continuous(limits = c(p.lb, p.ub))

  if (arrows) {
    annotations <- data.frame(
      metric = c("Misclassification Rate", "Misclassification Rate"), # Assuming 'Misclassification Rate' is in the first facet
      x = c(1 - 0.2, 1 - 0.2),
      y.start = c(p.ub - 0.3 * p.ub, p.lb + 0.3 * p.ub),
      y.end = c(p.ub, p.lb),
      label = p.label
    )
    annotations_bottom <- annotations %>%
      filter(label == p.label[1]) # Assuming this is your bottom annotation

    annotations_top <- annotations %>%
      filter(label == p.label[2]) # Assuming this is your top annotation

    p <- p +
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

  return(p)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "est_NA", "est_lb", "est_ub", "se_NA", "se_lb", "se_ub"
  ))
}