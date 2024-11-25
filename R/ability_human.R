#' Visualize Difference in Risk (Human+AI v. Human)
#'
#' Visualize the the difference in risk between human+AI and human decision makers using AIPW estimators.
#' Generate a plot based on the overall agreement and subgroup-specific agreement.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
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
#' @param p.label A vector of two labels for the annotations (character). Default c("PSA harms", "PSA helps").
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_diff_human_aipw(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   nuis_funcs = nuis_func,
#'   true.pscore = rep(0.5, nrow(NCAdata)),
#'   l01 = 1,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male"),
#'   p.title = NULL, p.lb = -0.3, p.ub = 0.3
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
plot_diff_human_aipw <- function(Y, D, Z, V = NULL, l01 = 1,
                                 nuis_funcs = NULL, true.pscore = NULL,
                                 subgroup1, subgroup2,
                                 label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                                 x.order = NULL,
                                 p.title = NULL, p.lb = -1, p.ub = 1, y.lab = "Impact of PSA", p.label = c("PSA harms", "PSA helps")) {
  if (is.null(nuis_funcs)) {
    nuis_funcs <- compute_nuisance_functions(Y = Y, D = D, Z = Z, V = V)
  }
  if (!is.null(true.pscore)) {
    nuis_funcs$pscore <- true.pscore
  }
  diff_hwa <- compute_stats_aipw(Y = Y, D = D, Z = Z, nuis_funcs = nuis_funcs, X = NULL, l01 = 1) |>
    mutate(cov = "Overall", X = "Overall") |>
    relocate(Z_focal, Z_compare, X) |>
    bind_rows(compute_stats_aipw(Y = Y, D = D, Z = Z, nuis_funcs = nuis_funcs, X = subgroup1, l01 = 1) |> mutate(cov = label.subgroup1) |> relocate(cov)) |>
    bind_rows(compute_stats_aipw(Y = Y, D = D, Z = Z, nuis_funcs = nuis_funcs, X = subgroup2, l01 = 1) |> mutate(cov = label.subgroup2) |> relocate(cov))
  diff_hwa_vis <- diff_hwa |>
    pivot_longer(cols = -c("X", "cov", "Z_focal", "Z_compare")) |>
    mutate(
      metric = case_when(
        str_detect(name, "loss_diff") ~ "Misclassification Rate",
        str_detect(name, "fn_diff") ~ "False Negative Proportion",
        str_detect(name, "fp_diff") ~ "False Positive Proportion",
        TRUE ~ ""
      ),
      se = ifelse(str_detect(name, "_se"), "se", "est")
    ) |>
    filter(metric != "") |>
    select(-name) |>
    pivot_wider(values_from = "value", names_from = "se") |>
    mutate(cov = factor(cov, levels = c("Overall", label.subgroup1, label.subgroup2)))
  vis_diff_human(diff_hwa_vis, 
                 label.subgroup1 = label.subgroup1,
                 label.subgroup2 = label.subgroup2,
                 x.order = x.order,
                 p.title = p.title, p.lb = p.lb, p.ub = p.ub, y.lab = y.lab, p.label = p.label)
}

#' Visualize Difference in Risk (Human+AI v. Human)
#'
#' Visualize the the difference in risk between human+AI and human decision makers using difference-in-means estimators.
#' Generate a plot based on the overall agreement and subgroup-specific agreement.
#'
#' @param Y An observed outcome (binary: numeric vector of 0 or 1).
#' @param D An observed decision (binary: numeric vector of 0 or 1).
#' @param Z A treatment indicator (binary: numeric vector of 0 or 1).
#' @param l01 Ratio of the loss between false positives and false negatives. Default 1.
#' @param subgroup1 A pretreatment covariate used for subgroup analysis (vector).
#' @param subgroup2 A pretreatment covariate used for subgroup analysis (vector).
#' @param label.subgroup1 A label for subgroup1 (character). Default "Subgroup 1".
#' @param label.subgroup2 A label for subgroup2 (character). Default "Subgroup 2".
#' @param x.order An order for the x-axis (character vector). Default NULL.
#' @param p.title A title for the plot (character). Default NULL.
#' @param p.lb A lower bound for the y-axis (numeric). Default -1.
#' @param p.ub An upper bound for the y-axis (numeric). Default 1.
#' @param y.lab A label for the y-axis (character). Default "Impact of PSA".
#' @param p.label A vector of two labels for the annotations (character). Default c("PSA harms", "PSA helps").
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_diff_human(
#'   Y = NCAdata$Y,
#'   D = ifelse(NCAdata$D == 0, 0, 1),
#'   Z = NCAdata$Z,
#'   l01 = 1,
#'   subgroup1 = ifelse(NCAdata$White == 1, "White", "Non-white"),
#'   subgroup2 = ifelse(NCAdata$Sex == 1, "Male", "Female"),
#'   label.subgroup1 = "Race",
#'   label.subgroup2 = "Gender",
#'   x.order = c("Overall", "Non-white", "White", "Female", "Male"),
#'   p.title = NULL, p.lb = -0.3, p.ub = 0.3
#' )
#'
#' @useDynLib aihuman, .registration = TRUE
#' @export
plot_diff_human <- function(Y, D, Z, l01 = 1,
                            subgroup1, subgroup2,
                            label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                            x.order = NULL,
                            p.title = NULL, p.lb = -1, p.ub = 1, y.lab = "Impact of PSA", p.label = c("PSA harms", "PSA helps")) {
  diff_hwa <- compute_stats(Y = Y, D = D, Z = Z, X = NULL, l01 = 1) |>
    mutate(cov = "Overall", X = "Overall") |>
    relocate(Z_focal, Z_compare, X) |>
    bind_rows(compute_stats(Y = Y, D = D, Z = Z, X = subgroup1, l01 = 1) |> mutate(cov = label.subgroup1) |> relocate(cov)) |>
    bind_rows(compute_stats(Y = Y, D = D, Z = Z, X = subgroup2, l01 = 1) |> mutate(cov = label.subgroup2) |> relocate(cov))
  diff_hwa_vis <- diff_hwa |>
    pivot_longer(cols = -c("X", "cov", "Z_focal", "Z_compare")) |>
    mutate(
      metric = case_when(
        str_detect(name, "loss_diff") ~ "Misclassification Rate",
        str_detect(name, "fn_diff") ~ "False Negative Proportion",
        str_detect(name, "fp_diff") ~ "False Positive Proportion",
        TRUE ~ ""
      ),
      se = ifelse(str_detect(name, "_se"), "se", "est")
    ) |>
    filter(metric != "") |>
    select(-name) |>
    pivot_wider(values_from = "value", names_from = "se") |>
    mutate(cov = factor(cov, levels = c("Overall", label.subgroup1, label.subgroup2)))
  vis_diff_human(diff_hwa_vis, 
                 label.subgroup1 = label.subgroup1,
                 label.subgroup2 = label.subgroup2,
                 x.order = x.order,
                 p.title = p.title, p.lb = p.lb, p.ub = p.ub, y.lab = y.lab, p.label = p.label)
}


#' Visualize Risk (Human+AI v. Human; internal)
#'
#' Internal function to visualize the difference in risk between human+AI and human decision makers.
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
#'
#' @return A ggplot object.
#'
#' @importFrom forcats fct_relevel
#'
#'
#' @keywords internal
#' @useDynLib aihuman, .registration = TRUE
#'
vis_diff_human <- function(df, 
                           label.subgroup1 = "Subgroup 1", label.subgroup2 = "Subgroup 2",
                           x.order = NULL,
                           p.title = NULL, p.lb = -0.2, p.ub = 0.2, y.lab = "Impact of PSA", p.label = c("PSA harms", "PSA helps")) {
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

  df %>%
    filter(
      metric %in% c("Misclassification Rate", "False Negative Proportion", "False Positive Proportion"),
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
    facet_grid(. ~ fct_relevel(metric, "Misclassification Rate"), space = "free_x", scales = "free") +
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
    "metric"
  ))
}
